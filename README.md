# Vanadis Solver

[![CMAS 2026](https://img.shields.io/badge/CMAS_2026-Oral_Presentation-blue)](https://www.cmascenter.org/conference/2026/agenda.cfm)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21352013.svg)](https://doi.org/10.5281/zenodo.21352013)

**Project website:** [marek-ac.meri.pl](https://marek-ac.meri.pl)

---

## Overview

**Vanadis 3D** is a three-dimensional finite-element model for atmospheric pollutant transport and dispersion. The numerical formulation is based on **Directional Residual (DR) stabilization** and uses 8-node hexahedral finite elements (HEX8).

Vanadis is built around an **Element-by-Element (EbE)** architecture. Element matrices are retained locally and used directly by the iterative solver, avoiding the construction of a conventional global sparse matrix. CPU matrix assembly is parallelized with **OpenMP**, while the linear solution stage can be executed on the CPU or accelerated on NVIDIA GPUs using **CUDA**.

The CUDA implementation supports both atomic Element-by-Element operations and a graph-coloring variant for conflict-free processing of independent element groups. The solver uses DBCG with lightweight diagonal preconditioning.

Vanadis **v2026.3.1** is the current development line. It extends the model toward **terrain-following meshes**, spatially and temporally variable physical fields, **27-point HEX8 volume integration**, **9-point Q4 boundary integration**, and full mass-balance diagnostics.

Vanadis **v2026.3.0** introduced a **nonlinear concentration-dependent reaction/decay coefficient P(S)**, solved by **Picard iteration**, together with generalized **time-dependent source handling Q(t)**.

Vanadis **v2026.2.1** remains the final reference release for the CMAS 2026 study.

---

## Key features

- three-dimensional atmospheric advection-diffusion-reaction modeling
- finite-element discretization with HEX8 elements
- Directional Residual (DR) stabilization
- Element-by-Element matrix operations without a conventional global sparse matrix
- OpenMP parallel CPU assembly
- CPU and CUDA iterative DBCG solution paths
- diagonal/Jacobi preconditioning
- CUDA atomic and graph-coloring Element-by-Element variants
- implicit `theta = 2/3` transient time integration
- Dirichlet and flux-type boundary conditions
- nonlinear concentration-dependent reaction/decay handling P(S)
- Picard iteration for the nonlinear problem
- generalized time-dependent source handling Q(t)
- locally evaluated physical fields `v(x,y,z,t)`, `K(x,y,z,t)`, and `alpha(x,y,z,t)`
- terrain-following HEX8 meshes
- terrain-consistent test velocity fields
- 27-point (`3 x 3 x 3`) Gauss integration for HEX8 volume terms
- 9-point (`3 x 3`) Gauss integration for Q4 boundary terms
- full mass-balance diagnostics
- buffered result output and post-processing support

---

## Numerical architecture at a glance

Vanadis separates the **physical description of the atmospheric problem**, the **local FEM/Directional Residual formulation**, and the **linear-algebra backend**. For problems that remain within the present three-dimensional advection-diffusion-reaction model, the numerical core is intended to remain unchanged while the case-specific geometry, boundary conditions, sources, reaction parameters, and transport fields are supplied through the corresponding data and field routines.

In compact form, the v2026.3.1 computation path is:

**physical data -> HEX8 / 27-point Gauss integration -> local Jacobian and transport fields -> Directional Residual stabilization -> `theta = 2/3` transient FEM -> OpenMP element assembly -> Element-by-Element DBCG/Jacobi -> CPU or CUDA solution -> mass-balance diagnostics**

### Physical quantities and their spatial resolution

| Quantity | Representation in v2026.3.1 |
| --- | --- |
| mesh geometry | 8-node hexahedral elements (HEX8) |
| volume integration | 27 Gauss points (`3 x 3 x 3`) per HEX8 element |
| boundary integration | 9 Gauss points (`3 x 3`) per Q4 face |
| velocity `v(x,y,z,t)` | evaluated locally at physical Gauss points |
| diffusion `K(x,y,z,t)` / diagonal `Lambda(x,y,z,t)` | evaluated locally at physical Gauss points |
| Robin/deposition field `alpha(x,y,z,t)` | evaluated locally at boundary Gauss points |
| source `Q(t)` | assigned at element level by `source_term`; the routine also receives the element number |
| nonlinear reaction/decay `P(S)` | evaluated from an element-representative concentration and updated by Picard iteration |
| Directional Residual parameter `tau_DR` | evaluated locally at each volume Gauss point |

The use of 27 Gauss points **does not mean that Vanadis uses a 27-node element**. The interpolation remains HEX8; the `3 x 3 x 3` rule provides higher-order numerical integration of geometry-dependent and spatially varying terms.

### Directional Residual stabilization

At every volume Gauss point, Vanadis evaluates the physical position, the local velocity and diffusion fields, and the HEX8 Jacobian

```text
J = d(x,y,z) / d(xi,eta,zeta).
```

For non-zero velocity,

```text
e = v / |v|
```

defines the local streamline direction. The characteristic element length used by DR is obtained from the actual element mapping:

```text
h_stream = 2 / || J^(-1) e ||.
```

For the diagonal diffusion tensor

```text
Lambda = diag(lambda_1, lambda_2, lambda_3),
```

the diffusion acting in the streamline direction is

```text
lambda_s = e^T Lambda e.
```

The local Peclet number and DR parameter are then

```text
Pe = |v| h_stream / (2 lambda_s)

tau_DR = h_stream / (2 |v|)
         * (coth(Pe) - 1/Pe).
```

The implementation includes numerically safe small- and large-Peclet limits and the pure-advection limit

```text
lambda_s -> 0  =>  tau_DR -> h_stream / (2 |v|).
```

The perturbed test weight is

```text
W_i = N_i + tau_DR * (v . grad N_i).
```

Consequently, the stabilization responds automatically to:

- the local velocity direction,
- the local size and distortion of the HEX8 element,
- the element orientation relative to the flow,
- anisotropic diagonal diffusion,
- spatially and temporally varying `v` and `K`.

For isotropic diffusion, `Lambda = K I` and therefore `lambda_s = K`.

### Weak-form structure

The diffusion part intentionally remains the standard Galerkin term

```text
grad(N_i)^T Lambda grad(N_j).
```

The DR perturbation enters the transient, advective, reaction/source, and applicable boundary weighting through `W_i`. In simplified volume form, the element contributions are of the form

```text
H_ij =
  integral [
      grad(N_i)^T Lambda grad(N_j)
    + W_i (v . grad(N_j))
    + W_i P N_j
  ] dOmega

M_ij =
  integral W_i N_j dOmega

F_i =
  integral W_i q_e dOmega,
```

with the corresponding boundary contributions added separately.

This distinction is deliberate: Vanadis does **not** replace the Galerkin diffusion gradient by `grad(W_i)`. The DR term is used as directional residual weighting, so derivatives of `W_i` or `tau_DR` are not required by the present formulation.

### Transient formulation

Vanadis v2026.3.1 uses an implicit `theta = 2/3` scheme. Because the velocity, diffusion, boundary fields, and therefore the DR-weighted operator may vary with time, the element operator `H` and DR-weighted mass matrix `M` are evaluated at both time levels.

The assembled time-discrete equation is

```text
( 2 H^(n+1) + 3/dt M^(n+1) ) C^(n+1)
 =
( 3/dt M^n - H^n ) C^n
+ 2 F^(n+1) + F^n.
```

For nonlinear `P(S)`, the new-time-level problem is solved by Picard iteration while `C^n` remains fixed.

### Parallel and GPU architecture

Element matrices are stored and processed in **Element-by-Element (EbE)** form rather than assembled into a conventional global sparse matrix. For HEX8 this gives regular local `8 x 8` blocks.

The computational work is divided as follows:

- **OpenMP / CPU** — parallel construction of element matrices and vectors,
- **EbE DBCG with diagonal/Jacobi preconditioning** — iterative linear solution,
- **CUDA / GPU** — accelerated EbE matrix-vector operations and DBCG solution,
- **atomic or graph-coloring CUDA variants** — alternative handling of shared-node updates.

When the EbE data and solver vectors fit in GPU memory, the linear iterations can be performed on a single GPU without repeatedly transferring a conventional global sparse matrix between CPU and GPU.

### Design consequence

For atmospheric transport problems covered by the present model, a new application normally changes the **data and physical-field definitions**, not the FEM/DR core. Typical case-specific inputs include:

- mesh and terrain geometry,
- boundary-condition assignment,
- element source histories,
- reaction/decay parameters,
- velocity fields,
- diffusion fields,
- deposition/Robin fields,
- time-step and simulation parameters.

The same local HEX8/DR formulation, transient assembly, OpenMP parallelization, and EbE CPU/CUDA solver can then be reused without redesigning the numerical method.

---

## Project structure

- `src_v2026_1_1/` — historical v2026.1.1 source snapshot, supporting files, sample results, and corrections specific to that version
- `src_v2026_2_1/` — source snapshot corresponding to the CMAS 2026 reference release
- `src_v2026_3_0/` — v2026.3.0 released development snapshot with nonlinear `P(S)` and generalized `Q(t)`
- `src_v2026_3_1/` — current v2026.3.1 development baseline, including terrain-following meshes, local physical fields, higher-order Gauss integration, full mass-balance diagnostics, documentation, and validation cases
- `docs/v2026_3_0/` — technical documentation and architectural comparison material for v2026.3.0
- `source_code.txt` — guide to key source files and release structure
- **Releases** — complete distributable archives for published Vanadis versions

The source directories in `main` are intended for inspection, development tracking, documentation, and validation. Complete packaged distributions for published versions are provided through GitHub Releases.

---

## Download / Installation

### Fastest way to run Vanadis

The easiest way to run Vanadis is to download a complete packaged distribution from the **Releases** section. Release packages contain the solver sources, configuration files, documentation, input datasets, and sample results.

For inspection of the latest development work, the repository also contains versioned source snapshots, including the current `src_v2026_3_1/` development baseline.

### CMAS 2026 reference release — v2026.2.1

**Directional_Residual_stabilization_Vanadis.zip**  
https://github.com/marek-ac/vanadis_solver-directional_residual/releases/tag/v2026.2.1

Vanadis v2026.2.1 is the final reference release for the CMAS 2026 study. Numerical results reported in the CMAS 2026 extended abstract and presentation were obtained using this and earlier compatible Vanadis releases. Later releases and development lines may contain additional model developments that were not part of the CMAS 2026 study.

### Released development snapshot — v2026.3.0

Vanadis **v2026.3.0** extends the model with a **nonlinear concentration-dependent reaction/decay coefficient P(S)** and generalized **time-dependent source handling Q(t)**. The nonlinear problem is solved using **Picard iteration** within each time step, while preserving the existing Directional Residual stabilization, Element-by-Element formulation, OpenMP CPU assembly, and CPU/CUDA iterative solution architecture.

The source-term formulation supports constant, switched, pulsed, and otherwise time-varying emission profiles.

**Release:**  
https://github.com/marek-ac/vanadis_solver-directional_residual/releases/tag/v2026.3.0

This development was introduced after the CMAS 2026 reference version and was **not part of the CMAS 2026 study**.

The release ZIP archive includes:

- complete directory structure with Vanadis source files
- example configuration files
- documentation
- input datasets
- sample results

### Current development baseline — v2026.3.1

The current **v2026.3.1 development baseline** is available in:

`src_v2026_3_1/`

The v3.1 line extends Vanadis toward:

- terrain-following meshes
- spatially and temporally variable physical fields
- terrain-consistent test velocity fields
- 27-point HEX8 volume integration
- 9-point Q4 boundary integration
- full mass-balance diagnostics
- validation on distorted terrain-following meshes

v2026.3.1 is currently a **development version**, not yet the final v3.1 release.

---

## Documentation

### v2026.3.0

Detailed technical documentation and architectural comparison materials for Vanadis 3D v2026.3.0 are available directly in the repository.

The technical description is based directly on the Fortran and CUDA source code of Vanadis v2026.3.0.

#### Technical description

- [English](docs/v2026_3_0/Vanadis_3D_v2026.3.0_full_technical_description_EN.pdf)

#### Model architecture comparison

- [English](docs/v2026_3_0/Vanadis_model_architecture_comparison_EN.pdf)

The architecture comparison discusses Vanadis in relation to Fluidity-Atmosphere, MFEM, CMAQ and FLEXPART. It is intended as a comparison of numerical and software architecture, not as a ranking of scientific accuracy or model capability.

### v2026.3.1 development documentation

Documentation for the current v2026.3.1 development baseline is stored in:

`src_v2026_3_1/DOC/`

It includes guides for:

- input data and model tailoring
- full mass-balance analysis
- terrain-following geometry
- terrain-consistent and externally supplied velocity fields

The current v3.1 baseline also includes reproducible validation material under:

`src_v2026_3_1/validation/`

---

## Community development and independent forks

Vanadis is intended as a **numerical transport engine** that research groups can use as the basis for independent domain-specific models and forks.

Users are encouraged to download a complete release, study the reference implementation, and develop specialized variants in their own repositories. Examples may include extensions for particular atmospheric processes, source types, deposition mechanisms, particle classes, chemistry, or application domains.

The main Vanadis repository is maintained as the **reference implementation, documentation source, and release archive**.

Vanadis follows a **maintainer-authored development model**. External source-code contributions and pull requests are not incorporated into the reference implementation. Community discussions, bug reports, validation results, scientific ideas, feature proposals, and interoperability work are welcome. Ideas selected for the reference Vanadis line are independently implemented and validated by the maintainer.

Independent forks are encouraged and may evolve separately from the reference Vanadis implementation.

GitHub Discussions may be used to share ideas, results, extensions, validation experience, and interoperability proposals with the wider Vanadis community.

---

## Historical development

Vanadis has a documented numerical-development history dating back to the 1990s, including published 3-D FEM atmospheric-transport work and preserved source code from 1997.

See [Historical publications and project history](docs/history/README.md).

---

## CMAS 2026 reference

Vanadis v2026.2.1 is the reference software release associated with the accepted oral presentation at the **25th Annual CMAS Conference, Chapel Hill, NC, October 2026**.

Conference agenda:  
https://www.cmascenter.org/conference/2026/agenda.cfm

Later Vanadis releases and development lines, including v2026.3.0 and v2026.3.1, document continued development after the CMAS 2026 reference version.

---

## Licensing

Vanadis and the software implementation of Directional-Residual Stabilization are distributed under a dual-licensing model:

### 1. Open-source license: GPL-3.0

The open-source edition of Vanadis and its implementation of Directional-Residual Stabilization is released under the  
GNU General Public License, Version 3 (GPL-3.0).

This license allows:

- academic and research use
- educational use
- commercial use under GPL-3.0
- independent open-source development and forks
- modification and redistribution subject to GPL-3.0 requirements

If software incorporating or linking to Vanadis or its implementation of Directional-Residual Stabilization is distributed, the resulting combined or derivative work must comply with the applicable GPL-3.0 requirements, including provision of the corresponding source code where required.

The GPL-3.0 license therefore preserves copyleft obligations for distributed derivative and combined works.

See: **[GNU GPL-3.0](https://www.gnu.org/licenses/gpl-3.0.html)**

---

### 2. Commercial license

A commercial license is available for organizations that wish to:

- use Vanadis or its implementation of Directional-Residual Stabilization in closed-source or proprietary software
- integrate Vanadis or its implementation of Directional-Residual Stabilization into commercial products without GPL-3.0 copyleft obligations
- distribute proprietary derivative or combined works
- obtain professional support or custom development

The commercial license grants, subject to the applicable commercial license agreement:

- rights to use Vanadis or its implementation of Directional-Residual Stabilization in proprietary applications
- the ability to keep derivative or combined works closed-source
- alternative licensing terms that do not impose GPL-3.0 copyleft requirements on the commercially licensed use
- optional support and consulting

Commercial licensing inquiries:  
**Marek Chodorski — marek_ac@wp.pl**

See: `LICENSE_COMMERCIAL.txt`

---

### 3. License files included in this repository

- `LICENSE_GPL3.txt` — GPL-3.0 notice and link to the full license text
- `LICENSE_COMMERCIAL.txt` — terms and contact information for commercial licensing

---

### 4. Summary

Vanadis and its implementation of Directional-Residual Stabilization may be used, modified, and distributed under GPL-3.0, including for commercial purposes, provided that the GPL-3.0 requirements are met.

Users who wish to incorporate Vanadis or its implementation of Directional-Residual Stabilization into proprietary or closed-source software without complying with GPL-3.0 copyleft requirements must obtain a commercial license.

This dual-licensing model protects the intellectual property embodied in Vanadis and its DR implementation while enabling both open scientific collaboration and commercial deployment.
