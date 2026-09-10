# Vanadis Solver

[![CMAS 2026](https://img.shields.io/badge/CMAS_2026-Oral_Presentation-blue)](https://www.cmascenter.org/conference/2026/agenda.cfm)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21352013.svg)](https://doi.org/10.5281/zenodo.21352013)

**Project website:** [marek-ac.meri.pl](https://marek-ac.meri.pl)

---

## Overview

**Vanadis 3D** is a three-dimensional finite-element model for atmospheric pollutant transport and dispersion. The numerical formulation is based on **Directional Residual (DR) stabilization** and uses 8-node hexahedral finite elements (HEX8).

Vanadis is built around an **Element-by-Element (EbE)** architecture. Element matrices are retained locally and used directly by the iterative solver, avoiding the construction of a conventional global sparse matrix. CPU matrix assembly is parallelized with **OpenMP**, while the linear solution stage can be executed on the CPU or accelerated on NVIDIA GPUs using **CUDA**.

The CUDA implementation supports both atomic Element-by-Element operations and a graph-coloring variant for conflict-free processing of independent element groups. The solver uses DBCG with lightweight diagonal preconditioning.

Vanadis **v2026.3.0** extends the model with a **nonlinear concentration-dependent reaction/decay coefficient P(S)**, solved by **Picard iteration**, together with generalized **time-dependent source handling Q(t)**.

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
- implicit transient time integration
- Dirichlet and flux-type boundary conditions
- nonlinear concentration-dependent reaction/decay handling P(S)
- Picard iteration for the nonlinear problem
- generalized time-dependent source handling Q(t)
- buffered result output and post-processing support

---

## Project structure

- `src_v2026_1_1/` — historical v2026.1.1 source snapshot, supporting files, sample results, and corrections specific to that version
- `src_v2026_2_1/` — source snapshot corresponding to the CMAS 2026 reference release
- `src_v2026_3_0/` — selected key source files for the current v2026.3.0 development release
- `docs/v2026_3_0/` — technical documentation and architectural comparison material for v2026.3.0
- `source_code.txt` — guide to key source files and release structure
- **Releases** — complete distributable archives for published Vanadis versions

The source directories in `main` are intended for inspection and documentation. Complete packaged distributions are provided through GitHub Releases.

---

## Download / Installation

To run the Vanadis solver, **cloning the repository is not sufficient**.

The fully packaged, ready-to-run version must be downloaded from the Releases section:

**Directional_Residual_stabilization_Vanadis.zip**  
https://github.com/marek-ac/vanadis_solver-directional_residual/releases/tag/v2026.2.1

Vanadis v2026.2.1 is the final reference release for the CMAS 2026 study. Numerical results reported in the CMAS 2026 extended abstract and presentation were obtained using this and earlier compatible Vanadis releases. Later releases may contain additional model developments that were not part of the CMAS 2026 study.

### Current development release — v2026.3.0

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

---

## Documentation — v2026.3.0

Detailed technical documentation and architectural comparison materials for Vanadis 3D v2026.3.0 are available directly in the repository.

The technical description is based directly on the Fortran and CUDA source code of Vanadis v2026.3.0.

### Technical description

- [English](docs/v2026_3_0/Vanadis_3D_v2026.3.0_full_technical_description_EN.pdf)

### Model architecture comparison

- [English](docs/v2026_3_0/Vanadis_model_architecture_comparison_EN.pdf)

The architecture comparison discusses Vanadis in relation to Fluidity-Atmosphere, MFEM, CMAQ and FLEXPART. It is intended as a comparison of numerical and software architecture, not as a ranking of scientific accuracy or model capability.

---

## Historical development

Vanadis has a documented numerical-development history dating back to the 1990s, including published 3-D FEM atmospheric-transport work and preserved source code from 1997.

See [Historical publications and project history](docs/history/README.md).

---
## CMAS 2026 reference

Vanadis v2026.2.1 is the reference software release associated with the accepted oral presentation at the **25th Annual CMAS Conference, Chapel Hill, NC, October 2026**.

Conference agenda:  
https://www.cmascenter.org/conference/2026/agenda.cfm

Later Vanadis releases, including v2026.3.0, document continued development after the CMAS 2026 reference version.

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
- open-source development and contributions
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
