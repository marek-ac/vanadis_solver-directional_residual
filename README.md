# Vanadis Solver

[![CMAS 2026](https://img.shields.io/badge/CMAS_2026-Oral_Presentation-blue)](https://www.cmascenter.org/conference/2026/agenda.cfm)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22016556.svg)](https://doi.org/10.5281/zenodo.22016556)

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
