# Historical Vanadis source code

This directory preserves selected source-code milestones from the early development of the Vanadis atmospheric-transport model.

The files are retained as **historical software artifacts**. They are not part of the current supported Vanadis release and should not be treated as production software. Wherever possible, the original source bytes should remain unchanged. Compiler-specific constructs, hard-coded research cases, comments, experimental branches and historical implementation limitations are part of the record.

If a historical program is restored for a modern compiler, the modified version should be stored separately (for example under a `restored/` subdirectory) while the historical source remains unchanged.

## Relationship to the 1997 publications

The preserved source code can now be compared with two distinct publications from 1997 documenting the atmospheric-transport FEM work.

The Chodorski-Pietrzyk paper:

> Marek Chodorski, Maciej Pietrzyk, *Propozycja zastosowania metody elementów skończonych do symulacji rozprzestrzeniania się zanieczyszczeń w atmosferze*, Ochrona Powietrza i Problemy Odpadów, 31(3), 90-93, 1997.

presented a two-dimensional steady FEM implementation and discussed the numerical difficulty of convection-dominated transport.

A separate POL-IMIS 1997 proceedings paper:

> Marek Chodorski, Zbigniew Malinowski, Stanisław Słupek, Andrzej Buczek, *Zastosowanie metody elementów skończonych do modelu rozprzestrzeniania się zanieczyszczeń w atmosferze*, in: *Ocena wielkości imisji zanieczyszczeń powietrza*, POL-IMIS 1997, Szklarska Poręba, 19-22 June 1997, pp. 13-22.

already documented a **three-dimensional atmospheric-dispersion FEM model using eight-node hexahedral elements**. The paper also presented calculation variants involving:

- a transient concentration field 100 s after the start of emission;
- lateral convection;
- an obstacle placed in the pollutant plume;
- a concentration-dependent removal coefficient of the form

$$P=a_1+a_2 S+a_3 S^2.$$

The POL-IMIS paper provides published evidence for these model capabilities. The archived source code in this directory provides complementary software evidence and, in some cases, exposes implementation details not stated in the proceedings paper, such as the Element-by-Element operator structure.

The proceedings contents are independently listed by PZITS:

https://www.pzits.not.pl/docs/ksiazki/pol_1997.html

---

## Development overview

| Period | Directory / file | Historical role | Main characteristics |
|---|---|---|---|
| 1995-1996 | `1995_2d_pascal/ART_SPEE.PAS` | early 2-D atmospheric-transport implementation | stationary FEM, triangular elements, convection-diffusion-reaction, Peclet-dependent upwind weighting, nonsymmetric band system |
| 1997 | `1997_3d_stationary_nonlinear/GR_P(S).FOR` | 3-D stationary nonlinear prototype | HEX8 FEM, concentration-dependent reaction/removal coefficient, fixed-point iteration, EBE operator, iterative use of `A` and `A^T` |
| 1997 | `1997_3d_transient_ebe/fm_t_ok.for` | early 3-D transient EBE branch | HEX8 FEM, transient mass matrix, implicit time integration, Element-by-Element operator, iterative nonsymmetric solve |
| 1997 | `1997_3d_transient_lu/LU_T_OK.FOR` | alternative transient implementation | HEX8 FEM, transient mass matrix, assembled global band matrix, LINPACK band-LU solution |
| 1990s | `ftn90_visualization/1.f90` | interactive concentration-field viewer | Salford FTN90, VGA graphics, bilinear interpolation, colour mapping, mouse readout |
| 1998 | `1998_defence_animation/` | transient-animation software used for the doctoral-defence presentation | PCX frame generation, frame preloading, VGA playback, simulation-time display |

The dates above describe the historical development period and archival context. Historical archive timestamps are treated as provenance evidence rather than cryptographic proof of authorship dates.

The archived copy corresponding to the transient EBE source has an internal archive timestamp of **9 August 1997**. This is useful provenance evidence, but it should not be interpreted as cryptographic proof of the authorship date.

---

## `1995_2d_pascal/ART_SPEE.PAS`

An early DOS/Pascal implementation of a two-dimensional stationary finite-element model for atmospheric transport.

The code includes:

- triangular finite elements;
- directional diffusion coefficients (`kx`, `ky`);
- advective velocity components (`vx`, `vy`);
- a volumetric source term (`qv`);
- a first-order reaction/removal coefficient (`pzanik`);
- experiments with first- and third-kind boundary conditions;
- compact band storage for a nonsymmetric global system;
- a dedicated asymmetric Gaussian-elimination routine;
- optional Peclet-dependent upwind weighting.

The upwind parameter uses a hyperbolic-function expression equivalent to the classical form

$$\coth(Pe)-\frac{1}{Pe}.$$

This code is a mathematical precursor of the later Vanadis stabilisation work, but it should **not** be described as the modern Directional Residual formulation.

### Historical significance

The program documents the 2-D stage of the Vanadis development line: a working convection-diffusion-reaction FEM implementation with explicit treatment of advection-dominated transport and nonsymmetric systems before the transition to 3-D HEX8 and Element-by-Element methods.

---

## `1997_3d_stationary_nonlinear/GR_P(S).FOR`

A three-dimensional stationary transport prototype on a structured HEX8 mesh.

The code contains:

- 8-node hexahedral finite elements;
- three directional diffusion coefficients;
- three velocity components;
- weighted/upwind test functions in the three element directions;
- local 8 x 8 element matrices stored consecutively in the `bb` array;
- Element-by-Element implementations of both `A x` and `A^T x`;
- an iterative nonsymmetric linear solver;
- a nonlinear concentration-dependent reaction/removal coefficient.

For the archived test case, the nonlinear coefficient is evaluated elementwise as

$$P(\bar S_e)=4\times10^{-6}+2\times10^{-4}|\bar S_e|+0.1|\bar S_e|^2.$$

The nonlinear problem is solved by repeated assembly, linear solution and concentration update until the change between successive nodal solutions falls below a prescribed tolerance. In modern terminology this is a **Picard/fixed-point-type iteration**.

### Published 1997 cross-check

The separate POL-IMIS 1997 proceedings paper independently documents a concentration-dependent removal experiment. Its fifth calculation variant uses the general relationship

$$P=a_1+a_2 S+a_3 S^2.$$

The preserved source therefore provides an implementation-level counterpart to a nonlinear model feature that was already described publicly in 1997. The exact numerical coefficients in this archived source belong to the saved research case and should not be assumed to be identical to every calculation reported in the proceedings paper unless the corresponding case data are independently matched.

### Historical significance

This source documents several ideas that later remained important in Vanadis: 3-D HEX8 transport, element-local matrix storage, EBE application of the global operator, use of both `A` and `A^T`, and concentration-dependent physics.

### Historical implementation notes

The original source is intentionally preserved unchanged. Two details should be noted when reproducing the program:

1. Expressions such as `3/2*ALFA...` use integer division in standard Fortran, so `3/2` evaluates to `1` rather than `1.5`. A restored version should use a real constant if `1.5` is intended.
2. The initial solution vector is controlled by run-specific commented code. A restored version should initialise the starting vector explicitly.

The original exponential implementation of the hyperbolic upwind parameter also predates the numerical safeguards used in modern Vanadis.

---

## `1997_3d_transient_ebe/fm_t_ok.for`

An early transient 3-D HEX8 implementation using Element-by-Element operator application instead of a conventional assembled global sparse matrix for the iterative solve.

Principal features include:

- 3-D structured HEX8 mesh;
- transient mass matrix `CMA_N`;
- implicit time integration with the historical 2/3-1/3 weighting;
- local 8 x 8 time-step matrices stored element-by-element in `bb`;
- explicit EBE implementations of `A x` and `A^T x`;
- iterative solution of a nonsymmetric system;
- test logic for an obstacle and time-dependent changes of source/flow conditions;
- output of successive concentration fields for post-processing and animation.

### Published 1997 cross-check

The POL-IMIS proceedings paper published earlier in 1997 independently documents that the atmospheric-dispersion model had already reached a **3-D HEX8** stage. It also includes a transient calculation showing the concentration field **100 s after the start of emission**, as well as variants with lateral convection and an obstacle.

The archived transient source provides the corresponding implementation evidence for the 3-D transient development branch and additionally documents the **Element-by-Element** architecture and explicit `A x` / `A^T x` operators.

The distinction is important: the proceedings paper supports the claims of 3-D HEX8 modelling and a transient calculation, while **the EBE implementation is established by the preserved source code rather than by an explicit statement in that paper**.

### Historical significance

This file is one of the clearest early examples of the architecture that later characterised Vanadis: 3-D HEX8 finite elements combined with Element-by-Element operator application and an iterative nonsymmetric-system strategy.

The archived copy corresponding to this source has an internal archive timestamp of **9 August 1997**, providing a dated provenance point for this implementation stage.

### Historical implementation note

The original declares `character*12 wynik_t(100)` but indexes it with the simulation-time counter, which can exceed 100 in the archived configuration. A restored version should use a scalar filename or consistently sized storage. The historical source should remain unchanged.

---

## `1997_3d_transient_lu/LU_T_OK.FOR`

A transient three-dimensional HEX8 implementation using an assembled global band matrix and LU factorisation.

The source contains the transient mass matrix `CMA_N` and an implicit time-step formulation algebraically equivalent to the historical 2/3-1/3 weighting used in the EBE branch.

The global matrix is stored in band form. Classic LINPACK-style band routines are used for factorisation and solution; the factorisation is reused while the matrix remains constant.

### Historical significance

This source documents an alternative implementation strategy developed alongside the EBE work. It provides a useful comparison between a conventional assembled band-matrix approach and the Element-by-Element architecture that later became central to Vanadis.

### External numerical-library code

The file contains LINPACK/BLAS-derived routines. Their original comments and attribution should be preserved and they should not be presented as Vanadis-authored algorithms.

A modern smoke test with GNU Fortran and legacy fixed-form compatibility options showed that this historical program can still be compiled and run through its programmed test stop without modifying the source.

---

## `ftn90_visualization/1.f90`

A historical post-processing utility written for the Salford FTN90 environment.

The program reads a two-dimensional concentration plane, constructs quadrilateral cell connectivity and interpolates nodal values onto a raster using four-node bilinear shape functions. It renders the concentration field with a VGA colour palette and supports interactive mouse-based coordinate/concentration readout.

### Historical significance

This program documents the visualisation side of the early Vanadis workflow: solver results were converted into interactive concentration maps using dedicated Fortran graphics software.

### Portability

The program depends on Salford FTN90/DOS graphics extensions such as `vga@`, `set_pixel@`, `draw_text@` and mouse routines. It is therefore not expected to compile unchanged with a standard modern Fortran compiler.

---

## `1998_defence_animation/`

This directory preserves the FTN90 software used to prepare and display animated transient Vanadis concentration fields, together with a later screen capture of the surviving program.

### `xb.f90` — frame generator

`xb.f90` reads concentration results for successive simulation times, interpolates the solution on a two-dimensional plane using four-node bilinear shape functions, assigns concentration ranges to a VGA colour palette and writes the resulting frames as PCX images.

The generated filenames follow the simulation time, for example:

```text
X1_0005.pcx
X1_0010.pcx
...
X1_1980.pcx
```

The program uses Salford FTN90 graphics routines including `screen_block_to_pcx@`.

### `run1b.f90` — animation player

`run1b.f90` loads the generated PCX frames into memory before playback and then restores them sequentially to the VGA screen. The simulation time is updated on screen for each displayed frame.

The program uses a `sleep@(0.1)` delay between frames. Preloading the image sequence avoided repeated disk I/O during playback and allowed the animation to run smoothly in the original DOS environment.

### `vanadis1.webm` — preserved modern screen capture

`vanadis1.webm` is **not an original video recording from 1998**. It is a later screen capture of the historical FTN90/DOS animation running inside a virtual machine.

The original animation was executed directly under DOS on contemporary hardware and was displayed smoothly. The WebM file is included to document the appearance and behaviour of the surviving historical software.

An animation of transient Vanadis results was demonstrated during Marek Chodorski's public doctoral defence in 1998. The surviving animation software, transient solver source and defence documentation together preserve this stage of the Vanadis development history.

The newly recovered POL-IMIS 1997 proceedings paper shows that a transient calculation had already been published the previous year. The 1998 defence therefore provides an additional dated public demonstration of the transient development line, including its animated visualisation, rather than the first evidence that transient calculations existed.

### Portability

These programs depend on Salford FTN90/DOS graphics functions such as `vga@`, `screen_block_to_pcx@`, `pcx_to_screen_block@`, `restore_screen_block@`, `draw_text@` and related routines. Porting to a modern compiler would require replacement of the graphics layer.

A minor historical typo (`erro_code` instead of `error_code`) exists in one error-handling branch of `xb.f90`; the original source should remain unchanged.

---

## Recommended repository layout

```text
docs/history/source/
├── README.md
├── SHA256SUMS.txt
├── 1995_2d_pascal/
│   └── ART_SPEE.PAS
├── 1997_3d_stationary_nonlinear/
│   └── GR_P(S).FOR
├── 1997_3d_transient_ebe/
│   └── fm_t_ok.for
├── 1997_3d_transient_lu/
│   └── LU_T_OK.FOR
├── ftn90_visualization/
│   └── 1.f90
└── 1998_defence_animation/
    ├── xb.f90
    ├── run1b.f90
    └── vanadis1.webm
```

The filenames shown above follow the current repository paths. Historical/original filenames and any portability-related renaming are documented in the README files of the corresponding subdirectories.

For any future reconstruction, a useful pattern is:

```text
<case>/
├── original/     # historical source preserved unchanged
├── restored/     # modernised/compilable derivative
└── README.md     # changes, compiler, inputs and reproduced results
```

## Integrity

`SHA256SUMS.txt` contains SHA-256 digests of the historical source files and the preserved animation capture. The checksums were calculated from the uploaded historical copies before any source modification.

## Interpretation

These files should not be used to claim that every feature of the current Vanadis formulation was already present in the earliest code. They document a **development sequence**: 2-D FEM transport and Peclet-dependent upwinding, followed by 3-D HEX8 formulations, nonlinear experiments, transient integration, Element-by-Element operator application, alternative LU solution, and dedicated DOS/FTN90 visualisation and animation tools.

The 1997 record is now supported by two complementary evidence types:

- **published evidence** — the POL-IMIS proceedings paper documents a 3-D HEX8 atmospheric-dispersion model, a transient calculation, lateral convection, an obstacle case and concentration-dependent removal;
- **software evidence** — the preserved 1997 Fortran sources document how the 3-D, transient, nonlinear and Element-by-Element branches were actually implemented.

The two should be kept distinct. In particular, the published POL-IMIS paper does not by itself establish that the reported calculations used the EBE branch; the EBE architecture is directly evidenced by the preserved source code.

Their combined value is historical and technical: they show how the Vanadis numerical architecture evolved through published scientific work and working research software rather than appearing as a single modern rewrite.
