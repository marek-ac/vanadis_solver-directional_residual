# Early 3D transient Element-by-Element Vanadis source

This directory preserves an early three-dimensional transient Fortran implementation from the development line that later evolved into the modern **Vanadis** solver.

The source is retained as a historical software artifact and is **not part of the current Vanadis production release**.

## Source file

- `FM_T_OK.FOR`

Original historical filename:

```text
FM_(T)OK.FOR
```

The repository filename has been simplified for portability. The file contents should remain unchanged.

## Historical context

The preserved source belongs to the 1997 stage of Vanadis development.

The archived copy corresponding to this source has an internal archive timestamp of **9 August 1997**. This timestamp is treated as archival provenance rather than cryptographic proof of the date of authorship.

The code is historically important because it already contains a working combination of:

- three-dimensional transport,
- eight-node hexahedral finite elements,
- transient time integration,
- Element-by-Element matrix storage,
- direct `A v` and `A^T v` operators,
- iterative solution of the nonsymmetric system,
- time-dependent model parameters,
- local obstacle removal,
- and generation of concentration fields for successive simulation times.

## Computational grid

The preserved test configuration uses:

```fortran
parameter (nax=7,nay=14,naz=23)
```

which corresponds to:

```text
2024 HEX8 elements
2254 nodes
```

The connectivity array stores eight node numbers per element, and the Element-by-Element operator data are held as 64 coefficients per element:

```fortran
dimension iSS((nax-1)*(nay-1)*(naz-1),8)
dimension bb((nax-1)*(nay-1)*(naz-1)*64)
```

## Transient formulation

The saved test uses a time increment:

```fortran
i_d_czas=5
```

and repeatedly calls the transient solution routine until the simulated time passes 2000 s.

The element formulation contains a separate transient matrix:

```fortran
CMA_N(8,8)
```

assembled from the finite-element interpolation functions:

```fortran
CMA_N(i,j)=CMA_N(i,j)+ro*cp*cg(i,j)*djac
```

The historical source explicitly labels this term:

```fortran
C DO STANU NIESTAC
```

The system used for the new time level is stored Element-by-Element through:

```fortran
bb(k)=2.0d0*hma(i,j)+3.0d0/i_d_czas*CMA_N(i,j)
```

while the previous solution contributes through:

```fortran
hma_t(i,j)=hma(i,j)-3.0d0/i_d_czas*CMA_N(i,j)
```

and the local right-hand-side assembly:

```fortran
fmv3(I_OT_ELEM(j))=
     fmv3(I_OT_ELEM(j))+3.0d0*fmv(j)-fmv_t(j)
```

This README intentionally documents the algebra as implemented rather than assigning a modern scheme name to the historical formulation.

## Time-dependent physical configuration

The routine:

```fortran
subroutine atmosfera(...)
```

changes the saved atmospheric configuration after:

```text
t = 1000 s
```

For the first stage, the velocity includes:

```fortran
v(1)=0.
v(2)=0.
v(3)=1.554128
```

while after 1000 s the saved case changes the second velocity component:

```fortran
v(2)=-0.50d+00
v(3)=1.554128
```

The diffusion and first-order removal parameters are also defined in this routine.

The source therefore demonstrates a transient calculation in which the physical configuration can change during the run.

## Source terms

The preserved test contains localized volumetric sources associated with selected elements.

For example:

```fortran
if (nr_el.eq.273) qv=8.0d+00*1d-09
```

and, in the later time interval:

```fortran
if (nr_el.eq.687) qv=1.111*1d-09
```

These values belong to the historical test configuration and are not general model restrictions.

## Element-by-Element architecture

The global transport operator is applied without assembling a conventional full global sparse matrix for the iterative solve.

The code contains:

```fortran
subroutine asub(...)
subroutine atsub(...)
```

implementing the actions

$$y=A x$$

and

$$y=A^T x.$$

Each routine traverses the HEX8 elements, applies the corresponding local `8 x 8` matrix, and accumulates contributions into the global vector.

This is a direct historical predecessor of the Element-by-Element operator architecture retained in modern Vanadis.

## Historical iterative solver

The routine:

```fortran
subroutine sparse(...)
```

solves the nonsymmetric system using repeated `A v` and `A^T v` operations.

The implementation uses residual-based search directions and both the forward and transpose operators.

It should be described as a **historical nonsymmetric iterative solver** rather than identified with the current production DBCG routine, which is a later implementation.

## Boundary conditions

The saved test applies first-kind / Dirichlet boundary conditions through:

```fortran
call wb_1(...)
```

on several domain faces.

The call for one face is commented out in the preserved configuration, reflecting the particular historical test case.

This should not be interpreted as a limitation of the general finite-element formulation.

## Obstacle representation

The source contains a simple geometric experiment in which three finite elements are omitted:

```fortran
C     USUWANIE 3 ELEMENTOW - przeszkoda

      IF (JK.EQ.739) GOTO 100
      IF (JK.EQ.740) GOTO 100
      IF (JK.EQ.741) GOTO 100
```

This represents an early obstacle / excluded-volume experiment within the structured HEX8 domain.

It is historically interesting because obstacle effects were later discussed in published Vanadis work.

## Output for successive time levels

The solver writes concentration fields for successive simulation times to files such as:

```text
t_0005.TXT
t_0010.TXT
...
```

Each record contains:

```text
X  Y  Z  concentration
```

The concentration output is scaled by `1e6` in the saved test.

These time-dependent result files formed the basis for the contemporary visualization and animation tools written in Salford FTN90.

## Connection to the 1998 defence animation

The transient fields generated by this development line were later visualized using dedicated FTN90 graphics programs.

The repository separately preserves:

```text
1998_defence_animation/
├── xb.f90
├── run1b.f90
└── vanadis1.webm
```

The surviving doctoral-defence documentation independently records that animated transient results were demonstrated during the public defence in 1998.

The preserved WebM is a modern screen capture of the historical visualization software running in a virtual machine; it is not an original video recording from 1998.

## Known historical implementation issue

The source declares:

```fortran
character*12 wynik_t(100)
```

but later indexes the array with:

```fortran
wynik_t(i_czas)
```

while the simulation time index extends beyond 100.

This can address elements outside the declared array bounds.

The intended use appears to require only one filename string at a time, so a restored implementation could replace this with a scalar character variable or a correctly sized array.

The historical source should remain unchanged.

## Modern compiler check

The preserved source can still be compiled with a modern GNU Fortran compiler using legacy compatibility options such as:

```bash
gfortran -std=legacy \
         -ffixed-line-length-none \
         -fallow-argument-mismatch \
         FM_T_OK.FOR -o fm_t_ok
```

Successful compilation does **not** imply that the program is a supported modern release; runtime use still depends on the original file workflow and historical test assumptions.

## Preservation policy

`FM_T_OK.FOR` should be retained verbatim.

If a corrected or modernized version is prepared later, keep it separate, for example:

```text
1997_3d_transient_ebe/
├── README.md
├── original/
│   └── FM_T_OK.FOR
└── restored/
    └── FM_T_OK_restored.for
```

This keeps the historical artifact distinct from subsequent reconstruction.

## Integrity

SHA-256 of the preserved source:

```text
da7987804c1577d0fe57f93274190d4d4143a148a257a0e59bd01ab86dacc1f7  FM_T_OK.FOR
```

File size:

```text
36129 bytes
```

## Significance

`FM_T_OK.FOR` is one of the most important preserved historical Vanadis sources.

It documents that the Vanadis development line had already reached, by the late 1990s, a technically recognizable architecture containing:

- 3D HEX8 finite elements,
- transient transport,
- Element-by-Element matrix storage,
- nonsymmetric iterative solution,
- explicit `A v` and `A^T v` operators,
- changing atmospheric conditions,
- localized sources,
- obstacle experiments,
- and time-resolved concentration output.

The source therefore provides a direct software link between the early Pascal/FEM work and the modern transient EBE/DR/OpenMP/CUDA Vanadis solver.
