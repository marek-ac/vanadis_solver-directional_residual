# Early 3D transient Vanadis — banded LU branch

This directory preserves a historical three-dimensional transient Fortran implementation from the early development of **Vanadis**.

The source is retained as a historical software artifact and is **not part of the current Vanadis production release**.

## Source file

- `LU_T_OK.FOR`

The historical filename is already portable and is preserved unchanged.

## Historical context

This source belongs to the 1997 transient-development stage of Vanadis.

It represents an alternative solution branch to the Element-by-Element iterative version preserved separately in:

```text
../1997_3d_transient_ebe/FM_T_OK.FOR
```

Both codes use the same broad 3D transient finite-element formulation, but this version assembles the global system in **banded storage** and solves it by LU factorization.

The file is therefore useful historically because it shows that two different linear-system strategies were being tested within the same Vanadis development line:

```text
3D transient FEM
├── Element-by-Element iterative solver
└── assembled band matrix + LU
```

## Computational grid

The preserved test configuration uses:

```fortran
parameter (nax=5,nay=10,naz=13)
```

which gives:

```text
432 HEX8 elements
650 nodes
```

The structured-grid bandwidth parameter is computed as:

```fortran
mband=nax*nay+nax+2
```

For this saved case:

```text
MBAND = 57
```

and the band-storage leading dimension used by the LU routines is:

```fortran
lda=3*MBAND-2
```

giving:

```text
LDA = 169
```

The global matrix is stored as:

```fortran
common /duza_tB/ hma2B(169,650),fmv2B(650)
```

## Transient finite-element formulation

The code contains the transient element matrix:

```fortran
CMA_N(8,8)
```

assembled from the finite-element interpolation functions:

```fortran
CMA_N(i,j)=CMA_N(i,j)+ro*cp*cg(i,j)*DJac
```

The historical source explicitly marks this contribution:

```fortran
C DO STANU NIESTAC
```

The saved time step is:

```fortran
i_d_czas=1
```

The transient system assembled into the band matrix contains:

```fortran
2.*hma(i,j)+3./i_d_czas*CMA_N(i,j)
```

while the old-solution contribution uses:

```fortran
hma_t(i,j)=hma(i,j)-3./i_d_czas*CMA_N(i,j)
```

As with the EBE historical source, this README documents the algebra exactly as implemented and does not retroactively assign a modern time-integration label unless independently verified.

## Banded global assembly

Unlike the EBE version, this branch assembles element contributions into a global band matrix:

```fortran
hma2B(I_OT_ELEM(I)-I_OT_ELEM(j)+M,I_OT_ELEM(j))
```

with:

```fortran
M=2*MBAND-1
```

The assembled system therefore uses a conventional band representation suitable for direct factorization.

This branch is valuable for comparison with the EBE architecture because it shows the memory/solver alternative that existed during early Vanadis development.

## LU solution

The code solves the banded system using routines derived from the LINPACK/BLAS numerical-software lineage.

The main solve sequence is:

```fortran
IF (I_CZAS.EQ.1)
     & call SGBCO(HMA2B,LDA,N,ML,MU,IPVT,RCOND,Z)

call dgbsl(HMA2B,lda,n,ml,mu,ipvt,FMV2B,job)
```

The LU factorization is performed only at the first time step in the saved test:

```fortran
IF (I_CZAS.EQ.1)
```

and the factorized band matrix is then reused for subsequent right-hand sides.

This is possible in this historical configuration because the system matrix is treated as unchanged over the simulated interval.

## Bundled numerical routines

The source includes historical LINPACK/BLAS-style routines directly in the same Fortran file.

The preserved comments identify the provenance, including text such as:

```text
linpack. this version dated 08/14/78 .
cleve moler, university of new mexico, argonne national lab.
```

The source also documents band factorization and solve operations such as:

```text
SGBCO
DGBSL
```

and associated BLAS routines.

These routines are **not original Vanadis algorithms** and should remain accompanied by their historical attribution comments.

Their presence is typical of scientific-software practice from the period, when numerical-library routines were often distributed directly with application source code.

## Physical test configuration

The routine:

```fortran
subroutine atmosfera(...)
```

defines the saved transport parameters.

Among the historical test values are:

```fortran
ab_pasquil=1.43
lambda(1)=8.15
lambda(2)=lambda(1)*ab_pasquil
lambda(3)=lambda(1)*ab_pasquil

pzanik=0.000693147

v(1)=0.
v(2)=0.
v(3)=0.5
```

A localized volumetric source is applied through:

```fortran
if (nr_el.eq.126) qv=8*1e-9
```

These values describe one historical verification/test case and are not fixed limits of the formulation.

## Output

The saved test runs until:

```text
t = 100 s
```

and then writes the solution to a text file containing:

```text
X  Y  Z  concentration
```

with concentration scaled by `1e6`.

The source terminates intentionally after producing the `t = 100 s` result:

```fortran
if (i_czas.eq.100) then
   ...
   stop
endif
```

## Relation to the EBE branch

The most important historical distinction is:

### `FM_T_OK.FOR`

```text
local HEX8 matrices
        ↓
Element-by-Element Av / A^T v
        ↓
iterative solver
```

### `LU_T_OK.FOR`

```text
local HEX8 matrices
        ↓
assembled global band matrix
        ↓
band LU factorization
```

The two preserved branches therefore document an early architectural comparison between:

- lower-memory operator-based iterative solution,
- and a conventional assembled direct-solver approach.

The later Vanadis line retained the Element-by-Element strategy, which became increasingly attractive as model size grew.

## Modern compiler test

The preserved file was tested with GNU Fortran using:

```bash
gfortran -std=legacy \
         -ffixed-line-length-none \
         -fallow-argument-mismatch \
         LU_T_OK.FOR -o lu_t_ok
```

The historical source compiles successfully without source modification.

The resulting executable was also run with the preserved built-in test configuration and reached the programmed termination at:

```text
czas = 100
```

This is a useful reproducibility result, but it should not be interpreted as certification of the historical model for present-day scientific use.

## Preservation policy

`LU_T_OK.FOR` should be kept unchanged.

If a restored or modernized implementation is prepared later, it should be stored separately, for example:

```text
1997_3d_transient_lu/
├── README.md
├── original/
│   └── LU_T_OK.FOR
└── restored/
    └── LU_T_OK_restored.for
```

This keeps the original historical artifact distinct from later reconstruction.

## Integrity

SHA-256 of the preserved source:

```text
b03335369fa56b2675b1c75be2b1a65af2d3caec913f8778838614e121293af9  LU_T_OK.FOR
```

File size:

```text
62559 bytes
```

## Significance

`LU_T_OK.FOR` documents an important alternative branch of early Vanadis development.

It shows that the transient 3D HEX8 formulation was not tied to a single linear-system architecture: the same class of finite-element problem was also implemented using a conventional assembled band matrix and direct LU factorization.

Together with `FM_T_OK.FOR`, this source provides unusually clear historical evidence of the solver-design choices that eventually led Vanadis toward the Element-by-Element architecture used in later versions.
