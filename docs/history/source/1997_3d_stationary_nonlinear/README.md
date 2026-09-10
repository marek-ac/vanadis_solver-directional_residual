# 3D stationary nonlinear Vanadis prototype

This directory preserves a historical Fortran implementation from the early three-dimensional development of **Vanadis**.

The source is retained as a historical software artifact and is **not part of the current Vanadis release**.

## Source file

- `GR_P_S.FOR`

Original historical filename:

```text
GR#P(S)(2).FOR
```

The repository filename has been simplified to avoid special characters. The file contents should remain unchanged.

## Historical context

This source belongs to the early 3D Fortran stage of Vanadis development and is preserved in the repository chronology as a **1997 stationary nonlinear prototype**.

The program predates the modern Directional Residual, OpenMP and CUDA implementations.

Its importance is that it already combines:

- a three-dimensional finite-element transport model,
- eight-node hexahedral elements,
- Element-by-Element storage and matrix-vector operations,
- an iterative solver operating through `A v` and `A^T v`,
- and a concentration-dependent reaction/removal coefficient updated by a fixed-point iteration.

## Computational mesh

The saved test configuration uses:

```fortran
nx1 = 6
nx2 = 13
nx3 = 22
```

which gives:

```text
1716 HEX8 elements
2254 nodes
```

Each element contributes a local `8 x 8` matrix. The Element-by-Element matrix data are stored in:

```fortran
bb(109824)
```

corresponding to:

```text
1716 x 64
```

stored element-matrix coefficients.

## Governing transport terms

The element routine includes:

- three-dimensional convection,
- anisotropic diffusion,
- volumetric source terms,
- reaction / removal,
- finite-element boundary terms.

The preserved test configuration contains, among other values:

```fortran
ab_pasquil = 1.43

lambdaX = 8.15
lambdaY = lambdaX*ab_pasquil
lambdaZ = lambdaX*ab_pasquil

vx = 0.
vy = 0.
vz = 1.554128
```

These are historical test parameters and should not be interpreted as fixed limitations of the formulation.

## Nonlinear concentration-dependent removal

The most distinctive feature of this source is the concentration-dependent coefficient `pzanik`.

For each element, an element concentration measure is obtained as the mean of the eight nodal solution values:

```fortran
tnw3(je)=tnw2(lnd2(j))/8+tnw3(je)
```

The local reaction/removal coefficient is then evaluated as:

```fortran
pzanik=4e-6
     & +dabs(tnw3(kel))*2e-4
     & +(dabs(tnw3(kel)))**2*0.1
```

or, mathematically,

$$P(\bar S_e) = 4\times10^{-6} + 2\times10^{-4}|\bar S_e| + 0.1|\bar S_e|^2.$$

This makes the stationary transport problem nonlinear.

## Fixed-point / Picard-type iteration

The nonlinear problem is solved by repeatedly:

1. computing the current element-average concentration,
2. updating `pzanik`,
3. rebuilding the Element-by-Element operator,
4. solving the resulting linear problem,
5. comparing the new nodal solution with the previous one.

The historical convergence measure is:

```fortran
z_kontrola=0.
do i=1,non
   z_kontrola=(dabs(tnw4(i))-dabs(tnw2(i)))**2+z_kontrola
enddo
```

and the nonlinear loop is repeated while:

```fortran
z_kontrola.gt.0.00000000000001
```

This is best described as a **Picard-type fixed-point iteration** for the concentration-dependent coefficient.

## Element-by-Element operator

The global matrix is not required as a conventional assembled sparse matrix for the iterative solution.

The source provides separate Element-by-Element operator routines:

```fortran
asub
atsub
```

implementing:

$$y = A x$$

and

$$y = A^T x.$$

`asub` accumulates the contributions of the local `8 x 8` matrices directly into the global output vector.

The explicit transpose operator in `atsub` is important historically because it anticipates the operator structure later used by nonsymmetric Krylov solvers in Vanadis.

## Historical iterative linear solver

The routine:

```fortran
subroutine sparse(...)
```

uses both `asub` and `atsub`.

It is a conjugate-gradient-style iterative minimization of the residual norm for the nonsymmetric Element-by-Element operator. It should not be confused with the modern Vanadis DBCG implementation.

The key operations include:

```fortran
call asub(x,xi,...)
call atsub(xi,g,...)
```

followed by iterative search-direction updates.

This historical implementation is preserved as part of the evolution of Vanadis linear-system solution methods.

## Boundary conditions in the saved test

The preserved configuration applies first-kind / Dirichlet conditions on boundary faces:

```text
1, 2, 4, 5, 6
```

while the call for face `3` is commented out.

This is a property of the saved historical test case, not a general restriction of the finite-element formulation.

## Restart / initial-vector workflow

The historical program uses the binary file:

```text
wyn.bin
```

to save and reload the solution vector between nonlinear iterations / runs.

The source also contains a commented initialization block:

```fortran
c     DO I=1,NON
c     TNW2(I)=0
c     ENDDO
```

with the note that it should be enabled for the first run.

This reflects the original development workflow and is intentionally preserved rather than modernized.

## Known historical implementation note

The source contains expressions such as:

```fortran
3/2*ALFA1
3/2*ALFA2
3/2*ALFA3
```

In Fortran, `3/2` is integer division and therefore evaluates to `1`, not `1.5`.

A restored version intended for renewed numerical use should replace these expressions with an explicit real value, for example:

```fortran
1.5d0*ALFA1
```

and regression-test the resulting formulation.

The historical source in this directory should remain unchanged.

## Software status

The code is fixed-form legacy Fortran and uses conventions typical of its original development environment.

It can be treated as:

```text
historical source
not production code
not a drop-in replacement for current Vanadis
```

If a modernized version is created, it should be stored separately from the original source.

## Integrity

SHA-256 of the preserved source:

```text
aeacccbbee1ddc7193241141f3cac89452d36ccb6e782f15a70bfe8dab5009a1  GR_P_S.FOR
```

File size:

```text
31765 bytes
```

## Significance

`GR_P_S.FOR` documents an important intermediate stage in the Vanadis development line.

It shows that the project had already moved beyond the early 2D Pascal implementation to a three-dimensional HEX8 finite-element architecture with:

- Element-by-Element matrix storage,
- direct `A v` and `A^T v` operators,
- iterative solution of a nonsymmetric transport system,
- and nonlinear concentration-dependent reaction/removal treated by fixed-point iteration.

These ideas form part of the numerical lineage that later evolved into the modern transient, stabilized, OpenMP/CUDA Vanadis solver.
