# Early 2D Pascal FEM atmospheric-transport model

This directory preserves an early Pascal implementation from the development line that later evolved into **Vanadis**.

The source is retained as a historical software artifact and is **not part of the current Vanadis release**.

## Source file

- `ART_SPEE.PAS`

The preserved source is a two-dimensional finite-element model for atmospheric transport. It uses triangular elements and solves a steady convection–diffusion–reaction problem with optional upwind weighting.

## Historical context

The Vanadis research line originated during doctoral work begun in 1994.

Related archived Pascal sources document a working 2D FEM atmospheric-transport code by 1995. The preserved copy of `ART_SPEE.PAS` belongs to this early development period; its archived internal timestamp is **4 February 1996**.

The file is therefore best treated as part of the **1995–1996 early 2D Pascal stage** of Vanadis development rather than as a modern reconstructed source.

## Numerical characteristics

The code contains:

- two-dimensional finite-element discretization,
- three-node triangular elements,
- convection/advection in two spatial directions,
- anisotropic diffusion coefficients `kx` and `ky`,
- a volumetric source term,
- a first-order reaction/decay coefficient,
- support for boundary-condition terms,
- banded global matrix storage,
- a custom Gaussian-elimination solver for the nonsymmetric band system,
- optional upwind weighting for convection-dominated cases.

The source parameters include, for example:

```pascal
kx=10;
ky=10;
vx=0;
vy=0;
pzanik=0.0000;
```

These values represent one historical test configuration and should not be interpreted as fixed limitations of the model.

## Upwind formulation

The code contains a switch:

```pascal
ab=1;      {1 - stosowany upwinding , 0 - niestosowany upwinding}
```

and evaluates a Peclet-dependent weighting factor through:

```pascal
gamma:=vx*dlugelx/kx;

alfa2:=(exp(gamma*0.5)+exp(-0.5*gamma))
       /(exp(gamma*0.5)-exp(-0.5*gamma))-2/gamma;
```

Mathematically, this corresponds to a weighting of the form

$$\alpha = \coth(Pe)-\frac{1}{Pe},$$

with the local Peclet number represented by `gamma/2`.

This is an early convection-stabilization mechanism and is historically related to the later development of stabilized Vanadis formulations. It should **not** be identified with the present Directional Residual (DR) formulation, which uses a different three-dimensional directional metric and directional diffusivity.

## Matrix solution

The global finite-element system is stored in banded form. The routine

```pascal
procedure gaussasy2;
```

implements a custom Gaussian-elimination procedure for the resulting nonsymmetric system.

This reflects the memory constraints and implementation environment of the original DOS-era code.

## Boundary-condition support

The source contains handling associated with several boundary-condition types, including parameters for:

- prescribed boundary values,
- flux terms,
- third-kind / Robin-type terms.

The exact configuration used in a given run is controlled by hard-coded parameters in the historical source.

## Software environment

The program was written for the DOS-era Pascal environment and uses:

```pascal
uses crt,dos;
```

together with compiler directives and memory settings typical of that platform.

The source has intentionally been preserved in its original style. It may require adaptation to compile with a modern Pascal compiler.

## Preservation policy

The historical file should be kept **unchanged**.

If a modernized or restored version is created in the future, it should be placed separately, for example:

```text
1995_2d_pascal/
├── README.md
├── original/
│   └── ART_SPEE.PAS
└── restored/
    └── ART_SPEE_restored.pas
```

This keeps the original historical artifact distinct from later maintenance work.

## Integrity

SHA-256 of the preserved source:

```text
e68ad2df5f841c2c174e206ac8dbc972fc9dec7572ee272371108b8a9bcdce33  ART_SPEE.PAS
```

File size:

```text
24745 bytes
```

## Significance

`ART_SPEE.PAS` documents the early stage of the Vanadis development line before the transition to the later three-dimensional Fortran implementations.

It shows that the project already included:

- a finite-element atmospheric-transport formulation,
- treatment of convection, diffusion and reaction,
- stabilization for convection-dominated transport,
- explicit control of boundary terms,
- and a dedicated numerical solution procedure.

The source is preserved primarily for **historical documentation, reproducibility of the project history, and study of the numerical evolution of Vanadis**.
