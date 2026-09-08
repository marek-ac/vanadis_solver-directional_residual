# Vanadis experiments

This directory contains experimental numerical methods, alternative solver
implementations and formulation variants that are being preserved for testing,
comparison and possible future integration into Vanadis.

Code in this directory is **not part of the current production release** unless
explicitly stated otherwise.

The material may include:

- alternative Krylov solvers,
- alternative preconditioners,
- numerical robustness experiments,
- Directional Residual formulation variants,
- research prototypes,
- code intended for benchmark or regression testing.

## Current experiments

### Linear solvers

Alternative Element-by-Element Krylov solver implementations, including
BiCGSTAB and restarted GMRES with local HEX8 LU preconditioning.

See [`linear_solvers`](linear_solvers/).

### Directional Residual

Experimental and robustness-oriented variants of the Directional Residual
formulation.

See [`directional_residual`](directional_residual/).