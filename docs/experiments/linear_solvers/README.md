# Experimental linear solvers

These routines are preserved as alternative research solvers for Vanadis.

They are not drop-in replacements for the current production solver and may
require adaptation to the current Dirichlet projection, nonlinear Picard loop
and memory architecture.

Included implementations:

- `dbicgstab_ebe_lu_pivot_fixed.for`
  - BiCGSTAB
  - Element-by-Element matrix-vector operations
  - local HEX8 LU preconditioning
  - partial pivoting

- `dgmres_ebe_pivoting_replacement.for`
  - restarted GMRES(20)
  - Element-by-Element matrix-vector operations
  - local HEX8 LU preconditioning
  - partial pivoting and reorthogonalization

The current Vanadis production solver remains DBCG with Jacobi
preconditioning.

See the accompanying solver comparison PDF for architectural and memory
considerations.