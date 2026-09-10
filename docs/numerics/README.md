### Alternative BDF2 time discretization

Vanadis currently uses the A-stable \(\theta\)-scheme with

$$
\theta=\frac{2}{3}.
$$

An experimental BDF2 implementation was also prepared in order to compare the two second-order time-integration approaches.

An interesting implementation property is that the existing \(\theta=2/3\) formulation and BDF2 lead to the **same scaled left-hand-side matrix**.

For the \(\theta=2/3\) scheme, the discrete system can be written as

$$
\left(2H+\frac{3}{\Delta t}C\right)u^{n+1}
=
3f-\left(H-\frac{3}{\Delta t}C\right)u^n,
$$

where \(H\) denotes the spatial transport operator and \(C\) is the transient (mass) matrix.

For BDF2,

$$
C\frac{3u^{n+1}-4u^n+u^{n-1}}{2\Delta t}
+Hu^{n+1}=f,
$$

which, after multiplication by 2, becomes

$$
\left(2H+\frac{3}{\Delta t}C\right)u^{n+1}
=
2f+\frac{4}{\Delta t}Cu^n
-\frac{1}{\Delta t}Cu^{n-1}.
$$

Therefore the BDF2 variant does not require changes to the EBE operator representation or to the iterative linear solver. The main implementation differences are:

* storage of one additional solution level, \(u^{n-1}\);
* use of the existing \(\theta=2/3\) scheme for the first timestep, since BDF2 requires two previous solution levels;
* a different transient right-hand-side assembly for subsequent timesteps;
* advancement of the solution history after each step.

The implementation difference is intentionally small and is provided as a reference patch:

[`docs/numerics/theta23_to_bdf2.diff`](docs/numerics/theta23_to_bdf2.diff)

The production Vanadis solver continues to use the \(\theta=2/3\) scheme. In the tests performed during development, differences between the two formulations were small, while the \(\theta=2/3\) formulation remained the preferred scheme for Vanadis.
