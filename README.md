### Directional Residual stabilization

At every volume Gauss point, Vanadis evaluates the physical position, the local velocity and diffusion fields, and the HEX8 Jacobian

$$
J=
\frac{\partial(x,y,z)}
{\partial(\xi,\eta,\zeta)}.
$$

For non-zero velocity,

$$
\mathbf{e}
=
\frac{\mathbf{v}}
{\lVert \mathbf{v}\rVert}
$$

defines the local streamline direction.

The characteristic element length used by DR is obtained from the actual element mapping:

$$
h_{\mathrm{stream}}
=
\frac{2}
{\lVert J^{-1}\mathbf{e}\rVert}.
$$

For the diagonal diffusion tensor

$$
\Lambda
=
\operatorname{diag}
\left(
\lambda_1,\lambda_2,\lambda_3
\right),
$$

the diffusion acting in the streamline direction is

$$
\lambda_s
=
\mathbf{e}^{T}
\Lambda
\mathbf{e}.
$$

The local Peclet number is

$$
Pe
=
\frac{
\lVert\mathbf{v}\rVert
h_{\mathrm{stream}}
}{
2\lambda_s
}.
$$

The Directional Residual parameter is then

$$
\tau_{DR}
=
\frac{
h_{\mathrm{stream}}
}{
2\lVert\mathbf{v}\rVert
}
\left(
\coth(Pe)
-
\frac{1}{Pe}
\right).
$$

The implementation includes numerically safe small- and large-Peclet limits and the pure-advection limit

$$
\lambda_s \rightarrow 0
\qquad\Longrightarrow\qquad
\tau_{DR}
\rightarrow
\frac{
h_{\mathrm{stream}}
}{
2\lVert\mathbf{v}\rVert
}.
$$

The perturbed test weight is

$$
W_i
=
N_i
+
\tau_{DR}
\left(
\mathbf{v}\cdot\nabla N_i
\right).
$$

Consequently, the stabilization responds automatically to:

- the local velocity direction,
- the local size and distortion of the HEX8 element,
- the element orientation relative to the flow,
- anisotropic diagonal diffusion,
- spatially and temporally varying `v` and `K`.

For isotropic diffusion,

$$
\Lambda = K I,
$$

and therefore

$$
\lambda_s = K.
$$


### Weak-form structure

The diffusion part intentionally remains the standard Galerkin term

$$
\nabla N_i^{T}
\Lambda
\nabla N_j.
$$

The DR perturbation enters the transient, advective, reaction/source, and applicable boundary weighting through \(W_i\).

In simplified volume form, the element operator is

$$
H_{ij}
=
\int_{\Omega_e}
\left[
\nabla N_i^{T}
\Lambda
\nabla N_j
+
W_i
\left(
\mathbf{v}\cdot\nabla N_j
\right)
+
W_i P N_j
\right]
\,d\Omega .
$$

The DR-weighted mass matrix is

$$
M_{ij}
=
\int_{\Omega_e}
W_i N_j
\,d\Omega .
$$

For the element volumetric source \(q_e\), the load-vector contribution is

$$
F_i
=
\int_{\Omega_e}
W_i q_e
\,d\Omega .
$$

The corresponding boundary contributions are added separately.

This distinction is deliberate: Vanadis does **not** replace the Galerkin diffusion gradient by

$$
\nabla W_i.
$$

The diffusion operator remains Galerkin, while the Directional Residual contribution enters through residual weighting by \(W_i\).

Consequently, derivatives of \(W_i\) or \(\tau_{DR}\) are not required by the present formulation.


### Transient formulation

Vanadis v2026.3.1 uses an implicit

$$
\theta=\frac{2}{3}
$$

time-integration scheme.

Because the velocity, diffusion, boundary fields, and therefore the DR-weighted operator may vary with time, both the element operator \(H\) and the DR-weighted mass matrix \(M\) are evaluated at the old and new time levels.

The assembled time-discrete equation is

$$
\left(
2H^{n+1}
+
\frac{3}{\Delta t}
M^{n+1}
\right)
C^{n+1}
=
\left(
\frac{3}{\Delta t}
M^n
-
H^n
\right)
C^n
+
2F^{n+1}
+
F^n .
$$

For nonlinear concentration-dependent reaction/decay,

$$
P=P(S),
$$

the new-time-level problem is solved by Picard iteration while \(C^n\) remains fixed during the nonlinear iterations.
