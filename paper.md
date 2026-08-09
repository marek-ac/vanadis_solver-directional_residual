Summary

Vanadis is an open-source finite element solver designed for large-scale 3D transport problems in atmospheric modelling. The code implements Directional Residual (DR) stabilization, a stabilization approach exhibiting good global mass-balance behavior in advection-dominated flows, and provides a parallel CUDA/GPU implementation suitable for high-performance environmental simulations. The CUDA implementation accelerates the solution of the linear system using a parallel Element-by-Element strategy combined with lightweight diagonal preconditioning.

Boundary conditions at the ground surface include either a Dirichlet constraint or a flux-type deposition condition, with the latter controlled by the parameter α, which controls pollutant transfer to the surface.

Time integration is performed using an implicit time-discrete formulation, in which element matrices and right-hand-side contributions are constructed locally and the resulting algebraic system is solved iteratively at each time step using an Element-by-Element representation. In the tested atmospheric dispersion cases, this approach demonstrated robust numerical stability and transient behavior for stiff advection–diffusion systems. The combination of DR spatial discretization, GPU-accelerated Element-by-Element computations, and implicit time integration results in a robust and efficient framework for large-scale atmospheric dispersion simulations.

Vanadis includes:

   - a stabilized FEM formulation based on Directional Residual stabilization,

   - a CUDA-accelerated Element-by-Element iterative solver,

   - support for large three-dimensional hexahedral meshes,

   - lightweight diagonal preconditioning,

   - open-source code and reproducible examples.

The project aims to provide a lightweight, efficient, and scientifically transparent tool for atmospheric and environmental modelling.

Statement of need

Accurate simulation of atmospheric pollutant transport requires numerical methods that remain stable and exhibit good global mass-balance behavior under strongly advection-dominated conditions. Traditional finite element discretizations may exhibit spurious oscillations or loss of numerical stability when applied to strongly advection-dominated transport problems, while stabilization techniques may introduce additional numerical diffusion. Vanadis addresses these challenges by combining Directional Residual stabilization with GPU-accelerated computation, enabling efficient and stable simulations on large meshes.

State of the field

Several established atmospheric transport models exist, including Eulerian grid-based solvers and finite volume approaches commonly used in air quality modelling. Stabilized finite element methods are also widely studied, with SUPG, GLS, and discontinuous Galerkin methods being the most common.

Vanadis contributes to this landscape by:

- implementing Directional Residual stabilization designed for local,
  streamline-oriented treatment of advection-dominated transport,

- providing an Element-by-Element CUDA solver with low memory requirements,

- combining the stabilization with lightweight diagonal preconditioning and iterative solution methods,

- focusing specifically on scalar advection–diffusion transport,

- offering a compact and transparent codebase suitable for experimentation
  with stabilization and time-integration schemes.

This positions Vanadis as a complementary tool to larger atmospheric modelling frameworks, particularly for research on numerical stabilization and GPU-accelerated transport solvers.

Functionality

Vanadis provides:

- stabilized finite element discretization for 3D advection–diffusion transport,
- CUDA-accelerated iterative solution using an Element-by-Element strategy,
- support for Dirichlet and flux-type deposition boundary conditions,
- implicit transient time integration,
- support for large nonuniform HEX8 meshes,
- example cases for atmospheric dispersion scenarios.

The solver is implemented in Fortran/C/CUDA and distributed under an open-source license. Example input files and scripts allow users to reproduce published results and adapt the solver to new transport problems.

Quality control

The codebase includes:

   - automated build scripts,
   - validation examples including comparison with a Pasquill-based Gaussian plume reference model,
   - benchmark cases demonstrating GPU performance,
   - reproducible atmospheric dispersion scenarios.

All examples are included in the repository and documented on the project homepage.

Acknowledgements

The initial concept behind the Vanadis solver was influenced by earlier work on convection-dominated heat transfer in FEM (Kryzhanowski & Pietrzyk, 1995).
The author expresses gratitude to Prof. M. Pietrzyk and Prof. Z. Malinowski for their supervision during the early stages of the author’s doctoral research (1995–1998),
which provided a valuable scientific background for later developments in finite-element modelling of atmospheric transport.
The Directional Residual (DR) method and the Vanadis solver were developed independently by the author at a later stage of this research trajectory.
The author acknowledges the use of AI-assisted programming tools (Microsoft Copilot) for code refactoring and automation of selected routines, as well as the use of an AI language model for text
editing and language refinement.


## Project links

- Homepage: [https://marek-ac.meri.pl/](https://marek-ac.meri.pl/)
- DOI: [https://doi.org/10.5281/zenodo.21352014](https://doi.org/10.5281/zenodo.21352014)
- GitHub: [https://github.com/marek-ac/vanadis_solver-directional_residual/releases/tag/v2026.1.1](https://github.com/marek-ac/vanadis_solver-directional_residual/releases/tag/v2026.1.1)
