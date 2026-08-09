Summary

Vanadis is an open-source finite element solver designed for large-scale 3D transport problems in atmospheric modelling. The code implements Directional Residual (DR) stabilization, a mass-conserving method for advection-dominated flows, and provides a fully parallel CUDA/GPU implementation suitable for high-performance environmental simulations. The CUDA implementation accelerates the solution of the linear system using a parallel element-by-element strategy, taking advantage of the diagonal structure preserved by DR stabilization.

Boundary conditions at the ground surface include either a Dirichlet constraint or a flux-type deposition condition, with the latter controlled by the parameter α, which governs pollutant penetration into the surface.

Time integration is performed using an implicit time-discrete formulation, in which the system matrix and right-hand-side vector are assembled element-wise and the resulting global system is solved iteratively at each time step. This approach provides strong numerical stability and robust transient behavior for stiff advection–diffusion systems. The combination of DR spatial discretization, GPU-accelerated element-wise computations, and the implicit time-integration scheme results in a robust and efficient framework for large-scale atmospheric dispersion simulations.

The combination of DR spatial discretization, GPU-accelerated element-wise solving, and the stabilized multistep time scheme results in a robust and efficient tool for large-scale atmospheric dispersion simulations.

Vanadis includes:

    a stabilized FEM formulation based on Directional Residuals,

    a GPU-accelerated linear solver using CUDA,

    support for large meshes (millions of nodes),

    open-source code and reproducible examples.

The project aims to provide a lightweight, efficient, and scientifically transparent tool for atmospheric and environmental modelling.

Statement of need

Accurate simulation of atmospheric pollutant transport requires numerical methods that remain stable and mass-conserving under strongly advection-dominated conditions. Traditional finite element approaches often suffer from numerical diffusion or instability when applied to large-scale 3D transport problems. Vanadis addresses these challenges by combining Directional Residual stabilization with GPU-accelerated computation, enabling efficient and stable simulations on large meshes.

Researchers working in atmospheric modelling, air quality, environmental engineering, and computational fluid dynamics can use Vanadis as a specialized tool for convection–diffusion transport. Its open-source nature and reproducible examples make it suitable for scientific studies, method development, and educational purposes.
State of the field

Several established atmospheric transport models exist, including Eulerian grid-based solvers and finite volume approaches commonly used in air quality modelling. Stabilized finite element methods are also widely studied, with SUPG, GLS, and discontinuous Galerkin methods being the most common.

Vanadis contributes to this landscape by:

    implementing Directional Residual stabilization, which preserves diagonal structure and enables efficient GPU parallelization,

    providing a fully element-wise CUDA solver,

    focusing specifically on scalar transport rather than full Navier–Stokes systems (The Directional Residual concept can be naturally extended to vector systems such as the Navier-Stokes equations),

    offering a compact, transparent codebase suitable for experimentation with stabilization and time-integration schemes.

This positions Vanadis as a complementary tool to larger atmospheric modelling frameworks, particularly for research on numerical stabilization and GPU-accelerated transport solvers.

Functionality

Vanadis provides:

    stabilized finite element discretization for 3D advection–diffusion transport,

    CUDA-accelerated linear system solving using an element-by-element strategy,

    support for Dirichlet and flux-type deposition boundary conditions,

    implicit second-order multistep time integration,

    mesh handling for large unstructured grids,

    example cases for atmospheric dispersion scenarios.

The solver is implemented in Fortran/C/CUDA and distributed under an open-source license. Example input files and scripts allow users to reproduce published results and adapt the solver to new transport problems.

Quality control

The codebase includes:

    automated build scripts,

    validation examples comparing numerical results with analytical solutions for advection–diffusion transport,

    benchmark cases demonstrating GPU speedup,

    reproducible atmospheric dispersion scenarios.

All examples are included in the repository and documented on the project homepage.

Acknowledgements

The initial concept behind the Vanadis solver was influenced by earlier work on convection-dominated heat transfer in FEM (Kryzhanowski & Pietrzyk, 1995).
The author expresses gratitude to Prof. M. Pietrzyk and Prof. Z. Malinowski for their
supervision during the early stages of the author’s doctoral research (1995–1998), which provided
a valuable scientific background for later developments in finite-element modelling of atmospheric
transport. The Directional Residual (DR) method and the Vanadis solver were developed indepen-
dently by the author at a later stage of this research trajectory.
The author acknowledges the use of AI-assisted programming tools (Microsoft Copilot) for code
refactoring and automation of selected routines, as well as the use of an AI language model for text
editing and language refinement.


Project links

    Homepage: https://marek-ac.meri.pl/

    DOI: https://doi.org/10.5281/zenodo.21352014

    GitHub: https://github.com/marek-ac/vanadis_solver-directional_residual/releases/tag/v2026.1.1
