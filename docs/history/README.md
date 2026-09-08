Historical publications and project history

Vanadis has a documented development history dating back to the 1990s.
The materials in this directory preserve the early publications and archival
documents that show how the model evolved from an initial two-dimensional FEM
formulation into the present three-dimensional transient solver.
Project history

    1997 — first published finite-element formulation for atmospheric pollutant
    dispersion. The published numerical implementation was two-dimensional,
    while the governing transport equation and the planned development already
    pointed toward a broader three-dimensional formulation.

    1998 — a working transient extension of the model was publicly demonstrated
    during the doctoral defence of Marek Chodorski on 29 June 1998. The defence
    protocol records the presentation of animated results showing the
    time-dependent evolution of pollutant transport.

    1999 — a three-dimensional FEM version was published and compared with the
    Pasquill dispersion model. The paper also demonstrated cases with a lateral
    wind component and an obstacle affecting the pollutant plume.

    Today — Vanadis is a 3D transient FEM solver with Directional Residual (DR)
    stabilization, an Element-by-Element (EbE) solution architecture, OpenMP CPU
    assembly and CUDA iterative solvers.

Historical publications
1997

First published finite-element formulation for atmospheric pollutant dispersion

The 1997 paper presented the use of the finite element method for atmospheric
pollutant dispersion as a methodological proposal.

The published implementation was two-dimensional and steady-state, but the
underlying transport formulation already included convection, diffusion,
emission and pollutant removal.

The paper also identified the numerical difficulty associated with
convection-dominated transport and discussed the need for a weighted-residual
approach with nonsymmetric weighting functions. This issue later became an
important part of the numerical development of Vanadis.
1998

Public demonstration of the transient extension

During the public doctoral defence of Marek Chodorski on 29 June 1998, a
working transient extension of the atmospheric transport model was demonstrated.

The defence protocol records the presentation of animated results showing the
time-dependent evolution of pollutant transport. This provides a dated archival
record that the transient implementation was already operational in 1998.
1999

Three-dimensional FEM model and comparison with the Pasquill model

The 1999 publication presented a three-dimensional finite-element atmospheric
dispersion model and compared its results with the Pasquill model.

For a case corresponding closely to the assumptions of the Pasquill approach,
the two models produced similar concentration distributions. Additional cases
demonstrated capabilities beyond the simple Gaussian formulation, including a
lateral wind component and an obstacle affecting the plume.

This publication represents the first clearly published three-dimensional stage
of the Vanadis development line.
History note

A bilingual Polish/English history note in this directory provides a short
overview of the development of Vanadis from the first published FEM formulation
in 1997, through the transient model demonstrated in 1998 and the published 3D
model in 1999, to the modern Vanadis architecture.

The historical record shows a continuous numerical-development line:

1997
2D FEM proof of concept
        |
        v
1998
working transient extension demonstrated publicly
        |
        v
1999
published 3D FEM model and Pasquill comparison
        |
        v
later development
EbE + nonsymmetric Krylov solvers + DR
        |
        v
current Vanadis
3D transient HEX8 + DR + EbE + OpenMP + CUDA

These documents are preserved as archival material documenting the origins and
long-term development of the Vanadis project.