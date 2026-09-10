# Historical publications and project history

Vanadis has a documented development history dating back to the 1990s.
The research line that later developed into Vanadis originated during doctoral
work begun in 1994. Archived source code demonstrates that a working
two-dimensional FEM atmospheric-transport model existed by 1995.

The materials in this directory preserve early publications, doctoral-defence
documentation and historical source code showing how the model evolved from an
initial two-dimensional FEM implementation into a three-dimensional transient
solver and, later, into the present Vanadis architecture.

An important point in the reconstructed chronology is that the transition to
three dimensions occurred earlier than the 1999 journal publication alone would
suggest. A POL-IMIS conference proceedings paper published in 1997 already
described a three-dimensional finite-element atmospheric-dispersion model based
on eight-node hexahedral elements and presented both stationary and transient
calculations.

## Project history

### 1994 — doctoral research begun

The research programme that eventually led to Vanadis began as doctoral work on
numerical modelling of atmospheric pollutant transport and dispersion.

### 1995–1996 — working 2-D FEM implementation

Archived Pascal source code documents a working two-dimensional stationary FEM
atmospheric-transport model.

The implementation already included:

- convection-diffusion-reaction transport,
- anisotropic diffusion coefficients,
- volumetric emission and first-order removal,
- Dirichlet and Robin/flux-type boundary terms,
- triangular finite elements,
- a nonsymmetric banded algebraic system,
- and Peclet-dependent upwind / Petrov-Galerkin-type weighting.

This stage predates the modern Directional Residual formulation and should not
be described as DR.

### 1997 — published 2-D and 3-D FEM stages

Two distinct publications from 1997 document different stages of the work.

#### 1997 journal paper — 2-D FEM proof of concept

Marek Chodorski and Maciej Pietrzyk published:

> **Propozycja zastosowania metody elementów skończonych do symulacji
> rozprzestrzeniania się zanieczyszczeń w atmosferze**

in *Ochrona Powietrza i Problemy Odpadów*, vol. 31, no. 3, pp. 90–93,
1997.

The paper presented a two-dimensional steady FEM implementation and discussed
the numerical difficulty of convection-dominated transport, including the need
for weighted-residual treatment with nonsymmetric weighting functions.

#### 1997 POL-IMIS proceedings paper — published 3-D HEX8 model

A separate paper by Marek Chodorski, Zbigniew Malinowski, Stanisław Słupek and
Andrzej Buczek was published in the proceedings of the II POL-IMIS conference:

> **Zastosowanie metody elementów skończonych do modelu rozprzestrzeniania się
> zanieczyszczeń w atmosferze**

Proceedings of **POL-IMIS 1997 — Ocena wielkości imisji zanieczyszczeń
powietrza**, Szklarska Poręba, 19–22 June 1997, PZITS publication no. 727,
pp. 13–22.

This proceedings paper already describes a three-dimensional atmospheric
transport model using eight-node hexahedral finite elements and three
directional convection and diffusion components.

The mathematical derivation shown in the paper is written for the stationary
case, but the numerical examples also include an explicitly transient variant:
a concentration field evaluated 100 s after the beginning of emission.

The paper further presents calculation variants with:

- lateral convection through a non-zero transverse wind component,
- an obstacle placed in the pollutant plume,
- concentration-dependent reaction/removal,
- volumetric emission,
- first-kind boundary conditions,
- and ground-surface absorption / exchange terms.

For the nonlinear variant the removal coefficient is written in the form

$$P = a_1 + a_2 S + a_3 S^2.$$

This publication therefore provides independent published evidence that the
research code had already reached a three-dimensional HEX8 stage in 1997 and
that transient and nonlinear experiments were being performed at that time.

### August 1997 — archived 3-D transient HEX8 / EBE source

Archived Fortran source from the same development period documents a working
three-dimensional transient implementation using eight-node hexahedral elements
and Element-by-Element operator application.

The preserved `FM_T_OK.FOR` source has an archival timestamp of 9 August 1997.
The timestamp is treated as provenance evidence rather than cryptographic proof
of authorship date.

The source contains:

- 3-D HEX8 finite elements,
- a transient mass matrix,
- implicit time integration,
- local 8 x 8 matrices stored element by element,
- direct Element-by-Element implementations of `A x` and `A^T x`,
- iterative treatment of a nonsymmetric system,
- time-dependent model parameters,
- localized volumetric sources,
- obstacle-related test logic,
- and output of successive transient concentration fields.

A separate 1997 source branch, `GR_P(S).FOR`, documents a stationary nonlinear
3-D HEX8 model with a concentration-dependent removal coefficient and a
Picard/fixed-point-type iteration.

These archived sources provide software evidence complementary to the 1997
POL-IMIS publication.

### 1998 — public doctoral-defence demonstration

A working transient extension of the model was publicly demonstrated during the
doctoral defence of Marek Chodorski on 29 June 1998.

The defence protocol records the presentation of animated results showing the
time-dependent evolution of pollutant transport. Historical FTN90 programs used
to prepare and display the animation are also preserved in this repository.

The 1998 defence is therefore an independent dated record of the transient model
being demonstrated publicly, although the newly recovered 1997 POL-IMIS paper
shows that transient calculations had already been documented earlier.

### 1999 — published 3-D FEM / Pasquill comparison

Marek Chodorski and Zbigniew Malinowski published:

> **Porównanie modelu Pasquilla z trójwymiarowym rozwiązaniem metodą elementów
> skończonych. Zagadnienia dyfuzji atmosferycznej**

in *Ochrona Powietrza i Problemy Odpadów*, vol. 33, no. 1, pp. 5–9, 1999.

The paper presented a three-dimensional finite-element atmospheric-dispersion
model and compared its results with the Pasquill dispersion model.

For a case corresponding closely to the assumptions of the Pasquill approach,
the two models produced similar concentration distributions. Additional cases
demonstrated the influence of a lateral wind component and an obstacle on the
pollutant plume.

This paper should now be viewed as a later published 3-D stage of the
development line and as the first preserved publication in the project history
that directly compares the 3-D FEM model with the Pasquill approach. It is no
longer accurate to describe it as the first published three-dimensional stage,
because the POL-IMIS proceedings paper already documented a 3-D HEX8 model in
1997.

### Later development

The model subsequently evolved through further numerical and implementation
developments, including:

- continued transient 3-D calculations,
- Element-by-Element solution architecture,
- nonsymmetric Krylov iterative solvers,
- Directional Residual stabilization,
- anisotropic directional Peclet treatment,
- OpenMP parallel CPU assembly,
- CUDA GPU iterative solution,
- nonlinear concentration-dependent removal,
- and generalized time-dependent source handling.

The modern Vanadis code is therefore best viewed not as a recent standalone
implementation, but as the continuation of a numerical-development line begun
in the 1990s.

## Historical timeline

```text
1994
doctoral research begun
        |
        v
1995–1996
archived working 2-D FEM source code
        |
        v
1997
2-D journal paper: FEM proof of concept and convection-dominated treatment
        |
        +---- POL-IMIS 1997 proceedings paper:
        |     published 3-D HEX8 FEM model
        |     + transient calculation
        |     + lateral convection
        |     + obstacle experiment
        |     + concentration-dependent P(S)
        |
        +---- archived August 1997 source:
              3-D transient HEX8 + Element-by-Element operator
        |
        v
1998
transient animation demonstrated during the public PhD defence
        |
        v
1999
published 3-D FEM / Pasquill comparison
        |
        v
later development
EbE + nonsymmetric Krylov solvers + Directional Residual stabilization
        |
        v
current Vanadis
3-D transient HEX8 + DR + EbE + OpenMP + CUDA + nonlinear P(S) + Q(t)
```

## Historical publications and documents

### 1997 — Chodorski and Pietrzyk

**M. Chodorski, M. Pietrzyk**  
*Propozycja zastosowania metody elementów skończonych do symulacji
rozprzestrzeniania się zanieczyszczeń w atmosferze*  
*Ochrona Powietrza i Problemy Odpadów*, vol. 31, no. 3, pp. 90–93, 1997.

Repository scan:

- `article_scan_nr 3_97.pdf`

### 1997 — POL-IMIS 3-D FEM paper

**M. Chodorski, Z. Malinowski, S. Słupek, A. Buczek**  
*Zastosowanie metody elementów skończonych do modelu rozprzestrzeniania się
zanieczyszczeń w atmosferze*  
Proceedings of the II POL-IMIS Conference, *Ocena wielkości imisji
zanieczyszczeń powietrza*, Szklarska Poręba, 19–22 June 1997,
PZITS publication no. 727, pp. 13–22.

Suggested repository filename for the preserved paper:

- `1997_POL_IMIS_3D_FEM.pdf`

If the surviving original word-processing file is also preserved, it should be
stored separately and identified explicitly as an archival source file rather
than as the published proceedings scan.

### 1998 — doctoral defence

**Public doctoral-defence protocol, 29 June 1998**  
Doctoral dissertation: *Modelowanie propagacji zanieczyszczeń w atmosferze*.

Repository document:

- `1998_doctoral_defence_protocol.pdf`

### 1999 — Chodorski and Malinowski

**M. Chodorski, Z. Malinowski**  
*Porównanie modelu Pasquilla z trójwymiarowym rozwiązaniem metodą elementów
skończonych. Zagadnienia dyfuzji atmosferycznej*  
*Ochrona Powietrza i Problemy Odpadów*, vol. 33, no. 1, pp. 5–9, 1999.

Repository scan:

- `article_scan_nr 1_99.pdf`

## Historical source code

Selected original source-code milestones are preserved under:

[`docs/history/source/`](source/)

They include:

- `1995_2d_pascal/` — early 2-D atmospheric-transport FEM implementation,
- `1997_3d_stationary_nonlinear/` — 3-D HEX8 nonlinear stationary prototype,
- `1997_3d_transient_ebe/` — early 3-D transient Element-by-Element branch,
- `1997_3d_transient_lu/` — alternative assembled band-LU transient branch,
- `ftn90_visualization/` — historical concentration-field visualization,
- `1998_defence_animation/` — software used for the transient animation workflow.

The historical source files are retained as archival software artifacts. They
should not be interpreted as current production code, and known historical
limitations are documented rather than silently corrected.

## Interpretation of the historical record

The surviving evidence should not be used to claim that every feature of the
current Vanadis formulation already existed in the earliest code.

In particular, the modern Directional Residual formulation, OpenMP and CUDA are
later developments.

What the historical publications and source code do document is a continuous
technical development sequence:

**2-D FEM atmospheric transport → 3-D HEX8 FEM → transient and nonlinear
experiments → Element-by-Element operator architecture → later nonsymmetric
Krylov solvers and DR → modern OpenMP/CUDA Vanadis.**

The distinction between published evidence and software evidence is also
important. The 1997 POL-IMIS paper independently documents the 3-D HEX8 model,
transient example, lateral convection, obstacle experiment and nonlinear
removal variant. The archived 1997 source code additionally documents the
Element-by-Element implementation and other internal numerical details not
fully described in the proceedings paper.

## History note

The bilingual Polish/English historical note in this directory was prepared
before the 1997 POL-IMIS proceedings paper was recovered and may therefore
contain an older, simplified chronology in which the 1999 paper appears as the
first clearly documented 3-D publication.

That interpretation is superseded by the evidence summarized in this README.
A future revision of the bilingual history note should incorporate the POL-IMIS
1997 publication and the corrected chronology.

## Preservation principle

Historical publications, source files and archival documents are preserved to
show the development of the numerical model as it actually occurred.

Original historical source files should remain unchanged wherever possible.
Corrected or modernized versions should be stored separately and clearly marked
as restored or derived versions.
