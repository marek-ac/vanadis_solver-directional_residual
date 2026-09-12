# Vanadis 3D — v2026.3.1 development line

Vanadis 3D is a three-dimensional finite-element model for transient atmospheric pollutant transport.

This directory contains the current **v2026.3.1 development baseline**.

The v3.1 line extends Vanadis beyond the regular-grid v3.0 configuration toward terrain-following meshes, spatially and temporally variable physical fields, improved numerical diagnostics, and more robust integration of distorted HEX8 elements.

## Main features of v2026.3.1

- 3D transient convection-diffusion-reaction transport
- HEX8 finite elements
- Directional Residual (DR) stabilization
- nonlinear decay coefficient `P(S)`
- time-dependent source term
- locally evaluated fields:
  - `v(x,y,z,t)`
  - `K(x,y,z,t)`
  - `alpha(x,y,z,t)`
- terrain-following mesh support
- terrain-consistent test velocity field
- 27-point (`3 x 3 x 3`) Gauss integration for HEX8 volume terms
- 9-point (`3 x 3`) Gauss integration for Q4 boundary terms
- full mass-balance diagnostics
- OpenMP element assembly
- CUDA DBCG linear solver
- GPU element coloring option

## Main files

### Solver

`VANADIS.FOR`

Main Fortran source code for Vanadis v2026.3.1.

### CUDA solver

`cudatest.cu`

Current default CUDA DBCG implementation.

Additional comparison versions:

- `cudatest_color.cu` — coloring-enabled version
- `cudatest_no_color.cu` — non-coloring reference version

### Configuration

`vanadis.cfg`

Runtime configuration file containing the principal physical parameters.

### Build

`Makefile`

Build configuration for the Fortran/CUDA version.

### Post-processing

`VANADIS_BMP.exe`

Utility used to generate bitmap concentration maps from Vanadis output.

## Documentation

See the `DOC/` directory.

### English

- `Vanadis_input_data_and_tailoring_guide_v3_1_EN.pdf`
- `Vanadis_how_to_analyze_full_mass_balance_v3_1_EN.pdf`
- `VANADIS_terrain_and_velocity_field_guide_EN.pdf`

### Polish

- `Vanadis_instrukcja_wprowadzania_danych_v3_1_PL.pdf`
- `Vanadis_jak_analizowac_pelny_bilans_masy_v3_1_PL.pdf`
- `VANADIS_instrukcja_teren_i_pole_predkosci_PL.pdf`

For terrain-following simulations, the recommended starting document is:

`VANADIS_terrain_and_velocity_field_guide_EN.pdf`

or its Polish counterpart:

`VANADIS_instrukcja_teren_i_pole_predkosci_PL.pdf`

## Terrain-following geometry

In Vanadis the first coordinate, `x`, is the vertical coordinate.

Terrain is therefore represented as

`x_g = x_g(y,z)`.

The terrain-following mesh is generated from nodal ground elevations stored in the horizontal `y-z` grid.

The current development baseline contains a controlled terrain validation case based on the nodal elevations:

`0, 1, 4, 7 m`.

The mesh geometry and the terrain-consistent test velocity field use the same terrain representation.

## Velocity fields over terrain

For the controlled terrain test, the velocity normal to the ground is corrected so that

`v . n = 0`

at the terrain surface.

For real applications, Vanadis may instead use externally supplied velocity fields, for example from:

- measurements,
- meteorological models,
- CFD calculations,
- interpolated observational datasets.

An externally supplied velocity field that already represents the flow over the actual terrain should normally be used directly rather than being artificially corrected by the terrain-slope test formula.

## Validation

Reference results for the current terrain baseline are stored in:

`validation/terrain_0_1_4_7/`

The directory contains:

- `MASS_BALANCE.TXT`
- `vanadis.log`
- `X_02000_001.bmp`

This case is intended as a reproducible validation checkpoint for the v2026.3.1 development line.

## Current status

v2026.3.1 is currently a **development version**, not yet the final v3.1 release.

The present baseline has verified:

- flat-terrain regression,
- terrain-following mesh generation,
- terrain-consistent velocity treatment,
- full mass balance,
- nonlinear `P(S)`,
- CUDA DBCG operation,
- GPU coloring.

Further validation is planned for stronger terrain deformation and comparison of 8-point and 27-point volume integration on warped meshes.

## Coordinate and unit convention

Vanadis uses SI units internally.

The first coordinate is vertical:

- `x` — vertical coordinate
- `y`, `z` — horizontal coordinates

Concentration is stored internally in:

`kg/m^3`

and converted to:

`mg/m^3`

only when results are written for output.

## Author

Marek Chodorski

Vanadis 3D — atmospheric pollutant transport FEM solver