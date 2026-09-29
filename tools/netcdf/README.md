# Vanadis NetCDF converter

`vanadis2netcdf.f90` converts a sequence of Vanadis `T_XXXXX.TXT` files into a
single compressed NetCDF-4 dataset.

The converter is a standalone postprocessor. No change to the Vanadis solver is
required.

## Purpose

Vanadis normally writes one ASCII file for each output time. Each record
contains:

```text
X  Y  Z  concentration
```

For a time-dependent calculation this produces files such as:

```text
T_00100.TXT
T_00200.TXT
T_00300.TXT
...
T_02000.TXT
```

`vanadis2netcdf` combines these files into one dataset containing spatial
coordinates, simulation time, and the complete concentration field.

The resulting file can be opened by standard NetCDF-aware scientific tools,
for example NASA Panoply, xarray, and other NetCDF software.

## Grid detection

The grid dimensions do **not** need to be supplied manually.

The converter reads the first available `T_XXXXX.TXT` file and infers
`nx`, `ny`, and `nz` from the coordinates.

The current converter expects the structured Vanadis tensor-product ordering:

```text
X varies fastest, then Y, then Z
```

The coordinate spacing does not need to be uniform.

Every subsequent time file is checked against the grid from the first file. If
the node ordering or coordinates differ, the converter stops instead of writing
an inconsistent dataset.

## Time-file selection

Usage:

```text
vanadis2netcdf [input_dir] [t_start] [t_end] [output.nc]
```

Defaults:

```text
input_dir  = OUTPUT
t_start    = 100
t_end      = 2000
output.nc  = vanadis_00100_02000.nc
```

The converter scans every integer second in the requested range but reads only
files that actually exist. Therefore, no output interval has to be supplied.

For example, if the directory contains only:

```text
T_00100.TXT
T_00200.TXT
...
T_02000.TXT
```

the output will contain 20 time levels.

Example:

```bash
./vanadis2netcdf.exe . 100 2000 vanadis.nc
```

## NetCDF structure

For the included example, the generated dataset has the logical structure:

```text
x(nx)
y(ny)
z(nz)
time(time)

concentration(time, z, y, x)
```

The coordinate variables use metres and `time` uses seconds. Concentration is
stored as:

```text
mg m-3
```

The concentration field is written as NetCDF `float` data with NetCDF-4
compression enabled (`shuffle` plus deflate level 2).

The `time` dimension is unlimited.

Negative concentration values present in the Vanadis output are preserved.
This program is a format converter, not a numerical filter.

### Negative-value diagnostic

Small negative concentration values caused by numerical oscillations are retained
in the original Vanadis output and are transferred unchanged to the NetCDF dataset.
The standard concentration bitmap uses a display scale intended for the physically
relevant positive concentration range and therefore does not emphasize these small
undershoots. A separate **small oscillations view** can be used to visualize them
explicitly. The same negative values remain visible in standard NetCDF viewers such
as NASA Panoply, which displays the actual data range rather than filtering or
clipping the concentration field.

The file records Vanadis model coordinates only. No geodetic coordinate
reference system (CRS) is currently encoded.

## Build

The program requires the NetCDF C library and the NetCDF-Fortran interface.

### Windows with MSYS2 UCRT64

Install the required packages:

```bash
pacman -S mingw-w64-ucrt-x86_64-gcc-fortran \
          mingw-w64-ucrt-x86_64-netcdf \
          mingw-w64-ucrt-x86_64-netcdf-fortran
```

Check that the tools are available:

```bash
which gfortran
which nf-config
```

Typical result:

```text
/ucrt64/bin/gfortran
/ucrt64/bin/nf-config
```

Compile:

```bash
gfortran -O2 vanadis2netcdf.f90 \
    $(nf-config --fflags) \
    $(nf-config --flibs) \
    -o vanadis2netcdf.exe
```

### Linux

With NetCDF-Fortran installed:

```bash
gfortran -O2 vanadis2netcdf.f90 \
    $(nf-config --fflags) \
    $(nf-config --flibs) \
    -o vanadis2netcdf
```

## Verify the output

The NetCDF command-line tools can be used to inspect the generated file:

```bash
ncdump -h vanadis.nc
```

To inspect the time coordinate:

```bash
ncdump -v time vanadis.nc
```

To also display storage/chunking/compression metadata:

```bash
ncdump -hs vanadis.nc
```

## Example files

The `example/` directory contains:

- `vanadis_T_files_0100_2000.zip` — the example Vanadis
  `T_XXXXX.TXT` files from `100 s` through `2000 s`;
- `vanadis.nc` — NetCDF-4 file generated from those example files;
- `Panoply_sample.PNG` — example visualization of `vanadis.nc` in NASA
  Panoply.

For the supplied test case the converter detects:

```text
nx = 7
ny = 20
nz = 20
```

with 20 output times from `100 s` to `2000 s`.

A ground-level view can be obtained in Panoply by selecting the appropriate
Vanadis X level (for the included example, `X = 0 m`) and plotting the remaining
Vanadis Y-Z plane.

## Example workflow

```text
Vanadis
   |
   +-- T_00100.TXT
   +-- T_00200.TXT
   +-- ...
   +-- T_02000.TXT
          |
          v
   vanadis2netcdf
          |
          v
      vanadis.nc
          |
          v
   Panoply / xarray / other NetCDF tools
```

This keeps the solver output simple while providing a standard exchange format
for downstream analysis and visualization.
