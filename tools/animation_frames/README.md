# Vanadis animation-frame generator

`vanadis_animation_frames.c` is a standalone batch postprocessor for generating
BMP visualization frames from Vanadis X/Y/Z plane-output files.

It is intended for time-dependent Vanadis calculations where many plane files
are produced and the resulting BMP sequence is later assembled into an
animation by an external tool.

The program does **not** modify the Vanadis solver or its output files.

## Supported input files

The program scans one directory for files named:

```text
X_?????_???.TXT
Y_?????_???.TXT
Z_?????_???.TXT
```

For example:

```text
X_00100_001.TXT
X_00200_001.TXT
Y_00100_010.TXT
Z_00100_010.TXT
```

The filename contains:

```text
<plane family>_<simulation time>_<plane index>.TXT
```

The program extracts the plane family (`X`, `Y`, or `Z`), simulation time, and
plane index automatically.

## Plane conventions

The current version follows the standard Vanadis test-grid ordering directly.

### X planes

An `X_?????_???.TXT` file represents a plane at constant Vanadis `X`.

Each record contains:

```text
Y  Z  concentration
```

The plane contains:

```text
20 x 20 nodes
```

`Y` varies fastest in the file.

These planes are rendered as horizontal Y-Z concentration maps.

### Y planes

A `Y_?????_???.TXT` file represents a plane at constant Vanadis `Y`.

Each record contains:

```text
Z  X  concentration
```

The plane contains:

```text
20 x 7 nodes
```

Vanadis `X` is the vertical coordinate and varies fastest.

### Z planes

A `Z_?????_???.TXT` file represents a plane at constant Vanadis `Z`.

Each record contains:

```text
Y  X  concentration
```

The plane contains:

```text
20 x 7 nodes
```

Vanadis `X` is again the vertical coordinate and varies fastest.

## Current test-grid geometry

The current implementation is tied to the standard Vanadis test grid:

```text
NAX = 7
NAY = 20
NAZ = 20
```

with the coordinate arrays embedded in the source code.

The Vanadis X levels are:

```text
0, 20, 35, 85, 135, 230, 420 m
```

The Y and Z coordinates use the standard 20-node test-grid coordinates from
`-4000 m` to `+4000 m`.

Therefore this version should be described as a batch renderer for the standard
Vanadis test cases, not as a fully general arbitrary-grid renderer.

## Rendering method

For every matching input file the program:

1. reads the nodal concentration field,
2. determines the appropriate color-scale maximum,
3. interpolates nodal values onto image pixels using bilinear Q4 shape
   functions,
4. clips negative concentrations to zero for visualization,
5. applies a small visualization cutoff,
6. maps concentration values to the Vanadis color palette,
7. draws coordinate axes and a concentration legend,
8. writes a 32-bit top-down BMP image.

The interpolation is performed from the nodal Vanadis field and is not simple
nearest-neighbour coloring.

### Visualization cutoff

For every individual input file, values below:

```text
0.1% of that file's positive maximum
```

are treated as zero when rendering the frame.

This cutoff is independent of the selected color-scale mode.

## Plot ranges

For X planes, the rendered horizontal view is:

```text
Y = -1500 ... +1500 m
Z = -1500 ... +1500 m
```

with a 500 x 500 pixel plot.

For Y and Z vertical sections, the rendered view is:

```text
horizontal coordinate = -600 ... +600 m
Vanadis X (height)     = 0 ... 420 m
```

The vertical plots preserve the same physical metres-per-pixel scale in both
directions, so the geometry is not stretched vertically.

## Color-scale modes

Three scale modes are available.

### `global`

```text
vanadis_animation_frames.exe OUTPUT global
```

One common concentration maximum is used for **all X, Y, and Z files at all
times**.

This is the strictest mode for comparing every generated frame directly.

### `family`

```text
vanadis_animation_frames.exe OUTPUT family
```

One common scale is used for the same plane family and plane index across all
times.

For example:

```text
X_*_001.TXT
```

uses one common scale for every simulation time, while:

```text
X_*_002.TXT
```

uses its own common scale.

This mode is especially useful for time animations of one fixed plane because
the same color represents the same concentration throughout the sequence.

### `local`

```text
vanadis_animation_frames.exe OUTPUT local
```

Each input file receives its own scale based on its own maximum.

This can show the internal structure of weak plumes more clearly, but colors
cannot be compared quantitatively between frames without reading the legend.

If the scale mode is omitted, the program uses:

```text
global
```

## Build

### GCC / MinGW / MSYS2

```bash
gcc -O3 -std=c11 -Wall -Wextra -pedantic \
    vanadis_animation_frames.c -lm \
    -o vanadis_animation_frames.exe
```

The original source can also be compiled under its historical filename:

```bash
gcc bmp_xyz_axes.c -O3 -o bmp_xyz.exe -lm
```

The program contains separate directory-scanning implementations for Windows
and POSIX systems.

## Usage

General form:

```text
vanadis_animation_frames DIRECTORY [global|family|local]
```

Examples:

```bash
vanadis_animation_frames.exe C:\tmp\output
```

```bash
vanadis_animation_frames.exe C:\tmp\output global
```

```bash
vanadis_animation_frames.exe C:\tmp\output family
```

```bash
vanadis_animation_frames.exe C:\tmp\output local
```

On Linux or in an MSYS2 shell, a directory can be supplied in the usual POSIX
form:

```bash
./vanadis_animation_frames ./OUTPUT family
```

## Output files

The BMP files are written into the same directory as the corresponding input
files.

The input filename is preserved and only the extension is changed:

```text
X_00100_001.TXT  ->  X_00100_001.bmp
X_00200_001.TXT  ->  X_00200_001.bmp

Y_00100_010.TXT  ->  Y_00100_010.bmp
Z_00100_010.TXT  ->  Z_00100_010.bmp
```

This naming scheme makes the resulting files naturally sortable by plane
family, simulation time, and plane index.

The generated BMP sequence can then be assembled into an animation with an
external video or GIF tool.

## Typical workflow

```text
Vanadis
   |
   +-- X_00100_001.TXT
   +-- X_00200_001.TXT
   +-- X_00300_001.TXT
   +-- ...
          |
          v
vanadis_animation_frames
          |
   +-- X_00100_001.bmp
   +-- X_00200_001.bmp
   +-- X_00300_001.bmp
   +-- ...
          |
          v
 external animation tool
```

The same workflow applies to Y and Z plane families.

## Notes

- Negative concentrations are clipped to zero for visualization.
- The original Vanadis text files are never modified.
- The 0.1% cutoff affects visualization only.
- The BMP output includes coordinate values and a concentration legend.
- `family` scaling is normally the most useful choice for animations of one
  fixed plane through time.
- The current implementation assumes the standard `7 x 20 x 20` Vanadis test
  grid and its embedded coordinate arrays.
