# Vanadis bitmap postprocessor

`vanadis_bitmap.c` is a small standalone postprocessor for a structured 2-D
Vanadis output plane. It converts nodal concentration values from a text file
into a 32-bit BMP image.

The program is intentionally separate from the Vanadis solver. It provides a
simple native visualization path for one selected output plane without adding
graphics or image-writing code to the solver itself.

## Input format

The input file contains one node per line:

```text
coord1  coord2  concentration
```

The nodes must form a rectangular structured grid. `coord1` must vary fastest,
followed by `coord2`.

For the included example, the input file is:

```text
example/X_02000_001.TXT
```

This is the Vanadis ground-level plane for the example calculation at
`t = 2000 s`.

## Method

The postprocessor:

1. reads all nodal coordinates and concentrations,
2. detects the 2-D grid dimensions automatically,
3. checks that the grid is rectangular and consistently ordered,
4. locates the containing grid cell for every output pixel,
5. evaluates bilinear Q4 shape functions inside that cell,
6. interpolates the concentration from the four nodal values,
7. clips negative interpolated concentrations to zero for visualization,
8. maps the result to the historical Vanadis color palette,
9. writes a 32-bit BMP image.

The interpolation is therefore performed from the Vanadis nodal field; the
bitmap is not obtained by simple nearest-neighbour coloring.

## Build

### GCC / MinGW / MSYS2

```bash
gcc -O2 -std=c11 -Wall -Wextra -pedantic \
    vanadis_bitmap.c -lm -o vanadis_bitmap.exe
```

On toolchains where the math library is linked automatically, `-lm` may be
omitted.

## Usage

With no arguments, the historical example settings are used:

```bash
./vanadis_bitmap.exe
```

Equivalent explicit command:

```bash
./vanadis_bitmap.exe \
    X_02000_001.TXT \
    vanadis_ground.bmp \
    -1500 1500 \
    -1500 1500 \
    500 500
```

General form:

```text
vanadis_bitmap INPUT OUTPUT XMIN XMAX YMIN YMAX WIDTH HEIGHT
```

The grid size itself does **not** need to be supplied. It is inferred from the
coordinates in the input file.

The requested rendering bounds must lie inside the input grid.

## Example files

The `example/` directory contains:

- `X_02000_001.TXT` — example Vanadis 2-D ground-level output at `t = 2000 s`;
- `vanadis_ground_raw.bmp` — direct bitmap output from the postprocessor;
- `vanadis_ground_with_scale.bmp` — presentation version of the same result
  with an added concentration scale.

The `with_scale` image is provided for interpretation and presentation. The
current `vanadis_bitmap.c` writes the raw raster image itself and does not add
text, axes, or a numerical legend.

## Notes

- Negative interpolated concentrations are clipped to zero **only for bitmap
  visualization**. The original Vanadis output file is not modified.
- The color scale uses the maximum positive nodal concentration found in the
  input file.
- The program supports non-uniform but strictly monotonic coordinate spacing.
- The bitmap renderer is intended as a lightweight Vanadis-specific
  postprocessor. For exchange with general scientific visualization software,
  see `../netcdf/`.
