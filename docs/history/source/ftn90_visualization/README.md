# Historical FTN90 concentration-field visualizer

This directory preserves a historical **Salford FTN90** graphics program used to visualize Vanadis concentration results under DOS.

The source is retained as a historical software artifact and is **not part of the current Vanadis production or visualization workflow**.

## Source file

- `1.f90`

The preserved program name is:

```fortran
program analiza_1
```

It is a standalone post-processing / visualization utility rather than a transport solver.

## Purpose

The program reads a two-dimensional concentration field extracted from a Vanadis three-dimensional calculation, interpolates the nodal values over the displayed domain, renders a colour concentration map in VGA graphics mode, and allows the user to inspect local coordinates and concentration values with the mouse.

In simplified form, the workflow is:

```text
Vanadis result plane
        ↓
read nodal coordinates and concentration values
        ↓
construct 2D cell connectivity
        ↓
interpolate concentration inside each cell
        ↓
map concentration to VGA colour levels
        ↓
display interactive concentration map
```

## Historical grid configuration

The saved source uses the same structured-grid dimensions associated with the historical 3D Vanadis case:

```fortran
nx1=6
nx2=13
nx3=22
```

The visualized plane contains:

```text
(nx2 + 1) x (nx3 + 1)
= 14 x 23
= 322 nodes
```

which is reflected by:

```fortran
dimension x1(322),x2(322),s(322)
```

The program builds:

```fortran
i_t(286,4)
```

representing the 286 four-node quadrilateral cells of the displayed plane.

## Input format

The input filename is requested interactively:

```fortran
call soua@('Podaj nazwe zbioru wejsciowego ZBIOR --->   ')
read(*,1)dane
```

The input file is then read as triples:

```text
coordinate 1
coordinate 2
concentration
```

using:

```fortran
read(9,3)x1(i),x2(i),s(i)
```

with the historical format:

```fortran
format (e10.4,2x,e10.4,2x,e10.4)
```

## Interpolation

The program does not simply display nodal values.

For every raster point it identifies the corresponding four-node cell and interpolates the concentration using bilinear Q4 shape functions:

```fortran
f_n(1)=0.25*(1-p_ksi)*(1-p_eta)
f_n(2)=0.25*(1+p_ksi)*(1-p_eta)
f_n(3)=0.25*(1+p_ksi)*(1+p_eta)
f_n(4)=0.25*(1-p_ksi)*(1+p_eta)
```

The concentration is then evaluated as:

```fortran
s=0
do i=1,4
   s=f_n(i)*s_w(i)+s
enddo
```

or:

\[
S(\xi,\eta)=\sum_{i=1}^{4}N_i(\xi,\eta)\,S_i.
\]

The auxiliary routine:

```fortran
wsp_ksi_eta
```

maps the physical raster coordinates into local element coordinates before interpolation.

## Rasterization

The displayed concentration field is rasterized onto:

```text
511 x 157
```

sample points:

```fortran
dimension stezenia(511,157)
```

The historical coordinate mapping is:

```fortran
x=-650+1300./156.*(...)
y=-250+4250./510.*(...)
```

which covers the saved physical cross-section of the historical test domain.

For each raster point, the interpolated concentration is assigned to a VGA palette index and drawn with:

```fortran
call set_pixel@(i+1,j+1,ij)
```

## Concentration colour scale

The program draws a discrete colour legend in:

```text
mg/m^3
```

with levels from approximately:

```text
0.000
0.002
0.003
0.004
0.005
0.006
0.007
0.008
0.010
0.012
0.014
0.016
0.018
0.020
```

The palette and thresholds are hard-coded for the historical visualization case.

This is appropriate for the preserved software artifact but should not be interpreted as the general plotting scale used by later Vanadis tools.

## Interactive mouse inspection

The program contains an interrupt-driven mouse handler:

```fortran
subroutine mouse_trap
```

and initializes the DOS/FTN90 mouse interface through:

```fortran
call initialise_mouse@
call set_trap@(mouse_trap,q,4)
call set_mouse_interrupt_mask@(11)
call set_mouse_bounds@(10,10,10+510,10+156)
call display_mouse_cursor@
```

When the left mouse button is pressed, the program obtains the current mouse position and converts it back to physical coordinates.

It then displays:

```text
X=
Y=
S=
```

on the screen.

The concentration value is taken from the already interpolated raster field:

```fortran
stezenia(ii-9,ii2-9)
```

The right mouse button exits the program and restores text mode.

This makes `1.f90` an **interactive scientific visualization tool**, not merely a static plot generator.

## FTN90 / DOS graphics environment

The program depends heavily on Salford FTN90 graphics and DOS extensions, including routines such as:

```text
vga@
set_video_dac_block@
set_video_dac@
set_pixel@
fill_rectangle@
draw_text@
initialise_mouse@
set_trap@
set_mouse_interrupt_mask@
set_mouse_bounds@
display_mouse_cursor@
get_mouse_position@
text_mode@
```

These routines are compiler/environment-specific.

The source is therefore not expected to compile with standard modern GNU Fortran without replacing the graphics and mouse layer.

## Relation to the defence-animation software

This program should be distinguished from the separate animation pipeline preserved in:

```text
../1998_defence_animation/
├── xb.f90
├── run1b.f90
└── vanadis1.webm
```

`1.f90` is an **interactive visualizer for a single concentration field**, whereas the defence-animation software processes and displays a sequence of transient fields.

The programs are nevertheless part of the same historical FTN90/DOS visualization environment developed around Vanadis.

## Known historical implementation notes

The program is working historical research software, not a modern plotting package.

Several implementation details reflect that context.

For example, the colour-threshold sequence contains:

```fortran
if (stezenia(...).le.0.001) ij=23
if (stezenia(...).ge.0.002) ij=24
```

which leaves the open interval:

```text
0.001 < S < 0.002
```

without an explicit new assignment in that sequence.

Also, the cell search uses upper bounds of the form:

```fortran
x.lt....
y.lt....
```

so the exact outermost upper boundary requires care.

These points should be regarded as historical implementation details. The original source should remain unchanged.

## Preservation policy

`1.f90` should be preserved verbatim.

If a modern equivalent is ever created, it should be placed separately, for example:

```text
ftn90_visualization/
├── README.md
├── original/
│   └── 1.f90
└── restored/
    └── concentration_viewer_modern.f90
```

or implemented in a modern visualization language while retaining the original source as the historical artifact.

## Integrity

SHA-256 of the preserved source:

```text
7e8b5a12813e112011a8023a5a4c8c22506267e3af361570beb639521cb7bb9f  1.f90
```

File size:

```text
8533 bytes
```

## Significance

`1.f90` documents that the historical Vanadis development included its own dedicated visualization software rather than relying only on raw numerical output.

The program combines:

- reading Vanadis concentration-plane results,
- four-node finite-element interpolation,
- rasterization of the continuous field,
- a discrete scientific colour scale,
- DOS/VGA graphics,
- and interactive mouse-based interrogation of the computed concentration field.

Together with the transient-animation programs, it shows that the early Vanadis environment already included a complete numerical workflow from simulation to graphical inspection of the results.
