# 1998 Vanadis doctoral-defence animation

This directory preserves the historical FTN90/DOS animation software associated with the transient Vanadis results demonstrated during Marek Chodorski's public doctoral defence in **1998**.

The material is retained as historical documentation and is **not part of the current Vanadis production or visualization workflow**.

## Files

```text
xb.f90
run1b.f90
vanadis1.webm
```

### `xb.f90`

Historical Salford FTN90 program used to generate animation frames from transient Vanadis concentration results.

### `run1b.f90`

Historical Salford FTN90 program used to preload and display the generated frames as an animation.

### `vanadis1.webm`

A modern screen capture of the preserved historical animation running inside a virtual machine.

It is **not an original recording from the 1998 doctoral defence**.

## Historical context

During the 1998 public defence of the doctoral dissertation on atmospheric-pollution propagation modelling, transient Vanadis results were demonstrated in animated form.

The surviving FTN90 source code documents how that visualization workflow operated.

The original programs were run directly under DOS on contemporary hardware. The animation was displayed smoothly in that environment.

The preserved `vanadis1.webm` file was recorded much later from the same class of historical software running in a virtual machine, so its playback smoothness and timing should not be treated as an exact reproduction of the original DOS presentation.

## Animation workflow

The preserved software implements a two-stage pipeline:

```text
Vanadis transient result files
        ↓
xb.f90
        ↓
interpolated VGA concentration fields
        ↓
PCX frame sequence
        ↓
run1b.f90
        ↓
preloaded frames in memory
        ↓
animated transient concentration display
```

This separation was important for smooth playback on the original DOS system.

## `xb.f90` — frame generator

The program reads concentration results for successive simulation times.

The preserved configuration uses:

```fortran
PARAMETER (Nx1=6,Nx2=13,nx3=22,iczas=1980)
```

The generated animation therefore covers a transient sequence up to approximately:

```text
1980 s
```

in increments of:

```text
5 s
```

The program generates filenames such as:

```text
X1_0005.pcx
X1_0010.pcx
X1_0015.pcx
...
X1_1980.pcx
```

## Interpolation of the FEM result

For each time level, `xb.f90` reads a two-dimensional concentration plane and interpolates the nodal solution over the display raster.

The program uses standard four-node bilinear shape functions:

```fortran
f_n(1)=0.25*(1-p_ksi)*(1-p_eta)
f_n(2)=0.25*(1+p_ksi)*(1-p_eta)
f_n(3)=0.25*(1+p_ksi)*(1+p_eta)
f_n(4)=0.25*(1-p_ksi)*(1+p_eta)
```

and evaluates the concentration from:

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

The display is therefore based on interpolation of the finite-element field rather than on nearest-node colouring.

## Raster and physical domain

The concentration field is rasterized to:

```text
511 x 157
```

display points.

The saved coordinate mapping covers approximately:

```text
-650 to 650
```

in one displayed direction and:

```text
-250 to 4000
```

in the other.

These values belong to the historical test case.

## Concentration colour mapping

The interpolated values are assigned to discrete VGA colour levels.

The preserved thresholds include approximately:

```text
0.001
0.003
0.004
0.005
0.006
0.007
0.008
0.010
0.012
0.016
0.020
0.024
0.028
```

The generated raster is drawn with:

```fortran
call set_pixel@(i+1,j+1,ij)
```

and captured into an in-memory graphics block before being written as PCX.

## PCX frame generation

After drawing each field, `xb.f90` captures the VGA screen region:

```fortran
call get_screen_block@(10,10,520,166,buffer)
```

and writes it to a PCX file using:

```fortran
call screen_block_to_pcx@(dane_pcx(i_iczas),buffer,error_code)
```

This produces the frame sequence used by the animation player.

## `run1b.f90` — animation player

The animation player first constructs the same PCX filename sequence and enters VGA mode.

Before playback, it loads all required frames:

```fortran
do i_czas=5,iczas,5
   call pcx_to_screen_block@(
     & dane_pcx(i_czas),
     & buffer(i_czas),
     & palette,
     & error_code(i_czas))
enddo
```

This means the animation is **preloaded into memory before playback**.

The design avoids repeated disk access during the displayed sequence and helps explain why the animation was smooth when run directly under DOS.

## Playback

During animation, each preloaded frame is restored with:

```fortran
call restore_screen_block@(10,10,buffer(i_czas),0,
     & error_code(i_czas))
```

The program then updates the displayed simulation time:

```text
czas: 0005s
czas: 0010s
...
```

using:

```fortran
call draw_text@(czas(i_czas),220,250,23)
```

The playback loop includes:

```fortran
call sleep@(0.1)
```

corresponding nominally to approximately 10 displayed frames per second.

## Display legend

The player also draws a concentration legend in:

```text
mg/m^3
```

with values extending up to approximately:

```text
0.028 mg/m^3
```

for the preserved case.

The legend, palette and concentration thresholds are hard-coded for this historical demonstration.

## FTN90 / DOS graphics environment

Both programs depend on Salford FTN90 graphics extensions, including routines such as:

```text
vga@
set_video_dac_block@
set_video_dac@
set_pixel@
get_screen_block@
screen_block_to_pcx@
pcx_to_screen_block@
restore_screen_block@
fill_rectangle@
draw_text@
sleep@
text_mode@
```

These are not standard Fortran procedures.

The programs are therefore not expected to compile with a standard modern Fortran compiler without replacement of the graphics layer.

## Preserved video

`vanadis1.webm` is a modern screen recording of the historical visualization software running inside a virtual machine.

The recording is useful because it shows the actual appearance and dynamics of the preserved visualization system, including:

- the developing concentration plume,
- the concentration colour scale,
- the simulation-time counter,
- and the general presentation style of the historical Vanadis animation.

However:

> **The WebM file is not an original 1998 recording.**

The original defence demonstration was run directly under DOS and was smoother than the later virtual-machine capture.

## Relationship to `1.f90`

The separate historical program:

```text
../ftn90_visualization/1.f90
```

is an interactive concentration-field viewer for a single field.

By contrast:

```text
xb.f90 + run1b.f90
```

form the transient-animation pipeline used to generate and play a sequence of concentration fields.

Together, these programs document the historical post-processing environment developed around Vanadis.

## Known historical implementation note

In `xb.f90`, the error branch contains:

```fortran
erro_code=1
```

instead of:

```fortran
error_code=1
```

This appears to be a typographical error in an error-handling path.

The original source should remain unchanged.

## Preservation policy

All files in this directory should be preserved as historical artifacts.

If modernized versions of the programs are ever created, keep them separately, for example:

```text
1998_defence_animation/
├── README.md
├── original/
│   ├── xb.f90
│   ├── run1b.f90
│   └── vanadis1.webm
└── restored/
    └── ...
```

The historical source should not be silently rewritten or reformatted.

## Integrity

SHA-256 values of the preserved files:

```text
25fe6e24ead58161bb0a5c0bb3b8bc3684de3ff5d742412b11c3140ef100b7c6  xb.f90
1de9f2b86d29b2671813e4d9d09a7444aa8400287f8fa1df054fa997478d9fe3  run1b.f90
f63e085743abd136e8d361827a8198a0f1f0e729cd020829508bafd6f77e8822  vanadis1.webm
```

## Significance

This directory preserves more than a plotting utility.

It documents a complete historical transient-visualization workflow:

- time-dependent Vanadis FEM results,
- interpolation onto a display raster,
- generation of successive PCX frames,
- preloading of frames into memory,
- smooth DOS/VGA playback,
- concentration legend,
- and simulation-time display.

Together with the surviving 1997 transient solver source and the 1998 doctoral-defence documentation, these files provide direct evidence that the early Vanadis environment included both transient simulation and dedicated dynamic visualization of the calculated concentration field.
