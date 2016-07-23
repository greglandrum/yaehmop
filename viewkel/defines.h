/*******************************************************

Copyright (C) 1995 Greg Landrum
All rights reserved

This file is part of yaehmop.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

********************************************************************/

/******************************

  defines.h

  This file contains the #defines used by everything

  Created by gL June 1995

*******************************/

/***
  Recent Edit History:
   10.05.98 gL:
     added support for arbitrary MO planes
   04.09.98 gL:
     added #define for ABS
   08.09.98 gL:
     #defines for atom and bond types
   08.09.98 gL:
     support for colored and/or shaded atoms
   24.09.98 gL:
     support for colored and/or shaded atoms upgraded
   18.10.98 gL:
     support for fatbands added.

***/

#ifndef _MY_DEFINES_
#define _MY_DEFINES_

#ifndef USE_BZERO
#define bzero(a, b) (memset((void *)(a), 0, (b)))
#define bcopy(a, b, c) (memcpy((void *)(b), (const void *)(a), (c)))
#endif

/* the number of gray's used */
#define NUM_COLORS 8

#define MAX_STR_LEN 5000
#define NORMAL_STR_LEN 240

#define NUM_TITLE_LINES 3

/* default tolerance for drawing bonds */
#define DEF_BOND_TOL 0.10

/* the font size used in Postscript printing */
#define TEXT_SIZE 12

#ifndef PI
#define PI 3.141592653589793
#endif
#ifndef M_PI
#define M_PI PI
#endif

/* the dimension of the matrices used ( 3 for 2D )*/
#define DIM 4

#define UP 0
#define DOWN 1
#define RIGHT 2
#define LEFT 3
#define BIGUP 10
#define BIGDOWN 11
#define BIGRIGHT 12
#define BIGLEFT 13

#define X_AX 0
#define Y_AX 1
#define Z_AX 2

/* the number of samples used in the parametric algorithms */
#define SAMPLES 8

/* mode options for main_mode */
#define NONE 0
#define ROT 1
#define SCALE 2
#define TRANS 3
#define CENT 4
#define CHOOSE 5

/* drawing modes for parametric surfaces */
#define ALL 1
#define FRONT 2
#define BACK 3

/*******
 default sizes for graphs
******/
#define DEF_GRAPH_X 200
#define DEF_GRAPH_Y 200

/*******
 default sizes for contour plots
******/
#define DEF_CONT_PLOT_X 200
#define DEF_CONT_PLOT_Y 200

/********
  we need some kind of default scaling thing so that graphics
   can be properly displayed on the Mac and other small
   screen displays, but still print the same size
   as if they came from an X display
*********/
#ifdef USING_THE_MAC
#define GRAPHICS_SCALE 0.75
#else
#define GRAPHICS_SCALE 1.0
#endif

/* default length of tic marks */
#define TIC_DIM 5

/* the dimensions of the various windows */
#define GWINWIDTH 700
#define GWINHEIGHT 600
#define ORTHOWIDTH 300
#define ORTHOHEIGHT 225
#define BUTWIDTH 220
#define BUTHEIGHT 600

/* the primitive types */
#define MOLECULE 0
#define GRAPH 2
#define PROP_GRAPH 3
#define BAND_GRAPH 4
#define WALSH_GRAPH 5
#define FMO_DIAGRAM 6
#define MO_SURF 7
#define LABEL 8
#define CONT_PLOT 9

/* types of contour plot */
#define CONT_GEN 0
#define CONT_FCO 1
#define CONT_MO 2

#define DOS_PROP 0
#define COOP_PROP 1

#define USE_MATRIX 0
#define STORE_MATRIX 1

/* number of initial particles / center */
#define INIT_PARTS_PER_CENTER 20

/**********

  defines for Tektronix emulation

**********/
#define TEK_XMAX 1024
#define TEK_YMAX 780
#define TEK_VCHAR 25
#define TEK_HCHAR 14
#define TEK_VTIC 11
#define TEK_HTIC 11

/*****
  defaults for FMO diagrams
*****/
#define FMO_DEGEN_TOL 1e-3
#define FMO_MIN_CONTRIB .10
#define FMO_LEVEL_WIDTH 20
#define FMO_LEVEL_THICKNESS 3
#define FMO_DEGEN_LEVEL_SKIP 5
#define FMO_ELECTRON_LENGTH 10

/* drawing modes for electron fillings in FMO diagrams */
#define FMO_FILL_NONE 0
#define FMO_FILL_HOMO 1
#define FMO_FILL_ALL 2

/* stuff to deal with delaunay polygonalizations */
#define NNSORT_POLYG 0
#define NNSORT_CONVEX_HULL 1
#define NNSORT_TOL 2

#ifdef INCLUDE_ADF_PLOTS
/* Max number of MOs */
#define MAX_MOS 512
/* the max number of AO's on a center. */
#define MAX_AOS 100
#else
/* Max number of MOs */
#define MAX_MOS 256
/* the max number of AO's on a center. */
#define MAX_AOS 16
#endif

/* default contraction coefficient for orbitals */
#define DEFAULT_CONTRACT 1.5
/* lookup table stuff */
#define DEFAULT_NUM_ENTRIES 10000
#define DEFAULT_LOOKUP_MIN 0.0
#define DEFAULT_LOOKUP_MAX 5.0

/**** more stuff for contour plots ****/
#define CONT_PLOT_PROPS_SCALE 0.5

#define INTERP_NOTHING 0 /* Kind of interpolations on contours. */
#define INTERP_CUBIC 1   /* Cubic spline interp. */
#define APPROX_BSPLINE 2 /* Bspline interpolation. */

#define LEVELS_AUTO 0             /* How contour levels are set */
#define LEVELS_INCREMENTAL 1      /* user specified start & incremnet */
#define LEVELS_DISCRETE 2         /* user specified discrete levels */
#define DEFAULT_NUM_OF_ZLEVELS 20 /* Some dflt values (setable via flags). */
#define DEFAULT_NUM_APPROX_PTS 5
#define DEFAULT_BSPLINE_ORDER 3

/********* stuff for MO contour plots ***********/
#define ORIENT_X 0
#define ORIENT_Y 1
#define ORIENT_Z 2
#define ORIENT_ARBITRARY 3

/***** toggles for atom and bond styles *********/
#define ATOM_PLAIN_FILL 0
#define ATOM_SHADE_FILL 1
#define ATOM_COLOR_FILL 2
#define ATOM_CSHADE_FILL 3

#define BOND_PLAIN 0
#define BOND_SHADE 1
#define BOND_CSHADE 2

#ifdef SUPPORT_COLOR_X
#define NUM_X_SHADES 10
#endif

#ifdef SUPPORT_FATBANDS
#define FATBANDS_SHADE 0
#define FATBANDS_LINE 1
#endif

/********

  stuff used by the buttons that's probably useful to have
  around elsewhere.

********/
#define MAX_MODES 10 /* any more than this is too many */
#define MAX_LINES 10 /* the maximum number of lines in a multi-line title */
#define MAX_ARGS 10  /* the maximum number of arguments a function can take */

#define TEXTOFF 10
#define BUTTONOFF 10

#define CURVE_SEG_LEN 50 /* length (in pixels) of sample linestyles shown */

/* this is stolen from steve gifford */
#ifndef PROTO
#if defined(_NO_PROTO) || defined(_alpha) || defined(MIPSEL)
#define PROTO(x) ()
#else /* defined(_NO_PROTO) || defined(__alpha) || defined(MIPSEL) */
#define PROTO(x) x
#endif /* defined(_NO_PROTO) || defined(__alpha) || defined(MIPSEL) */
#endif /* PROTO */

/********

  rint seems to do weird things on some systems,
  so I made my own....

*********/
#define rint(a) (floor((a) + 0.5))

#define ABS(__a__) ((__a__) > 0 ? (__a__) : -(__a__))

#endif
