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

/***
  Recent Edit History:
   26.04.98 gL:
     added near plane clipping toggle extern;
   30.05.98 gL:
     added polyhedron outline toggle extern;
   08.09.98 gL:
     support for colored and/or shaded atoms
   18.01.99 gL:
     added dump_grid toggle extern
   29.01.99 gL:
     added global_read_on toggle extern
   29.06.1999 gL:
     added #define for single_float in order to be able to get
     real floats when that is required
***/

/******************** INCLUDES *************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef X_GRAPHICS
/* These are the XWindows include files */
#include <X11/Xlib.h>
#include <X11/Xos.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>
#else

#ifdef SUPPORT_VOLUMES
#include <fcntl.h>
#endif

/* define some things that might be useful */
typedef struct { short x, y; } XPoint;

typedef int Window;
#endif

#define FATAL_BUG(a) fatal_bug((a), __FILE__, __LINE__)

#ifndef MEM_DEBUG
#define D_MALLOC(_a_) malloc(_a_)
#define D_CALLOC(_a_, _b_) calloc(_a_, _b_)
#define D_REALLOC(_a_, _b_) realloc(_a_, _b_)
#define D_FREE(_a_) free(_a_)
#else
#define D_MALLOC(_a_) d_malloc(_a_, __FILE__, __LINE__)
#define D_CALLOC(_a_, _b_) d_calloc(_a_, _b_, __FILE__, __LINE__)
#define D_REALLOC(_a_, _b_) d_realloc(_a_, _b_, __FILE__, __LINE__)
#define D_FREE(_a_) d_free(_a_, __FILE__, __LINE__)
#endif

#ifndef USE_FLOATS
#define float double
#endif

#include "defines.h"
#include "buttons.h"
#include "basic_objects.h"
#include "graph_objects.h"
#include "3D_objects.h"
#include "ps_stuff.h"
#include "matrix_ops.h"
#include "contour.h"
#include "labels.h"
#include "polyhed.h"

#ifdef USING_THE_MAC
#include "Mac_Fopen.h"
#include <unistd.h>
#include <string.h>

#endif
/******************** TYPE DECLARATIONS ****************************/

/******

  this is used to allow individual objects to have their own
  button windows (in X) or menu bars (on the Mac)
  to allow options and such things.

********/
typedef struct button_win_type_def button_win_type;
struct button_win_type_def {
  button_type *button_list;
  int xdim, ydim;
#ifdef X_GRAPHICS
  Window which_win;
#endif
#ifdef MAC_GRAPHICS
  MenuHandle theMenu;
  short menuID;
#endif
  button_win_type *child;
  button_win_type *next;
};

/*************

  the primitive

**************/
typedef struct {
  int which;
  molec_type *molec;
  graph_type *graph;
  prop_graph_type *prop_graph;
  band_graph_type *band_graph;
  walsh_graph_type *walsh_graph;
  param_surf_type *p_surf;
  MO_surface_type *MO_surf;
  FMO_diagram_type *FMO_diagram;
  contour_plot_type *cont_plot;
  label_type *label;
  /* this allows option windows for each primitive */
  button_win_type *but_win;
} prim_type;

/*********

  the key frame node

**********/
typedef struct snap_type_def {
  float xc, yc, zc, xt, yt, zt;
  float xs, ys, zs, xr, yr, zr;
  struct snap_type_def *next;
} snap_type;

/**********

  this is the keyframe forest node

**********/
typedef struct shead_type_def {
  point_type lf, la, vup; /* the camera parameters */
  float foclength;        /* hsize and yaspect do not need to be interpolated */
  int which;              /* what the keyframe number is */
  snap_type *snap;
  struct shead_type_def *next;
} shead_type;

/**********

  here are the object nodes

***********/
typedef struct object_type_def {
  point_type cent;       /* the center of the object */
  point_type trans;      /* the translations */
  point_type scale;      /* the scale */
  point_type rot;        /* the rotations */
  point_type bmin, bmax; /* the bounding box */
  prim_type *prim;       /* the primitive at this node */
  struct object_type_def *children[25];
} object_type;

/***********

  the forest nodes

************/
typedef struct head_type_def {
  object_type *obj;

  struct head_type_def *next;
} head_type;

/**********

  the camera

***********/
typedef struct camera_type_def {
  point_type lf, la, vup;
  float foclength, hsize, yaspect;
} camera_type;

/************* GLOBAL VARIABLES EXTERN DECLARATIONS ***********/
extern matrix_type *stack, *orthoscalemat, *mainortho, *ident;
extern head_type *head, *parenthead;
extern int mainmode;
extern char non_blocking_event_loop;
extern char animating_molecule;
extern char quit, orthoviewson, projviewon, drawcameraon;
extern char redo_projection;
extern object_type *whichobj, *parentobj, *newparent;
extern point_type globalmin, globalmax;
extern camera_type *camera;
extern shead_type *shead;
extern button_win_type *button_wins;
extern point_type localmin, localmax, globalmin, globalmax;

extern int num_selected;

extern char fill_projections;
extern char near_plane_clipping_on, outline_polyhed_on;
extern char dump_grids_on;
extern PS_options_type PS_options;

/* Graphics globals */
extern unsigned long fcolor, bcolor;
extern int g_xmax, g_ymax;
extern int g_color;
extern int g_error;
extern float mouse_rotn_scale;

#ifdef X_GRAPHICS
extern Display *disp;
extern GC graphgc, bigtextgc, smalltextgc, blackgc, colorgc, colorgc2;
extern GC gray1gc, gray2gc;
extern GC gray3gc, graygc[NUM_COLORS], widegc;
extern int screen_depth;
extern Window root, graph_win, gwin, butwin;
extern Pixmap gpix;
extern XFontStruct *big_font, *small_font;
extern int screen;
extern char doing_X;
extern char doing_tek;

extern char useButtons;
#endif

extern char doing_ps;
extern FILE *psfile;
extern FILE *Tek_file;

#ifdef MAC_GRAPHICS
extern char doing_Mac;
#endif

/******

  these are globals that really shouldn't be global, but are
  for the implicit polygonalization thing.

*******/
extern MO_info_type *gMO_info;
extern MO_center_list_type *gcenters;
extern int gnum_centers;
extern MO_surface_type *gsurf;
extern char global_read_on;
extern FILE *rayfile;

#ifdef X_GRAPHICS
extern char refresh_all_colormaps;
#endif

#include "prototypes.h"
