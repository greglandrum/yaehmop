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
     added near plane clipping toggle
   30.05.98 gL:
     added polyhedron outline toggle
   08.09.98 gL:
     support for colored and/or shaded atoms
   18.01.99 gL:
     added dump_grid toggle
   29.01.99 gL:
     added global_read_on toggle

***/

/****************************************************************************
*
*                   globals.c
*
* These are the global variable definitions
*
*****************************************************************************/
#include "viewkel.h"

int mainmode;
char quit, orthoviewson, projviewon, drawcameraon;
char non_blocking_event_loop=0;
char animating_molecule=0;
char tgifon,redo_projection=1;
head_type *head, *parenthead;
matrix_type *stack, *orthoscalemat, *mainortho, *ident;
object_type *whichobj, *parentobj, *newparent;
point_type globalmin, globalmax;
camera_type *camera;
shead_type *shead;
button_win_type *button_wins;
point_type localmin, localmax, globalmin, globalmax;

int num_selected;

char fill_projections;
char near_plane_clipping_on,outline_polyhed_on;
char dump_grids_on;
PS_options_type PS_options;

#ifdef X_GRAPHICS
/* XWindows globals */
Display *disp;
GC graphgc, bigtextgc, smalltextgc, blackgc, colorgc,colorgc2,gray1gc,gray2gc;
GC  gray3gc,graygc[NUM_COLORS],widegc;
int screen_depth;
Window graph_win, gwin, butwin,root;
Pixmap gpix;
XFontStruct *big_font, *small_font;
int screen;
unsigned long fcolor, bcolor;

char useButtons=1;
#endif
int g_xmax, g_ymax;
int but_xmax,but_ymax;
int g_color;
int g_error;
float mouse_rotn_scale=10;

FILE *psfile,*Tek_file=0;
char doing_ps=0;
#ifdef X_GRAPHICS
char doing_X=0;
char doing_tek=0;
#endif

#ifdef MAC_GRAPHICS
char doing_Mac=0;
#endif
float radius_multiple=1;


MO_info_type *gMO_info;
MO_center_list_type *gcenters;
int gnum_centers;
MO_surface_type *gsurf;
char global_read_on;
FILE *rayfile;

#ifdef X_GRAPHICS
char refresh_all_colormaps;
#endif

