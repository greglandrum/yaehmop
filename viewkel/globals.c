/*******************************************************
*      Copyright (C) 1995, 1998, 1999 Greg Landrum
*
*  This file is part of yaehmop.
*
*   This is free software.
* 
*  Permission is granted to modify, or otherwise fold, spindle, and mutilate this
*    code provided all copyright notices are left intact.
*
*  This code may be distributed to your heart's content, in whatever form,
*    provided no fee is charged for the distribution, all copyright notices are
*    left intact, and the source is distributed (without fee) along with any
*    binaries to anyone who requests it.
*
*  There are, of course, no warranties at all on this program.
*
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

