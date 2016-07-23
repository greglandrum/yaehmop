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

  graph_objects.h

  This file contains the definitions required to draw 2D graphs

  Created by gL June 1995

*******************************/

/***
  Recent Edit History:
   26.04.98 gL:
     added occupations field to FMO_fragment_type and FMO_diagram_type
  18.10.98 gL:
     support for fatbands added.

  06.03.99 gL:
     support for printing curve names added.

***/

#ifndef _GRAPH_OBJECTS_
#define _GRAPH_OBJECTS_

#ifndef _MY_DEFINES_
#include "defines.h"
#endif

#ifndef _BASIC_OBJECTS_
#include "basic_objects.h"
#endif

/****

  this structure is used to display the energy levels in FMO
  diagrams nicely (so that degeneracies are represented correctly)

****/
typedef struct FMO_level_type_def {
  char highest;
  int num_degen;
  int number;
  float num_electrons;
  float energy;

  float xmin, xmax, yloc;

  struct FMO_level_type_def *next;
} FMO_level_type;

/***********

  a fragment (for FMO)

***********/
typedef struct {
  char label[NORMAL_STR_LEN];
  char show_label;

  int num_orbs;
  float num_electrons;

  float *raw_energies;
  /* sometimes we want to specify occupations in the .FMO file */
  float *occups;

  FMO_level_type *levels;
} FMO_fragment_type;

/***********

  this is used to connect the levels in FMO diagrams

************/
typedef struct FMO_connect_type_def {
  int which_frag;
  char linestyle;
  float contrib;
  int fragment_level, main_level;
  struct FMO_connect_type_def *next;
} FMO_connect_type;

/*************

  used to display FMO diagrams

**************/
typedef struct {
  char filename[NORMAL_STR_LEN];
  char xlegend[NORMAL_STR_LEN], ylegend[NORMAL_STR_LEN];
  char title[NUM_TITLE_LINES][NORMAL_STR_LEN];
  char label[NORMAL_STR_LEN];
  char do_title;
  char show_data;
  char do_y_tics;
  char show_box;
  int electron_filling_mode;

  int num_frags;
  int tot_num_orbs;
  float tot_num_electrons;

  /* the width of the lines used to draw energy levels (in pixels) */
  int level_width;

  /* the thickness of the lines used to draw energy levels (in pixels) */
  int level_thickness;

  /* the length of lines used to indicate electron filling (in pixels) */
  int electron_length;

  int left_fragment, right_fragment;

  float max_y, min_y;
  float old_max_y, old_min_y;
  float tic_sep_y, num_tics_y, tic_start_y;

  /* used to determine what is considered a degeneracy */
  float degen_tol;

  FMO_fragment_type *frags;
  float *raw_energies;
  FMO_level_type *levels;
  float *chg_mat;

  /* sometimes we want to specify occupations in the .FMO file */
  float *occups;

  float min_contrib;
  char show_connects;
  FMO_connect_type *connects;

} FMO_diagram_type;

/********

  a graph

*********/
typedef struct {
  char filename[NORMAL_STR_LEN];
  int num_p;
  int num_curves;
  char xlegend[NORMAL_STR_LEN], ylegend[NORMAL_STR_LEN];
  char title[NUM_TITLE_LINES][NORMAL_STR_LEN];
  char do_x_tics, do_y_tics, do_title;
  char *curves_to_display;
  char *styles, *fills;
  char *curve_names;
  point_type2D *data, *raw_data;

  float min_x, max_x;
  float max_y, min_y;
  float old_min_x, old_max_x;
  float old_max_y, old_min_y;
  float tic_sep_x, num_tics_x, tic_start_x;
  float tic_sep_y, num_tics_y, tic_start_y;
} graph_type;

/*********

  for displaying properties calculations results

**********/
typedef struct {
  char filename[NORMAL_STR_LEN];
  graph_type *the_data;
  graph_type *the_integration;
  float max_x, max_y, min_x, min_y;
  float old_max_x, old_max_y, old_min_x, old_min_y;
  char type;
  char show_fermi;
  char integs_for_tics;
  float Fermi_E;
} prop_graph_type;

/**********

  special points (for band structures)

**********/
typedef struct {
  point_type loc;
  char name[80];
} special_point_type;

/************

  band structures

*************/
typedef struct {
  char filename[NORMAL_STR_LEN];
  int num_special_points;
  int points_per_line;
  int num_orbs;
  char show_fermi;
  float Fermi_E;
  special_point_type *special_points;
  graph_type *the_data;
#ifdef SUPPORT_FATBANDS
  int fatband_fill;
  int num_fatbands;
  char fatbands_on;
  float fatband_scale;
#endif

} band_graph_type;

/************

  walsh diagrams

*************/
typedef struct {
  char filename[NORMAL_STR_LEN];
  int num_orbs, num_p;
  graph_type *the_data;
  graph_type *total_E;
} walsh_graph_type;

/********

  a contour plot

*********/
#ifndef _CONTOUR_
#include "contour.h"
#endif

typedef struct {
  char filename[NORMAL_STR_LEN];
  int type;
  int num_p, disp_num_p;
  int num_curves;
  char xlegend[NORMAL_STR_LEN], ylegend[NORMAL_STR_LEN];
  char title[NUM_TITLE_LINES][NORMAL_STR_LEN];
  char do_x_tics, do_y_tics, do_title;
  char total_DOS_on_y;
  char *curves_to_display;
  char *styles;
  iso_curve_type *data, *raw_data;
  point_type2D *data2D, *raw_data2D;
  gnuplot_contour_type **contours;

  /* contour details */
  int num_levels, num_approx_pts, interp_kind, order, levels_kind;

  float min_x, max_x, step_x;
  float max_y, min_y, step_y;
  float max_DOS;
  int raw_num_x, raw_num_y;
  int num_x, num_y;
  float old_min_x, old_max_x;
  float old_max_y, old_min_y;
  float tic_sep_x, num_tics_x, tic_start_x;
  float tic_sep_y, num_tics_y, tic_start_y;
  float min_z, max_z;
  float old_min_z, old_max_z;
  float begin_conts;
  int num_conts;
  double *levels_list;

} contour_plot_type;

#endif
