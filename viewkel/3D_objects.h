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

  3D_objects.h

  This file contains the definitions required to draw 3D
  objects (molecules and surfaces)

  Created by gL June 1995

*******************************/

/***
  Recent Edit History:
   28.02.98 gL:
     modifications to reduce memory usage for MO storage.
   10.05.98 gL:
     additions to support contouring of arbitrary planes.
   18.05.98 gL:
     additions to molec_type to store polyedra centers and cut_offs
   04.09.98 gL:
     changed is_selected field of atom_type to an int.  This was
     necessary with the new click-drag selection scheme.  I think it may
     clear up a few other mysterious bugs as well.
   08.09.98 gL:
     support for colored and/or shaded atoms
   26.09.98 gL:
     support for use of symmetry in MOs
   29.01.99 gL:
     updated form of lattice vectors
     added box stuff and num_along field to molec_type
***/

#ifndef _3D_OBJECTS_
#define _3D_OBJECTS_

#ifndef _MY_DEFINES_
#include "defines.h"
#endif

#ifndef _BASIC_OBJECTS_
#include "basic_objects.h"
#endif

/****
  used to store the characters of MOs w.r.t. symmetry operations

****/
typedef struct { int planes[3]; } MO_character_type;

/***********

  parametric surfaces

************/
typedef struct {
  int samples1, samples2;
  point_type *points;
  int *colors;
} param_surf_type;

/*******

  used to store information about MO gradients and values at
  a given point

********/
typedef struct {
  point_type loc;
  point_type grad;
  float val;
} MO_info_type;

/********

  used to construct a lookup table to speed function
   evaluations

********/
typedef struct {
  char filled;
  int num_entries;
  float min_val, max_val;
  float step;
  float *values;
} lookup_table_type;

/* these macros are used to pull things out of lookup tables */
#define READ_FROM_LOOKUP_TBL(tbl, x)                                           \
  ((x < tbl->max_val && x > tbl->min_val)                                      \
       ? tbl->values[(int)floor((x - tbl->min_val) / tbl->step + 0.5)]         \
       : 0.0)
#define INSERT_INTO_LOOKUP_TBL(tbl, x, val)                                    \
  ((x < tbl->max_val && x > tbl->min_val)                                      \
       ? tbl->values[(int)floor((x - tbl->min_val) / tbl->step + 0.5)] = val   \
       : 0)

/***

  this is used to construct a list of AO's contributing to
  a given center

***/
typedef struct AO_list_type_def {
  /* the type of this AO (used as a tab into the radial lookup table) */
  int type;

  /* the coefficient of this orbital in the MO */
  float coeff[MAX_MOS], coeffI[MAX_MOS];
  /*  float *coeff,*coeffI;*/
  float zeta1, zeta2, C1, C2;

#ifdef INCLUDE_ADF_PLOTS
  /* these are used for the ADF plotting */
  int kx, ky, kz, kr;
  float norm_fact;
#endif

  void(*ang_func) PROTO((point_type *, float, float, float *, point_type *));
  void(*rad_func) PROTO((point_type *, float, float, float, float, float, float,
                         float *, point_type *));

  lookup_table_type *rad_lookup_tbl;

  struct AO_list_type_def *next;
} AO_list_type;

/****

  this structure is the linked list of centers (atoms) contributing
  to an MO isosurface

*****/
typedef struct MO_center_list_type_def {
  /* a toggle used to indicate if this center is excluded from the calculation
   */
  char exclude;

  /* this center's type */
  int atom_type;
  int color;
  char *type;

  /* the location of this center */
  point_type *loc;

  /* the list of AO components */
  int num_AOs;
  /*  AO_list_type AO_list[MAX_AOS];*/
  AO_list_type *AO_list;

} MO_center_list_type;

/***********

  an atom

***********/
typedef struct {
  char exclude;
  char type[4];
  char color;
  /* stuff for custom atom displays */
  char custom, crosses_on, outlines_on, shading_on;

  int num;
  int *linesto;
  int num_lines_out;
  char lineto;
  point_type loc;
#ifdef INCLUDE_ADF_PLOTS
  point_type *displacements;
#endif
  float rad;
  param_surf_type *p_surf;

  /* used for click-selecting atoms */
  int is_selected;
  point_type2D screen_loc;
  int screen_rad, draw_order;

  /* used for coordination polyhedra */
  int num_triangles_in;
  int *triangles_in;

  /* information for bond length->bond_valence plots */
  float ri, ci;

  /* stuff for shading/color of atoms */
  float atom_shade;
  float atom_color[3];
  long int Cpixel_val[NUM_X_SHADES], Gpixel_val[NUM_X_SHADES];
} atom_type;

/***********

  a molecule

***********/
typedef struct {
  /* some booleans that will be used */
  char filename[240];
  char numbers_on, draw_connectors, outlines_on, shading_on, hydrogens_on;
  char fancy_lines, symbols_on, axes_on, dummies_on;
  char breaking_lines, tubes_on, crosses_on, trapezoidal_bonds;
  char draw_lattice;
  char draw_box;
  char draw_polyhed;
  int num_atoms;
  int num_frames, current_frame;
  atom_type *atoms;

  /* the "bonds" drawn with the molecule */
  float bond_tol;
  float old_bond_tol;
  float rad_mult;
  int *num_lines;
  line_type *lines;
  int line_width;
  float bond_rad;

  /* this is stuff for dealing with solids */
  int num_dim, num_atoms_in_cell;
  point_type lattice_vect[4], orig_lattice[4];
  point_type cell_box[8];
  int num_along[3];

  /* for doing coordination polyhedra */
  int num_triangles, num_polyhed_verts;
  vertex_type *polyhed_verts;
  triangle_type *triangles;
  int *polyhed_centers;
  float *polyhed_rads;

#ifdef USE_LASSP_ROTATE
  int rotate_tool_on;
  matrix_type rot_matrix;
#endif

#ifdef INCLUDE_ADF_PLOTS
  int num_vibrations, active_vibn;
  float vibration_scale;
#endif

#ifdef INCLUDE_BOND_VALENCE
  char valence_for_bonds;
  float bond_tol2, old_bond_tol2;
#endif
} molec_type;

#ifndef _CONTOUR_
#include "contour.h"
#endif

typedef struct MO_contours_def {
  int num_pts;
  double value;

  /********

    these are used in the hidden line removal

  *********/
  /* the list of hidden points */
  char *hidden_points;
  /* 1/slopes and intercepts */
  float *inv_slope, *intercept;
  /* the orientation */
  char orientation;
  /* the coefficients of the equation of the plane */
  float A_C, B_C, D_C;
  /* the max and min values */
  point_type min_vals, max_vals;

  /* the actual data */
  point_type *coords;

  struct MO_contours_def *next;
} MO_contours_type;

typedef struct {

  /* contour details */
  int num_levels, num_approx_pts, interp_kind, order, levels_kind;

  int num_conts;
  int num_a, num_b;
  point_type left_top_corner, right_bottom_corner;
  int orientation;
  iso_curve_type *data;
  float begin_conts;

  /* used for arbitrary planes */
  point_type Bas1, Bas2, Origin;

  double *levels_list;
  MO_contours_type *contours;
  gnuplot_contour_type *gnu_contours;
} MO_contour_plot_type;

/***********

  an MO isosurface

************/
typedef struct {
  char filename[240];

  int orig_num_centers;
  int num_centers, num_centers_in_cell;

  /* keep track of which MO we're looking at */
  int num_MOs, active_MO;
  int *MO_numbers;

  char display_molec, display_surf;
  char do_shading, do_lines;

  molec_type *molec;

  /* stuff for the radial lookup tables */
  int num_lookup_entries;
  float lookup_min, lookup_max;

  /* the unique center list (for radial lookup tables) */
  int num_unique;
  MO_center_list_type *unique_centers;

  /* the list of centers contributing to the surface */
  MO_center_list_type *MO_centers;
  MO_center_list_type *raw_MO_centers;

  /* the k points (used for extended systems) */
  point_type *kpoints;

  /* the contraction coefficient for exponents */
  float contraction;

  float surface_value;
  float surface_tolerance;

  /* stuff for doing gridded surfaces */
  float slop, search_radius;
  point_type bmin, bmax;
  int num_steps;
  float voxel_size;

  /* for drawing triangulated surfaces */
  int num_triangles;
  triangle_type *triangles;
  explicit_triangle_type *explicit_triangles;
  vertex_type *triangle_vertices;
  int num_vertices;

  /* for doing contours */
  char display_conts;
  char sort_conts;
  char invert_phase;
  int levels_kind;
  char plane_dirs[80];
  char hidden_surf;
  MO_contour_plot_type *MO_contours;

  char adf_plot;

  MO_character_type *characters;

} MO_surface_type;

typedef struct {
  point_type *loc;
  point_type *next;
  float val;
} contour_point_type;

/**************

  This is stuff for generic 3D objects, it's used to allow depth sorting
   of all 3D objects (atoms, surface polygons, etc.) and to generally
   clean up the 3D code.

***************/

#define ATOM_3D 0
#define CONT_POINT_3D 1
#define TRIANGLE_3D 2
#define LINE_3D 3

typedef struct { point_type *end1, *end2; } line_3D_type;

typedef union {
  contour_point_type contour_point;
  atom_type *atom;
  triangle_type *triangle;
  line_3D_type line;
} storage_3D;

typedef struct {
  char type;
  storage_3D object;
} generic_3D_object;

typedef struct { point_type points[4]; } axis_type;

#endif
