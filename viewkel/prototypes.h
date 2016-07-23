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
/********

  the prototypes for all functions

  Created by greg Landrum November 1995

********/

/* this is stolen from steve gifford */
#ifndef PROTO
#if defined(_NO_PROTO) || defined(_alpha) || defined(MIPSEL)
#define PROTO(x) ()
#else /* defined(_NO_PROTO) || defined(__alpha) || defined(MIPSEL) */
#define PROTO(x) x
#endif /* defined(_NO_PROTO) || defined(__alpha) || defined(MIPSEL) */
#endif /* PROTO */

/* main.c */
extern void main PROTO((int, char **));

/* 3D_objects.c */
extern int find_numbered_atoms_in_objects
PROTO((generic_3D_object *, int, int));
extern void draw_axes PROTO((axis_type *));
extern void draw_atom
PROTO((atom_type * atom, int atom_num, generic_3D_object *objects,
       int num_objects, molec_type *molec));
extern void draw_triangle
PROTO((triangle_type *, vertex_type *, MO_surface_type *));
extern int object_zcompare PROTO((const void *, const void *));
extern void draw_3D_objects PROTO((prim_type *, object_type *));
extern line_type *find_the_line PROTO((int num1, int num2, molec_type *molec));

/* FMO_diags.c */
extern void find_FMO_tic_sep PROTO((FMO_diagram_type *));
extern void fill_levels_with_electrons
PROTO((FMO_level_type *, float, float *));
extern FMO_level_type *map_orbital_to_level PROTO((FMO_level_type *, int));
extern void preprocess_FMO_data PROTO((FMO_diagram_type *));
extern void change_FMO_degen_tol PROTO((int, char **));
extern void read_FMO_data PROTO((FILE *, FMO_diagram_type *));
extern void new_FMO_diagram PROTO((char *));
extern void FMO_draw_electrons PROTO((FMO_diagram_type *, FMO_level_type *));
extern void draw_FMO_diagram PROTO((prim_type *, object_type *));
extern line_type *find_the_line PROTO((int num1, int num2, molec_type *molec));

/* band_graphs.c */
extern void preprocess_band_graph_data PROTO((band_graph_type * band_graph));
extern void read_band_data PROTO((FILE *, band_graph_type *));
extern void new_band_graph PROTO((char *));
extern void draw_band_graph PROTO((prim_type *, object_type *));

/* buttons.c */
#ifdef X_GRAPHICS
extern void g_floatbutton PROTO((Window, button_type *, int, int));
extern void g_intbutton PROTO((Window, button_type *, int, int));
extern void g_stringbutton PROTO((Window, button_type *, int, int));
extern void g_multistringbutton PROTO((Window, button_type *, int, int));
extern void g_functionbutton PROTO((Window, button_type *, int, int));
extern void g_curvebutton PROTO((Window, button_type *, int, int));
extern void g_togglebutton PROTO((Window, button_type *, int, int));
extern void g_modebutton PROTO((Window, button_type *, int, int));
extern void g_drawbuttons PROTO((button_win_type *));
extern void g_draw_all_buttons PROTO((button_win_type *));
#endif
extern button_type *newbutton PROTO((void));
extern void free_but_win PROTO((button_win_type *, button_win_type **));
extern void free_child_but_win PROTO((button_win_type *));
extern void handlebutton PROTO((button_win_type *, int, int, int, int));
extern void find_button_win
PROTO((button_win_type *, Window, int, int, int, int));
extern button_win_type *new_button_win PROTO((button_win_type *, char *));
extern void set_button_pixmaps PROTO((button_win_type *));
extern void build_main_button_window PROTO((button_win_type **));
extern void build_PS_options_button_window PROTO((button_win_type **));
extern void build_label_button_window PROTO((button_win_type **, label_type *));
extern void build_molec_button_window PROTO((button_win_type **, molec_type *));
extern void build_graph_button_window PROTO((button_win_type **, graph_type *));
extern void build_walsh_button_window
PROTO((button_win_type **, walsh_graph_type *));
extern void build_band_button_window
PROTO((button_win_type **, band_graph_type *));
extern void build_cont_plot_button_window
PROTO((button_win_type **, contour_plot_type *));
extern void build_prop_graph_button_window
PROTO((button_win_type **, prop_graph_type *));
extern void build_FMO_button_window
PROTO((button_win_type **, FMO_diagram_type *));
extern void build_MO_surf_button_window
PROTO((button_win_type **, MO_surface_type *));

/* enhpost.c */
extern char *ENHPS_RememberFont PROTO((char *));
extern void ENHPS_init PROTO((void));
extern char *ENHPS_recurse PROTO((char *, char, char *, float, float, char));
extern int ENHPS_put_text PROTO((unsigned int, unsigned int, char *, char));
extern void open_ps_button_window PROTO((void));
extern void do_ps_output PROTO((void));

/* fileio.c */
extern void general_read PROTO(());
extern void save_MO_surf PROTO((FILE *, prim_type *, object_type *));
extern void save_FMO PROTO((FILE *, prim_type *, object_type *));
extern void save_molecule PROTO((FILE *, prim_type *, object_type *));
extern void save_graph PROTO((FILE *, prim_type *, object_type *));
extern void save_prop_graph PROTO((FILE *, prim_type *, object_type *));
extern void save_band_graph PROTO((FILE *, prim_type *, object_type *));
extern void save_prim PROTO((FILE *, prim_type *, object_type *));
extern void save_all PROTO((void));
extern int read_MO_surf_object PROTO((FILE *, char *));
extern int read_FMO PROTO((FILE *, char *));
extern int read_molecule_object PROTO((FILE *, char *));
extern int read_graph_object PROTO((FILE *, char *));
extern int read_band_object PROTO((FILE *, char *));
extern int read_prop_graph_object PROTO((FILE *, char *));
extern int read_object PROTO((FILE *));
extern void read_all PROTO((char *));
extern void dump_3D_objects PROTO((prim_type *, object_type *));
extern void dump_VRML PROTO((prim_type *, object_type *));

/* fit_orbs.c */
extern int choose_MO PROTO((MO_surface_type *));
extern void exclude_atoms PROTO((int, char **));
extern void include_atoms PROTO((int, char **));
extern void change_active_MO PROTO((int, char **));
extern int read_MO_surface PROTO((FILE *, MO_surface_type *));
extern void new_MO_surface PROTO((char *filename));
#ifdef INCLUDE_ADF_PLOTS
extern void new_ADF_MO_surface PROTO((char *filename));
extern int read_ADF_MO_surface PROTO((FILE *, MO_surface_type *));
#endif

/* genutil.c */
extern void fatal_bug PROTO((char *, char *, int));
extern void fatal PROTO((char *));
extern void error PROTO((char *));
extern void printmat PROTO((matrix_type *));
extern void display PROTO((char *));
extern void readfloatparm PROTO((char *, float *));
extern void readintparm PROTO((char *, int *));
extern void readcharparm PROTO((char *, char *));
extern void readstringparm PROTO((char *, char **));
extern void readmultistringparm PROTO((char *, int, char **));
extern char *safe_strcpy PROTO((char *, char *));
extern int skipcomments PROTO((FILE *, char *));
extern void upcase PROTO((char *));
extern void dump_molecular_coords PROTO((molec_type *));
extern float my_drand PROTO((float));
extern void parse_integer_string PROTO((char *, int **, int *));
extern void hide_selected_atoms PROTO((int, object_type *));
extern void show_all_atoms PROTO((object_type *));
extern void show_selected_data PROTO((int, object_type *, int, int));
extern void unselect_all_atoms(int, object_type *);
extern void invert_selected_atoms(int *, object_type *);
#ifdef INCLUDE_ADF_PLOTS
extern void map_AO_number_to_center
PROTO((int AO_num, MO_surface_type *MO_surf, MO_center_list_type **, int *));
#endif
extern void angles_from_bond_vect
PROTO((point_type * p1, point_type *p2, float *theta_y, float *theta_z,
       float *len));
;

/* graphics.c */
extern void g_line PROTO((float, float, float, float));
extern void g_lines PROTO((XPoint *, point_type2D *, int, char));
extern void g_center_text PROTO((float, float, char *));
extern void g_label PROTO((float, float, label_type *));
extern void g_bot_center_label PROTO((float, float, char *));
extern void g_cent_center_label PROTO((float, float, char *));
extern void g_right_text PROTO((float, float, char *));
extern void g_left_text PROTO((float, float, char *));
extern void g_change_linestyle PROTO((int));
extern void g_change_linewidth PROTO((float));
extern void g_draw_breaking_line PROTO((float, float, float, float, int));
extern void g_draw_stop_line
PROTO((float, float, float, float, int, float, float));
extern void g_draw_tube_line PROTO((float, float, float, float, int));
extern void g_draw_stop_tube
PROTO((float, float, float, float, int, float, float));
extern void g_draw_dashed_breaking_line
PROTO((float, float, float, float, int, int));
extern void g_draw_dashed_tube_line
PROTO((float, float, float, float, int, int));
extern void g_draw_dashed_stop_line
PROTO((float, float, float, float, int, int, float, float));
extern void g_draw_dashed_stop_tube
PROTO((float, float, float, float, int, int, float, float));
extern void g_rectangle PROTO((float, float, float, float));
extern void g_change_color PROTO((int));
extern void g_open_circle PROTO((float, float, float));
extern void g_crossed_circle PROTO((float, float, float));
extern void g_filled_circle
PROTO((float, float, float, float, float *, long int *, long int *));
extern void g_white_circle PROTO((float, float, float));
extern void g_filled_polygon PROTO((XPoint *, int));
extern void g_shaded_polygon PROTO((XPoint *, int, point_type *));
extern void g_open_polygon PROTO((XPoint *, int));
extern void g_xlegend PROTO((float, float, char *));
extern void g_ylegend PROTO((float, float, char *));
extern void g_title
PROTO((float x, float y, char title[NUM_TITLE_LINES][NORMAL_STR_LEN]));
extern void g_initgraphics PROTO((int, int));
extern void g_clear_screen PROTO((void));
extern void g_switch_buffers PROTO((void));
extern void g_insert_comment PROTO((char *));

/* graphs.c */
extern void find_tic_sep PROTO((graph_type *));
extern void count_curves PROTO((char *, int *));
extern void read_graph_data PROTO((FILE *, graph_type *, char, char));
extern void preprocess_graph_data PROTO((graph_type *));
extern void new_graph PROTO((char *));
extern void draw_graph PROTO((prim_type *, object_type *));

/* help.c */
extern void show_help PROTO((void));

/* hierarchy.c */
extern void makenewobject PROTO((void));
extern void addchild PROTO((object_type *, object_type *));
extern void free_obj PROTO((object_type *));
extern void clearall PROTO((void));
extern void clear_labels PROTO((void));
extern void select_object PROTO((int, int));
extern int select_atom PROTO((molec_type *, int, int));
extern int select_atoms_in_region PROTO((molec_type *, int, int, int, int));

/* implicit_polyg.c (just the externally visible routine) */
extern void construct_MO_isosurface PROTO((int, char **));

/* interface.c */
#ifdef X_GRAPHICS
extern void do_keypress PROTO((XEvent event));
extern void do_button PROTO((XEvent event));
extern void do_events PROTO((void));
#endif
#ifdef MAC_GRAPHICS
extern void do_keypress PROTO((EventRecord * event));
#endif

extern void parse_commands PROTO((void));
extern void show_keys PROTO((void));

/* manipulate.c */
extern void addprojectmat PROTO((void));
extern void drawit PROTO((prim_type * prim, object_type *obj));
extern void instantiate PROTO((object_type * object, char));
extern void redrawgraph PROTO((void));
extern void redraw PROTO((void));

/* matrix_ops.c */
extern double V3SquaredLength PROTO((Vector3 * a));
extern double V3Length PROTO((Vector3 * a));
extern Vector3 *V3Negate PROTO((Vector3 * v));
extern Vector3 *V3Normalize PROTO((Vector3 * v));
extern Vector3 *V3ConstantScale PROTO((Vector3 * v, double newlen));
extern Vector3 *V3Scale PROTO((Vector3 * v, double newlen));
extern Vector3 *V3Add PROTO((Vector3 * a, Vector3 *b, Vector3 *c));
extern Vector3 *V3Sub PROTO((Vector3 * a, Vector3 *b, Vector3 *c));
extern double V3Dot PROTO((Vector3 * a, Vector3 *b));
extern Vector3 *V3Lerp
PROTO((Vector3 * lo, Vector3 *hi, double alpha, Vector3 *result));
extern Vector3 *V3Combine
PROTO((Vector3 * a, Vector3 *b, Vector3 *result, double ascl, double bscl));
extern Vector3 *V3Mul PROTO((Vector3 * a, Vector3 *b, Vector3 *result));
extern double V3DistanceBetween2Points PROTO((Point3 * a, Point3 *b));
extern double V3AngleBetween3Points PROTO((Point3 * a, Point3 *b, Point3 *c));
extern double V3DihedralAngle
PROTO((Point3 * a, Point3 *b, Point3 *c, Point3 *d));
extern Vector3 *V3Cross PROTO((Vector3 * a, Vector3 *b, Vector3 *c));
extern void *V3MulPointByProjMatrix PROTO((Point3 * pin, Point3 *pout));

/* vibrations.c */
extern void new_vibration PROTO((char *filename));
extern void read_vibration_data PROTO((FILE * infile, molec_type *molec));

/* molecule.c */
extern void show_coord_env PROTO((int, object_type *));
extern void adjust_style PROTO((int, object_type *));
extern void adjust_color PROTO((int, object_type *, int));
extern int zcompare PROTO((atom_type * e1, atom_type *e2));
extern int find_numbered_atom
PROTO((atom_type * array, int length, int number));
extern void center_molecule PROTO((int num_args, char **molec_ptr));
extern void determine_mol_bounds
PROTO((molec_type * molec, point_type *bmin, point_type *bmax));
extern void hide_atoms PROTO((int num_args, char **molec_ptr));
extern void show_atoms PROTO((int num_args, char **molec_ptr));
extern void determine_connections PROTO((molec_type * molec));
extern void fill_atomic_info PROTO((int num_atoms, atom_type *atoms));
extern void read_molecule_data PROTO((FILE * infile, molec_type *molec));
extern void new_molecule PROTO((char *filename));

/* orbitals.c  externally visible functions only */
extern void calc_MO_value
PROTO((int which_MO, MO_info_type *MO_info, MO_center_list_type *centers,
       int num_centers, char ADF_plot));
extern void build_radial_lookup_table PROTO((MO_surface_type * MO_surf));

/* prop_graphs.c */
extern void new_prop_graph PROTO((char *filename));
extern object_type *traverseobj
PROTO((object_type * obj, int cor1, int cor2, Window win));
extern void bbsearch PROTO((int cor1, int cor2, Window win));
extern snap_type *readobj PROTO((snap_type * snapptr, object_type *obj));
extern void readsnap PROTO((int which, shead_type *shead));
extern snap_type *snapobj PROTO((object_type * obj));
extern void takesnap PROTO((int which));
extern void animateit PROTO((int from, int to));
extern void doanimation PROTO((void));
extern void getkeyframe PROTO((void));
extern void takekeyframe PROTO((void));
extern void grow_solid PROTO((int num_args, char *solid_p[MAX_ARGS]));
extern void grow_solid_with_surface
PROTO((int num_args, char *solid_p[MAX_ARGS]));

/* stack.c */
extern void loadmatrix PROTO((matrix_type * mat));
extern void pushmatrix PROTO((matrix_type * mat));
extern void popmatrix PROTO((void));
extern void multmatrix PROTO((matrix_type * mat));
extern void premultmatrix PROTO((matrix_type * mat));
extern void transform PROTO((point_type * point));
extern void transform_norm PROTO((point_type * point));
extern void getmatrix PROTO((matrix_type * mat));
extern void xrot PROTO((float theta));
extern void yrot PROTO((float theta));
extern void zrot PROTO((float theta));
extern void translate PROTO((float xtrans, float ytrans, float ztrans));
extern void scale PROTO((float xscale, float yscale, float zscale));
extern void shear PROTO((float sh1, float sh2, int which));
extern void ortho2 PROTO((float left, float right, float bottom, float top));

/* tek_lib.c */
extern void TEK_init PROTO((void));
extern void TEK_reset PROTO((void));
extern void TEK_graphics PROTO((void));
extern void TEK_text PROTO((void));
extern void TEK_move PROTO((unsigned int x, unsigned int y));
extern void TEK_vector PROTO((unsigned int x, unsigned int y));
extern void TEK_linetype PROTO((int linetype));
extern void TEK_put_text PROTO((unsigned int x, unsigned int y, char *str));
extern void TEK_center_text PROTO((unsigned int x, unsigned int y, char *str));
extern void TEK_right_text PROTO((unsigned int x, unsigned int y, char *str));

/* triangles.c */
extern void convert_explicit_triangles
PROTO((explicit_triangle_type * explicit, int num_explicit,
       triangle_type **triangles, int *num_tri, point_type **points,
       int *num_points));
extern void remove_degen_triangles PROTO((MO_surface_type * MO_surf));
extern void calc_triangle_centers PROTO((MO_surface_type * MO_surf));
extern void save_triangle_locs PROTO((int num_args, char *MO_surf_p[MAX_ARGS]));
extern void read_triangle_locs PROTO((int num_args, char *MO_surf_p[MAX_ARGS]));

/* walsh_graph.c */
extern void read_walsh_data
PROTO((FILE * infile, graph_type *graph, graph_type *tot_E));
extern void new_walsh_graph PROTO((char *));
extern void draw_walsh_graph PROTO((prim_type * prim, object_type *obj));

/* xstuff.c */
extern void X_initgraphics PROTO((int xsize, int ysize));

/* cont_plots.c */
extern void cont_find_tic_sep PROTO((contour_plot_type * cont_plot));
extern void read_3D_data PROTO((FILE * infile, contour_plot_type *cont_plot));
extern void preprocess_cont_plot_data PROTO((contour_plot_type * cont_plot));
extern void new_cont_plot PROTO((char *filename));
extern void draw_cont_plot PROTO((prim_type * prim, object_type *obj));
extern void contour_data_set PROTO((int num_args, char **cont_plot_ptr));

/* contour.c */
extern struct gnuplot_contours *contour_data
PROTO((int, int, iso_curve_type *, int, int, int, int, int, double *));

/* labels.c */
extern void new_label PROTO((char *));
extern void draw_label PROTO((prim_type *, object_type *));
extern void adjust_label PROTO((label_type *));

/* MO_conts.c */
extern void eval_MO_plane PROTO((MO_surface_type * surf));
extern void construct_MO_contours PROTO((int num_args, char **MO_surf_ptr));
extern void MO_contour_surf PROTO((int num_args, char **MO_surf_ptr));
extern void draw_contours
PROTO((MO_surface_type * surf, MO_contours_type *contours));
extern void find_contour_plane PROTO((MO_contours_type * contour));

#ifdef INCLUDE_BOND_VALENCE
/* valence.c */
extern void bond_length_to_bond_valence
PROTO((atom_type * atom1, atom_type *atom2, float length, float *RO_val,
       float *valence));
extern float calc_R0 PROTO((float c1, float r1, float c2, float r2));
#endif

/* polyhed.c */
extern void gen_coord_polyhed PROTO((molec_type *));
extern void find_coord_polyhed PROTO((molec_type *, int, float, int));

/* chull.c */
extern void gen_chull PROTO((atom_type *, int, triangle_list_type **));

#ifdef USING_THE_MAC
extern FILE *choose_mac_file(char *, char);
#endif

#ifdef USE_READLINE
extern char *readline PROTO((char *));
extern void add_history PROTO((char *));
#endif

#ifdef MEM_DEBUG
extern char *d_malloc PROTO((int size, char *filename, int linenum));
extern char *d_calloc PROTO((int num, int size, char *filename, int linenum));
extern void d_free PROTO((char *address, char *filename, int linenum));
extern char *d_realloc
PROTO((char *address, int size, char *filename, int linenum));
extern void d_check_core PROTO(());
#endif
