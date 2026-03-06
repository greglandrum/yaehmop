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

  Created by greg Landrum October 1995

********/

/***
  Edit History:

  11.04.98 gL:
   prototypes for fatal_bug and nonfatal_bug modified

  August '98: WG
    - intercell_COOP_check added

***/

/* this is stolen from steve gifford */
#ifndef PROTO
#if defined(_NO_PROTO) || defined(_alpha) || defined(MIPSEL)
#define PROTO(x) ()
#else /* defined(_NO_PROTO) || defined(__alpha) || defined(MIPSEL) */
#define PROTO(x) x
#endif /* defined(_NO_PROTO) || defined(__alpha) || defined(MIPSEL) */
#endif /* PROTO */

extern real eval_COOP PROTO((COOP_type *, detail_type *, cell_type *, int,
                             avg_prop_info_type *, hermetian_matrix_type,
                             K_orb_ptr_type *, int *));
extern void gen_COOP PROTO((detail_type *, cell_type *, int,
                            avg_prop_info_type *, hermetian_matrix_type,
                            K_orb_ptr_type *, int *));
extern void gen_avg_COOPs PROTO((detail_type *, cell_type *, int,
                                 avg_prop_info_type *, hermetian_matrix_type,
                                 K_orb_ptr_type *, int *));
extern void intercell_COOP_check PROTO((COOP_type *));
extern void gen_total_DOS PROTO((detail_type *, cell_type *, int,
                                 avg_prop_info_type *, K_orb_ptr_type *));
extern real orb_contribution PROTO((int, int, int, avg_prop_info_type *, int));
extern real FMO_contribution PROTO((int, int, int, avg_prop_info_type *, int));
extern real atom_contribution PROTO((int, int, int, int, avg_prop_info_type *,
                                     int, int *));
extern void gen_projected_DOS PROTO((detail_type *, cell_type *, int,
                                     avg_prop_info_type *, K_orb_ptr_type *,
                                     int *));
extern void build_k_hamil_FAT PROTO((cell_type *, hermetian_matrix_type,
                                     hermetian_matrix_type,
                                     hermetian_matrix_type, int));
extern void build_k_hamil_THIN PROTO((cell_type *, hermetian_matrix_type,
                                      hermetian_matrix_type,
                                      hermetian_matrix_type, int));
extern void build_k_overlap_FAT PROTO((cell_type *, k_point_type *,
                                       hermetian_matrix_type,
                                       hermetian_matrix_type, int));
extern void build_k_overlap_THIN PROTO((cell_type *, detail_type *,
                                        k_point_type *, hermetian_matrix_type,
                                        hermetian_matrix_type, int));
extern void build_all_K_overlaps PROTO((cell_type *, detail_type *,
                                        hermetian_matrix_type,
                                        hermetian_matrix_type, int, int,
                                        int *));

extern void R_space_Hamiltonian PROTO((cell_type *, detail_type *,
                                       hermetian_matrix_type,
                                       hermetian_matrix_type, int, int *));
extern void full_R_space_Hamiltonian PROTO((cell_type *, detail_type *,
                                            hermetian_matrix_type,
                                            hermetian_matrix_type, int, int *,
                                            char));
extern void zero_overlaps PROTO((real *, detail_type *, int, int, char, int *));
extern void calc_R_overlap PROTO((real *, cell_type *, detail_type *, int,
                                  point_type, char, int *));
extern void R_space_overlap_matrix PROTO((cell_type *, detail_type *,
                                          hermetian_matrix_type, int, int,
                                          int *, int));
extern int find_atom PROTO((atom_type *, int, int));
extern void eval_Zmat_locs PROTO((atom_type *, int, int, char));
extern void calc_avg_occups PROTO((detail_type *, cell_type *, int,
                                   K_orb_ptr_type *, avg_prop_info_type *,
                                   prop_type *, real *));
extern void print_avg_occups PROTO((detail_type *, cell_type *, int,
                                    K_orb_ptr_type *, avg_prop_info_type *,
                                    prop_type, real *));
extern void calc_avg_FMO_occups PROTO((detail_type *, int, K_orb_ptr_type *,
                                       avg_prop_info_type *, real *));
extern void calc_avg_OP PROTO((detail_type *, cell_type *, int,
                               K_orb_ptr_type *, avg_prop_info_type *,
                               hermetian_matrix_type, prop_type));
extern void find_crystal_occupations PROTO((detail_type *, real, int,
                                            K_orb_ptr_type *, real *));
extern void store_avg_prop_info PROTO((detail_type *, int, eigenset_type,
                                       hermetian_matrix_type, int, real *,
                                       avg_prop_info_type *));
extern int sort_energies_helper PROTO((const void *, const void *));
extern void sort_avg_prop_info PROTO((detail_type *, int, avg_prop_info_type *,
                                      K_orb_ptr_type *));
extern void gen_symm_lines PROTO((band_info_type *));
extern void construct_band_structure PROTO(
    (cell_type *, detail_type *, hermetian_matrix_type, hermetian_matrix_type,
     hermetian_matrix_type, hermetian_matrix_type, complex *, complex *,
     eigenset_type, real *, real *, real *, complex *, int, int *));
extern void eval_charge_matrix PROTO((cell_type *, eigenset_type,
                                      hermetian_matrix_type, int, int *, real *,
                                      real *));
extern void reduced_charge_matrix PROTO((int, int, int *, real *, real *));
extern void check_a_cell PROTO((atom_type *, point_type, int, real, char *));
extern void check_nn_contacts PROTO((cell_type *, detail_type *details));
extern void build_distance_matrix PROTO((cell_type *, detail_type *details));
extern void dump_distance_mats PROTO((cell_type *, detail_type *details));
extern void fill_distance_matrix PROTO((cell_type *, int, float *,
                                        point_type *));
extern void display_lattice_parms PROTO((cell_type *));
extern int factorial PROTO((int));
extern void free_atomic_energy PROTO((cell_type *, real *, real *));
extern void AO_occupations PROTO((cell_type *, int, real *, int *, real *));
extern void eval_electrostatics PROTO((cell_type *, int, eigenset_type, real *,
                                       real *, int *, real *, real *, real *,
                                       real *, real *));
extern void read_geom_frag PROTO((FILE *, geom_frag_type *));
extern void write_atom_parms PROTO((detail_type *, atom_type *, int, char));
extern void write_atom_coords PROTO((atom_type *, int, char, char));
extern void fill_atomic_parms PROTO((atom_type *, int, FILE *, char *));
extern void parse_printing_options PROTO((FILE *, detail_type *, cell_type *));
extern void read_inputfile PROTO((cell_type *, detail_type *, char *, int *,
                                  int **, FILE *, char *));
extern void fatal PROTO((char *));
extern void fatal_bug PROTO((char *, char *, int));
extern void nonfatal_bug PROTO((char *, char *, int));
extern void error PROTO((char *));
extern void upcase PROTO((char *));
extern char *safe_strcpy PROTO((char *, char *));
extern void center_text_string PROTO((char *src, char *dest, int dest_len));
extern void left_just_text_string PROTO((char *src, char *dest, int dest_len));
extern void map_orb_num_to_name PROTO((char *, int, int *, int, atom_type *,
                                       int));
extern void debugmat PROTO((real *, int, int, real));
extern void printmat PROTO((real *, int, int, FILE *, real, char, int));
extern void print_labelled_mat PROTO((real *, int, int, FILE *, real,
                                      atom_type *, int, int *, int, char, char,
                                      int));
extern void print_sym_mat PROTO((real *, int, int, FILE *, char *, char *,
                                 int));
extern int skipcomments PROTO((FILE *, char *, char));
extern void find_atoms_orbs PROTO((int, int, int, int *, int *, int *));
extern int overlap_tab_from_vect PROTO((point_type *, cell_type *));
extern void check_for_errors PROTO((cell_type *, detail_type *, int));
extern void build_orbital_lookup_table PROTO((cell_type *, int *, int **));
extern void print_MOs PROTO((detail_type *, int, eigenset_type, int,
                             atom_type *, int, int, int *));
extern double d_sign PROTO((double, double));
extern void parse_integer_string PROTO((char *, int **, int *));
extern void dump_hermetian_mat PROTO((int, real *, int));
extern void dump_sparse_mat PROTO((FILE *, real *, int, real));
extern void loop_over_k_points
    PROTO((cell_type *, detail_type *, hermetian_matrix_type,
           hermetian_matrix_type, hermetian_matrix_type, hermetian_matrix_type,
           complex *, complex *, eigenset_type, real *, real *, real *,
           complex *, prop_type *, avg_prop_info_type *, int, int *));

extern void sparsify_hermetian_matrix PROTO((real, hermetian_matrix_type, int));
extern void sparsify_matrix PROTO((real, real *, real *, int));
extern void allocate_matrices
    PROTO((cell_type *, detail_type *, hermetian_matrix_type *,
           hermetian_matrix_type *, hermetian_matrix_type *,
           hermetian_matrix_type *, complex **, complex **, eigenset_type *,
           real **, real **, real **, complex **, prop_type *,
           avg_prop_info_type **, int, int *, int *, K_orb_ptr_type **));
extern void cleanup_memory PROTO(());
extern void mov PROTO((real *, real *, real *, real *, int, int, real, int, int,
                       int, int, atom_type *));
extern void calc_occupations PROTO((detail_type *, real, int, real *,
                                    eigenset_type));
extern void reduced_mulliken PROTO((int, int, int *, real *, real *));
extern void FMO_reduced_mulliken PROTO((detail_type *, int, int, real *,
                                        real *));
extern void eval_mulliken PROTO((cell_type *, eigenset_type,
                                 hermetian_matrix_type, int, real *, int *,
                                 real *, real *, real *));
extern void modified_mulliken PROTO((cell_type *, eigenset_type,
                                     hermetian_matrix_type, int, real *, int *,
                                     real *, real *, real *, real *));
extern void read_NEW3file PROTO((cell_type *, detail_type *, FILE *, char *));
extern void find_princ_axes PROTO((atom_type *, point_type *, real[3][3],
                                   real[3], int));

extern void vector_diff PROTO((point_type *, point_type *, point_type *));
extern void normalize_vector PROTO((point_type *, point_type *));
extern real dot_prod PROTO((point_type *, point_type *));
extern void scale_vector PROTO((point_type *, real));
extern void cross_prod PROTO((point_type *, point_type *, point_type *));
extern void mult_matrices PROTO((real *, real *, real *, int));
extern void auto_walsh PROTO((real *, int, real, real));
extern void walsh_update PROTO((cell_type *, detail_type *, int, char));
extern void walsh_output PROTO((detail_type *, cell_type *, int, eigenset_type,
                                hermetian_matrix_type, hermetian_matrix_type,
                                prop_type, int *, int));
extern void eval_xtal_coord_locs PROTO((cell_type *, char));
extern void update_zetas PROTO((cell_type *, real *, real, int *, char));
extern void init_FMO_file PROTO((detail_type *, int, real));
extern void build_FMO_overlap PROTO((detail_type *, int, int,
                                     hermetian_matrix_type, int *));
extern void build_FMO_hamil PROTO((detail_type *, int, int,
                                   hermetian_matrix_type, int *));
extern void diagonalize_FMO PROTO((detail_type *, real *, real *, real *,
                                   complex *, complex *, complex *));
extern void gen_FMO_tform_matrices PROTO((detail_type *));
extern void tform_wavefuncs_to_FMO_basis PROTO((detail_type *, int, int,
                                                eigenset_type, int *));
extern void tform_matrix_to_FMO_basis PROTO((detail_type *, int, int, real *,
                                             real *, real *, real *,
                                             complex_matrix_type, int *));
extern void tform_hermetian_matrix_to_FMO_basis
    PROTO((detail_type *, int, int, hermetian_matrix_type, real *, real *,
           hermetian_matrix_type, int *));

extern void charge_to_num_electrons PROTO((cell_type *));
extern void update_chg_it_parms PROTO((detail_type *, cell_type *, real *,
                                       int *, int, int *));
extern void fill_chg_it_parms PROTO((atom_type *, int, int, FILE *));
extern void parse_charge_iteration PROTO((FILE *, detail_type *, cell_type *));
extern void update_muller_it_parms PROTO((detail_type *, cell_type *, real *,
                                          int *, int, int *));
extern void calc_muller_init_parms PROTO((atom_type *));
extern void calc_muller_parms PROTO((atom_type * atom, real s_occup,
                                     real p_occup, real d_occup, real *s_E,
                                     real *s_Z, real *p_E, real *p_Z, real *d_E,
                                     real *d_Z));

extern void parse_muller_parms PROTO((FILE *, detail_type *, cell_type *));
extern void parse_equiv_atoms PROTO((FILE *, detail_type *, cell_type *));

extern void handle_sigint PROTO(());

/*********
  the fortran stuff
**********/
extern int abfns PROTO((real *, real *, real *, real *, real *, int *, int *,
                        int *, int *, int *, int *));
extern void lovlap PROTO((real *, real *, real *, real *, real *, real *, int *,
                          int *, int *, int *, int *, int *));
extern void cboris PROTO((int *, int *, real *, real *, real *, real *, real *,
                          real *, int *));

#ifndef SYM_OPS_DEFINED
#include "symmetry.h"
#endif
extern void name_sym_element PROTO((sym_op_type *, FILE *, int, int));
extern sym_op_type *make_new_sym_op PROTO((void));
extern void gen_sym_ops PROTO((sym_op_type **, int *));
extern void find_off_axis_sym_ops PROTO((detail_type *, cell_type *,
                                         point_type *, point_type *,
                                         point_type[3], sym_op_type *, int *));
extern void find_sym_ops PROTO((detail_type *, cell_type *));
extern void find_walsh_sym_ops PROTO((cell_type *, detail_type *));
extern void find_MO_symmetries PROTO((int, detail_type *, cell_type *,
                                      eigenset_type, hermetian_matrix_type,
                                      int *));
extern void compare_molecules PROTO((atom_type *, point_type *, point_type *,
                                     int, int *, char *, real));
extern void transform_p_orbs PROTO((real *, real[T_MAT_DIM][T_MAT_DIM]));
extern void transform_d_orbs PROTO((real *, real[D_T_MAT_DIM][D_T_MAT_DIM]));
extern void transform_orbitals PROTO((atom_type *, real *, sym_op_type *));
extern void full_transform PROTO((atom_type *, point_type, real[3][3], int));
extern void transform_atoms PROTO((atom_type *, real[T_MAT_DIM][T_MAT_DIM],
                                   int));
extern void transform_one_point PROTO((point_type *,
                                       real[T_MAT_DIM][T_MAT_DIM]));
extern void translate_atoms PROTO((point_type *, point_type, int));
extern void transform_atomic_locs PROTO((point_type *,
                                         real[T_MAT_DIM][T_MAT_DIM], int));
extern void transform_3x3_transpose PROTO((point_type *, real[3][3], int));

extern void postprocess_FMO PROTO((cell_type *, detail_type *,
                                   hermetian_matrix_type, hermetian_matrix_type,
                                   hermetian_matrix_type, hermetian_matrix_type,
                                   complex *, complex *, eigenset_type, real *,
                                   real *, real *, complex *, prop_type *,
                                   avg_prop_info_type *, int, int *));

extern void postprocess_FCO PROTO((cell_type *, detail_type *,
                                   hermetian_matrix_type, hermetian_matrix_type,
                                   hermetian_matrix_type, hermetian_matrix_type,
                                   complex *, complex *, eigenset_type, real *,
                                   real *, real *, complex *, prop_type *,
                                   avg_prop_info_type *, int, int *));

extern void postprocess_results
    PROTO((cell_type *, detail_type *, hermetian_matrix_type,
           hermetian_matrix_type, hermetian_matrix_type, hermetian_matrix_type,
           complex *, complex *, eigenset_type, real *, real *, real *,
           complex *, prop_type *, avg_prop_info_type *, int, int *));
extern void process_geom_frags PROTO((cell_type *));
extern void insert_geom_frags PROTO((cell_type *));

extern void compare_crystal_basis PROTO((cell_type *, point_type *,
                                         point_type *, point_type *, int, int *,
                                         char *, real));
extern void compare_crystal_lattice PROTO((cell_type *, point_type *,
                                           point_type *, char *, real, real *));
extern void reduce_kpoints PROTO((detail_type * details, cell_type *cell,
                                  point_type *raw_points, int num_raw_points,
                                  real *multiplicity, int *num_points,
                                  int orthogonal_axes));

extern int check_for_orthogonal_basis PROTO((point_type vects[3], int dim,
                                             real tol));
extern int atoms_are_equiv PROTO((cell_type * cell, point_type *loc1,
                                  point_type *loc2, int *which_cell,
                                  real symm_tol, point_type *));

extern void calc_reciprocal_lattice PROTO((cell_type * cell));
extern void gen_k_point_mesh PROTO((point_type * *points, int num_per_vect[3],
                                    int *num_generated,
                                    char include_high_symm_p, int dim,
                                    real offset));
extern void automagic_k_points PROTO((detail_type * details, cell_type *cell));

extern void set_details_defaults PROTO((detail_type *));
extern void set_cell_defaults PROTO((cell_type *));
extern void run_bind PROTO((char *, bool, char *));
extern void run_eht PROTO((FILE *));

extern int *my_malloc PROTO((long));
extern int *my_calloc PROTO((int, int));
extern int *my_realloc PROTO((int *, int));

#ifdef INCLUDE_NETCDF_SUPPORT
extern void netCDF_handle_error PROTO((int));
extern void netCDF_init_file PROTO((detail_type * details, cell_type *cell,
                                    int num_orbs, int));
extern void netCDF_close_file PROTO((detail_type * details));
extern void netCDF_write_Es PROTO((detail_type * details, int num_orbs,
                                   avg_prop_info_type *prop_info));

#endif

#ifdef USE_LAPACK
#ifndef F2C_INCLUDE
typedef long int integer;
typedef double doublereal;
#endif
extern int zhegv_ PROTO((integer * itype, char *jobz, char *uplo, integer *n,
                         doublecomplex *a, integer *lda, doublecomplex *b,
                         integer *ldb, doublereal *w, doublecomplex *work,
                         integer *lwork, doublereal *rwork, integer *info));
extern int zheev_ PROTO((char *jobz, char *uplo, integer *n, doublecomplex *a,
                         integer *lda, doublereal *w, doublecomplex *work,
                         integer *lwork, doublereal *rwork, integer *info));
#endif
