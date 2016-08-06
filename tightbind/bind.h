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

/****************************************************************************
*
*     this is the include file for the tightbinding program
*
*   contained herein are:
*
*       #include's and #define's that everyone needs
*       type definitions
*       extern definitions of global variables
*
*  created:  greg landrum  August 1993
*
*****************************************************************************/

/***
  Edit History:

  March 1998, Wingfield Glassey

  - modifications for f orbitals, abfns.c and COHP's

  11.04.98 gL:
   #defines for FATAL_BUG and NONFATAL_BUG added
   version_string changed to 2.5a, for the hell of it... it probably
   should be something like 3.0a, but I can't be bothered.
  07.05.98 gL:
   added netCDF support
   changed atom symbol length from hardcoded 10 to ATOM_SYMB_LEN
  18.08.98 gL:
   check for definition of PI before redefining it please.
  26.09.98 gL:
   added some symmetry storage to detail_type
  12.04.99 gL:
   version string updated to 3.0
***/

/******  includes for everyone ******/
#include <stdio.h>
#include <fcntl.h>
#include <math.h>
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#ifndef USING_THE_MAC
#include <sys/types.h>
#include <sys/stat.h>
#else
#include <SIOUX.h>
#endif

#ifdef INCLUDE_NETCDF_SUPPORT
#include <netcdf.h>
#endif

#ifdef ANAL_ABOUT_PROTOTYPES
/**********

  this is for development stuff, to make sure that everything is being called
  properly.  In general the user doesn't have to worry about this.

**********/
#ifndef USING_THE_MAC
#include <malloc.h>
#else
#include <unix.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#endif

/******  defines  ********/

#define FATAL_BUG(__a__) fatal_bug((__a__), __FILE__, __LINE__)
#define NONFATAL_BUG(__a__) nonfatal_bug((__a__), __FILE__, __LINE__)

#define VERSION_STRING "3.0"

/********
  The program will compile to use doubles instead of reals to represent
everything.
  If for some reason you want to use reals, then you must put -DUSE_FLOATS in
the
  CFLAGS field of the makefile
*********/
#ifndef USE_FLOATS
#define real double
#else
#define real float
#endif

#define MAX_STR_LEN 2048

/* fundamental constants */
#define BOHR .52918

/* some numbers to save work later. */
#define SQRT3 1.732050807568877
#define SQRT6 2.449489742783178
#define SQRT10 3.16227766016838
#define SQRT15 3.872983346207417
#define AUI 1.889644746
#define TWOPI 6.283185307179586

/******
  this is a little trick to convert T's and F's into 0's and 1's
******/
#define BOOL_CHAR_TO_VALUE(a) a = (a == 'T' ? 1 : 0)

#ifndef USING_THE_MAC
#define TRUE 1
#define FALSE 0
#endif

/*****
  used as tabs into the matrices
******/
#define BEGIN_S 0
#define BEGIN_P 1
#define END_P 3
#define BEGIN_D 4
#define END_D 8
#define BEGIN_F 9
#define END_F 15

#ifndef PI
#define PI 3.141592653589793
#endif
#define TWOPI 6.283185307179586

/* THE constant! */
#define THE_CONST 1.750000

/******
  These are used by the fileio routines.
*******/
#define FATAL 0
#define ERROR 1
#define IGNORE 2

/******
  Modes the program can run in... these determine memory usage and execution
   speed.
******/
#define FAT 6
#define THIN 1
#define MOLECULAR 27

/* used as generic indicators */
#define NORMAL 0
#define RESET 44

/* this is used to determine whether or not energy levels are degenerate */
#define DEGEN_TOL .000001

/* used for grouping peaks in the DOS diagram */
#define DOS_DEGEN_TOL .001

/* used to terminate the self consistent zeta variation */
#define ZETA_TOL .0001

/*******
  this is the multiple for the change in zeta with charge...
  .35 is Slater's value for shielding between electrons in the same
    shell.
*******/
#define ZETA_SCALE .35
#define ZETA_SCALE_1S .30
#define NEXT_SHELL_IN .15

/*******
  these are used to internally label printing options

  if you want to add printing options along walsh diagrams, then you
  need to add a new #define to this list for your printing option.

*******/
#define PRT_OP 0
#define PRT_ROP 1
#define PRT_OVERLAP 2
#define PRT_NET_CHG 3
#define PRT_CHG_MAT 4
#define PRT_RCHG_MAT 5
#define PRT_WAVE_FUNC 6
#define PRT_ELECTROSTAT 7
#define PRT_FERMI 8
#define PRT_DIST 9
#define PRT_HAMIL 10
#define PRT_AVG_E 11
#define PRT_AVG_OP 12
#define PRT_AVG_ROP 13
#define PRT_ENERGIES 14
#define PRT_ORB_ENERGY 15
#define PRT_ORB_COEFF 16
#define PRT_LEVELS 17
#define PRT_TRANSPOSE_FLAG 64

/************
  which parts of a matrix get labelled
************/
#define LABEL_ROWS 1
#define LABEL_COLS 2
#define LABEL_BOTH 3

/*******
  types of projected DOS
********/
#define P_DOS_ATOM 3
#define P_DOS_ORB 2
#define P_DOS_FMO 1

/******
  This is used to convert reals to integers
*******/
#define MULTIPLIER 1

/******
  This is the default longest distance between atoms in nearest neighbor
  cells that will be reported in the output file.
******/
#define NN_DISTANCE 2.5

/* offset from high symmetry points for automatic k-point generation */
#define K_OFFSET 0.01

/* max length of an atom symbol */
#define ATOM_SYMB_LEN 10

/*******

  used for systems that affix underscores to the names of
  fortran object files.

********/
#ifdef UNDERSCORE_FORTRAN
#define abfns abfns_
#define lovlap lovlap_
#define cboris cboris_
#define zhegv zhegv_
#define zheev zheev_
#endif

#define ABS(a) ((a) > 0 ? (a) : -(a))

#ifndef USE_BZERO
#define bzero(a, b) (memset((void *)(a), 0, (b)))
#define bcopy(a, b, c) (memcpy((void *)(b), (const void *)(a), (c)))
#endif

/* default values for the muller iteration procedure */
#define MULLER_MIX_DEF 0.20
#define MULLER_E_TOL_DEF 0.01
#define MULLER_Z_TOL_DEF 0.001

/****************************  type definitions *******************/

typedef char BOOLEAN;

#ifdef INCLUDE_NETCDF_SUPPORT
typedef struct {
  int file_ID;
  int num_orbs_ID, num_KPOINTS_ID, tot_num_orbs_ID, num_atoms_ID;
  int num_dim_ID, space_dim_ID, name_len_ID, param_len_ID;
} netCDF_info_type;
#endif

/*****

  a 3-D point

*******/
typedef struct { real x, y, z; } point_type;

/******

  a Z matrix element

*******/
typedef struct {
  int ref1, ref2, ref3;
  real bond_length;
  real alpha; /* the bond angle */
  real beta;  /* the dihedral */
} Z_mat_type;

/*********

  an atom

**********/
typedef struct {
  char symb[ATOM_SYMB_LEN];
  char chg_it_vary;
  int which_atom;
  int at_number;
  int num_valence;
  point_type loc;

  Z_mat_type Zmat_loc;

  int ns, np, nd, nf;                  /* principle quantum numbers */
  real coul_s, coul_p, coul_d, coul_f; /* coulomb integrals */
  real exp_s, exp_p, exp_d, exp_f, exp_d2, exp_f2; /* exponents */
  real coeff_d1, coeff_d2, coeff_f1,
      coeff_f2; /* coeffs for double zeta functions */
  /* charge iteration parms */
  real s_A, s_B, s_C;
  real p_A, p_B, p_C;
  real d_A, d_B, d_C;

  /******

    parameters for Edgar Muller's iterative scheme

  ******/
  real muller_s_E[7], muller_p_E[7], muller_d_E[7];
  real muller_s_Z[4], muller_p_Z[4], muller_d_Z[4];
  real init_s_occup, init_p_occup, init_d_occup;

} atom_type;

/* used to allow geometrical fragments */
typedef struct geom_frag_def {
  int which;
  char using_Z_mat;
  int num_atoms;
  atom_type *atoms;
  struct geom_frag_def *next;
} geom_frag_type;

/* used to allow specification with crystal coordinates */
typedef struct {
  real axis_lengths[3];
  real angles[3];
} xtal_defn_type;

/* used to define translation vectors */
typedef struct { int begin, end; } Tvect_type;

/* used to keep track of equivalent atoms for average properties */
typedef struct equiv_atom_def equiv_atom_type;
struct equiv_atom_def {
  int num_equiv;
  /* the list of equivalent atoms */
  int *equiv_atoms;

  /* tabs into the array of present symmetry operations */
  int *sym_ops;
  equiv_atom_type *next;
};

/**********

  a unit cell

************/
typedef struct {
  int dim; /* dimension of the cell */
  int num_atoms, num_raw_atoms;
  atom_type *atoms;

  geom_frag_type *geom_frags;

  char using_Zmat, using_xtal_coords;

  /******
    the distance matrix for the unit cell is stored as a symmetric matrix
      (see notes.outl for the format)
  *******/
  real *distance_mat;

  real num_electrons;
  real charge;

  char *sym_elems;

  Tvect_type tvects[3];      /* the translation vectors */
  point_type recip_vects[3]; /* the reciprocal space vectors */

  xtal_defn_type xtal_defn; /* allows specification with xtal coords */

  /* the number of overlaps in each direction */
  int overlaps[3];

  /* the location of the center of mass */
  point_type COM;

  /* the principle axes */
  real princ_axes[3][3];

  /* the list of equivalent atoms */
  equiv_atom_type *equiv_atoms;
} cell_type;

/***********

  k points

************/
typedef struct {
  point_type loc;
  real weight;
  real num_filled_bands; /* this needs to be a real to deal with degeneracies */
} k_point_type;

/**********

  special points (for band structures)

**********/
typedef struct {
  point_type loc;
  char name[80];
} special_point_type;

/***********

  This structure is used for dealing with performing band structure
   calculations

************/
typedef struct {
  int num_special_points;
  int points_per_line;

  /* this is the list of special points */
  special_point_type *special_points;

  /* one array will be used to hold all the kpoints which are gonna be used */
  k_point_type *lines;
} band_info_type;

/*************

  this structure is used to allow almost anything to be printed
  at each Walsh step (or each k point, or anywhere else)

**************/
typedef struct printing_info_def printing_info_type;
struct printing_info_def {
  char which_to_print;

  char type; /* this is either atom or orbital, or whatever */

  int contrib1, contrib2; /* used to indicate orbitals or atoms or FMOs */

  /* implemented as a linked list to make things easier */
  printing_info_type *next;
};

/*********

  stuff for a Walsh diagram

**********/
typedef struct {
  int num_steps;
  int num_vars;
  real *values;

  printing_info_type *things_to_print;
} walsh_details_type;

/*******

  This structure is used to keep track of which overlaps the user has
  specified should be turned off.

*********/
typedef struct {
  int type;
  int inter_cell;
  int which1, which2;
} overlap_cancel_type;

/*******************

  This is the data type used to keep track of which
   projected DOS's are being calculated.

********************/
typedef struct {
  int type; /* what kind of DOS it will be */

  int num_contributions;
  int *contributions;
  real *weights;
  real weight_sum;
} p_DOS_type;

/********************

  This is the data type for COOP's

   these are kept in a doubly linked list so that averaging is
    easy.

    The structure is:

    Types of COOP -------->

 things to   [ type 1 ]  ->  [ type 2] ->  [ type 3]
 average        |                                  |
   |            |                              |
   |                V                              V
   |         [ type 1 ]                    [ type 3]
   V

*********************/
typedef struct COOP_def_type COOP_type;

struct COOP_def_type {
  int type; /* which type of COOP */

  int energy_weight;

  int which;

  point_type cell; /* allows COOP's outside the unit cell */

  int contrib1, contrib2;

  real avg_value;

  COOP_type *next_type;
  COOP_type *next_to_avg;
};

/****************
  These structures (prop_type and avg_prop_info_type) contain pointers
   to the matrices/vectors/numbers/whatever that are used for
   properties and the data used to determine average properties.

  This is a way to isolate the details of which properties
    are calculated from most of the program, allowing easier
    addition of new properties.

  If changes are made here (i.e. properties added or deleted),
    the memory allocation routine must be changed, as well as
    any code that makes reference to particular members of the
    structures (obviously).
******************/
typedef struct {
  real *OP_mat;   /* mulliken overlap population */
  real *ROP_mat;  /* reduced mulliken overlap population */
  real *chg_mat;  /* charge matrix */
  real *Rchg_mat; /* reduced charge matrix */
  real *net_chgs; /* net charges */

  real *mod_OP_mat;   /* modified mulliken overlap population */
  real *mod_ROP_mat;  /* reduced modified mulliken overlap population */
  real *mod_net_chgs; /* net charges from modified mulliken analysis */

  real electrostat_E;
  real Fermi_E;
  real total_E;
} prop_type;

/********
  this is the information that will be used to calculate average
  properties in extended structure calculations.
*********/
typedef struct {
  float *orbs, *orbsI;         /* the wavefunctions */
  float *S;                    /* S(k) (the overlap matrix) */
  float *chg_mat;              /* the charge matrix */
  float *FMO_orbs, *FMO_orbsI; /* the wavefunctions in the FMO basis*/
  float *FMO_chg_mat;          /* the charge matrix in the FMO basis*/
  float *energies;
} avg_prop_info_type;

/*****
  This macro is used to pick out elements of the overlap matrix
  stored in an avg_prop_info_type structure
******/
#define S_ELEMENT_R(mat, dim, i, j)                                            \
  (j > i ? (real)mat[i * dim + j] / (real)MULTIPLIER                           \
         : (real)mat[j * dim + i] / (real)MULTIPLIER)
#define S_ELEMENT_I(mat, dim, i, j)                                            \
  (j > i ? (real)mat[j * dim + i] / (real)MULTIPLIER                           \
         : (j == i ? 0 : -(real)mat[i * dim + j] / (real)MULTIPLIER))

/*******
  this is used to order the orbitals when doing extended properties
  calculations
********/
typedef struct {
  real occup; /* used to store the occupation of this orbital */
  int Kpoint;
  int MO;
  float *energy;
} K_orb_ptr_type;

/*******************
  these data structures are used to conceal the specifics of the
  data storage from the programmer.

  i.e. they allow one to change the representation of the data
    without having to change the code that uses that data.

  NOTE: each data type should have a macro allowing access to the
    individual elements... for the sake of clarity those macros are
    defined immediately after the data types.

  Some parts of the code will still need to be changed if the representation,
  is changed, for example anywhere memory is allocated

  Also, for now, anywhere values are put into the matrices will need to be
  changed.
  This will remain true until I get around to fixing up all that stuff and
  rewriting the diagonalization code, which could very well be never.
*******************/

/* eigenvectors/values */
typedef struct {
  int dim;
  real *vectR, *vectI;
  real *val;
} eigenset_type;

#define EIGENVECT_R(set, row, column) set.vectR[row * set.dim + column]
#define EIGENVECT_I(set, row, column) set.vectI[row * set.dim + column]
#define EIGENVAL(set, which) set.val[which]

/* used for hermetian matrices */
typedef struct {
  int dim;
  real *mat;
} hermetian_matrix_type;

#define HERMETIAN_R(matrix, i, j)                                              \
  (j > i ? matrix.mat[i * matrix.dim + j] : matrix.mat[j * matrix.dim + i])
#define HERMETIAN_I(matrix, i, j)                                              \
  (j > i ? matrix.mat[j * matrix.dim + i] : matrix.mat[i * matrix.dim + j])

/* A generic complex matrix */
typedef struct {
  int dim;
  real *matR, *matI;
} complex_matrix_type;

/**********

  the structure used to store fragments for FMO analysis.

  each fragment will keep track of it's own eigenvectors, overlap
    matrices and properties.

***********/
typedef struct {
  /* FMO printing options */
  BOOLEAN distance_mat_PRT, OP_mat_PRT, ROP_mat_PRT;
  BOOLEAN chg_mat_PRT, Rchg_mat_PRT, wave_fn_PRT;
  BOOLEAN net_chg_PRT, overlap_mat_PRT, electrostat_PRT, fermi_PRT;
  BOOLEAN hamil_PRT, energies_PRT;
  cell_type *frag_cell;
  int num_atoms;
  int *atoms_in_frag;
  int num_orbs;
  real num_electrons;
  eigenset_type eigenset;
  hermetian_matrix_type overlap_R, overlap_K;
  hermetian_matrix_type hamil_R, hamil_K;
  complex_matrix_type tform_matrix;
  avg_prop_info_type *avg_prop_info;
  prop_type *properties;
  int *orbital_lookup_table;
} FMO_frag_type;

/******

  used to store properties information in the FMO basis

******/
typedef struct {
  eigenset_type eigenset;
  hermetian_matrix_type hamil, overlap;
  real *OP_mat;   /* mulliken overlap population */
  real *ROP_mat;  /* reduced mulliken overlap population */
  real *chg_mat;  /* charge matrix */
  real *Rchg_mat; /* reduced charge matrix */
  real *net_chgs; /* net charges on each fragment */
} FMO_prop_type;

/*********

  charge iteration parameters

**********/
typedef struct {
  char variable_step;
  int max_it;
  real adjust;
  real lambda;
  real damp1, damp2, damp3, lampri;
  real tolerance;
} chg_it_parm_type;

/********

  used to specify orbital occupations

*********/
typedef struct {
  int orb;
  real occup;
} orbital_occup_type;

/*******

  "experimental" details

********/
typedef struct {
  char title[240];
  char filename[240];
  int Execution_Mode;

  /* some toggles */
  BOOLEAN just_geom, avg_props, gradients, save_energies;
  BOOLEAN use_symmetry, find_princ_axes;
  BOOLEAN vary_zeta;
  BOOLEAN eval_electrostat;
  BOOLEAN weighted_Hij;

  /* for dumping binary files containing the matrices */
  BOOLEAN dump_overlap, dump_hamil;

  /* for dumping sparse matrix files */
  BOOLEAN dump_sparse_mats;

  /* for dumping the distance matrix (for find_coops) */
  BOOLEAN dump_dist_mat;

  /*******
    printing options
  ********/

  /* limits on which levels are printed */
  int upper_level_PRT, lower_level_PRT;
  real max_dist_PRT; /* maximum distance that will be printed */
  BOOLEAN distance_mat_PRT, OP_mat_PRT, ROP_mat_PRT;
  BOOLEAN chg_mat_PRT, Rchg_mat_PRT, wave_fn_PRT;
  BOOLEAN net_chg_PRT, overlap_mat_PRT, electrostat_PRT, fermi_PRT;
  BOOLEAN hamil_PRT, energies_PRT, levels_PRT;
  BOOLEAN avg_OP_mat_PRT, avg_ROP_mat_PRT;
  BOOLEAN no_total_DOS_PRT;
  BOOLEAN just_avgE, just_matrices;
  BOOLEAN mod_OP_mat_PRT, mod_ROP_mat_PRT, mod_net_chg_PRT;
  BOOLEAN orbital_mapping_PRT;

  /* the line width for printed output */
  int line_width;

  char diag_wo_overlap;
  char store_R_overlaps;

  /******
    these are printing options that will be done at each walsh
     step/whatever
  ******/
  printing_info_type *step_print_options;

  /* stuff for printing MO coefficients */
  int num_MOs_to_print;
  int *MOs_to_print;

  /* whatever the hell this is */
  real rho;

  /* goodies for charge iteration */
  char do_chg_it;
  chg_it_parm_type chg_it_parms;

  /* this is the cutoff for printing of close nearest neighbor contacts */
  real close_nn_contact;

  /* for dealing with multiple occupations */
  int num_occup_AVG;
  real occup_AVG_step;

  /* used to occupy specific orbitals */
  int num_orbital_occups;
  orbital_occup_type *orbital_occups;

  /* the k points */
  int num_KPOINTS;
  k_point_type *K_POINTS;

  /* setup for automagic k-points */
  char use_automatic_kpoints, use_high_symm_p;
  int points_per_axis[3];
  real k_offset;

  /* multiple occupations at each k point */
  int num_occup_KPOINTS;
  real *occup_KPOINTS;

  /* overlap population stuff */
  int num_bonds_OOP;

  /* FMO things */
  int num_FMO_frags;
  int num_FCO_frags;
  FMO_frag_type *FMO_frags;
  FMO_prop_type *FMO_props;

  /* DOS details */
  int num_proj_DOS;
  p_DOS_type *proj_DOS;

  /* COOP stuff */
  COOP_type *the_COOPS;

  /* do a moments analysis? */
  char do_moments;
  int num_moments;
  real *moments;

  /* info for the walsh diagram */
  walsh_details_type walsh_details;

  /* info for band structures. */
  band_info_type *band_info;

  /* used to allow the user to switch off overlaps */
  int num_overlaps_off;
  overlap_cancel_type *overlaps_off;

  /* this is K, the constant used for the Hij's */
  real the_const;

  real sparsify_value;

  /*******
    the tolerance for atoms being considered equivalent in the
    symmetry analysis
  *******/
  real symm_tol;
  /****
    stuff for symmetry analysis
  ****/
  int num_sym_ops;
  real *characters;

  /*******
    stuff for doing Muller iteration
  ********/
  char do_muller_it;
  int *atoms_to_vary;
  real muller_mix, muller_E_tol, muller_Z_tol;

#ifdef INCLUDE_NETCDF_SUPPORT
  /*******
    stuff for writing netCDF files for post-processing by
    other applications
  *********/
  char do_netCDF;
  /****
    define this as a pointer so we don't have to do a make clean
    every time something new gets added to the netCDF structure
  ****/
  netCDF_info_type *netCDF_info;
#endif

} detail_type;

/* these are for faking the Fortran complex and complex16 types */
typedef struct { double r, i; } doublecomplex;

typedef struct { float r, i; } floatcomplex;

#ifdef USE_FLOATS
typedef floatcomplex complex;
#else
typedef doublecomplex complex;
#endif

/****** globals *******/
extern FILE *status_file, *output_file, *walsh_file, *band_file, *FMO_file;
extern FILE *MO_file;
extern int temp_file;
extern cell_type *unit_cell;
extern detail_type *details;
extern eigenset_type eigenset;
extern hermetian_matrix_type Hamil_R, Hamil_K;
extern hermetian_matrix_type Overlap_R, Overlap_K;
extern complex *cmplx_hamil, *cmplx_overlap, *cmplx_work;
extern real *work1, *work2, *work3;
extern real *OP_mat, *net_chgs;
extern int num_orbs, tot_overlaps;
extern int *orbital_lookup_table;
extern atom_type *unique_atoms;
extern int num_unique_atoms;
extern prop_type properties;
extern avg_prop_info_type *avg_prop_info;
extern K_orb_ptr_type *orbital_ordering;

extern real electrostatic_term, eHMO_term, total_energy;

#include "prototypes.h"
