/*******************************************************
*      Copyright (C) 1995 Greg Landrum
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



/****************************************************************************
*
*     this has all the global variables for the program bind
*
*  created:  greg landrum  August 1993
*
*****************************************************************************/

#include "bind.h"
#include "symmetry.h"

FILE *status_file,*output_file,*walsh_file,*band_file, *FMO_file;
FILE *MO_file;
int temp_file;

cell_type *unit_cell;
detail_type *details;

/* the matrices */
eigenset_type eigenset;
hermetian_matrix_type Hamil_R,Hamil_K;
hermetian_matrix_type Overlap_R,Overlap_K;

complex *cmplx_hamil,*cmplx_overlap,*cmplx_work;

real *work1,*work2,*work3;
real *OP_mat,*net_chgs;

/* dimensions */
int num_orbs,tot_overlaps;

int *orbital_lookup_table;

atom_type *unique_atoms;
int num_unique_atoms;

real electrostatic_term,eHMO_term,total_energy;

sym_op_type *sym_ops_present=0;

prop_type properties;
avg_prop_info_type *avg_prop_info;
K_orb_ptr_type *orbital_ordering;
