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

bool print_progress = false;
