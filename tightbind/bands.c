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
*   These are the things needed for doing a band structure calculation
*
*   Please note that this is not a "smart" band structure routine,
*    the bands are not generated using the symmetry of the wavefunctions.
*    Everything is handled by using a large number of points along
*    each symmetry line and assuming that that will take care of
*    avoided crossings.
*
*  created:  greg landrum  June 1994
*
*****************************************************************************/

/******
  Recent Edit History:

  28.01.99 gL:  print "#BAND_DATA" to first line of .band file.

******/
#include "bind.h"





/****************************************************************************
 *
 *                   Procedure gen_symm_lines
 *
 * Arguments:  bands: pointer to band_info_type
 *
 * Returns: none
 *
 * Action:  Automagically generates the kpoints along each symmetry line.
 *
 ****************************************************************************/
void gen_symm_lines(bands)
  band_info_type *bands;
{
  int i,j;
  int points_per_line;
  point_type spacing,curr_loc;

  points_per_line = bands->points_per_line;

  /* loop over the special points */
  for( i=0; i<bands->num_special_points-1; i++){

    /****
      determine the spacing between k points in each direction
      along this line
      ****/
    spacing.x = (bands->special_points[i+1].loc.x -
                 bands->special_points[i].loc.x) / points_per_line;
    spacing.y = (bands->special_points[i+1].loc.y -
                 bands->special_points[i].loc.y) / points_per_line;
    spacing.z = (bands->special_points[i+1].loc.z -
                 bands->special_points[i].loc.z) / points_per_line;

    curr_loc.x = bands->special_points[i].loc.x;
    curr_loc.y = bands->special_points[i].loc.y;
    curr_loc.z = bands->special_points[i].loc.z;

    for(j=0;j<points_per_line;j++){
      bands->lines[i*points_per_line+j].loc.x = curr_loc.x;
      bands->lines[i*points_per_line+j].loc.y = curr_loc.y;
      bands->lines[i*points_per_line+j].loc.z = curr_loc.z;

      curr_loc.x += spacing.x;
      curr_loc.y += spacing.y;
      curr_loc.z += spacing.z;
    }

    /* now put in the last special point */
    bands->lines[(bands->num_special_points-1)*points_per_line].loc.x =
      bands->special_points[bands->num_special_points-1].loc.x;
    bands->lines[(bands->num_special_points-1)*points_per_line].loc.y =
      bands->special_points[bands->num_special_points-1].loc.y;
    bands->lines[(bands->num_special_points-1)*points_per_line].loc.z =
      bands->special_points[bands->num_special_points-1].loc.z;
  }
  /* that's all she wrote folks. */
}

/****************************************************************************
 *
 *                   Procedure construct_band_structure
 *
 * Arguments:  cell: pointer to cell type
 *          details: pointer to detail type
 *          overlapR: hermetian_matrix_type
 *            hamilR: hermetian_matrix_type
 *          overlapK: hermetian_matrix_type
 *            hamilK: hermetian_matrix_type
 *   cmplx_hamil, cmplx_overlap: pointers to complex
 *          eigenset: eigenset_type
 * work1,work2,work3: pointers to reals
 *        cmplx_work: pointer to complex
 *          num_orbs: int
 * orbital_lookup_table: pointer to int.
 *
 * Returns: none
 *
 * Action:  Loops through the k points which make up the symmetry
 *  lines which were specified for this calculation.
 *  At each point the following things must be done:
 *    - Generate the K space overlap and hamiltonian matrices.
 *    - solve the eigenvalue eqn:
 *      H(k) * Y = S(k) * E * Y
 *    - write any required info to the output file
 *
 *
 *   The work arrays are used as temporary memory in the various functions called
 *    by this one.  The dimensions should be:
 *      work1,work2: num_orbs;
 *            work3: num_orbs*num_orbs;
 *
 *
 *  NOTE:  Since the only things which are used from these calculations are
 *   the eigenvalues, it's really un-necessary and inefficient to fully
 *   diagonalize the matrices.  It would be much better to just find the
 *   eigenvalues.... Mabye I'll get to this later.
 *
 *
 *  If the LAPACK diagonalization routine is used, only the eigenvalues
 *     are generated.  Efficiency!
 *
 ****************************************************************************/
void construct_band_structure(cell,details,overlapR,hamilR,overlapK,hamilK,
                              cmplx_hamil,cmplx_overlap,
                              eigenset,work1,work2,work3,cmplx_work,
                              num_orbs,orbital_lookup_table)
  cell_type *cell;
  detail_type *details;
  hermetian_matrix_type overlapR,hamilR;
  hermetian_matrix_type overlapK,hamilK;
  complex *cmplx_hamil,*cmplx_overlap;
  eigenset_type eigenset;
  real *work1,*work2,*work3;
  complex *cmplx_work;
  int num_orbs;
  int *orbital_lookup_table;
{
  static char (*label)[4]=0;
  k_point_type *kpoint;
  int i,j,k,l,m;
  int jtab,ktab,ltab,mtab;
  int diag_error;
  real temp;
  real total_energy,tot_chg;
  real *mat_save;
  int electrons_so_far;
  real *occupations;
  int num_KPOINTS;
  band_info_type *bands;
#ifdef USE_LAPACK
  char jobz, uplo;
  int info, itype;
  int num_orbs2;
#endif

  if( details->Execution_Mode == FAT && !details->store_R_overlaps )
    mat_save = overlapK.mat;

  bands = details->band_info;

  num_KPOINTS = (bands->num_special_points-1) *
                  bands->points_per_line + 1;

  /* write out some status information */
  fprintf(status_file,"Generating band structure.\n");

  /*****
    put the important information about this calculation into the output file.
  ******/
  fprintf(band_file,"#BAND_DATA\n");
  fprintf(band_file,";Band structure calculation results\n");
  fprintf(band_file,"%d Special points were done with\n",
          bands->num_special_points);
  fprintf(band_file,"%d k points connecting them. There are\n",
          bands->points_per_line);

  fprintf(band_file,"%d orbitals in the unit cell.\n",num_orbs);

  for(i=0;i<bands->num_special_points;i++){
    fprintf(band_file,"%s %lf %lf %lf\n",
            bands->special_points[i].name,bands->special_points[i].loc.x,
            bands->special_points[i].loc.y,bands->special_points[i].loc.z);
  }

  fprintf(band_file,"; Begin band data.\n");

  /********

    here's the loop over the k point set.

  ********/
  for(i=0;i<num_KPOINTS;i++){
    /* get a pointer to the k point we're working on */
    kpoint = &(bands->lines[i]);


    /*****

      build the overlap matrix and hamiltonian

    *****/
    switch(details->Execution_Mode){
    case FAT:
      if( details->store_R_overlaps ){
        build_k_overlap_FAT(cell,kpoint,overlapR,overlapK,num_orbs);
      } else{
        fatal("construct_band_structure called while store_R_overlaps was 0.\n \
\tThis is wrong.");
      }
      build_k_hamil_FAT(cell,hamilR,hamilK,overlapK,num_orbs);
      break;
    case THIN:
      build_k_overlap_THIN(cell,details,kpoint,overlapR,overlapK,num_orbs);
      build_k_hamil_THIN(cell,hamilR,hamilK,overlapK,num_orbs);
      break;
    default:
      FATAL_BUG("Somehow a bogus execution mode got passed to generate_band_structure.");
    }
//fprintf(stderr,">");

    /******
      The matrix diagonalization routine destroys the overlap matrix,
       but that's okay, because we don't need it again later when doing band
       structures.  skipping a bcopy will save some time.
    *******/

#ifndef USE_LAPACK
    /*******

      now diagonalize that beast by calling the FORTRAN subroutine used
       to diagonalize stuff in new3 and CACAO.

      THIS REALLY SHOULD BE REPLACED with a routine written in C, so if you
       happen to have some time on your hands....

    ********/
    cboris(&(num_orbs),&(num_orbs),hamilK.mat,overlapK.mat,eigenset.vectI,
           eigenset.val,work1,
           work2,&diag_error);

#else
    for(j=0;j<num_orbs;j++){
      jtab = j*num_orbs;
      for(k=j+1;k<num_orbs;k++){
        ktab = k*num_orbs;
        cmplx_hamil[jtab+k].r = hamilK.mat[jtab+k];
        cmplx_hamil[jtab+k].i = hamilK.mat[ktab+j];
        cmplx_overlap[jtab+k].r = overlapK.mat[jtab+k];
        cmplx_overlap[jtab+k].i = overlapK.mat[ktab+j];
        cmplx_hamil[ktab+j].r = 0.0;
        cmplx_hamil[ktab+j].i = 0.0;
        cmplx_overlap[ktab+j].r = 0.0;
        cmplx_overlap[ktab+j].i = 0.0;
      }
      cmplx_hamil[jtab+j].r = hamilK.mat[jtab+j];
      cmplx_hamil[jtab+j].i = 0.0;
      cmplx_overlap[jtab+j].r = overlapK.mat[jtab+j];
      cmplx_overlap[jtab+j].i = 0.0;

    }

    itype = 1;
    jobz = 'N';
    uplo = 'L';
    num_orbs2 = num_orbs*num_orbs;
    fprintf(stderr,"{");
    if(!details->diag_wo_overlap){
      zhegv((long *)&(itype),&jobz,&uplo,(long *)&num_orbs,cmplx_hamil,
                              (long *)&num_orbs,cmplx_overlap,
                  (long *)&num_orbs,eigenset.val,cmplx_work,
                  (long *)&num_orbs2,work3,(long *)&diag_error);
    }else{
      zheev(&jobz,&uplo,(long *)&num_orbs,cmplx_hamil,(long *)&num_orbs,
                                    eigenset.val,cmplx_work,(long *)&num_orbs2,work3,(long *)&diag_error);
    }
    fprintf(stderr,"}");
#endif

    //fprintf(stderr,"<\n");
    fprintf(status_file,"Error value from Diagonalization (0 is good): %d\n",diag_error);


    /*******
      write just the energies to the output file.
    ********/
    fprintf(band_file,"; K point: %lf %lf %lf\n",
            kpoint->loc.x,kpoint->loc.y,kpoint->loc.z);
    for(j=0;j<num_orbs;j++){
      fprintf(band_file,"%10.8lg\n",EIGENVAL(eigenset,j));
    }

  } /* end of k point loop */
  if( details->Execution_Mode == FAT && !details->store_R_overlaps ){
    overlapK.mat = mat_save;
  }

  // Indicate that we have finished the band data
  fprintf(band_file, "#END_BAND_DATA\n");
}




