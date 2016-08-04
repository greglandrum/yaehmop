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
*     this file has the routines for looping over kpoints and dealing with the
*      results.
*
*  created:  greg landrum  September 1993
*
*****************************************************************************/
#include "bind.h"

/****
  Recent edit history

  26.09.98 gL:
   swapped order of postprocess_results and print_MOs
****/


/****************************************************************************
 *
 *                   Procedure loop_over_k_points
 *
 * Arguments:  cell: pointer to cell type
 *          details: pointer to detail type
 *          overlapR: hermetian_matrix_type
 *            hamilR: hermetian_matrix_type
 *          overlapK: hermetian_matrix_type
 *            hamilK: hermetian_matrix_type
 *   cmplx_hamil, cmplx_overlap: pointers to complex
 *          eigenset: eigenset_type
 *        properties: pointer to prop_type
 *     avg_prop_info: pointer to avg_prop_info_type
 * work1,work2,work3: pointers to reals
 *        cmplx_work: pointer to complex
 *          num_orbs: int
 * orbital_lookup_table: pointer to int.
 *
 * Returns: none
 *
 * Action:  Loops over all the k-points for this run.  At each k-point the following
 *     things are done:
 *    - Generate the K space overlap and hamiltonian matrices.
 *    - solve the eigenvalue eqn:
 *      H(k) * Y = S(k) * E * Y
 *    - do any required calculations for this k-point (properties, etc.)
 *    - write any required info to the output file
 *
 *   Some helpful definitions:
 *
 *
 *              N                     |              N
 *            -----                   |            -----
 *             \                      |             \
 *    H(k) :=   )   exp(-i k . r) H(r)|    S(k) :=   )   exp(-i k . r) S(r)
 *             /                      |             /
 *            -----                   |            -----
 *            r = 1                   |            r = 1
 *                                    |
 *                                    |
 *
 *   The work arrays are used as temporary memory in the various functions called
 *    by this one.  The dimensions should be:
 *         work1,work2: num_orbs;
 *    work3,cmplx_work: num_orbs*num_orbs;
 *
 *
 *
 ****************************************************************************/
void loop_over_k_points(cell,details,overlapR,hamilR,overlapK,hamilK,
                        cmplx_hamil,cmplx_overlap,
                        eigenset,work1,work2,work3,cmplx_work,
                        properties,avg_prop_info,
                        num_orbs,orbital_lookup_table)
  cell_type *cell;
  detail_type *details;
  hermetian_matrix_type overlapR,hamilR;
  hermetian_matrix_type overlapK,hamilK;
  complex *cmplx_hamil,*cmplx_overlap;
  eigenset_type eigenset;
  real *work1,*work2,*work3;
  complex *cmplx_work;
  prop_type *properties;
  avg_prop_info_type *avg_prop_info;
  int num_orbs;
  int *orbital_lookup_table;
{
  static char tempfilename[240];
  static char *label=0;
  static FILE *sparse_OVfile,*sparse_HAMfile;
  k_point_type *kpoint;
  real *mat_save;
  int i,j,k,l,m;
  int itab,jtab,ktab;
  int ltab,mtab;
  int diag_error;
  real temp;
  real total_energy,tot_chg;
  int electrons_so_far;
  real *occupations;
  real *chg_mat;
  int num_KPOINTS;
  int overlap_file,hamil_file;
#ifdef USE_LAPACK
  char jobz, uplo;
  int info, itype;
  int num_orbs2;
#endif

  if( details->Execution_Mode == FAT && !details->store_R_overlaps )
    mat_save = overlapK.mat;

  if( !label ){
    /* get space for the array storing the atomic labels... */
    label = (char *)calloc(4*cell->num_atoms,sizeof(char));
    if( !label )fatal("Can't allocate label array in loop_over_k_points");

    /* put each atom in the label array */
    for(i=0;i<cell->num_atoms;i++){
      bcopy(cell->atoms[i].symb,&(label[4*i]),4*sizeof(char));
    }
  }

  /* make sure that we loop once for a molecular calculation */
  if( details->Execution_Mode != MOLECULAR ){
    num_KPOINTS = details->num_KPOINTS;
  }
  else{
    num_KPOINTS = 1;
  }
  /********

    here's the loop over the k point set.

  ********/
#ifdef IRIX_MP
#pragma parallel
#pragma local (kpoint, overlapK, hamilK, work3, eigenset, work1, work2, diag_error, occupations, i, total_energy, j, tot_chg, k, jtab, chg_mat)
#pragma shared (properties, avg_prop_info, details, cell,  hamilR,num_KPOINTS, num_orbs,orbital_lookup_table,overlapR )
#pragma pfor iterate (i=0;num_KPOINTS;1)
#endif
  for(i=0;i<num_KPOINTS;i++){
    /* get a pointer to the k point we're working on */
    kpoint = &(details->K_POINTS[i]);

    /* print some status information */
    if( cell->dim > 0){
      fprintf(status_file,"Kpoint: %d\n",i+1);
      fprintf(output_file,";***& Kpoint: %d (%6.4lf %6.4lf %6.4lf) Weight: %lf\n",i+1,
              kpoint->loc.x,kpoint->loc.y,kpoint->loc.z,kpoint->weight);
    }

    /*****

      build the overlap matrix and hamiltonian

      *****/
    switch(details->Execution_Mode){
    case FAT:
      if( details->store_R_overlaps ){
        build_k_overlap_FAT(cell,kpoint,overlapR,overlapK,num_orbs);
      } else{
        overlapK.mat = &mat_save[i*num_orbs*num_orbs];
      }

      if( details->sparsify_value > 0.0 ){
        fprintf(stderr,"Overlap Sparsification\n");
        sparsify_hermetian_matrix(details->sparsify_value,
                                  overlapK,num_orbs);
      }
      build_k_hamil_FAT(cell,hamilR,hamilK,overlapK,num_orbs);
      break;
    case THIN:
      build_k_overlap_THIN(cell,details,kpoint,overlapR,overlapK,num_orbs);
      build_k_hamil_THIN(cell,hamilR,hamilK,overlapK,num_orbs);
      break;
    case MOLECULAR:
      /******

        for molecular calculations, we don't need to evaluate anything in k space,
        so just set the pointers here....

        *******/
      overlapK = overlapR;
      if( details->sparsify_value > 0.0 ){
        fprintf(stderr,"Overlap Sparsification\n");
        sparsify_hermetian_matrix(details->sparsify_value,
                                  overlapK,num_orbs);
      }
      hamilK = hamilR;
      break;
    default:
      FATAL_BUG("Somehow a bogus execution mode got passed to loop_over_kpoints.");
    }

    /* do we need to print out the overlap matrix? */
    if( details->overlap_mat_PRT ){
      fprintf(output_file,
              ";\t\t --- Overlap Matrix ");
      if( details->Execution_Mode == MOLECULAR ){
        fprintf(output_file,"S(R) ---\n");
      }
      else{
        fprintf(output_file,"S(K) ---\n");
      }
      print_labelled_mat(overlapK.mat,num_orbs,num_orbs,output_file,1e-4,
                         cell->atoms,cell->num_atoms,orbital_lookup_table,
                         num_orbs,details->overlap_mat_PRT & PRT_TRANSPOSE_FLAG,
                         LABEL_BOTH,details->line_width);
    }

    /* What about the hamiltonian? */
    if( details->hamil_PRT ){
      fprintf(output_file,
              ";\t\t --- Hamiltonian ");
      if( details->Execution_Mode == MOLECULAR ){
        fprintf(output_file,"H(R) ---\n");
      }
      else{
        fprintf(output_file,"H(K) ---\n");
      }
      print_labelled_mat(hamilK.mat,num_orbs,num_orbs,output_file,1e-4,
                         cell->atoms,cell->num_atoms,orbital_lookup_table,
                         num_orbs,details->hamil_PRT & PRT_TRANSPOSE_FLAG,
                         LABEL_BOTH,details->line_width);
    }

    /* do we need to do binary dumps of the matrices? */
    if( details->dump_overlap ){
      /* if this is the first call, the open the file */
      if(i==0){
        sprintf(tempfilename,"%s.OV",details->filename);

#ifndef USING_THE_MAC
        overlap_file = open(tempfilename,
                            O_RDWR|O_TRUNC|O_APPEND|O_CREAT,S_IRUSR|S_IWUSR);
#else
        overlap_file = open(tempfilename,O_RDWR|O_TRUNC|O_APPEND|O_CREAT);
#endif

        if( overlap_file == -1 ){
          fatal("Can't open .OV file for binary I/O");
        }
        write(overlap_file,(const char *)&num_KPOINTS,sizeof(int));
        write(overlap_file,(const char *)&num_orbs,sizeof(int));

      }
      dump_hermetian_mat(overlap_file,overlapK.mat,num_orbs);
    }

    if( details->dump_sparse_mats ){
      /* do we need to open the sparse matrix files? */
      if(i==0){
        sprintf(tempfilename,"%s.SPARSE.OV",details->filename);

        sparse_OVfile = fopen(tempfilename,"w+");
        if( !sparse_OVfile ) fatal("Can't open .SPARSE.OV file!\n");
        sprintf(tempfilename,"%s.SPARSE.HAM",details->filename);
        sparse_HAMfile = fopen(tempfilename,"w+");
        if( !sparse_HAMfile ) fatal("Can't open .SPARSE.HAM file!\n");
      }
      dump_sparse_mat(sparse_OVfile,overlapK.mat,num_orbs,1e-08);
      dump_sparse_mat(sparse_HAMfile,hamilK.mat,num_orbs,1e-08);
    }

    if( details->dump_hamil ){
      /* if this is the first call, the open the file */
      if(i==0){
        sprintf(tempfilename,"%s.HAM",details->filename);
#ifndef USING_THE_MAC
        hamil_file = open(tempfilename,
                          O_RDWR|O_APPEND|O_CREAT,S_IRUSR|S_IWUSR);
#else
        hamil_file = open(tempfilename,O_RDWR|O_APPEND|O_CREAT);
#endif
        if( hamil_file == -1 ){
          fatal("Can't open .HAM file for binary I/O");
        }
        write(hamil_file,(const char *)&num_KPOINTS,sizeof(int));
        write(hamil_file,(const char *)&num_orbs,sizeof(int));

      }
      dump_hermetian_mat(hamil_file,hamilK.mat,num_orbs);
    }

    if( !details->just_matrices ){

      /*****

        do the FCO analysis if it's needed

        if we're doing a molecule, the appropriate matrices have
        already been built, so we don't have to worry about those

      ******/
      if(details->num_FCO_frags && details->Execution_Mode != MOLECULAR){
        /* first build the matrices */
        build_FMO_overlap(details,num_orbs,unit_cell->num_atoms,Overlap_K,
                          orbital_lookup_table);
        build_FMO_hamil(details,num_orbs,unit_cell->num_atoms,Hamil_K,
                        orbital_lookup_table);
        /* now diagonalize them */
        diagonalize_FMO(details,work1,work2,work3,cmplx_hamil,cmplx_overlap,cmplx_work);

        /* generate the transform matrices */
        gen_FMO_tform_matrices(details);
      }
      //fprintf(stdout,"%d >",i+1);

#ifndef USE_LAPACK
      /******
        The matrix diagonalization routine destroys the overlap and hamiltonian
        matrices, so if we need to (i.e. we are printing elements of them)
        we make a copy of the overlap matrix in work3 and the
        hamiltonian matrix in eigenset.vectR.  We'll move things around
        later to get everything straightened out.
        *******/
      if(!details->diag_wo_overlap){
        bcopy((char *)overlapK.mat,(char *)work3,num_orbs*num_orbs*sizeof(real));
      } else {
        bzero((char *)work3,num_orbs*num_orbs*sizeof(real));
        for(j=0;j<num_orbs;j++) work3[j*num_orbs+j] = 1.0;
      }
      if( details->hamil_PRT ){
        bcopy((char *)hamilK.mat,(char *)eigenset.vectR,num_orbs*num_orbs*sizeof(real));
      }

      /*******

        now diagonalize that beast by calling the FORTRAN subroutine used
        to diagonalize stuff in new3 and CACAO.

        THIS REALLY SHOULD BE REPLACED with a routine written in C, so if you
        happen to have some time on your hands....

        ********/
      cboris(&(num_orbs),&(num_orbs),hamilK.mat,work3,eigenset.vectI,eigenset.val,work1,
             work2,&diag_error);

      /********

        This is some comic relief aimed at members of the Hoffmann group.
        If you want to do something similar for your site, uncomment this
        section of code and change the uid's (you can find these in the
        file /etc/passwd) and messages.

        ********/
#if 0
      switch(getuid()){
      case 1426: fprintf(stderr,"Jahn-Teller is REAL!"); break;
      case 1501: fprintf(stderr,"Ultimate Man!"); break;
      case 1649: fprintf(stderr,"Done Fishing?"); break;
      case 1559: fprintf(stderr,"Damn texan!"); break;
      case 1622: fprintf(stderr,"Back to the library!"); break;
      case 1645: fprintf(stderr,"More Helices?"); break;
      }
#endif

      /*********

        at this point, hamilK.mat contains the real part of the eigenvectors,
        eigenset.vectI contains the imaginary part,
        eigenset.val has the energies,
        and eigenset.vectR contains the hamiltonian matrix.

        rearrange things so that eigenset.vectR and hamilK.mat store the
        proper information.


        **********/
      if( details->hamil_PRT ){
        bcopy((char *)eigenset.vectR,(char *)work3,num_orbs*num_orbs*sizeof(real));
        bcopy((char *)hamilK.mat,(char *)eigenset.vectR,num_orbs*num_orbs*sizeof(real));
        bcopy((char *)work3,(char *)hamilK.mat,num_orbs*num_orbs*sizeof(real));
      } else{
        bcopy((char *)hamilK.mat,(char *)eigenset.vectR,num_orbs*num_orbs*sizeof(real));
      }

#else
      /**********

        we're using LAPACK to diagonalize and we need to copy the matrices into those
        used by the LAPACK diagonalizer

        **********/
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
      if( details->just_avgE ){
        jobz = 'N';
        //fprintf(stdout,".");
      } else{
        jobz = 'V';
      }
      uplo = 'L';
      num_orbs2 = num_orbs*num_orbs;
      fprintf(stdout,"{");
      if(!details->diag_wo_overlap){
        zhegv((long *)&(itype),&jobz,&uplo,(long *)&num_orbs,cmplx_hamil,
                                (long *)&num_orbs,cmplx_overlap,
              (long *)&num_orbs,eigenset.val,cmplx_work,
              (long *)&num_orbs2,work3,(long *)&diag_error);
      }else{
        zheev(&jobz,&uplo,(long *)&num_orbs,cmplx_hamil,(long *)&num_orbs,
              eigenset.val,cmplx_work,(long *)&num_orbs2,work3,
              (long *)&diag_error);
      }
      fprintf(stdout,"}");

      /* now copy stuff back out of the results */
      if( !details->just_avgE ){
        for(j=0;j<num_orbs;j++){
          jtab = j*num_orbs;
          for(k=0;k<num_orbs;k++){
            ktab = k*num_orbs;
            eigenset.vectR[jtab+k] = cmplx_hamil[jtab+k].r;
            eigenset.vectI[jtab+k] = cmplx_hamil[jtab+k].i;
          }
        }
      }
#endif

      //fprintf(stdout,"<\n");
      fprintf(status_file,"Error value from Diagonalization (0 is good): %d\n",
              diag_error);
      fflush(status_file);
      if( diag_error != 0 ){
        error("Problems in the diagonalization, try more overlaps.");
      }

      if( !details->just_avgE ){

        postprocess_results(cell,details,overlapR,hamilR,overlapK,hamilK,
                            cmplx_hamil,cmplx_overlap,
                            eigenset,work1,work2,work3,cmplx_work,
                            properties,avg_prop_info,
                            num_orbs,orbital_lookup_table);

        if(details->num_MOs_to_print > 0) print_MOs(details,num_orbs,eigenset,i,
                                                    unique_atoms,num_unique_atoms,
                                                    cell->num_atoms,orbital_lookup_table);


        if(details->num_FMO_frags != 0 ){
          postprocess_FMO(cell,details,overlapR,hamilR,overlapK,hamilK,
                          cmplx_hamil,cmplx_overlap,
                          eigenset,work1,work2,work3,cmplx_work,
                          properties,avg_prop_info,
                          num_orbs,orbital_lookup_table);

        }

        if(details->num_FCO_frags != 0 ){
          postprocess_FCO(cell,details,overlapR,hamilR,overlapK,hamilK,
                          cmplx_hamil,cmplx_overlap,
                          eigenset,work1,work2,work3,cmplx_work,
                          properties,avg_prop_info,
                          num_orbs,orbital_lookup_table);

        }

      }
      /*********

        if we are doing an extended system calculation and are evaluating
        average properties, then store the information required for the
        properties calculation which will be done later in the execution.

        *********/
      if( details->avg_props ){
        store_avg_prop_info(details,i,eigenset,overlapK,num_orbs,
                            properties->chg_mat,avg_prop_info);
      }
    } /* end of if(!details->just_matrices) */
  } /* end of k point loop */

  if( details->Execution_Mode == FAT && !details->store_R_overlaps ){
    overlapK.mat = mat_save;
  }

  if( details->dump_hamil ) close(hamil_file);
  if( details->dump_overlap) close(overlap_file);

}



