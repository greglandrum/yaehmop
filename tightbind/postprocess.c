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
*     this file has the routines for dealing with postprocessing the
*      wavefunctions, etc. that are generated at each k point.
*
*  created:  greg landrum  June 1996
*
*****************************************************************************/
#include "bind.h"

/****************************************************************************
 *
 *                   Procedure postprocess_FMO
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
 ****************************************************************************/
void postprocess_FMO(cell,details,overlapR,hamilR,overlapK,hamilK,
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
  real *occupations;
  real *chg_mat;

  int i,j,k;
  int itab,jtab,ktab;
  real tot_chg;


  /********
    use the work2 array to store the occupation numbers for later use
    (in case they become needed later....)

    do this by setting a pointer to work2 to make things a little more
    readable.
    *********/
  occupations = work2;

  /**************

    do the FMO transformation and print out the results

    ***************/
  if( details->num_FMO_frags != 0 ){
    tform_wavefuncs_to_FMO_basis(details,num_orbs,cell->num_atoms,eigenset,
                                 orbital_lookup_table);

    if( details->wave_fn_PRT ){
      fprintf(output_file,"\n;\t Wavefunctions in FMO basis \n");
      printmat(details->FMO_props->eigenset.vectR,num_orbs,num_orbs,output_file,
               1e-5,details->wave_fn_PRT & PRT_TRANSPOSE_FLAG,details->line_width);

      if( details->Execution_Mode != MOLECULAR ){
        printmat(details->FMO_props->eigenset.vectI,num_orbs,num_orbs,
                 output_file,1e-5,details->wave_fn_PRT & PRT_TRANSPOSE_FLAG,
                 details->line_width);
      }
      fprintf(output_file,"\n");
    }

    tform_hermetian_matrix_to_FMO_basis(details,num_orbs,cell->num_atoms,
                                        overlapK,
                                        work3,details->FMO_props->chg_mat,
                                        details->FMO_props->overlap,
                                        orbital_lookup_table);

    /* do we need to print out the overlap matrix? */
    if( details->overlap_mat_PRT ){
      fprintf(output_file,
              ";\t\t --- Overlap Matrix in FMO basis");
      if( details->Execution_Mode == MOLECULAR ){
        fprintf(output_file,"FMO_S(R) ---\n");
      }
      else{
        fprintf(output_file,"FMO_S(K) ---\n");
      }
      printmat(details->FMO_props->overlap.mat,num_orbs,num_orbs,output_file,
               1e-4, details->overlap_mat_PRT & PRT_TRANSPOSE_FLAG,details->line_width);
    }
    /* what about the hamiltonian? */
    if( details->hamil_PRT ){
      tform_hermetian_matrix_to_FMO_basis(details,num_orbs,cell->num_atoms,
                                          hamilK,
                                          work3,details->FMO_props->chg_mat,
                                          details->FMO_props->hamil,
                                          orbital_lookup_table);
      fprintf(output_file,
              ";\t\t --- Hamiltonian Matrix in FMO basis");
      if( details->Execution_Mode == MOLECULAR ){
        fprintf(output_file,"FMO_H(R) ---\n");
      }
      else{
        fprintf(output_file,"FMO_H(K) ---\n");
      }
      printmat(details->FMO_props->hamil.mat,num_orbs,num_orbs,output_file,
               1e-4, details->overlap_mat_PRT & PRT_TRANSPOSE_FLAG,details->line_width);
    }

    if( details->OP_mat_PRT || details->ROP_mat_PRT || details->net_chg_PRT
       || details->vary_zeta ){
      eval_mulliken(cell,details->FMO_props->eigenset,
                    details->FMO_props->overlap,num_orbs,
                    occupations,orbital_lookup_table,
                    details->FMO_props->OP_mat,
                    work3,work1);
    }

    /* this stores the reduced overlap population matrix */
    if( details->ROP_mat_PRT ){
      FMO_reduced_mulliken(details,cell->num_atoms,num_orbs,
                           details->FMO_props->OP_mat,
                           details->FMO_props->ROP_mat);
    }
    if( details->OP_mat_PRT ){
      fprintf(output_file,
              "\n\n; \tq-q-q-q-q-q-q-q  Mulliken Analysis in FMO basis q-q-q-q-q-q-q-q\n");

      fprintf(output_file,
              "\n;    FMO Mulliken Overlap Population Matrix for %6.3lf electrons:\n",
              cell->num_electrons);
      printmat(details->FMO_props->OP_mat,num_orbs,num_orbs,output_file,
               1e-05,details->OP_mat_PRT & PRT_TRANSPOSE_FLAG,details->line_width);

    }
    if( details->ROP_mat_PRT ){
      fprintf(output_file,
              "\n;   Reduced FMO Mulliken Overlap Population Matrix for %6.3lf electrons:\n",
              cell->num_electrons);
      print_sym_mat(details->FMO_props->ROP_mat,details->num_FMO_frags,
                    details->num_FMO_frags,
                    output_file,(char *)0,(char *)0,details->line_width);
    }

    /* find the charge matrix in the FMO basis */
    eval_charge_matrix(cell,details->FMO_props->eigenset,
                       details->FMO_props->overlap,
                       num_orbs,orbital_lookup_table,
                       details->FMO_props->chg_mat,work1);

    if( details->chg_mat_PRT || details->Rchg_mat_PRT ){
      fprintf(output_file,"\n\n# \tq-F-q-F-q-F-q-F  Charge Matrix in FMO basis F-q-F-q-F-q-F-q\n");
    }
    if( details->chg_mat_PRT ){
      fprintf(output_file,";  Complete FMO Charge Matrix Independant of Occupation \
<normalized to 2 electrons>\n");
      printmat(details->FMO_props->chg_mat,num_orbs,num_orbs,output_file,
               1e-05,details->chg_mat_PRT & PRT_TRANSPOSE_FLAG,details->line_width);
    }

    /*******

      if we're not doing an extended calculation, write some information into
      the FMO output file

      *******/

    /* first the molecular orbital energies need to be written */
    if( details->Execution_Mode == MOLECULAR ){
      fprintf(FMO_file,"; Molecular orbital energies\n");
      for(j=0;j<num_orbs;j++){
        fprintf(FMO_file,"%lg\n",EIGENVAL(eigenset,j));
      }
      /* now write out the elements of the charge matrix for each MO */
      fprintf(FMO_file,"; charge matrix\n");
      chg_mat = details->FMO_props->chg_mat;
      for(j=0;j<num_orbs;j++){
        jtab = j*num_orbs;
        for(k=0;k<num_orbs;k++){
          if( fabs(chg_mat[jtab+k]) >= 1e-06 ){
            fprintf(FMO_file,"%6.4lg ",chg_mat[jtab+k]);
          }
          else{
            fprintf(FMO_file,"%6.4lg ",0.0);
          }
        }
        fprintf(FMO_file,"\n");
      }
    }
  } else FATAL_BUG("postprocess_FMO called with no FMO fragments");
}


/****************************************************************************
 *
 *                   Procedure postprocess_FCO
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
 ****************************************************************************/
void postprocess_FCO(cell,details,overlapR,hamilR,overlapK,hamilK,
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

  static int FCO_file=0;
  char FCO_filename[240];

  real tot_K_weight;
  real *occupations;
  real *chg_mat;
  int i,j,k;
  int itab,jtab,ktab;
  int test_int;
  real tot_chg;

  /********
    use the work2 array to store the occupation numbers for later use
    (in case they become needed later....)

    do this by setting a pointer to work2 to make things a little more
    readable.
    *********/
  occupations = work2;

#if 0
  if( details->Execution_Mode != MOLECULAR &&
     details->num_FCO_frags != 0 ){
#endif
    if( details->num_FCO_frags != 0 ){
    /******

      on the first call, we need to open the FCO output file
      and write some header information.

    *******/
    if(!FCO_file){
      sprintf(FCO_filename,"%s.FCO",details->filename);
#ifndef USING_THE_MAC
      FCO_file = open(FCO_filename,O_RDWR|O_TRUNC|O_CREAT,S_IRUSR|S_IWUSR);
#else
      FCO_file = open(FCO_filename,O_RDWR|O_APPEND|O_CREAT);
#endif
      if( FCO_file == -1 ){
        fatal("Can't open .FCO file for binary I/O");
      }

      /* it's open, now write the header... */

      /*************

        we start with an int of known value to allow us to
        check at read time that things are fine (i.e. it's a similar
        machine type).  Otherwise fit_FCO will read garbage.

      **************/
      test_int = 4231;
      write(FCO_file,(const char *)&(test_int),sizeof(int));
      write(FCO_file,(const char *)&(details->num_KPOINTS),sizeof(int));
      tot_K_weight = 0.0;
      for(i=0;i<details->num_KPOINTS;i++){
        tot_K_weight += details->K_POINTS[i].weight;
      }
      write(FCO_file,(const char *)&(tot_K_weight),sizeof(real));
      write(FCO_file,(const char *)&(num_orbs),sizeof(int));
      write(FCO_file,(const char *)&(details->num_FCO_frags),sizeof(int));
      for(i=0;i<details->num_FCO_frags;i++){
        write(FCO_file,(const char *)&(details->FMO_frags[i].num_orbs),sizeof(int));
      }
    }


    /**************

      do the FCO transformation and print out the results

      ***************/

    tform_wavefuncs_to_FMO_basis(details,num_orbs,cell->num_atoms,eigenset,
                                 orbital_lookup_table);

    if( details->wave_fn_PRT ){
      fprintf(output_file,"\n;\t Wavefunctions in FCO basis \n");
      printmat(details->FMO_props->eigenset.vectR,num_orbs,num_orbs,output_file,
               1e-5,details->wave_fn_PRT & PRT_TRANSPOSE_FLAG,details->line_width);

      if( details->Execution_Mode != MOLECULAR ){
        printmat(details->FMO_props->eigenset.vectI,num_orbs,num_orbs,
                 output_file,1e-5,details->wave_fn_PRT & PRT_TRANSPOSE_FLAG,details->line_width);
      }
      fprintf(output_file,"\n");
    }

    tform_hermetian_matrix_to_FMO_basis(details,num_orbs,cell->num_atoms,
                                        overlapK,
                                        work3,details->FMO_props->chg_mat,
                                        details->FMO_props->overlap,
                                        orbital_lookup_table);

    /* do we need to print out the overlap matrix? */
    if( details->overlap_mat_PRT ){
      fprintf(output_file,
              ";\t\t --- Overlap Matrix in FCO basis");
      if( details->Execution_Mode == MOLECULAR ){
        fprintf(output_file,"FCO_S(R) ---\n");
      }
      else{
        fprintf(output_file,"FCO_S(K) ---\n");
      }
      printmat(details->FMO_props->overlap.mat,num_orbs,num_orbs,output_file,
               1e-4, details->overlap_mat_PRT & PRT_TRANSPOSE_FLAG,details->line_width);
    }
    /* what about the hamiltonian? */
    if( details->hamil_PRT ){
      tform_hermetian_matrix_to_FMO_basis(details,num_orbs,cell->num_atoms,
                                          hamilK,
                                          work3,details->FMO_props->chg_mat,
                                          details->FMO_props->hamil,
                                          orbital_lookup_table);
      fprintf(output_file,
              ";\t\t --- Hamiltonian Matrix in FCO basis");
      if( details->Execution_Mode == MOLECULAR ){
        fprintf(output_file,"FCO_H(R) ---\n");
      }
      else{
        fprintf(output_file,"FCO_H(K) ---\n");
      }
      printmat(details->FMO_props->hamil.mat,num_orbs,num_orbs,output_file,
               1e-4, details->overlap_mat_PRT & PRT_TRANSPOSE_FLAG,details->line_width);
    }

    if( details->OP_mat_PRT || details->ROP_mat_PRT || details->net_chg_PRT
       || details->vary_zeta ){
      eval_mulliken(cell,details->FMO_props->eigenset,
                    details->FMO_props->overlap,num_orbs,
                    occupations,orbital_lookup_table,
                    details->FMO_props->OP_mat,
                    work3,work1);
    }

    /* this stores the reduced overlap population matrix */
    if( details->ROP_mat_PRT ){
      FMO_reduced_mulliken(details,cell->num_atoms,num_orbs,
                           details->FMO_props->OP_mat,
                           details->FMO_props->ROP_mat);
    }
    if( details->OP_mat_PRT ){
      fprintf(output_file,
              "\n\n; \tq-q-q-q-q-q-q-q  Mulliken Analysis in FCO basis q-q-q-q-q-q-q-q\n");

      fprintf(output_file,
              "\n;    FCO Mulliken Overlap Population Matrix for %6.3lf electrons:\n",
              cell->num_electrons);
      printmat(details->FMO_props->OP_mat,num_orbs,num_orbs,output_file,
               1e-05,details->OP_mat_PRT & PRT_TRANSPOSE_FLAG,details->line_width);

    }
    if( details->ROP_mat_PRT ){
      fprintf(output_file,
              "\n;   Reduced FCO Mulliken Overlap Population Matrix for %6.3lf electrons:\n",
              cell->num_electrons);
      print_sym_mat(details->FMO_props->ROP_mat,details->num_FMO_frags,
                    details->num_FMO_frags,
                    output_file,(char *)0,(char *)0,details->line_width);
    }

    /* find the charge matrix in the FMO basis */
    eval_charge_matrix(cell,details->FMO_props->eigenset,
                       details->FMO_props->overlap,
                       num_orbs,orbital_lookup_table,
                       details->FMO_props->chg_mat,work1);

    if( details->chg_mat_PRT || details->Rchg_mat_PRT ){
      fprintf(output_file,"\n\n# \tq-F-q-F-q-F-q-F  Charge Matrix in FCO basis F-q-F-q-F-q-F-q\n");
    }
    if( details->chg_mat_PRT ){
      fprintf(output_file,";  Complete FCO Charge Matrix Independant of Occupation \
<normalized to 2 electrons>\n");
      printmat(details->FMO_props->chg_mat,num_orbs,num_orbs,output_file,
               1e-05,details->chg_mat_PRT & PRT_TRANSPOSE_FLAG,details->line_width);
    }

    /*******

      write some information into
      the FCO output file

    *******/

    /* start out with the energies */
    write(FCO_file,(const char *)(eigenset.val),num_orbs*sizeof(real));
    for(i=0;i<details->num_FCO_frags;i++){
      write(FCO_file,(const char *)(details->FMO_frags[i].eigenset.val),
            details->FMO_frags[i].num_orbs*sizeof(real));
    }
    /* write the charge matrix */
    write(FCO_file,(const char *)details->FMO_props->chg_mat,
          num_orbs*num_orbs*sizeof(real));


  } else FATAL_BUG("postprocess_FCO called with no FMO fragments");
}




/****************************************************************************
 *
 *                   Procedure postprocess_results
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
 * Action:  this is basically a place to collect all the crap that
 *   needs to be done after a diagonalization is complete.
 *   putting it all in here makes the k point loop a lot cleaner
 *   and easier to deal with.
 *
 ****************************************************************************/
void postprocess_results(cell,details,overlapR,hamilR,overlapK,hamilK,
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

  real *occupations;
  real *chg_mat;
  int i,j,k;
  int itab,jtab,ktab;
  real tot_chg;


  /* do we need to print out the wave functions? */
  if( details->wave_fn_PRT ){
    if( !(details->wave_fn_PRT & PRT_TRANSPOSE_FLAG) ){
      fprintf(output_file,
              ";\t*** Wavefunctions *** (MO's in rows, AO's in columns, energies \
INCREASE down)\n");
    }else{
      fprintf(output_file,
              ";\t*** Wavefunctions ***  (AO's in rows, MO's in columns, energies \
INCREASE to the right)\n");
    }

    fprintf(output_file,";\t***> REAL:\n");
    print_labelled_mat(eigenset.vectR,num_orbs,num_orbs,output_file,1e-4,
                       cell->atoms,cell->num_atoms,orbital_lookup_table,
                       num_orbs,details->wave_fn_PRT & PRT_TRANSPOSE_FLAG,
                       LABEL_COLS,details->line_width);

    /********

      there's no need to print out imaginary components of the
      wavefunctions for molecular calculations.

      ********/
    if( details->Execution_Mode != MOLECULAR ){
      fprintf(output_file,";\t***> IMAGINARY:\n");

      print_labelled_mat(eigenset.vectI,num_orbs,num_orbs,output_file,1e-4,
                         cell->atoms,cell->num_atoms,orbital_lookup_table,
                         num_orbs,details->wave_fn_PRT & PRT_TRANSPOSE_FLAG,
                         LABEL_COLS,details->line_width);
    }
  }

  /* print out the energies and the total energy */
  total_energy = 0;

  /********
    use the work2 array to store the occupation numbers for later use
    (in case they become needed later....)

    do this by setting a pointer to work2 to make things a little more
    readable.
    *********/
  occupations = work2;

  /* zero out the occupations array */
  bzero((char *)occupations,num_orbs*sizeof(real));

  calc_occupations(details,cell->num_electrons,num_orbs,occupations,eigenset);

  if( details->levels_PRT || details->Execution_Mode == MOLECULAR){
    fprintf(output_file,"\n#\t******* Energies (in eV)  and Occupation Numbers *******\n");
    for(j=0;j<num_orbs;j++){
      fprintf(output_file,"%d:--->  %8.6lg  [%4.3lf Electrons]\n",j+1,
              EIGENVAL(eigenset,j), occupations[j]);
      total_energy += occupations[j]*EIGENVAL(eigenset,j);
    }
    fprintf(output_file,"Total_Energy: %8.6lg\n",total_energy);
    properties->total_E = total_energy;
  }

  /******************

    Determine the symmetry of the wavefunctions if that is necessary

    *******************/
  if( details->use_symmetry && details->Execution_Mode == MOLECULAR){
    find_MO_symmetries(num_orbs,details,cell,eigenset,overlapK,
                       orbital_lookup_table);
  }

  /******************

    Do the properties calculations

    *******************/

  /* only do mulliken population analysis if we need to */
  if( details->OP_mat_PRT || details->ROP_mat_PRT || details->net_chg_PRT
     || details->vary_zeta || details->mod_OP_mat_PRT ||
     details->mod_ROP_mat_PRT || details->mod_net_chg_PRT ){
    /* this stores the overlap population matrix and the net charges */
    eval_mulliken(cell,eigenset,overlapK,num_orbs,
                  occupations,orbital_lookup_table,properties->OP_mat,
                  properties->net_chgs,work1);

    /* this stores the reduced overlap population matrix */
    if( details->ROP_mat_PRT ){
      reduced_mulliken(cell->num_atoms,num_orbs,orbital_lookup_table,
                       properties->OP_mat,properties->ROP_mat);
    }

    if( details->mod_OP_mat_PRT || details->mod_ROP_mat_PRT ){
      modified_mulliken(cell,eigenset,overlapK,num_orbs,
                        occupations,orbital_lookup_table,properties->OP_mat,
                        properties->mod_OP_mat,
                        properties->mod_net_chgs,work1);

    }
    /*************

      for modified mulliken, we can use the standard
      reduced_mulliken routine, because the structure
      of the OP matrix is the same, just the way
      in which it is evaluated has been changed

      **************/
    if(details->mod_ROP_mat_PRT ){
      reduced_mulliken(cell->num_atoms,num_orbs,orbital_lookup_table,
                       properties->mod_OP_mat,properties->mod_ROP_mat);
    }
  }
  if( details->OP_mat_PRT ){
    fprintf(output_file,
            "\n\n; \t\tq-q-q-q-q-q-q-q  Mulliken Analysis q-q-q-q-q-q-q-q\n");

    fprintf(output_file,
            "\n;    Mulliken Overlap Population Matrix for %6.3lf electrons:\n",
            cell->num_electrons);
    print_labelled_mat(properties->OP_mat,num_orbs,num_orbs,output_file,1e-5,
                       cell->atoms,cell->num_atoms,orbital_lookup_table,
                       num_orbs,details->OP_mat_PRT & PRT_TRANSPOSE_FLAG,LABEL_BOTH,
                       details->line_width);
  }
  if( details->ROP_mat_PRT ){
    fprintf(output_file,
            "\n;   Reduced Mulliken Overlap Population Matrix for %6.3lf electrons:\n",
            cell->num_electrons);
    print_sym_mat(properties->ROP_mat,cell->num_atoms,cell->num_atoms,output_file,(char *)0,
                  (char *)0,details->line_width);
  }

  if( details->net_chg_PRT ){
    fprintf(output_file,"\n;    Net Atomic Charges for %6.3lf electrons:\n",
            cell->num_electrons);
    tot_chg = 0.0;
    for(j=0;j<cell->num_atoms;j++){
      tot_chg += properties->net_chgs[j];
      fprintf(output_file,"%d %s: %8.6lf\n", j+1, cell->atoms[j].symb,
              properties->net_chgs[j]);
    }
    fprintf(output_file,";      Total Charge is: %8.6lf\n",tot_chg);
  }



  if( details->mod_OP_mat_PRT ){
    fprintf(output_file,
            "\n\n; \t\tq-q-q-q-q-q  Modified Mulliken Analysis q-q-q-q-q-q\n");

    fprintf(output_file,
            "\n;    Modified Mulliken Overlap Population Matrix for %6.3lf electrons:\n",
            cell->num_electrons);
    print_labelled_mat(properties->mod_OP_mat,num_orbs,num_orbs,output_file,1e-5,
                       cell->atoms,cell->num_atoms,orbital_lookup_table,
                       num_orbs,details->mod_OP_mat_PRT & PRT_TRANSPOSE_FLAG,
                       LABEL_BOTH,details->line_width);
  }
  if( details->mod_ROP_mat_PRT ){
    fprintf(output_file,
            "\n;   Reduced Modified Mulliken Overlap Population Matrix for %6.3lf electrons:\n",
            cell->num_electrons);
    print_sym_mat(properties->mod_ROP_mat,cell->num_atoms,cell->num_atoms,
                  output_file,(char *)0,(char *)0,details->line_width);
  }

  if( details->mod_net_chg_PRT ){
    fprintf(output_file,"\n;   Modified Mulliken Net Atomic Charges for %6.3lf electrons:\n",
            cell->num_electrons);
    tot_chg = 0.0;
    for(j=0;j<cell->num_atoms;j++){
      tot_chg += properties->mod_net_chgs[j];
      fprintf(output_file,"%d %s: %8.6lf\n", j+1, cell->atoms[j].symb,
              properties->mod_net_chgs[j]);
    }
    fprintf(output_file,";      Total Charge is: %8.6lf\n",tot_chg);
  }


  /********

    The charge matrix needs to be evaluated not only when it is being
    printed out, but also when average properties calculations
    are being done.  The charge matrix is used in these instances
    to determine projected densities of states.

    *********/
  if( details->chg_mat_PRT || details->Rchg_mat_PRT ||
     (details->avg_props && !details->no_total_DOS_PRT)){
    /* find the charge matrix */
    eval_charge_matrix(cell,eigenset,overlapK,num_orbs,orbital_lookup_table,
                       properties->chg_mat,work1);

    if( details->sparsify_value > 0.0 ){
      fprintf(stderr,"Charge matrix sparsification\n");
      sparsify_matrix(details->sparsify_value,properties->chg_mat,0,
                      num_orbs);
    }


    if( details->chg_mat_PRT || details->Rchg_mat_PRT ){
      fprintf(output_file,"\n\n; \t\tq-q-q-q-q-q-q-q  Charge Matrix q-q-q-q-q-q-q-q\n");
    }
    if( details->chg_mat_PRT ){
      fprintf(output_file,
              ";  Complete Charge Matrix Independant of Occupation\n");

      print_labelled_mat(properties->chg_mat,num_orbs,num_orbs,output_file,1e-5,
                         cell->atoms,cell->num_atoms,orbital_lookup_table,
                         num_orbs,details->chg_mat_PRT & PRT_TRANSPOSE_FLAG,
                         LABEL_COLS,details->line_width);

    }
    /******
      put in reduced charge matrix stuff here
      *******/
  }

}
