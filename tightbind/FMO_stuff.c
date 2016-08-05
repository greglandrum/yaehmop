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
*   These are the things needed for doing FMO analysis
*
*
*  created:  greg landrum  October 1994
*
*****************************************************************************/
#include "bind.h"

/****************************************************************************
*
*                   Procedure init_FMO_file
*
* Arguments:    details: pointer to detail type
*              num_orbs: int
*         num_electrons: real
*
* Returns: none
*
* Action: writes the header to the FMO output file.
*
****************************************************************************/
void init_FMO_file(details,num_orbs,num_electrons)
  detail_type *details;
  int num_orbs;
  real num_electrons;
{
  int i;

  /* some error checking */
  if( !details->num_FMO_frags || !details->FMO_frags )
    FATAL_BUG("init_FMO_file called when FMO analysis isn't being done.");

  fprintf(FMO_file,"# FMO Results\n");
  fprintf(FMO_file,"#Total_number_of_orbitals: %d\n",num_orbs);
  fprintf(FMO_file,"#Number_of_electrons: %lf\n",num_electrons);
  fprintf(FMO_file,"#Number_of_fragments: %d\n",details->num_FMO_frags);
  fprintf(FMO_file,"; number of orbitals and electrons in each fragment.\n");

  for(i=0;i<details->num_FMO_frags;i++){
    fprintf(FMO_file,"%d %lf\n",details->FMO_frags[i].num_orbs,
            details->FMO_frags[i].num_electrons);
  }
}

/****************************************************************************
*
*                   Procedure build_FMO_overlap
*
* Arguments:    details: pointer to detail type
*              num_orbs: int
*             num_atoms: int
*               overlap: hermetian_matrix_type
*  orbital_lookup_table: pointer to int
*
* Returns: none
*
* Action: Copies the relevant pieces of the overlap matrix passed in
*    to the correct portions of the FMO structures contained within
*    'details.
*
****************************************************************************/
void build_FMO_overlap(details,num_orbs,num_atoms,overlap,orbital_lookup_table)
  detail_type *details;
  int num_orbs;
  int num_atoms;
  hermetian_matrix_type overlap;
  int *orbital_lookup_table;
{
  FMO_frag_type *FMO_frag;
  int i,j;
  int itab,jtab,FMO_itab,FMO_jtab;
  int begin1,end1,begin2,end2;
  int FMO_begin1,FMO_end1,FMO_begin2,FMO_end2;
  int frag,atom1,atom2;
  int num_frags;

  if( !details->num_FMO_frags ) num_frags = details->num_FCO_frags;
  else num_frags = details->num_FMO_frags;

  if( !details->FMO_frags ){
    FATAL_BUG("Bad FMO_frag_type array passed to build_FMO_overlap.");
  }

  /* loop over fragments */
  for(frag=0;frag<num_frags;frag++){
    FMO_frag = &(details->FMO_frags[frag]);

    /* loop over the atoms which contribute to this fragment */
    for( atom1=0;atom1<FMO_frag->num_atoms;atom1++){
      /******

        the location of the orbitals belonging to this atom must be found
        in BOTH the main overlap matrix and the FMO overlap matrix.

      *******/
      find_atoms_orbs(num_orbs,num_atoms,FMO_frag->atoms_in_frag[atom1],
                      orbital_lookup_table,&begin1,&end1);
      find_atoms_orbs(FMO_frag->num_orbs,FMO_frag->num_atoms,atom1,
                      FMO_frag->orbital_lookup_table,&FMO_begin1,&FMO_end1);
      /* ignore dummy atoms */
      if( begin1 >= 0 && FMO_begin1 >= 0){

        for( atom2=0;atom2<FMO_frag->num_atoms;atom2++){
          /* find the orbitals of the other atom */
          find_atoms_orbs(num_orbs,num_atoms,FMO_frag->atoms_in_frag[atom2],
                          orbital_lookup_table,&begin2,&end2);
          find_atoms_orbs(FMO_frag->num_orbs,FMO_frag->num_atoms,atom2,
                          FMO_frag->orbital_lookup_table,&FMO_begin2,&FMO_end2);

          if( begin2 >= 0 && FMO_begin2 >= 0){
            /* now loop over and copy the orbitals */
            for( i=0; i<(end1-begin1); i++){
              itab = (i+begin1)*num_orbs;
              FMO_itab = (i+FMO_begin1)*FMO_frag->num_orbs;
              for( j=0; j<(end2-begin2); j++){
                jtab = j+begin2;
                FMO_jtab = j+FMO_begin2;

                FMO_frag->overlap_K.mat[FMO_itab+FMO_jtab] =
                  overlap.mat[itab+jtab];
              }
            }
          }
        }
      }
    }
  }
}


/****************************************************************************
*
*                   Procedure build_FMO_hamil
*
* Arguments:    details: pointer to detail type
*              num_orbs: int
*             num_atoms: int
*                 hamil: hermetian_matrix_type
*  orbital_lookup_table: pointer to int
*
* Returns: none
*
* Action: Copies the relevant pieces of the hamiltonian matrix passed in
*    to the correct portions of the FMO structures contained within
*    'details.
*
****************************************************************************/
void build_FMO_hamil(details,num_orbs,num_atoms,hamil,orbital_lookup_table)
  detail_type *details;
  int num_orbs;
  int num_atoms;
  hermetian_matrix_type hamil;
  int *orbital_lookup_table;
{
  FMO_frag_type *FMO_frag;
  int i,j;
  int itab,jtab,FMO_itab,FMO_jtab;
  int begin1,end1,begin2,end2;
  int FMO_begin1,FMO_end1,FMO_begin2,FMO_end2;
  int frag,atom1,atom2;
  int num_frags;

  if( !details->num_FMO_frags ) num_frags = details->num_FCO_frags;
  else num_frags = details->num_FMO_frags;


  if( !details->FMO_frags ){
    FATAL_BUG("Bad FMO_frag_type array passed to build_FMO_hamil.");
  }


  /******

    this is the basically the same code used to build the overlap
    matrices

  *******/

  /* loop over fragments */
  for(frag=0;frag<num_frags;frag++){
    FMO_frag = &(details->FMO_frags[frag]);

    /* loop over the atoms which contribute to this fragment */
    for( atom1=0;atom1<FMO_frag->num_atoms;atom1++){
      /******

        the location of the orbitals belonging to this atom must be found
        in BOTH the main hamiltonian matrix and the FMO hamiltonian.

      *******/
      find_atoms_orbs(num_orbs,num_atoms,FMO_frag->atoms_in_frag[atom1],
                      orbital_lookup_table,&begin1,&end1);
      find_atoms_orbs(FMO_frag->num_orbs,FMO_frag->num_atoms,atom1,
                      FMO_frag->orbital_lookup_table,&FMO_begin1,&FMO_end1);
      /* ignore dummy atoms */
      if( begin1 >= 0 && FMO_begin1 >= 0){

        for( atom2=0;atom2<FMO_frag->num_atoms;atom2++){
          find_atoms_orbs(num_orbs,num_atoms,FMO_frag->atoms_in_frag[atom2],
                          orbital_lookup_table,&begin2,&end2);
          find_atoms_orbs(FMO_frag->num_orbs,FMO_frag->num_atoms,atom2,
                          FMO_frag->orbital_lookup_table,&FMO_begin2,&FMO_end2);

          if( begin2 >= 0 && FMO_begin2 >= 0){
            for( i=0; i<(end1-begin1); i++){
              itab = (i+begin1)*num_orbs;
              FMO_itab = (i+FMO_begin1)*FMO_frag->num_orbs;
              for( j=0; j<(end2-begin2); j++){
                jtab = j+begin2;
                FMO_jtab = j+FMO_begin2;

                FMO_frag->hamil_K.mat[FMO_itab+FMO_jtab] =
                  hamil.mat[itab+jtab];
              }
            }
          }
        }
      }
    }
  }
}


/****************************************************************************
*
*                   Procedure diagonalize_FMO
*
* Arguments:    details: pointer to detail type
*     work1,work2,work3: pointers to reals
*
*
*
* Returns: none
*
* Action: Copies the relevant pieces of the hamiltonian matrix passed in
*    to the correct portions of the FMO structures contained within
*    'details.
*
*   The work arrays are used as temporary memory in the various functions called
*    by this one.  The dimensions should be:
*      work1,work2: num_orbs;
*            work3: num_orbs*num_orbs;
*   where num_orbs is the total number of orbitals in the molecule.
*
*  Note that these are the same work arrays used later to diagonalize
*   the matrices for the full system.  It is not safe to assume that
*   they will continue to hold useful information after this function
*   returns (they will not).
*
****************************************************************************/
void diagonalize_FMO(details,work1,work2,work3,cmplx_hamil,cmplx_overlap,cmplx_work)
  detail_type *details;
  real *work1,*work2,*work3;
  complex *cmplx_hamil,*cmplx_overlap,*cmplx_work;
{
  FMO_frag_type *FMO_frag;
  int i,j,k,itab,jtab,ktab;
  int num_orbs,diag_error;
  real *occupations;
  int num_frags;
#ifdef USE_LAPACK
  char jobz, uplo;
  int info, itype;
  int num_orbs2;
#endif

  if( !details->num_FMO_frags ) num_frags = details->num_FCO_frags;
  else num_frags = details->num_FMO_frags;

  if( !details->FMO_frags ){
    FATAL_BUG("Bad FMO_frag_type array passed to diagonalize_FMO.");
  }

  fprintf(output_file,";------------------ FMO Analysis --------------\n");
  fprintf(output_file,"#NUM_FRAGMENTS: %d\n",num_frags);

  /* loop over each fragment and do it's orbitals individually */
  for( i=0; i<num_frags; i++){
    FMO_frag = &(details->FMO_frags[i]);
    num_orbs = FMO_frag->num_orbs;


    fprintf(output_file,"\n#Fragment %d <*><*><*><*><*><*><*><*><*><*><*><*>\n",i+1);

    if( details->overlap_mat_PRT ){
      fprintf(output_file,
              ";\t\t --- Overlap Matrix ");
      if( details->Execution_Mode == MOLECULAR ){
        fprintf(output_file,"S(R) ---\n");
      }
      else{
        fprintf(output_file,"S(K) ---\n");
      }
      printmat(FMO_frag->overlap_K.mat,num_orbs,num_orbs,output_file,1e-4,
               details->overlap_mat_PRT & PRT_TRANSPOSE_FLAG,details->line_width);
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
      printmat(FMO_frag->hamil_K.mat,num_orbs,num_orbs,output_file,1e-4,
               details->hamil_PRT & PRT_TRANSPOSE_FLAG,details->line_width);
    }

    fprintf(stdout,"]");

#ifndef USE_LAPACK
    /******
      The matrix diagonalization routine destroys the overlap matrix, so we
      make a copy of it in work3.
    *******/
    bcopy((char *)FMO_frag->overlap_K.mat,(char *)work3,num_orbs*num_orbs*sizeof(real));

    /*******

      now diagonalize that beast by calling the FORTRAN subroutine used
       to diagonalize stuff in new3 and CACAO.

    ********/
    cboris(&(num_orbs),&(num_orbs),FMO_frag->hamil_K.mat,
           work3,FMO_frag->eigenset.vectI,FMO_frag->eigenset.val,work1,
           work2,&diag_error);
    fprintf(status_file,"Error value from FMO diagonalization (fragment %d): %d\n",
            i,diag_error);
    fflush(status_file);
    if( diag_error != 0 ){
      error("Problems in the FMO diagonalization, try more overlaps.");
    }

    /*********

      at this point, hamilK.mat contains the real part of the eigenvectors,
      eigenset.vectI contains the imaginary part, and
      eigenset.val has the energies.

      copy the information from hamilK.mat into eigenset.vectR

    **********/
    bcopy((char *)FMO_frag->hamil_K.mat,FMO_frag->eigenset.vectR,
          num_orbs*num_orbs*sizeof(real));
#else
      /**********

        we're using LAPACK to diagonalize and we need to copy the matrices into those
        used by the LAPACK diagonalizer

        **********/
      for(j=0;j<num_orbs;j++){
        jtab = j*num_orbs;
        for(k=j+1;k<num_orbs;k++){
          ktab = k*num_orbs;
          cmplx_hamil[jtab+k].r = FMO_frag->hamil_K.mat[jtab+k];
          cmplx_hamil[jtab+k].i = FMO_frag->hamil_K.mat[ktab+j];
          cmplx_overlap[jtab+k].r = FMO_frag->overlap_K.mat[jtab+k];
          cmplx_overlap[jtab+k].i = FMO_frag->overlap_K.mat[ktab+j];
          cmplx_hamil[ktab+j].r = 0.0;
          cmplx_hamil[ktab+j].i = 0.0;
          cmplx_overlap[ktab+j].r = 0.0;
          cmplx_overlap[ktab+j].i = 0.0;
        }
        cmplx_hamil[jtab+j].r = FMO_frag->hamil_K.mat[jtab+j];
        cmplx_hamil[jtab+j].i = 0.0;
        cmplx_overlap[jtab+j].r = FMO_frag->overlap_K.mat[jtab+j];
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
      zhegv((long *)&(itype),&jobz,&uplo,(long *)&num_orbs,cmplx_hamil,
                              (long *)&num_orbs,cmplx_overlap,
                                           (long *)&num_orbs,FMO_frag->eigenset.val,cmplx_work,
                                           (long *)&num_orbs2,work3,(long *)&diag_error);
      fprintf(stdout,"}");

      /* now copy stuff back out of the results */
      if( !details->just_avgE ){
        for(j=0;j<num_orbs;j++){
          jtab = j*num_orbs;
          for(k=0;k<num_orbs;k++){
            ktab = k*num_orbs;
            FMO_frag->eigenset.vectR[jtab+k] = cmplx_hamil[jtab+k].r;
            FMO_frag->eigenset.vectI[jtab+k] = cmplx_hamil[jtab+k].i;
          }
        }
      }

#endif

    fprintf(stdout,"[ ");
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
      printmat(FMO_frag->eigenset.vectR,num_orbs,num_orbs,output_file,1e-4,
               details->wave_fn_PRT & PRT_TRANSPOSE_FLAG,details->line_width);
      /********

        there's no need to print out imaginary components of the
        wavefunctions for molecular calculations.

      ********/
      if( details->Execution_Mode != MOLECULAR ){
        fprintf(output_file,";\t***> IMAGINARY:\n");
        printmat(FMO_frag->eigenset.vectI,num_orbs,num_orbs,output_file,1e-4,
                 details->wave_fn_PRT & PRT_TRANSPOSE_FLAG,details->line_width);
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

    calc_occupations(0,FMO_frag->num_electrons,num_orbs,occupations,FMO_frag->eigenset);

    fprintf(output_file,"\n#\t****Fragment Energies (in eV) and Occupation Numbers ****\n");
    for(j=0;j<num_orbs;j++){
      fprintf(output_file,"%d:--->  %8.6lg  [%4.3lf Electrons]\n",j+1,
              EIGENVAL(FMO_frag->eigenset,j), occupations[j]);
      total_energy += occupations[j]*EIGENVAL(FMO_frag->eigenset,j);
    }
    fprintf(output_file,"Total_Fragment_Energy: %8.6lg\n",total_energy);

    /******

      FMO properties spoo goes HERE

    ******/

    /* write the relevant information into the FMO results file */
    if( details->num_FMO_frags ){
      fprintf(FMO_file,"; Fragment %d orbital energies\n",i+1);
      for(j=0;j<num_orbs;j++){
        fprintf(FMO_file,"%lg\n",EIGENVAL(FMO_frag->eigenset,j));
      }
    }
  }
  fprintf(stdout,"\n");
  fprintf(output_file,";-o-o-o-o-o-o-o-o- END FMO -o-o-o-o-o-o-\n\n");
}


/****************************************************************************
*
*                   Procedure gen_FMO_tform_matrices
*
* Arguments:    details: pointer to detail type
*
* Returns: none
*
* Action:  This forms the matrix Ct*S for each fragment.
*   (Ct is the transpose of the matrix of coefficients, S is the overlap
*     matrix of the fragment within the AO basis)
*
*  These are the blocks of T, the matrix used to transform the full system's
*   orbitals from the AO to the FMO basis.
*
****************************************************************************/
void gen_FMO_tform_matrices(details)
  detail_type *details;
{
  FMO_frag_type *FMO_frag;
  real *matR, *matI;
  int frag,i,j,k;
  int itab,jtab;
  int num_frags;

  if( !details->num_FMO_frags ) num_frags = details->num_FCO_frags;
  else num_frags = details->num_FMO_frags;

  /* some error checking */
  if( !details->FMO_frags ){
    FATAL_BUG("Bad FMO_frag_type array passed to gen_FMO_tform_matrices.");
  }

  /* loop over the individual fragments */
  for(frag=0;frag<num_frags;frag++){
    FMO_frag = &(details->FMO_frags[frag]);

    /* get pointers to some matrices to save pointer math later */
    matR = FMO_frag->tform_matrix.matR;
    matI = FMO_frag->tform_matrix.matI;

    /* now loop over orbitals within the fragments */
    for(i=0;i<FMO_frag->num_orbs;i++){
      itab = i*FMO_frag->num_orbs;
      for(j=0;j<FMO_frag->num_orbs;j++){
        jtab = j*FMO_frag->num_orbs;

        /******

          No games are played here with complex components, they are
          multiplied normally.  I'm pretty sure that averaging over K points
          should not be done here, it's taken care of later when these
          things are actually used.

        *******/
        matR[itab+j] = 0.0;
        matI[itab+j] = 0.0;
        for(k=0;k<FMO_frag->num_orbs;k++){
          matR[itab+j] += EIGENVECT_R(FMO_frag->eigenset,j,k)*
            S_ELEMENT_R(FMO_frag->overlap_K.mat,FMO_frag->num_orbs,k,i);

          /* no need to worry about imaginary parts if we're doing a molecule */
          if( details->Execution_Mode != MOLECULAR ){

            /*******

              there are almost certainly some problems here with
              the signs of contributions... I'm sure that there should
              be a complex conjugate somewhere that changes some of the
              signs, but I haven't had a chance to derive it yet.

              Thus far, after months of use everything seems fine, but I'm wary.

              BEWARE!!! :-)

            ********/
            matR[itab+j] -= EIGENVECT_I(FMO_frag->eigenset,j,k)*
            S_ELEMENT_I(FMO_frag->overlap_K.mat,FMO_frag->num_orbs,k,i);
            matI[itab+j] -= EIGENVECT_R(FMO_frag->eigenset,j,k)*
            S_ELEMENT_I(FMO_frag->overlap_K.mat,FMO_frag->num_orbs,k,i);
            matI[itab+j] -= EIGENVECT_I(FMO_frag->eigenset,j,k)*
            S_ELEMENT_R(FMO_frag->overlap_K.mat,FMO_frag->num_orbs,k,i);
          }
        }
      }
    }

#ifdef DEBUG
fprintf(output_file,"\n\n\n\t FMO transform matrix (fragment %d)\n",frag);
printmat(matR,FMO_frag->num_orbs,FMO_frag->num_orbs,output_file,1e-5,0,details->line_width);
printmat(matI,FMO_frag->num_orbs,FMO_frag->num_orbs,output_file,1e-5,0,details->line_width);
fprintf(output_file,"\n\n\n");
#endif
  }
}


/****************************************************************************
*
*                   Procedure tform_wavefuncs_to_FMO_basis
*
* Arguments:    details: pointer to detail type
*    num_orbs,num_atoms: integer
*              eigenset: eigenset_type
*  orbital_lookup_table: pointer to int
*
* Returns: none
*
* Action:   This is what does the work of transforming the matrix
*   of MO coefficients from the AO basis to the FMO basis.
*
*  The procedure is quite simple, the transformation matrix for the
*    FMOs is already constructed (in block form), so we just need
*    to multiply that by the matrix of coefficents in the AO basis (these
*    are stored in 'eigenset.  The new coefficients are stored in the
*    appropriate place within 'details->FMO_props.
*
*  The only complication is that the matrix multiplies are being
*    handled using sparse matrices.  This is because, as the number
*    of orbitals increases the tform matrices take up more and more space
*    and more and more of them is filled with zeroes (all terms between
*    fragments are zero).  A significant time and space savings can be
*    achieved using the sparse techniques.  This is documented more fully
*    below.
*
****************************************************************************/
void tform_wavefuncs_to_FMO_basis(details,num_orbs,num_atoms,eigenset,orbital_lookup_table)
  detail_type *details;
  int num_orbs;
  int num_atoms;
  eigenset_type eigenset;
  int *orbital_lookup_table;
{
  FMO_frag_type *FMO_frag;

  int num_orbs_this_frag;
  int frag,atom1,atom2,i,j,k;
  int itab,FMO_itab;
  int jtab,FMO_jtab;
  int ktab,FMO_ktab;
  int begin1,end1,FMO_begin1,FMO_end1;
  real *results_matR,*results_matI;
  real *T_matR,*T_matI;
  real accumR, accumI;


  /* start by getting pointers to and zeroing out the results matrices */
  results_matR = details->FMO_props->eigenset.vectR;
  results_matI = details->FMO_props->eigenset.vectI;
  bzero(results_matR,num_orbs*num_orbs*sizeof(real));
  bzero(results_matI,num_orbs*num_orbs*sizeof(real));



  for(i=0;i<num_orbs;i++){
    /* check to see if we need to move onto the next fragment */

    itab = i*num_orbs;

    /* now loop over columns (FMO's) */
    num_orbs_this_frag = 0;
    frag = 0;
    FMO_frag = &(details->FMO_frags[frag]);
    T_matR = FMO_frag->tform_matrix.matR;
    T_matI = FMO_frag->tform_matrix.matI;


    /* first loop over rows (MO's) in the new matrix */

    for(j=0;j<num_orbs;j++){
      jtab = j*num_orbs;

      if( num_orbs_this_frag == FMO_frag->num_orbs ){
        frag++;
        FMO_frag = &(details->FMO_frags[frag]);
        T_matR = FMO_frag->tform_matrix.matR;
        T_matI = FMO_frag->tform_matrix.matI;
        num_orbs_this_frag = 0;
      }

      FMO_jtab = num_orbs_this_frag;
      /***********

        to do the multiplication, loop over the atoms in the fragment and do the
        orbitals of each atom.

      ************/
      accumR = 0.0;
      accumI = 0.0;
      for(atom1=0; atom1<FMO_frag->num_atoms; atom1++){
        /******

          the location of the orbitals belonging to this atom must be found
          in BOTH the main and FMO coefficient matrices.

        *******/
        find_atoms_orbs(num_orbs,num_atoms,FMO_frag->atoms_in_frag[atom1],
                        orbital_lookup_table,&begin1,&end1);
        find_atoms_orbs(FMO_frag->num_orbs,FMO_frag->num_atoms,atom1,
                        FMO_frag->orbital_lookup_table,&FMO_begin1,&FMO_end1);

        if( begin1 >= 0 && FMO_begin1 >= 0 ){

          /* loop over the atom's orbitals */
          for(k=0;k<(end1-begin1);k++){
            ktab = k+begin1;
            FMO_ktab = (k+FMO_begin1)*FMO_frag->num_orbs;

            accumR += eigenset.vectR[ktab+itab]*T_matR[FMO_jtab + FMO_ktab];
            accumR -= eigenset.vectI[ktab+itab]*T_matI[FMO_jtab + FMO_ktab];
            accumI += eigenset.vectI[ktab+itab]*T_matR[FMO_jtab + FMO_ktab];
            accumI += eigenset.vectR[ktab+itab]*T_matI[FMO_jtab + FMO_ktab];
          }
        }

      }

      /* set the matrix element */
      results_matR[itab+j] = accumR;
      results_matI[itab+j] = accumI;
      num_orbs_this_frag++;
    }

  }
#ifdef DEBUG
fprintf(output_file,"\n\n\n\t FMO matrix \n");
printmat(results_matR,num_orbs,num_orbs,output_file,1e-5,0,details->line_width);
printmat(results_matI,num_orbs,num_orbs,output_file,1e-5,0,details->line_width);
fprintf(output_file,"\n\n\n");
#endif
}


/****************************************************************************
*
*                   Procedure tform_matrix_to_FMO_basis
*
* Arguments:    details: pointer to detail type
*    num_orbs,num_atoms: integer
*       AO_matR,AO_matI: pointers to reals
*   temp_matR,temp_matI: pointers to reals
*             cmplx_mat: complex_matrix_type
*  orbital_lookup_table: pointer to int
*
* Returns: none
*
* Action:   This is what does the work of transforming a generic matrix
*   from the AO basis to the FMO basis.
*
*   To do this, a similarity transformation is performed upon the AO basis
*    matrix ('AO_matR, and 'AO_matI), using the FMO coefficient matrices as
*    the transformation matrix.
*
*   If the matrix to be transformed is pure real, then 'AO_matI should just
*    be set to zero.
*
*   'AO_matR and 'AO_matI are assumed to be 'num_orbs x 'num_orbs.
*   'temp_matR and 'temp_matI, which are used to store intermediate results,
*    should also be at least 'num_orbs x 'num_orbs.
*
*    The results are placed in 'cmplx_mat.
*
****************************************************************************/
void tform_matrix_to_FMO_basis(details,num_orbs,num_atoms,AO_matR,AO_matI,
                               temp_matR,temp_matI,cmplx_mat,orbital_lookup_table)
  detail_type *details;
  int num_orbs;
  int num_atoms;
  real *AO_matR,*AO_matI,*temp_matR,*temp_matI;
  complex_matrix_type cmplx_mat;
  int *orbital_lookup_table;
{
  FMO_frag_type *FMO_frag;

  int num_orbs_this_frag;
  int frag,atom1,atom2,i,j,k;
  int itab,FMO_itab;
  int jtab,FMO_jtab;
  int ktab,FMO_ktab;
  int begin1,end1,FMO_begin1,FMO_end1;
  real *results_matR,*results_matI;
  real *T_matR,*T_matI;
  real accumR, accumI;

  /* start by getting pointers to and zeroing out the results matrix */
  results_matR = cmplx_mat.matR;
  results_matI = cmplx_mat.matI;
  bzero(results_matR,num_orbs*num_orbs*sizeof(real));
  bzero(results_matI,num_orbs*num_orbs*sizeof(real));

  /********

    Do the first matrix multiply and fill the temporary matrix

  *********/

  /* first loop over rows in the temporary matrix */
  for(i=0;i<num_orbs;i++){

    itab = i*num_orbs;

    frag = 0;
    FMO_frag = &(details->FMO_frags[frag]);
    T_matR = FMO_frag->eigenset.vectR;
    T_matI = FMO_frag->eigenset.vectI;
    num_orbs_this_frag = 0;

    /* now loop over columns */
    for(j=0;j<num_orbs;j++){

      /* check to see if we need to move onto the next fragment */
      if( num_orbs_this_frag == FMO_frag->num_orbs ){
        frag++;
        FMO_frag = &(details->FMO_frags[frag]);
        T_matR = FMO_frag->eigenset.vectR;
        T_matI = FMO_frag->eigenset.vectI;
        num_orbs_this_frag = 0;
      }

      jtab = j*num_orbs;
      FMO_jtab = num_orbs_this_frag*FMO_frag->num_orbs;

      /***********

        to do the multiplication, loop over the atoms in the fragment and do the
        orbitals of each atom.

      ************/
      accumR = 0.0;
      for(atom1=0; atom1<FMO_frag->num_atoms; atom1++){
        /******

          the location of the orbitals belonging to this atom must be found
          in BOTH the main and FMO coefficient matrices.

        *******/
        find_atoms_orbs(num_orbs,num_atoms,FMO_frag->atoms_in_frag[atom1],
                        orbital_lookup_table,&begin1,&end1);
        find_atoms_orbs(FMO_frag->num_orbs,FMO_frag->num_atoms,atom1,
                        FMO_frag->orbital_lookup_table,&FMO_begin1,&FMO_end1);

        if( begin1 >= 0 && FMO_begin1 >= 0 ){
          /* loop over the atom's orbitals */
          for(k=0;k<(end1-begin1);k++){
            FMO_ktab = k+FMO_begin1;
            ktab = (k+begin1);

            accumR += AO_matR[itab + ktab]*T_matR[FMO_ktab+FMO_jtab];
          }
        }
      }

      num_orbs_this_frag++;
      /* set the matrix element */
      temp_matR[itab+j] = accumR;
    }
  }

  /***********

    okay, now the temporary matrix is built, do the second matrix multiplication

  ************/
  frag = 0;
  FMO_frag = &(details->FMO_frags[frag]);
  T_matR = FMO_frag->eigenset.vectR;
  T_matI = FMO_frag->eigenset.vectR;
  num_orbs_this_frag = 0;

  /* first loop over rows (MO's) in the new matrix */
  for(i=0;i<num_orbs;i++){
    /* check to see if we need to move onto the next fragment */
    if( num_orbs_this_frag == FMO_frag->num_orbs ){
      frag++;
      FMO_frag = &(details->FMO_frags[frag]);
      T_matR = FMO_frag->tform_matrix.matR;
      T_matI = FMO_frag->tform_matrix.matI;
      num_orbs_this_frag = 0;
    }

    itab = i*num_orbs;

    FMO_itab = num_orbs_this_frag*FMO_frag->num_orbs;

    /* now loop over columns (FMO's) */
    for(j=0;j<num_orbs;j++){
      jtab = j;
      /***********

        to do the multiplication, loop over the atoms in the fragment and do the
        orbitals of each atom.

      ************/
      accumR = 0.0;
      for(atom1=0; atom1<FMO_frag->num_atoms; atom1++){
        /******

          the location of the orbitals belonging to this atom must be found
          in BOTH the main and FMO coefficient matrices.

        *******/
        find_atoms_orbs(num_orbs,num_atoms,FMO_frag->atoms_in_frag[atom1],
                        orbital_lookup_table,&begin1,&end1);
        find_atoms_orbs(FMO_frag->num_orbs,FMO_frag->num_atoms,atom1,
                        FMO_frag->orbital_lookup_table,&FMO_begin1,&FMO_end1);

        if( begin1 >= 0 && FMO_begin1 >= 0 ){
          /* loop over the atom's orbitals */
          for(k=0;k<(end1-begin1);k++){
            FMO_ktab = k+FMO_begin1;

            ktab = (k+begin1)*num_orbs;
            accumR += T_matR[FMO_itab + FMO_ktab]*temp_matR[ktab+jtab];
          }
        }
      }

      /* set the matrix element */
      results_matR[itab+j] = accumR;
    }
    num_orbs_this_frag++;
  }

#ifdef DEBUG
fprintf(output_file,"\n\n\n\t FMO matrix \n");
printmat(results_matR,num_orbs,num_orbs,output_file,1e-5,0,details->line_width);
fprintf(output_file,"\n\n\n");
#endif
}


/****************************************************************************
*
*                   Procedure tform_hermetian_matrix_to_FMO_basis
*
* Arguments:    details: pointer to detail type
*    num_orbs,num_atoms: integer
*              herm_mat: hermetian_matrix_type
*   temp_matR,temp_matI: pointers to reals
*               results: hermetian_matrix_type
*  orbital_lookup_table: pointer to int
*
* Returns: none
*
* Action:   This is what does the work of transforming a hermetian matrix
*   from the AO basis to the FMO basis.
*
*   To do this, a similarity transformation is performed upon the AO basis
*    matrix ('herm_mat), using the FMO coefficient matrices as
*    the transformation matrix.
*
*   'AO_matR and 'AO_matI are assumed to be 'num_orbs x 'num_orbs.
*   'temp_matR and 'temp_matI, which are used to store intermediate results,
*    should also be at least 'num_orbs x 'num_orbs.
*
*   The results of the multiplication are put into the matrix 'results
*
****************************************************************************/
void tform_hermetian_matrix_to_FMO_basis(details,num_orbs,num_atoms,herm_mat,
                               temp_matR,temp_matI,results,orbital_lookup_table)
  detail_type *details;
  int num_orbs;
  int num_atoms;
  hermetian_matrix_type herm_mat;
  real *temp_matR,*temp_matI;
  hermetian_matrix_type results;
  int *orbital_lookup_table;
{
  FMO_frag_type *FMO_frag;

  int num_orbs_this_frag;
  int frag,atom1,atom2,i,j,k;
  int itab,FMO_itab;
  int jtab,FMO_jtab;
  int ktab,FMO_ktab;
  int begin1,end1,FMO_begin1,FMO_end1;
  real *results_mat;
  real *T_matR,*T_matI;
  real accumR, accumI;


  /* start by zeroing out the results matrix */
  bzero(results.mat,num_orbs*num_orbs*sizeof(real));


  /********

    Do the first matrix multiply and fill the temporary matrix

  *********/

  /* first loop over rows in the temporary matrix */
  for(i=0;i<num_orbs;i++){

    itab = i*num_orbs;

    frag = 0;
    FMO_frag = &(details->FMO_frags[frag]);
    T_matR = FMO_frag->eigenset.vectR;
    T_matI = FMO_frag->eigenset.vectI;
    num_orbs_this_frag = 0;

    /* now loop over columns */
    for(j=0;j<num_orbs;j++){

      /* check to see if we need to move onto the next fragment */
      if( num_orbs_this_frag == FMO_frag->num_orbs ){
        frag++;
        FMO_frag = &(details->FMO_frags[frag]);
        T_matR = FMO_frag->eigenset.vectR;
        T_matI = FMO_frag->eigenset.vectI;
        num_orbs_this_frag = 0;
      }

      jtab = j*num_orbs;
      FMO_jtab = num_orbs_this_frag*FMO_frag->num_orbs;

      /***********

        to do the multiplication, loop over the atoms in the fragment and do the
        orbitals of each atom.

      ************/
      accumR = 0.0;
      accumI = 0.0;
      for(atom1=0; atom1<FMO_frag->num_atoms; atom1++){
        /******

          the location of the orbitals belonging to this atom must be found
          in BOTH the main and FMO coefficient matrices.

        *******/
        find_atoms_orbs(num_orbs,num_atoms,FMO_frag->atoms_in_frag[atom1],
                        orbital_lookup_table,&begin1,&end1);
        find_atoms_orbs(FMO_frag->num_orbs,FMO_frag->num_atoms,atom1,
                        FMO_frag->orbital_lookup_table,&FMO_begin1,&FMO_end1);

        if( begin1 >= 0 && FMO_begin1 >= 0 ){

          /* loop over the atom's orbitals */
          for(k=0;k<(end1-begin1);k++){
            FMO_ktab = k+FMO_begin1;
            ktab = (k+begin1);

            accumR += S_ELEMENT_R(herm_mat.mat,num_orbs,i,ktab)*
              T_matR[FMO_ktab+FMO_jtab];
            accumR += S_ELEMENT_I(herm_mat.mat,num_orbs,i,ktab)*
              T_matI[FMO_ktab+FMO_jtab];
            accumI += S_ELEMENT_I(herm_mat.mat,num_orbs,i,ktab)*
              T_matR[FMO_ktab+FMO_jtab];
            accumI -= S_ELEMENT_R(herm_mat.mat,num_orbs,i,ktab)*
              T_matI[FMO_ktab+FMO_jtab];
          }
        }
      }

      num_orbs_this_frag++;
      /* set the matrix element */
      temp_matR[itab+j] = accumR;
      temp_matI[itab+j] = accumI;
    }
  }

#ifdef DEBUG
fprintf(output_file,"\n\n\n\t FMO temporary matrix REAL part\n");
printmat(temp_matR,num_orbs,num_orbs,output_file,1e-5,0,details->line_width);
printmat(temp_matI,num_orbs,num_orbs,output_file,1e-5,0,details->line_width);
fprintf(output_file,"\n\n\n");
#endif

  /***********

    okay, now the temporary matrix is built, do the second matrix multiplication

  ************/
  frag = 0;
  FMO_frag = &(details->FMO_frags[frag]);
  T_matR = FMO_frag->eigenset.vectR;
  T_matI = FMO_frag->eigenset.vectI;
  num_orbs_this_frag = 0;

  /* first loop over rows (MO's) in the new matrix */
  for(i=0;i<num_orbs;i++){
    /* check to see if we need to move onto the next fragment */
    if( num_orbs_this_frag == FMO_frag->num_orbs ){
      frag++;
      FMO_frag = &(details->FMO_frags[frag]);
      T_matR = FMO_frag->eigenset.vectR;
      T_matI = FMO_frag->eigenset.vectI;
      num_orbs_this_frag = 0;
    }

    itab = i*num_orbs;

    FMO_itab = num_orbs_this_frag*FMO_frag->num_orbs;

    /* now loop over columns (FMO's) */
    for(j=i;j<num_orbs;j++){
      jtab = j;
      /***********

        to do the multiplication, loop over the atoms in the fragment and do the
        orbitals of each atom.

      ************/
      accumR = 0.0;
      accumI = 0.0;
      for(atom1=0; atom1<FMO_frag->num_atoms; atom1++){
        /******

          the location of the orbitals belonging to this atom must be found
          in BOTH the main and FMO coefficient matrices.

        *******/
        find_atoms_orbs(num_orbs,num_atoms,FMO_frag->atoms_in_frag[atom1],
                        orbital_lookup_table,&begin1,&end1);
        find_atoms_orbs(FMO_frag->num_orbs,FMO_frag->num_atoms,atom1,
                        FMO_frag->orbital_lookup_table,&FMO_begin1,&FMO_end1);

        if( begin1 >= 0 && FMO_begin1 >= 0 ){

          /* loop over the atom's orbitals */
          for(k=0;k<(end1-begin1);k++){
            FMO_ktab = k+FMO_begin1;
            ktab = (k+begin1)*num_orbs;
            accumR += T_matR[FMO_itab + FMO_ktab]*temp_matR[ktab+jtab];
            accumR -= T_matI[FMO_itab + FMO_ktab]*temp_matI[ktab+jtab];
            accumI += T_matR[FMO_itab + FMO_ktab]*temp_matI[ktab+jtab];
            accumI += T_matI[FMO_itab + FMO_ktab]*temp_matR[ktab+jtab];
          }
        }
      }

      /* set the matrix element */
      results.mat[itab+j] = accumR;
      if( j != i ){
        results.mat[j*num_orbs+i] = accumI;
      }
    }
    num_orbs_this_frag++;
  }
#ifdef DEBUG
fprintf(output_file,"\n\n\n\t FMO overlap matrix \n");
printmat(results.mat,num_orbs,num_orbs,output_file,1e-5,details->line_width);
fprintf(output_file,"\n\n\n");
#endif
}




