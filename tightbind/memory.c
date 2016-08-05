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
*     This file has the routines to deal with allocating the major arrays
*      used by the program.
*
*  created:  greg landrum  August 1993
*
*****************************************************************************/
#include "bind.h"




long tot_usage = 0;

/****************************************************************************
*
*                   Procedure my_malloc
*
* Arguments: size: long
*
* Returns: pointer to int
*
* Action: calls malloc and dumps some info into the status
*    file if the malloc fails
*
*****************************************************************************/
int *my_malloc(long size)
{
  int *tptr;

  tptr = (int *)malloc(size);
  if( !tptr ){
    fprintf(stderr,"Malloc failed getting %6.2f K.  Total allocated: %6.2f Meg\n",
            (float)size/1024,(float)tot_usage/(1024*1024));
    fprintf(status_file,
            "Malloc failed getting %6.2f K.  Total allocated: %6.2f Meg\n",
            (float)size/1024,(float)tot_usage/(1024*1024));
  }else{
    tot_usage += size;
  }
  return tptr;
}

/****************************************************************************
*
*                   Procedure my_calloc
*
* Arguments: num: int
*           size: int
*
* Returns: pointer to int
*
* Action: calls calloc and dumps some info into the status
*    file if the calloc fails
*
*****************************************************************************/
int *my_calloc(int num, int size)
{
  int *tptr;

  tptr = (int *)calloc(num,size);
  if( !tptr ){
    fprintf(stderr,"Calloc failed getting %6.2f K.  Total allocated: %6.2f Meg\n",
            (float)(num*size)/1024,(float)tot_usage/(1024*1024));
    fprintf(status_file,
            "Calloc failed getting %6.2f K.  Total allocated: %6.2f Meg\n",
            (float)(num*size)/1024,(float)tot_usage/(1024*1024));
  }else{
    tot_usage += num*size;
  }
  return tptr;
}

/****************************************************************************
*
*                   Procedure my_realloc
*
* Arguments: ptr: int *
*           size: int
*
* Returns: pointer to int
*
* Action: reallocates ptr.  dumps information into the status file
*    file if we fail
*
*****************************************************************************/
int *my_realloc(int *ptr, int size)
{
  int *tptr;

  tptr = (int *)calloc(1,size);
  if( !tptr ){
    fprintf(stderr,"Calloc failed getting %6.2f K in realloc.  Total allocated: %6.2f Meg\n",
            (float)(size)/1024,(float)tot_usage/(1024*1024));
    fprintf(status_file,
            "Calloc failed getting %6.2f K in realloc.  Total allocated: %6.2f Meg\n",
            (float)(size)/1024,(float)tot_usage/(1024*1024));
  }else{
    tot_usage += size;
  }

  /* we've got the new memory, copy in the old data */
  bcopy(ptr,tptr,size);
  /* free up the old pointer */
  free(ptr);
  return tptr;
}




/****************************************************************************
*
*                   Procedure allocate_matrices
*
* Arguments:         cell: pointer to cell type
*                 details: pointer to detail type
*                 H_R,S_R: pointer to hermetian_matrix_type
*                 H_K,S_K: pointer to hermetian_matrix_type
*           cmplx_hamil, cmplx_overlap: pointer to pointer to complex
*                eigenset: pointer to eigenset_type
*             work1,work2: pointers to pointers to reals
*                   work3: pointers to pointers to reals
*                   cmplx_work: pointer to pointer to complex
*              properties: pointer to prop_type
*           avg_prop_info: pointer to pointer to avg_prop_info_type
*                num_orbs: int
*            tot_overlaps: pointer to int
*    orbital_lookup_table: pointer to int
*        orbital_ordering: pointer to pointer to K_orb_ptr_type
*
*
* Returns: none
*
* Action:
*      This allocates the memory necessary for the Hamiltonian and overlap
*      matrices.  It also sets the variables indicating the dimensions
*      of these arrays.
*
*    The amount of memory used depends on the mode in which the program is
*     running.
*
* For Fat mode:
*    The overlap matrix in R space is num_orbs*num_orbs*tot_overlaps in size
*     where num_orbs is the total number of orbitals in the unit cell, and
*     tot_overlaps is the total number of symmetry distinct (in the general
*     case) overlaps possible.
*
*    The Hamiltonian matrix in R space is only num_orbs*num_orbs in
*      size.
*
*
* For Thin mode:
*    The overlap matrix is only num_orbs*num_orbs in size (only one is stored
*     at a time).
*    The Hamiltonian matrix in R space is num_orbs long (only diagonal elements).
*
*
* In both cases:
*    Each matrix in K space is num_orbs*num_orbs in size
*
*
*
*  In addition to the hamiltonian and overlap matrices, space is allocated
*   for the eigenvectors (num_orbs^2) and values(num_orbs), as well as to
*   provide working room for the diagonalization routines (2 num_orbs arrays),
*   the work3 array will be used to store the various matrices along the way,
*   it's num_orbs^2 in size.
*
*  For extended system calculations in FAT mode memory is also nabbed to
*   store the wavefunctions and overlap matrices in integer format as well
*   as the space needed for the energies at each kpoint.
*    The total memory usage here is 3*(num_orbs^2)*number of k points integers
*     and num_orbs*number of k points reals.
*
*  Space is also allocated for the overlap population matrix (num_orbs X num_orbs)
*   and the net charges array (cell->num_atoms)
*
*  Due to the fact that the matrices are hermitian, the real elements can be
*       stored in one triangle and the imaginary in the other.
*
*****************************************************************************/
void allocate_matrices(cell,details,H_R,S_R,
                       H_K,S_K,cmplx_hamil,cmplx_overlap,eigenset,work1,work2,
                       work3,cmplx_work,properties,avg_prop_info,num_orbs,tot_overlaps,
                       orbital_lookup_table,orbital_ordering)
  cell_type *cell;
  detail_type *details;
  hermetian_matrix_type *H_R,*S_R;
  hermetian_matrix_type *H_K,*S_K;
  complex **cmplx_hamil,**cmplx_overlap;
  eigenset_type *eigenset;
  real **work1,**work2,**work3;
  complex **cmplx_work;
  prop_type *properties;
  avg_prop_info_type **avg_prop_info;
  int num_orbs;
  int *tot_overlaps;
  int *orbital_lookup_table;
  K_orb_ptr_type **orbital_ordering;
{

  FMO_frag_type *FMO_frag;
  int i,j;
  int begin,end;
  int num_so_far;
  long mem_for_avg_props;
  long mem_per_overlapR,mem_per_hamR;
  long mem_per_overlapK,mem_per_hamK;
  real estimated_usage;
  long tot_usage=0;
  real *temp_mat;
  int num_frags;

  /* set the dimensionalities of the various matrices */
  H_R->dim = H_K->dim = S_R->dim = S_K->dim = eigenset->dim = num_orbs;

  /******

    figure out how many orbitals there are in each FMO fragment

  *******/
  if( (details->num_FMO_frags != 0 || details->num_FCO_frags != 0) &&
     details->FMO_frags ){
    if( !details->num_FMO_frags ) num_frags = details->num_FCO_frags;
    else num_frags = details->num_FMO_frags;
    fprintf(status_file,"There are %d fragments.\n",num_frags);
    for(i=0;i<num_frags;i++){
      FMO_frag = &(details->FMO_frags[i]);
      FMO_frag->num_orbs = 0;

      FMO_frag->orbital_lookup_table = (int *)my_calloc(FMO_frag->num_atoms,
                                                     sizeof(int));
      if( !FMO_frag->orbital_lookup_table )
        fatal("Can't get space for an FMO orbital_lookup_table.");

      for( j=0;j<FMO_frag->num_atoms;j++ ){
        find_atoms_orbs(num_orbs,cell->num_atoms,FMO_frag->atoms_in_frag[j],
                        orbital_lookup_table,&begin,&end);
        if( begin >= 0 ){
          FMO_frag->orbital_lookup_table[j] = FMO_frag->num_orbs;
          FMO_frag->num_orbs += end - begin;
        }else{
          FMO_frag->orbital_lookup_table[j] = -1;
        }
      }
      fprintf(status_file,"Fragment %d has %d atoms and %d orbitals.\n",
              i+1,FMO_frag->num_atoms,FMO_frag->num_orbs);

      /* set the dimensionalities of the matrices needed */
      FMO_frag->eigenset.dim = FMO_frag->overlap_R.dim =
        FMO_frag->overlap_K.dim = FMO_frag->hamil_R.dim =
          FMO_frag->hamil_K.dim = FMO_frag->num_orbs;
    }
  }


  /* here's some stuff that is specific to extended systems */
  if( cell->dim > 0){
    /* figure out how many distinct overlaps there will be */
    *tot_overlaps = (2*cell->overlaps[0]+1)*(2*cell->overlaps[1]+1)*
      cell->overlaps[2] + cell->overlaps[0] + 1 +
        (2*cell->overlaps[0]+1)*cell->overlaps[1];

    fprintf(status_file,"%d overlaps will be considered.\n",*tot_overlaps);

    /* figure out about how much memory is gonna be needed */
    if( details->Execution_Mode == FAT ){

      if( !details->num_KPOINTS || *tot_overlaps <= details->num_KPOINTS
         || details->the_COOPS
         || details->num_FMO_frags
         || details->num_FCO_frags
         || details->band_info ){
        mem_per_overlapR = num_orbs*(num_orbs)*(*tot_overlaps);
        mem_per_overlapK = num_orbs*(num_orbs);
        details->store_R_overlaps = 1;
      } else{
        mem_per_overlapR = num_orbs*(num_orbs);
        mem_per_overlapK = num_orbs*(num_orbs)*details->num_KPOINTS;
        details->store_R_overlaps = 0;
      }
#if 0
      mem_per_overlapR = num_orbs*(num_orbs)*(*tot_overlaps);
      mem_per_overlapK = num_orbs*(num_orbs);
      details->store_R_overlaps = 1;
#endif
      mem_per_hamR = num_orbs*(num_orbs);
      mem_per_hamK = num_orbs*(num_orbs);
    }
    else if( details->Execution_Mode == THIN){
      mem_per_overlapR = num_orbs*(num_orbs);
      mem_per_overlapK = num_orbs*(num_orbs);
      mem_per_hamR = num_orbs;
      mem_per_hamK = num_orbs*(num_orbs);
    }
  }
  else{
    mem_per_overlapR = num_orbs*(num_orbs);
    mem_per_hamR = num_orbs*(num_orbs);
    mem_per_overlapK = 0;
    mem_per_hamK =0;
  }

  /* this is the total amount needed for the average properties */
  mem_for_avg_props = num_orbs*(num_orbs)*details->num_KPOINTS*3 +
    (num_orbs)*details->num_KPOINTS;

  /* this is an _approximate_ measure of the amount of memory required */
  estimated_usage =
    mem_per_hamR + mem_per_hamK + mem_per_overlapR + mem_per_overlapK
      + 3*(num_orbs)*(num_orbs)+3*(num_orbs) + cell->num_atoms;


  /******

    it's gonna take about twice as much memory if we're doing FMO

    This is an *extremely* approximate measure.

  ****/
  if( details->num_FMO_frags || details->num_FCO_frags){
    mem_for_avg_props *= 2;
    estimated_usage += mem_per_hamK + mem_per_overlapK + num_orbs;
  }

  /***********
    Since the total amount of memory for the avg_props array will
    eventually be needed (to actually calculate the properties),
    we need to make sure that there will actually be enough space at
    that point to do so.

    If estimated_usage > (mem_for_avg_props + num_orbs^2) then we'll be
    fine.  Otherwise we need to check to see if we'll be able to
    get the memory.
  *************/
  if( details->Execution_Mode == THIN && details->avg_props ){

    if( estimated_usage < (mem_for_avg_props + (num_orbs*(num_orbs))) ){
      fprintf(status_file,"Checking to see if the total amount of memory needed\
is present.\n");

      /* try and allocate the memory required. */
      temp_mat = (real *)my_malloc((mem_for_avg_props + (num_orbs*(num_orbs)))*
                                 sizeof(real));
      if( !temp_mat ){
        fprintf(status_file,"Whoops! not enough memory.\n");
        fatal("Can't allocate required memory.\n");
      }

      /* We're safe, free that temporary array */
      free(temp_mat);
    }
  }
  else{
    if( details->avg_props ){
      estimated_usage += mem_for_avg_props;
    }
  }

  estimated_usage *= sizeof(real); /* convert to bytes */
  estimated_usage /= 1024; /* then to Kbytes */

  fprintf(status_file,"Allocating approximately %8.2lf Kbytes (%8.2lf Meg).\n",
          estimated_usage,estimated_usage/1024);

  fprintf(status_file,"Each matrix is: %8.2lf meg.\n",
          (real)mem_per_hamK*sizeof(real)/(real)(1024*1024));

  if( !details->just_geom ){
    /*********

      now actually get the memory
      (and make sure that we got it)

      ***********/
    H_R->mat = (real *)my_malloc(mem_per_hamR*sizeof(real));
    S_R->mat = (real *)my_malloc(mem_per_overlapR*sizeof(real));
    if( !(S_R->mat) ){
      fatal("Can't allocate space for overlap or hamiltonian matrices.");
    }
    if( cell->dim != 0 ){
      H_K->mat = (real *)my_malloc(mem_per_hamK*sizeof(real));
      S_K->mat = (real *)my_malloc(mem_per_overlapK*sizeof(real));
      if( !(S_K->mat) ){
        fatal("Can't allocate space for K space overlap or hamiltonian matrices.");
      }
    }
    if( !details->just_matrices ){
#ifdef USE_LAPACK
      /* allocate storage for the complex arrays used by the LAPACK routines */
      *cmplx_hamil = (complex *)my_malloc(num_orbs*num_orbs*sizeof(complex));
      *cmplx_overlap = (complex *)my_malloc(num_orbs*num_orbs*sizeof(complex));
      if( !(*cmplx_hamil) || !(*cmplx_overlap) ){
        fatal("Can't allocate complex matrices");
      }
#endif

      eigenset->vectR = (real *)my_calloc(num_orbs*(num_orbs),sizeof(real));
      eigenset->vectI = (real *)my_calloc(num_orbs*(num_orbs),sizeof(real));
      eigenset->val = (real *)my_calloc(num_orbs,sizeof(real));
      if( !(eigenset->vectR) || !(eigenset->val) )
        fatal("Can't allocate space for the eigenset storage.");

    /* only get space for mulliken population analysis if we need to */
      if( details->OP_mat_PRT || details->ROP_mat_PRT || details->net_chg_PRT
          || details->vary_zeta || details->avg_props){
        properties->OP_mat = (real *)my_calloc(num_orbs*(num_orbs),sizeof(real));
        properties->net_chgs = (real *)my_calloc(cell->num_atoms,sizeof(real));
        if( !(properties->net_chgs) || !(properties->OP_mat) )
          fatal("Can't get space for OP matrix.");
      }
      if( details->ROP_mat_PRT || details->avg_props ){
        properties->ROP_mat = (real *)my_calloc(cell->num_atoms*cell->num_atoms,sizeof(real));
        if( !(properties->ROP_mat) )
          fatal("Can't get space for reduced OP matrix.");
      }
    /* only get space for modified mulliken population analysis if we need to */
      if( details->mod_OP_mat_PRT || details->mod_ROP_mat_PRT ||
         details->mod_net_chg_PRT){
        properties->mod_OP_mat = (real *)my_calloc(num_orbs*(num_orbs),sizeof(real));
        properties->mod_net_chgs = (real *)my_calloc(cell->num_atoms,sizeof(real));
        if( !(properties->mod_net_chgs) || !(properties->mod_OP_mat) )
          fatal("Can't get space for modified OP matrix.");
      }
      if( details->mod_ROP_mat_PRT ){
        properties->mod_ROP_mat = (real *)my_calloc(cell->num_atoms*cell->num_atoms,sizeof(real));
        if( !(properties->mod_ROP_mat) )
          fatal("Can't get space for reduced modified OP matrix.");
      }

      if( details->chg_mat_PRT || details->Rchg_mat_PRT || details->avg_props ){
        properties->chg_mat = (real *)my_calloc(num_orbs*(num_orbs),sizeof(real));
        if( !(properties->chg_mat) )
          fatal("Can't get space for the charge matrix.");
        if( details->Rchg_mat_PRT ){
          properties->Rchg_mat = (real *)my_calloc(cell->num_atoms*cell->num_atoms,sizeof(real));
          if( !(properties->Rchg_mat) )
            fatal("Can't get space for the reduced charge matrix.");
        }
      }

      /* the temporary storage arrays */
      *work1 = (real *)my_calloc(num_orbs,sizeof(real));
      *work2 = (real *)my_calloc(num_orbs,sizeof(real));
      *work3 = (real *)my_calloc(num_orbs*(num_orbs),sizeof(real));
      /* check to see if we got the memory */
      if( !(*work2) || !(*work3) )
        fatal("Can't allocate space for the matrices.");
#ifdef USE_LAPACK
      *cmplx_work = (complex *)my_malloc(num_orbs*num_orbs*sizeof(complex));
      if( !(*cmplx_work) )
        fatal("Can't get space for complex work array");
#endif
    }
    /********

      get space for the average properties data
      ( if we need it and aren't running in THIN mode )

      *********/
    if( details->Execution_Mode != THIN && details->avg_props &&
        !details->just_matrices){
      *avg_prop_info = (avg_prop_info_type *)my_calloc(details->num_KPOINTS,
                                                    sizeof(avg_prop_info_type));
      if( !(*avg_prop_info) )
        fatal("Can't get info for average properties info array.");

      /******
        loop through and get the memory for each element of the
        avg_prop_info array
        *******/
      for(i=0;i<details->num_KPOINTS;i++){
        if( !details->just_avgE ){
          (*avg_prop_info)[i].orbs = (float *)my_calloc(num_orbs*(num_orbs),sizeof(float));
          (*avg_prop_info)[i].orbsI = (float *)my_calloc(num_orbs*(num_orbs),sizeof(float));
#ifdef KEEP_OVERLAP_MATS
          (*avg_prop_info)[i].S = (float *)my_calloc(num_orbs*(num_orbs),sizeof(float));
#endif
          (*avg_prop_info)[i].chg_mat = (float *)my_calloc(num_orbs*(num_orbs),
                                                        sizeof(float));

          if( !(*avg_prop_info)[i].chg_mat ){
            fatal("Can't get memory for an array within the avg_prop_info array.");
          }
        }
        (*avg_prop_info)[i].energies = (float *)my_calloc(num_orbs,sizeof(float));
        if( !(*avg_prop_info)[i].energies){
          fatal("Can't get memory for avg_prop_info energies\n");
        }
        /* get memory for the FMO properties stuff (if we need them) */
        if( details->num_FMO_frags || details->num_FCO_frags){
          (*avg_prop_info)[i].FMO_orbs =
            (float *)my_calloc(num_orbs*(num_orbs),sizeof(float));
          (*avg_prop_info)[i].FMO_orbsI =
            (float *)my_calloc(num_orbs*(num_orbs),sizeof(float));
          (*avg_prop_info)[i].FMO_chg_mat =
            (float *)my_calloc(num_orbs*(num_orbs),
                            sizeof(float));
          if( !(*avg_prop_info)[i].FMO_chg_mat )
            fatal("Can't get memory for an FMO array within the avg_prop_info array.");
        }
      }

      *orbital_ordering = (K_orb_ptr_type *)my_calloc(num_orbs*
                                                   details->num_KPOINTS,
                                                   sizeof(K_orb_ptr_type));
      if(!(*orbital_ordering))
        fatal("Can't get memory for orbital ordering array.");

    }

    /********

      allocate the memory needed for the FMO analysis

    ********/
    if( details->num_FMO_frags || details->num_FCO_frags ){

      details->FMO_props = (FMO_prop_type *)my_calloc(1,sizeof(FMO_prop_type));
      if( !details->FMO_props ) fatal("Can't get memory for FMO_props.");

      details->FMO_props->eigenset.vectR = (real *)my_calloc((num_orbs)*(num_orbs),
                                                     sizeof(real));
      details->FMO_props->eigenset.vectI = (real *)my_calloc((num_orbs)*(num_orbs),
                                                     sizeof(real));
      if( !details->FMO_props->eigenset.vectI )
        fatal("Can't get memory for eigenvectors in FMO basis.");

      details->FMO_props->eigenset.dim = num_orbs;
      /*****
        we don't need to get memory to store energies... that's already
        taken care of elsewhere
      *****/

      details->FMO_props->chg_mat = (real *)my_calloc((num_orbs)*(num_orbs),
                                                    sizeof(real));
      if(!details->FMO_props->chg_mat)
        fatal("Can't get memory for charge matrix in FMO basis.");

      details->FMO_props->OP_mat = (real *)my_calloc((num_orbs)*(num_orbs),
                                                    sizeof(real));
      if(!details->FMO_props->OP_mat)
        fatal("Can't get memory for overlap population matrix in FMO basis.");

      details->FMO_props->ROP_mat = (real *)my_calloc((num_frags)*
                                                    (num_frags),
                                                    sizeof(real));
      if(!details->FMO_props->ROP_mat)
        fatal("Can't get memory for reduced overlap population matrix in FMO basis.");

      details->FMO_props->overlap.mat = (real *)my_calloc((num_orbs)*(num_orbs),
                                                         sizeof(real));
      if(!details->FMO_props->overlap.mat)
        fatal("Can't get memory for overlap matrix in FMO basis.");

      details->FMO_props->overlap.dim = num_orbs;

      details->FMO_props->hamil.mat = (real *)my_calloc((num_orbs)*(num_orbs),
                                                         sizeof(real));
      if(!details->FMO_props->hamil.mat)
        fatal("Can't get memory for hamiltonian matrix in FMO basis.");

      details->FMO_props->hamil.dim = num_orbs;


      details->FMO_props->net_chgs = (real *)my_calloc(num_frags,
                                                     sizeof(real));
      if( !details->FMO_props->net_chgs ) fatal("Can't get space for FMO charges.");

      for(i=0; i<num_frags; i++){
        FMO_frag = &(details->FMO_frags[i]);


        if( details->Execution_Mode == FAT ){
          mem_per_overlapR = FMO_frag->num_orbs*(FMO_frag->num_orbs);
          mem_per_overlapK = FMO_frag->num_orbs*(FMO_frag->num_orbs);
          mem_per_hamR = FMO_frag->num_orbs*(FMO_frag->num_orbs);
          mem_per_hamK = FMO_frag->num_orbs*(FMO_frag->num_orbs);
        }

        FMO_frag->hamil_R.mat =
          (real *)my_calloc(mem_per_hamR,sizeof(real));
        FMO_frag->overlap_R.mat =
          (real *)my_calloc(mem_per_overlapR,sizeof(real));

        if( !(FMO_frag->overlap_R.mat) )
          fatal("Can't allocate space for FMO overlap or hamiltonian matrices.");

        FMO_frag->tform_matrix.dim = FMO_frag->num_orbs;
        FMO_frag->tform_matrix.matR =
          (real *)my_calloc(FMO_frag->num_orbs*FMO_frag->num_orbs,sizeof(real));
        FMO_frag->tform_matrix.matI =
          (real *)my_calloc(FMO_frag->num_orbs*FMO_frag->num_orbs,sizeof(real));
        if( !FMO_frag->tform_matrix.matI )
          fatal("Can't get space for FMO transform matrix.");

        if( cell->dim != 0 ){
          FMO_frag->hamil_K.mat =
            (real *)my_calloc(mem_per_hamK,sizeof(real));
          FMO_frag->overlap_K.mat =
            (real *)my_calloc(mem_per_overlapK,sizeof(real));
        }
        else{
          FMO_frag->hamil_K.mat =
            (real *)my_calloc(mem_per_hamR,sizeof(real));
          FMO_frag->overlap_K.mat =
            (real *)my_calloc(mem_per_overlapR,sizeof(real));
        }

        if( !(FMO_frag->overlap_K.mat) )
          fatal("Can't allocate space for K space FMO overlap or hamiltonian matrices.");

        FMO_frag->eigenset.vectR = (real *)my_calloc(FMO_frag->num_orbs*FMO_frag->num_orbs,
                                                   sizeof(real));
        FMO_frag->eigenset.vectI = (real *)my_calloc(FMO_frag->num_orbs*FMO_frag->num_orbs,
                                                   sizeof(real));
        FMO_frag->eigenset.val = (real *)my_calloc(FMO_frag->num_orbs,sizeof(real));

        if( !(FMO_frag->eigenset.vectR) || !(FMO_frag->eigenset.val) )
          fatal("Can't allocate space for the FMO eigenset storage.");

      }
    }
    /* that's that, everything is set.... */
  }
}










