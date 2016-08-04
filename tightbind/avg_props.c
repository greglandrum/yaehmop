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

/***
  Edit History:

  March '98: WG
  - print f AO occupations [print_avg_occups()]

***/


/****************************************************************************
*
*   This is the stuff for doing average properties
*
*  created:  greg landrum  March 1994
*
*****************************************************************************/
#include "bind.h"



/****************************************************************************
 *
 *                   Procedure calc_avg_occups
 *
 * Arguments:  details: pointer to detail_type
 *             cell: pointer to unit_cell_type
 *         num_orbs: int
 * orbital_ordering: pointer to K_orb_ptr_type
 *    avg_prop_info: pointer to avg_prop_info
 *       properties: pointer to prop_type
 *        AO_occups: pointer to real
 *
 *
 * Returns: none
 *
 * Action: Figures out the average occupation of all the atomic orbitals, and the
 *  net charges on all the atoms, and dumps the information into the output file.
 *
 *  The AO occupations are returned in AO_occups.
 *
 ****************************************************************************/
void calc_avg_occups(details,cell,num_orbs,orbital_ordering,avg_prop_info,
                     properties,AO_occups)
  detail_type *details;
  cell_type *cell;
  int num_orbs;
  K_orb_ptr_type *orbital_ordering;
  avg_prop_info_type *avg_prop_info;
  prop_type *properties;
  real *AO_occups;
{
  int i,j,k;
  real tot_num_K=0.0;
  float *MO_ptr,*MOI_ptr,*chg_mat_ptr;
  real accum;
  real avg_E_accum,tot_K_weight;
  int kpoint,MO;
  int begin,end;
  real num_occup_bands;
  atom_type *atom;
  real norm_fact;
  real contrib,total_electrons;

  /* zero out the occupations array. */
  bzero(AO_occups,num_orbs*sizeof(real));

  /* loop over all the occupied crystal orbitals */
  total_electrons = 0.0;


  accum = 0.0;
  avg_E_accum = 0.0;

  tot_num_K = 0.0;
  i = 0;
  while( orbital_ordering[i].occup > .0001 ){
    /* some pointers to make things a little more efficient */
    kpoint = orbital_ordering[i].Kpoint;
    avg_E_accum += ((real)*(orbital_ordering[i].energy))*orbital_ordering[i].occup*
      details->K_POINTS[kpoint].weight;

    accum += orbital_ordering[i].occup * details->K_POINTS[kpoint].weight;

    /* now figure out the contributions to the occupations */
    if( !details->just_avgE ){
      contrib = 0.0;
      MO = orbital_ordering[i].MO;
      chg_mat_ptr = &(avg_prop_info[kpoint].chg_mat[MO*num_orbs]);
      for(j=0;j<num_orbs;j++){
        AO_occups[j] += (real)chg_mat_ptr[j] * orbital_ordering[i].occup *
          details->K_POINTS[kpoint].weight;
      }
    }
    i++;
  }



  tot_num_K = 0.0;
  num_occup_bands = 0.0;
  tot_K_weight = 0.0;
  for( i=0; i<details->num_KPOINTS; i++){

    tot_K_weight += details->K_POINTS[i].weight;
    tot_num_K += details->K_POINTS[i].weight*
      details->K_POINTS[i].num_filled_bands;
    num_occup_bands += details->K_POINTS[i].num_filled_bands;
  }


  /******
    okay, now divide all the AO_occups by the total number of k points to get
    the average value.
  ******/
  if( !details->just_avgE ){
    for(i=0;i<num_orbs;i++){
      AO_occups[i] = AO_occups[i]*num_occup_bands/(tot_num_K*details->num_KPOINTS);
    }
  }
#ifdef DEBUG
  fprintf(output_file,"Num Electrons: %lf\n",accum);

  fprintf(output_file,"tot_num_K: %lf tot_K_weight: %lf\n",
          tot_num_K,tot_K_weight);
#endif
  accum = accum*num_occup_bands/(tot_num_K*details->num_KPOINTS);

#ifdef DEBUG
  fprintf(output_file,"num_filled: ");
  for(i=0;i<num_orbs;i++) fprintf(output_file,"%d ",
                                  (int)details->K_POINTS[i].num_filled_bands);
  fprintf(output_file,"\n");
  fprintf(output_file,"Num Electrons (post): %lf\n",accum);
#endif

  /* figure out the average energy */
  avg_E_accum /= tot_K_weight;
  properties->total_E = avg_E_accum;


}

/****************************************************************************
 *
 *                   Procedure print_avg_occups
 *
 * Arguments:  details: pointer to detail_type
 *             cell: pointer to unit_cell_type
 *         num_orbs: int
 * orbital_ordering: pointer to K_orb_ptr_type
 *    avg_prop_info: pointer to avg_prop_info
 *       properties: prop_type
 *        AO_occups: pointer to real
 *
 *
 * Returns: none
 *
 * Action: dumps the average occupation of all the atomic orbitals, and the
 *  net charges on all the atoms, to the output file.
 *
 ****************************************************************************/
void print_avg_occups(details,cell,num_orbs,orbital_ordering,avg_prop_info,
                     properties,AO_occups)
  detail_type *details;
  cell_type *cell;
  int num_orbs;
  K_orb_ptr_type *orbital_ordering;
  avg_prop_info_type *avg_prop_info;
  prop_type properties;
  real *AO_occups;
{
  int i,j,k,f_occup_print=0;
  real tot_num_K=0.0;
  float *MO_ptr,*MOI_ptr,*chg_mat_ptr;
  real accum;
  real avg_E_accum,tot_K_weight;
  int kpoint,MO;
  int begin,end;
  real num_occup_bands;
  atom_type *atom;
  real norm_fact;
  real contrib,total_electrons;

  total_electrons = 0.0;

  fprintf(output_file,"#Average Energy:  %lf eV\n\n",properties.total_E);
  if( !details->just_avgE ){
    fprintf(output_file,"# Atomic Orbital Occupations\n");
    fprintf(output_file,";        s      px      py      pz    dx2-y2    dz2     dxy     dxz     dyz\n");

    for( i=0; i<cell->num_atoms; i++){
      find_atoms_orbs(num_orbs,cell->num_atoms,i,orbital_lookup_table,&begin,&end);

      /*******
        if this isn't a dummy atom, then loop over its orbitals and print out their
        occupations
        *******/
      if( begin != -1 && end != -1 ){
        atom = &(cell->atoms[i]);
        j = begin;
        fprintf(output_file,"%d %s ",i+1,atom->symb);

        if( atom->ns != 0 ){
          fprintf(output_file,"% -6.4lf ",AO_occups[j]);
          j++;
        }
        else fprintf(output_file,"% -6.4lf ",0.0);
        if( atom->np != 0 ){
          fprintf(output_file,"% -6.4lf ",AO_occups[j++]);
          fprintf(output_file,"% -6.4lf ",AO_occups[j++]);
          fprintf(output_file,"% -6.4lf ",AO_occups[j++]);
        }
        else{
          fprintf(output_file,"% -6.4lf ",0.0);
          fprintf(output_file,"% -6.4lf ",0.0);
          fprintf(output_file,"% -6.4lf ",0.0);
        }
        if( atom->nd != 0 ){
          fprintf(output_file,"% -6.4lf ",AO_occups[j++]);
          fprintf(output_file,"% -6.4lf ",AO_occups[j++]);
          fprintf(output_file,"% -6.4lf ",AO_occups[j++]);
          fprintf(output_file,"% -6.4lf ",AO_occups[j++]);
          fprintf(output_file,"% -6.4lf ",AO_occups[j++]);
        }
        else{
          fprintf(output_file,"% -6.4lf ",0.0);
          fprintf(output_file,"% -6.4lf ",0.0);
          fprintf(output_file,"% -6.4lf ",0.0);
          fprintf(output_file,"% -6.4lf ",0.0);
          fprintf(output_file,"% -6.4lf ",0.0);
        }
        fprintf(output_file,"\n");

        /* figure out the net charge */
        properties.net_chgs[i] = atom->num_valence;
        for(j=begin;j<end;j++){
          properties.net_chgs[i] -= AO_occups[j];
          total_electrons += AO_occups[j];
        }
      }
    }

    /* now report f orbital occupations */

    for(i=0; i<cell->num_atoms; i++){
      if(cell->atoms[i].nf){
        f_occup_print =1;
      }
    }
    if( !details->just_avgE && f_occup_print ){
      fprintf(output_file,";       fz3    fxz2    fyz2    fxyz   fz(x2-y2)  fx(x2-3y2)  fy(3x2-y2)\n");

      /* loop over atoms and find atoms with f orbitals */

      for( i=0; i<cell->num_atoms; i++){
        if(cell->atoms[i].nf){
          find_atoms_orbs(num_orbs,cell->num_atoms,i,orbital_lookup_table,&begin,&end);

          /* check for dummy atom */

          if(begin != -1 && end != -1){
            atom = &(cell->atoms[i]);
            j=begin+BEGIN_F;
            fprintf(output_file,"%d %s ", i+1, atom->symb);
            fprintf(output_file,"% -6.4lf ",AO_occups[j++]);
            fprintf(output_file,"% -6.4lf ",AO_occups[j++]);
            fprintf(output_file,"% -6.4lf ",AO_occups[j++]);
            fprintf(output_file,"% -8.4lf ",AO_occups[j++]);
            fprintf(output_file,"% -11.4lf ",AO_occups[j++]);
            fprintf(output_file,"% -11.4lf ",AO_occups[j++]);
            fprintf(output_file,"% -12.4lf ",AO_occups[j++]);
            fprintf(output_file,"\n");

            /* no need to work out any more net atomic charges !!! */

          }
        }
      }
    }

    /* report the net charges */
    fprintf(output_file,"\n# AVERAGE NET CHARGES\n");

    for(i=0;i<cell->num_atoms;i++){
      fprintf(output_file,"%d %s %lf\n",i+1,cell->atoms[i].symb,properties.net_chgs[i]);

    }

    fprintf(output_file,"There are a total of %lf electrons.\n",total_electrons);
  }
}

/****************************************************************************
 *
 *                   Procedure calc_avg_FMO_occups
 *
 * Arguments:  details: pointer to detail_type
 *         num_orbs: int
 * orbital_ordering: pointer to K_orb_ptr_type
 *    avg_prop_info: pointer to avg_prop_info
 *        FMO_occups: pointer to real
 *
 *
 * Returns: none
 *
 * Action: Figures out the average occupation of all the fragment orbitals, and
 *  the net charges on all the framgents, and dumps the information into
 *  the output file.
 *
 *  The FMO occupations are returned in FMO_occups.
 *
 ****************************************************************************/
void calc_avg_FMO_occups(details,num_orbs,orbital_ordering,avg_prop_info,
                         FMO_occups)
  detail_type *details;
  int num_orbs;
  K_orb_ptr_type *orbital_ordering;
  avg_prop_info_type *avg_prop_info;
  real *FMO_occups;
{
  int i,j,k;
  real tot_num_K=0.0;
  float *MO_ptr,*MOI_ptr,*chg_mat_ptr;
  real accum;
  real avg_E_accum,tot_K_weight;
  int kpoint,MO;
  int begin,end;
  int tot_num_orbs;
  real num_occup_bands;
  atom_type *atom;
  real norm_fact;
  real contrib,total_electrons;

  /* error checking */
  if( !details->FMO_frags )
    FATAL_BUG("Null FMO_frags pointer in calc_avg_FMO_occups.");

  /* zero out the occupations array. */
  bzero(FMO_occups,num_orbs*sizeof(real));

  /* loop over all the occupied crystal orbitals */
  total_electrons = 0.0;

  accum = 0.0;
  avg_E_accum = 0.0;

  tot_num_K = 0.0;
  i = 0;
  while( orbital_ordering[i].occup > .0001 ){
    /* some pointers to make things a little more efficient */
    kpoint = orbital_ordering[i].Kpoint;
    MO = orbital_ordering[i].MO;
    chg_mat_ptr = &(avg_prop_info[kpoint].FMO_chg_mat[MO*num_orbs]);

    accum += orbital_ordering[i].occup * details->K_POINTS[kpoint].weight;

    /* now figure out the contributions to the occupations */
    contrib = 0.0;
    for(j=0;j<num_orbs;j++){
      FMO_occups[j] += (real)chg_mat_ptr[j] * orbital_ordering[i].occup *
        details->K_POINTS[kpoint].weight;
    }
    i++;
  }



  tot_num_K = 0.0;
  num_occup_bands = 0.0;
  tot_K_weight = 0.0;
  for( i=0; i<details->num_KPOINTS; i++){

    tot_K_weight += details->K_POINTS[i].weight;
    tot_num_K += details->K_POINTS[i].weight*
      details->K_POINTS[i].num_filled_bands;
    num_occup_bands += details->K_POINTS[i].num_filled_bands;
  }


  /******
    okay, now divide all the FMO_occups by the total number of k points to get
    the average value.
  ******/
  for(i=0;i<num_orbs;i++){
    FMO_occups[i] = FMO_occups[i]*num_occup_bands/(tot_num_K*details->num_KPOINTS);
  }
  accum = accum*num_occup_bands/(tot_num_K*details->num_KPOINTS);

  fprintf(output_file,"# Fragment MO Occupations\n");

  tot_num_orbs = 0;
  for( i=0; i<details->num_FMO_frags; i++){
    details->FMO_props->net_chgs[i] = 0.0;
    fprintf(output_file,";Fragment: %d\n",i+1);
    for(j=0;j<details->FMO_frags[i].num_orbs;j++){
      fprintf(output_file,"%d % -6.4lf\n",j+1,FMO_occups[tot_num_orbs]);
      details->FMO_props->net_chgs[i] += FMO_occups[tot_num_orbs];
      tot_num_orbs++;
    }
  }

  /* report the net charges */
  fprintf(output_file,"\n# AVERAGE FRAGMENT NET CHARGES (Based upon user supplied number of electrons)\n");

  for(i=0;i<details->num_FMO_frags;i++){
    fprintf(output_file,"%d % -6.4lf\n",i+1,
            details->FMO_frags[i].num_electrons-details->FMO_props->net_chgs[i]);

  }

}



/****************************************************************************
 *
 *                   Procedure calc_avg_OP
 *
 * Arguments:  details: pointer to detail_type
 *             cell: pointer to unit_cell_type
 *         num_orbs: int
 * orbital_ordering: pointer to K_orb_ptr_type
 *    avg_prop_info: pointer to avg_prop_info
 *         overlapR: hermetian_matrix_type
 *       properties: prop_type
 *
 *
 * Returns: none
 *
 * Action: Figures out the average overlap population and reduced overlap
 *   population matrices within the unit cell
 *
 *  The average overlap population and reduced avg OP matrices are
 *   stored in the appropriate places in 'properties.
 *
 ****************************************************************************/
void calc_avg_OP(details,cell,num_orbs,orbital_ordering,avg_prop_info,
                     overlapR,properties)
  detail_type *details;
  cell_type *cell;
  int num_orbs;
  K_orb_ptr_type *orbital_ordering;
  avg_prop_info_type *avg_prop_info;
  hermetian_matrix_type overlapR;
  prop_type properties;
{
  int i,j,k,l;
  real tot_num_K=0.0;
  float *MO_ptr,*MO_ptrI,*chg_mat_ptr;
  hermetian_matrix_type overlap;
  real accum,accumI;
  int kpoint,MO;
  atom_type *atom;
  real norm_fact;
  real contrib,total_electrons;
  int begin1,begin2,end1,end2;
  real num_electrons;
  int num_elements,num_occup_bands,num_occup_orbs;
  COOP_type *COOP_ptr;

  tot_num_K = 0.0;

  overlap.dim = num_orbs;

  /*******

    First determine the average OP and ROP matrices within the unit cell

  ********/

  /* zero them both out first */
  bzero(properties.OP_mat,num_orbs*num_orbs*sizeof(real));
  bzero(properties.ROP_mat,num_orbs*num_orbs*sizeof(real));

  overlap.mat = overlapR.mat;
  i = 0;
  while( orbital_ordering[i].occup > .0001 ){
    /* some pointers to make things a little more efficient */
    kpoint = orbital_ordering[i].Kpoint;
    MO = orbital_ordering[i].MO;

    /* get pointers to the orbital information */
    MO_ptr = &(avg_prop_info[kpoint].orbs[MO*num_orbs]);
    MO_ptrI = &(avg_prop_info[kpoint].orbsI[MO*num_orbs]);

    for(j=0;j<num_orbs;j++){
      for( k=j;k<num_orbs;k++){
        /*****
          within the unit cell there's no need to accumulate an imaginary
          contribution to the overlap population,
          *****/
        accum = ((real)MO_ptr[j] * MO_ptr[k] +
                 (real)MO_ptrI[j] * MO_ptrI[k]) /
                   ((real)MULTIPLIER*(real)MULTIPLIER);
        if(j != k){
          properties.OP_mat[j*num_orbs+k] +=
            (2.0*orbital_ordering[i].occup*details->K_POINTS[kpoint].weight*
              accum*HERMETIAN_R(overlap,j,k));

        }
        else{
          properties.OP_mat[j*num_orbs+k] +=
            (orbital_ordering[i].occup*details->K_POINTS[kpoint].weight*
              accum*HERMETIAN_R(overlap,j,k));
        }
      }
    }
    i++;
  }

  num_occup_bands = i;

  tot_num_K = 0.0;
  for( i=0; i<details->num_KPOINTS; i++){
    tot_num_K += details->K_POINTS[i].weight*details->K_POINTS[i].num_filled_bands;
  }

  /******
    okay, now divide all the elements of the
    overlap population matrix by the total number of k points to get
    the average value, and copy the elements across the diagonal.
  ******/
  num_electrons = 0.0;
  for(i=0;i<num_orbs;i++){
    for(j=i;j<num_orbs;j++){
      properties.OP_mat[i*num_orbs+j] = properties.OP_mat[i*num_orbs+j] *
        num_occup_bands / (tot_num_K*details->num_KPOINTS);
      properties.OP_mat[j*num_orbs+i] = properties.OP_mat[i*num_orbs+j];
      num_electrons += properties.OP_mat[i*num_orbs+j];
    }
  }

  /* figure out the average reduced overlap population matrix */
  num_elements = 0;
  for( i=0; i<cell->num_atoms; i++){
    find_atoms_orbs(num_orbs,cell->num_atoms,i,orbital_lookup_table,&begin1,&end1);

    /*******
      if this isn't a dummy atom, then loop over its orbitals and add
        up their contributions to the reduced overlap population.
    *******/
    if( begin1 != -1 && end1 != -1 ){

      /* do the cross terms first */
      for( j=0; j<i; j++){
        find_atoms_orbs(num_orbs,cell->num_atoms,j,orbital_lookup_table,&begin2,&end2);
        if( begin2 != -1 && end2 != -1 ){
          accum = 0.0;
          for( k=begin1; k<end1; k++ ){
            for( l=begin2; l<end2; l++){
              accum += properties.OP_mat[k*num_orbs+l];
            }
          }
          /*********
            store the element....
            treat the matrix as a symmetrical matrix (see notes.outl for the format)
          **********/
          properties.ROP_mat[num_elements++] = accum;
        }
      }

      /* now add up the diagonal elements which arise from this atom */
      accum = 0.0;
      for(j=begin1; j<end1; j++ ){
        accum += properties.OP_mat[j*num_orbs+j];
      }
      properties.ROP_mat[num_elements++] = accum;
    }
  }

  /* print out the matrices */
  if( details->avg_OP_mat_PRT ){
    fprintf(output_file,
            "\n; Average Mulliken Overlap Population Matrix W/in the Unit Cell.\n");
    printmat(properties.OP_mat,num_orbs,num_orbs,output_file,1e-05,0,details->line_width);
  }

  if( details->avg_ROP_mat_PRT ){
  fprintf(output_file,
          "\n; Average Reduced Mulliken Overlap Population Matrix W/in the Unit Cell.\n");
  print_sym_mat(properties.ROP_mat,cell->num_atoms,cell->num_atoms,output_file,
                (char *)0,(char *)0,details->line_width);
  }
  fprintf(output_file,"\n");


}


/****************************************************************************
 *
 *                   Procedure find_crystal_occupations
 *
 * Arguments:  details: pointer to detail_type
 *      electrons_per_cell: real
 *         num_orbs: int
 * orbital_ordering: pointer to K_orb_ptr_type
 *          Fermi_E: pointer to real
 *
 * Returns: none
 *
 * Action: this loops through the 'orbital_ordering array and populates
 *   the first 'electrons_per_cell * 'details->num_KPOINTS / 2.0 orbitals.
 *
 *  The Fermi Energy is stored in the the variable 'Fermi_E.
 *
 ****************************************************************************/
void find_crystal_occupations(details,electrons_per_cell,num_orbs,
                              orbital_ordering,Fermi_E)
  detail_type *details;
  real electrons_per_cell;
  int num_orbs;
  K_orb_ptr_type *orbital_ordering;
  real *Fermi_E;
{
  int i;
  int tot_orbs;
  int last_occup,begin_degen,end_degen,num_degen_levels;
  real num_degen_electrons,electrons_per_level;

  real electrons_left;
  real num_here;
  real accum;
  k_point_type *temp_kpoint;

  tot_orbs = num_orbs * details->num_KPOINTS;
  electrons_left = electrons_per_cell * (real)details->num_KPOINTS;

  /* zero out the orbital ordering array */
  for(i=0;i<tot_orbs;i++){
    orbital_ordering[i].occup = 0.0;
  }


  /* loop until we either run out of orbitals or electrons */
  accum = 0.0;
  for(i=0;i<tot_orbs && electrons_left > 0.0; i++){
    temp_kpoint = &(details->K_POINTS[orbital_ordering[i].Kpoint]);

    if( electrons_left >= 2.0 ){
      num_here = 2.0;

      orbital_ordering[i].occup = num_here;
      accum += num_here;

      electrons_left -= 2.0;
    }else if(electrons_left > 0.0 ){
      num_here = electrons_left;
      orbital_ordering[i].occup = num_here;
      electrons_left = 0.0;
      accum += num_here;
    }
  }

#ifdef DEBUG
  fprintf(output_file,"%%Tot num electrons in: %lf\n",accum);
#endif

  /*********

    check for degeneracies

  **********/
  last_occup = i-1;
  end_degen = i;
  while( end_degen < tot_orbs && fabs((real)*(orbital_ordering[last_occup].energy) -
              (real)*(orbital_ordering[end_degen].energy)) < DEGEN_TOL ){
    end_degen++;
  }

  begin_degen = last_occup - 1;
  while( begin_degen > -1 && fabs((real)*(orbital_ordering[last_occup].energy) -
              (real)*(orbital_ordering[begin_degen].energy)) < DEGEN_TOL){
    begin_degen--;
  }
  /* we went one step too far, so increment begin_degen */
  begin_degen++;

  num_degen_levels = end_degen-begin_degen;
  if( num_degen_levels > 1){
    fprintf(output_file,"; >>>>> The HOCO was found to be %d-fold degenerate.\n",
            num_degen_levels);

    /********
      now adjust the occupation numbers for the degenerate orbitals
        do this by finding the number of electrons in degenerate levels and
        then distributing them equally amongst those levels.
    ******/
    num_degen_electrons = 0;
    for(i=begin_degen;i<end_degen;i++){
      temp_kpoint = &(details->K_POINTS[orbital_ordering[i].Kpoint]);
      num_degen_electrons += orbital_ordering[i].occup;
    }
    electrons_per_level = num_degen_electrons / (real)num_degen_levels;
    for(i=begin_degen;i<end_degen;i++){
      temp_kpoint = &(details->K_POINTS[orbital_ordering[i].Kpoint]);
      orbital_ordering[i].occup = electrons_per_level;
    }
  }

  /******

    figure out the number of filled bands at each k point

  ******/
  /* first zero the array */
  for(i=0;i<details->num_KPOINTS;i++){
    details->K_POINTS[i].num_filled_bands = 0.0;
  }
  /* now fill it up */
  for(i=0;i<= last_occup;i++){
    temp_kpoint = &(details->K_POINTS[orbital_ordering[i].Kpoint]);
    temp_kpoint->num_filled_bands += orbital_ordering[i].occup / 2.0;
  }

  /* store the Fermi level */
  *Fermi_E = (real)*(orbital_ordering[last_occup].energy);
}





/****************************************************************************
 *
 *                   Procedure store_avg_prop_info
 *
 * Arguments:  details: pointer to detail_type
 *          which_k: int
 *         eigenset: eigenset_type
 *          overlap: hermetian_matrix_type
 *         num_orbs: int
 *          chg_mat: pointer to real
 *    avg_prop_info: pointer to avg_prop_info_type
 *
 *
 * Returns: none
 *
 * Action: This takes the wavefunctions and overlap and charge matrices
 *    and stores them either in the avg_prop_info element passed in or
 *    in the file indicated by the element (depending on execution mode)
 *
 ****************************************************************************/
void store_avg_prop_info(details,which_k,eigenset,overlap,num_orbs,
                         chg_mat,avg_prop_info)
  detail_type *details;
  int which_k;
  eigenset_type eigenset;
  hermetian_matrix_type overlap;
  int num_orbs;
  real *chg_mat;
  avg_prop_info_type *avg_prop_info;
{
  int i,j;
  int itab,jtab;
  avg_prop_info_type *temp_ptr;

  /* use temp_ptr to save us from TOO much pointer math */
  temp_ptr = &(avg_prop_info[which_k]);

  /* some error checking */
  if( (!details->just_avgE && (!temp_ptr->orbs || !temp_ptr->orbsI)) ||
     !temp_ptr->energies ){
    FATAL_BUG("Bogus avg_prop_info struct passed to store_avg_prop_info.");
  }

  /**********
    Loop over MO's, then AO's
  **********/

  for( i=0; i<num_orbs; i++){
    if( !details->just_avgE &&
       (details->num_proj_DOS || details->the_COOPS ||
        !details->no_total_DOS_PRT || details->num_FMO_frags) ){
      itab = i*num_orbs;
      for( j=0; j<num_orbs; j++){
        if( details->Execution_Mode != THIN ){
          temp_ptr->orbs[itab+j] = (float)EIGENVECT_R(eigenset,i,j);
          temp_ptr->orbsI[itab+j] = (float)EIGENVECT_I(eigenset,i,j);
          temp_ptr->chg_mat[itab+j] = (float)chg_mat[itab+j];

          /* do the FMO stuff if we need to */
          if( details->num_FMO_frags ){
            temp_ptr->FMO_orbs[itab+j] = (float)EIGENVECT_R(details->FMO_props->eigenset,i,j);
            temp_ptr->FMO_orbsI[itab+j] = (float)EIGENVECT_I(details->FMO_props->eigenset,i,j);
            temp_ptr->FMO_chg_mat[itab+j] = (float)details->FMO_props->chg_mat[itab+j];
          }
        }
        else{
          fatal("THIN mode properties calculations not implemented yet.");
        }
      }
    }
    temp_ptr->energies[i] = (float)EIGENVAL(eigenset,i);
  }
}


/*****

  this is the helper function used by qsort to sort the
    orbital energies.

*******/
int sort_energies_helper(orb1,orb2)
  const void *orb1, *orb2;
{
  real diff;

  diff = (real)*(((K_orb_ptr_type *)orb1)->energy) -
    (real)*(((K_orb_ptr_type *)orb2)->energy);
  if( diff > 0 ) return(1);
  else if( diff < 0 ) return(-1);
  else return( 0 );
}





/****************************************************************************
 *
 *                   Procedure sort_avg_prop_info
 *
 * Arguments:  details: pointer to detail_type
 *         num_orbs: int
 *    avg_prop_info: pointer to avg_prop_info_type
 * orbital_ordering: pointer to K_orb_ptr_type
 *
 *
 * Returns: none
 *
 * Action: This takes the details->num_KPOINTS long 'avg_prop_info
 *  array and "sorts" it in order of increasing energy level.
 *
 *  In order to make it easier to calculate the COOP (and whatever other
 *   average properties might be desired) later, the actual orbital ordering
 *   isn't rearranged.  Instead, an array of K_orb_ptr_type structures which
 *   point to their energies within the avg_prop_info array is sorted.
 *
 ****************************************************************************/
void sort_avg_prop_info(details,num_orbs,avg_prop_info,orbital_ordering)
  detail_type *details;
  int num_orbs;
  avg_prop_info_type *avg_prop_info;
  K_orb_ptr_type *orbital_ordering;
{
  int i,j;
  int itab;
  int num_so_far,num_elements;

  num_elements = details->num_KPOINTS * num_orbs;

  fprintf(status_file," Sorting %d crystal orbitals.\n",num_elements);

  /********

    fill the array

  *********/
  num_so_far = 0;

  /* loop over K points */
  for(i=0;i<details->num_KPOINTS;i++){
    itab = i;

    /* loop over MO's */
    for(j=0;j<num_orbs;j++){
      orbital_ordering[num_so_far].Kpoint = i;
      orbital_ordering[num_so_far].MO = j;
      if( details->Execution_Mode != THIN ){
        orbital_ordering[num_so_far].energy =
          &(avg_prop_info[i].energies[j]);
      }
      else{
        fatal("THIN mode properties calculations aren't implemented yet.");
      }
      num_so_far++;
    }
  }

  /********

    Now sort it using C's built in quick sort function.

  *********/
  qsort((void *)orbital_ordering,num_elements,sizeof(K_orb_ptr_type),
        sort_energies_helper);

  /* that's all we have to do here. */
}








