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
*   These are the things needed for doing a Density of States
*
*   (DOS here has nothing to do with MS-DOS, thank god)
*
*  created:  greg landrum  March 1994
*
*****************************************************************************/

/*******
  Recent edit history:

  10.07.1999 gL:
    started to add some Fermi surface functionality... this isn't
    finished (I don't have time) so it's commented out.

******/
#include "bind.h"


/****************************************************************************
 *
 *                   Procedure gen_total_DOS
 *
 * Arguments: details: pointer to detail_type
 *              cell: pointer to cell type
 *          num_orbs: int
 *     avg_prop_info: pointer to avg_prop_info_type
 *  orbital_ordering: pointer to K_orb_ptr_type
 *
 *
 * Returns: none
 *
 * Action:  Generates the total DOS and dumps it to the output file.
 *
 *   If a moments analysis is being done, this will calculate the moments
 *    up through 'num_moments;
 *      The basis of the moments code is from something written
 *        by Grigori Vagenine
 *
 *
 ****************************************************************************/
void gen_total_DOS(details,cell,num_orbs,avg_prop_info,orbital_ordering)
  detail_type *details;
  cell_type *cell;
  int num_orbs;
  avg_prop_info_type *avg_prop_info;
  K_orb_ptr_type *orbital_ordering;
{
  int i,j;
  int tot_num_orbs;
  real num_at_this_E;
  real diff;
  real this_E;
  real accum;
  real tot_K_weight;
  int num_entries=0;
  real *energies=0;
  real *weights=0;


  tot_K_weight = 0.0;
  for(i=0;i<details->num_KPOINTS;i++){
    tot_K_weight += details->K_POINTS[i].weight;
  }

  tot_num_orbs = num_orbs * details->num_KPOINTS;

  i=0;
  accum = 0.0;

  if( details->do_moments ){
    /* get memory to store the energies and weights */
    num_entries = 0;
    energies = (real *)calloc(tot_num_orbs,sizeof(real));
    weights = (real *)calloc(tot_num_orbs,sizeof(real));
    if(details->moments) free(details->moments);
    details->moments = (real *)calloc(details->num_moments,sizeof(real));
    if( !energies || !weights || !details->moments){
      error("Can't allocate memory for moments, not doing them.");
      if(details->moments) free(details->moments);
      details->moments = 0;
      if(energies)free(energies);
      energies = 0;
      if(weights)free(weights);
      weights=0;
    }
  }

  fprintf(output_file,"\n### TOTAL DENSITY OF STATES \n");
  fprintf(output_file,"%d states are present.\n",num_orbs*2);
  fprintf(output_file,"#BEGIN CURVE\n");
  while(i<tot_num_orbs){
    num_at_this_E = details->K_POINTS[orbital_ordering[i].Kpoint].weight;
    this_E = (real)*(orbital_ordering[i].energy);
    j=i+1;

    /*******
      just add up all of the contributions to the total DOS that fall
      within DOS_DEGEN_TOL of each other...
    ********/
    if( j < tot_num_orbs ){
      diff = (real)*(orbital_ordering[i].energy) - (real)*(orbital_ordering[j].energy);
      while(fabs(diff) < DOS_DEGEN_TOL && j<tot_num_orbs){
        num_at_this_E += details->K_POINTS[orbital_ordering[j].Kpoint].weight;
        j++;

        if( j < tot_num_orbs ){
          diff = (real)*(orbital_ordering[i].energy) - (real)*(orbital_ordering[j].energy);
        }
      }
    }

    i = j;
    /* write out the result */
#ifdef FS_TEST
    fprintf(output_file,"%lf %lf  ; %d\n",(real)num_at_this_E/tot_K_weight,
            this_E,orbital_ordering[i].Kpoint);
#else
    fprintf(output_file,"%lf %lf\n",(real)num_at_this_E/tot_K_weight,
            this_E);
#endif
    if( weights && energies ){
      weights[num_entries] = num_at_this_E/tot_K_weight;
      energies[num_entries] = this_E;
      num_entries++;
      if(num_entries >= tot_num_orbs) FATAL_BUG("num_entries too big in moments analysis.");
    }
  }
  fprintf(output_file,"#END CURVE\n");

  /****

    okay... if we're doing moments analysis, go ahead and do it.
    the 0th moment is calculated normally.
    higher moments are normalized by the 0th and referenced to the
    first.

  ***/
  if( weights && energies ){
    /* 0th moment */
    for(i=0;i<num_entries;i++)
      details->moments[0] += weights[i];

    /* 1st moment */
    for(i=0;i<num_entries;i++)
      details->moments[1] += weights[i]*energies[i];
    /* normalize it */
    details->moments[1] /= details->moments[0];

    /* higher moments */
    for(i=0;i<num_entries;i++){
      /* reference to the first moment */
      diff = energies[i]-details->moments[1];
      for(j=2;j<=details->num_moments;j++){
        details->moments[j] += weights[i]*pow(diff,(double)j);
      }
    }
    /* normalize them */
    for(j=2;j<=details->num_moments;j++){
      details->moments[j] /= details->moments[0];
    }

    /* we're done */
    free(weights);
    weights = 0;
    free(energies);
    energies = 0;
  }
}

/****************************************************************************
 *
 *                   Function orb_contribution
 *
 * Arguments: num_orbs: int
 *         kpoint,MO: ints
 *     avg_prop_info: pointer to avg_prop_info_type
 *          which_AO: int
 *
 * Returns: real
 *
 * Action:  Determines the contribution of 'which_AO to 'MO.
 *  This is just the charge matrix element for 'which_AO.
 *
 ****************************************************************************/
real orb_contribution(num_orbs,kpoint,MO,avg_prop_info,which_AO)
  int num_orbs,kpoint,MO;
  avg_prop_info_type *avg_prop_info;
  int which_AO;
{
  real contrib;

  /* just return the appropriate charge matrix element */
  contrib = (real)avg_prop_info[kpoint].chg_mat[MO*num_orbs+which_AO] /
    (real)MULTIPLIER;

  return(contrib);
}


/****************************************************************************
 *
 *                   Function FMO_contribution
 *
 * Arguments: num_orbs: int
 *         kpoint,MO: ints
 *     avg_prop_info: pointer to avg_prop_info_type
 *          which_FMO: int
 *
 * Returns: real
 *
 * Action:  Determines the contribution of 'which_FMO to 'MO.
 *  This is just the charge matrix element for 'which_FMO.
 *
 ****************************************************************************/
real FMO_contribution(num_orbs,kpoint,MO,avg_prop_info,which_FMO)
  int num_orbs,kpoint,MO;
  avg_prop_info_type *avg_prop_info;
  int which_FMO;
{
  real contrib;

  /* just return the appropriate charge matrix element */
  contrib = (real)avg_prop_info[kpoint].FMO_chg_mat[MO*num_orbs+which_FMO] /
    (real)MULTIPLIER;

  return(contrib);
}


/****************************************************************************
 *
 *                   Function atom_contribution
 *
 * Arguments: num_orbs: int
 *          num_atoms: int
 *         kpoint,MO: ints
 *     avg_prop_info: pointer to avg_prop_info_type
 *          which_atom: int
 * orbital_lookup_table: pointer to int;
 *
 * Returns: real
 *
 * Action:  sums up the charge matrix for 'MO and determines the contribution of
 *   'which_atom to this MO.
 *
 ****************************************************************************/
real atom_contribution(num_orbs,num_atoms,kpoint,MO,avg_prop_info,which_atom,
                        orbital_lookup_table)
  int num_orbs,num_atoms,kpoint,MO;
  avg_prop_info_type *avg_prop_info;
  int which_atom,*orbital_lookup_table;
{
  int i,atom_begin,atom_end;
  float *chg_mat_ptr;
  int AO;
  real contrib;

  chg_mat_ptr = &(avg_prop_info[kpoint].chg_mat[MO*num_orbs]);

  /* find this atom's orbitals */
  find_atoms_orbs(num_orbs,num_atoms,which_atom,orbital_lookup_table,
                  &atom_begin,&atom_end);
  contrib = 0.0;

  if( atom_begin >= 0 ){
    for(i=atom_begin;i<atom_end;i++){
      contrib += (real)chg_mat_ptr[i]/(real)MULTIPLIER;
    }
  }

  return(contrib);
}


/****************************************************************************
 *
 *                   Procedure gen_projected_DOS
 *
 * Arguments: details: pointer to detail_type
 *              cell: pointer to cell type
 *          num_orbs: int
 *     avg_prop_info: pointer to avg_prop_info_type
 *  orbital_ordering: pointer to K_orb_ptr_type
 * orbital_lookup_table: pointer to int
 *
 *
 * Returns: none
 *
 * Action:  Generates all of the projected DOS curves.
 *
 ****************************************************************************/
void gen_projected_DOS(details,cell,num_orbs,avg_prop_info,orbital_ordering,
                       orbital_lookup_table)
  detail_type *details;
  cell_type *cell;
  int num_orbs;
  avg_prop_info_type *avg_prop_info;
  K_orb_ptr_type *orbital_ordering;
  int *orbital_lookup_table;
{
  int i,j,k,l;
  int tot_num_orbs,num_states,begin,end;
  real num_at_this_E;
  real diff;
  real accum;
  real this_E;
  real tot_K_weight;

  tot_K_weight = 0.0;
  for(i=0;i<details->num_KPOINTS;i++){
    tot_K_weight += details->K_POINTS[i].weight;
  }

  tot_num_orbs = num_orbs * details->num_KPOINTS;


  for (k=0;k<details->num_proj_DOS;k++){
    i=0;
    /* write out some identifying information */
    fprintf(output_file,"\n### PROJECTED DENSITY OF STATES \n");
    switch(details->proj_DOS[k].type){
    case P_DOS_ORB:
      fprintf(output_file,"%lf states are present.\n",
              details->proj_DOS[k].weight_sum * 2.0);
      fprintf(output_file,"ORBITAL ");
      for(j = 0; j<details->proj_DOS[k].num_contributions; j++){
        fprintf(output_file,"%d %lf,",details->proj_DOS[k].contributions[j]+1,
                details->proj_DOS[k].weights[j]);
      }
      break;
    case P_DOS_ATOM:
      num_states = 0;
      for (j=0; j<details->proj_DOS[k].num_contributions; j++){
        find_atoms_orbs(num_orbs,cell->num_atoms,
                        details->proj_DOS[k].contributions[j],
                        orbital_lookup_table,&begin,&end);
        num_states += details->proj_DOS[k].weights[j]*(end - begin);
      }
      fprintf(output_file,"%lf states are present.\n",(real)num_states * 2.0);
      fprintf(output_file,"ATOM ");
      for(j = 0; j<details->proj_DOS[k].num_contributions; j++){
        fprintf(output_file,"%d %lf,",details->proj_DOS[k].contributions[j]+1,
                details->proj_DOS[k].weights[j]);
      }
      break;
    case P_DOS_FMO:
      num_states=1;
      fprintf(output_file,"%lf states are present.\n",
              details->proj_DOS[k].weight_sum * num_states * 2);
      fprintf(output_file,"FMO ");
      for(j = 0; j<details->proj_DOS[k].num_contributions; j++){
        fprintf(output_file,"%d %lf,",details->proj_DOS[k].contributions[j]+1,
                details->proj_DOS[k].weights[j]);
      }
      break;
    }
    fprintf(output_file,"\n");

    fprintf(output_file,"#BEGIN CURVE\n");
    while(i<tot_num_orbs){

      /* now figure out the contribution of whatever is being projected out */
      num_at_this_E = 0.0;
      for( l=0; l<details->proj_DOS[k].num_contributions; l++){
        switch(details->proj_DOS[k].type){
        case P_DOS_ORB:
          num_at_this_E += details->proj_DOS[k].weights[l] *
            orb_contribution(num_orbs,orbital_ordering[i].Kpoint,
                             orbital_ordering[i].MO,
                             avg_prop_info,
                             details->proj_DOS[k].contributions[l]);
          break;
        case P_DOS_ATOM:
          num_at_this_E += details->proj_DOS[k].weights[l] *
            atom_contribution(num_orbs,cell->num_atoms,
                              orbital_ordering[i].Kpoint,
                              orbital_ordering[i].MO,
                              avg_prop_info,
                              details->proj_DOS[k].contributions[l],
                              orbital_lookup_table);
          break;
        case P_DOS_FMO:
          num_at_this_E += details->proj_DOS[k].weights[l] *
            FMO_contribution(num_orbs,
                              orbital_ordering[i].Kpoint,
                              orbital_ordering[i].MO,
                              avg_prop_info,
                              details->proj_DOS[k].contributions[l]);
          break;
        default:
          FATAL_BUG("Invalid projection type in gen_proj_DOS.");
        }
      }

      /* multiply by the total DOS weighting */
      num_at_this_E *= (real)details->K_POINTS[orbital_ordering[i].Kpoint].weight;

      this_E = (real)*(orbital_ordering[i].energy);
      j=i+1;

      /*******
        just add up all of the contributions to the total DOS that fall
        within DOS_DEGEN_TOL of each other...
        ********/
      if( j < tot_num_orbs ){
        diff = (real)*(orbital_ordering[i].energy) - (real)*(orbital_ordering[j].energy);
        while(fabs(diff) < DOS_DEGEN_TOL && j<tot_num_orbs){

          accum = 0.0;
          for( l=0; l<details->proj_DOS[k].num_contributions; l++){
            switch(details->proj_DOS[k].type){
            case P_DOS_ORB:
              accum += details->proj_DOS[k].weights[l] *
                orb_contribution(num_orbs,orbital_ordering[j].Kpoint,
                                 orbital_ordering[j].MO,
                                 avg_prop_info,
                                 details->proj_DOS[k].contributions[l]);
              break;
            case P_DOS_ATOM:
              accum += details->proj_DOS[k].weights[l] *
                atom_contribution(num_orbs,cell->num_atoms,
                                  orbital_ordering[j].Kpoint,
                                  orbital_ordering[j].MO,
                                  avg_prop_info,
                                  details->proj_DOS[k].contributions[l],
                                  orbital_lookup_table);
              break;
            case P_DOS_FMO:
              accum += details->proj_DOS[k].weights[l] *
                FMO_contribution(num_orbs,
                                 orbital_ordering[j].Kpoint,
                                 orbital_ordering[j].MO,
                                 avg_prop_info,
                                 details->proj_DOS[k].contributions[l]);
              break;

            default:
              FATAL_BUG("Invalid projection type in gen_proj_DOS.");
            }
          }
          j++;

          if( j < tot_num_orbs ){
            diff = (real)*(orbital_ordering[i].energy) - (real)*(orbital_ordering[j].energy);
            accum *= (real)details->K_POINTS[orbital_ordering[j].Kpoint].weight;
          }

          num_at_this_E += accum;
        }
      }
      i = j;
      /* write out the result */
      fprintf(output_file,"%lf %lf\n",num_at_this_E/tot_K_weight,this_E);
    }
    fprintf(output_file,"#END CURVE\n");
    fprintf(output_file," \n");
  }
}
