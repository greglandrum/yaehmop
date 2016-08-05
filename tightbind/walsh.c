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
*     this file has the routines for handling walsh diagrams.
*
*  created:  greg landrum  January 1994
*
*****************************************************************************/
#include "bind.h"


/****
  Recent edit history

  30.10.98 gL:
   fixed problem with printing of distance matrix into .walsh file.

****/

/****************************************************************************
*
*                   Procedure auto_walsh
*
* Arguments:  values: pointer to type real
*          num_steps: int
*          begin,end: reals
*
* Returns: none
*
* Action:  Fills in all the steps for a walsh diagram running from 'begin to
*   'end in 'num_steps steps....
*
*  NOTE:  I'll give you a fence post problem....
*
*****************************************************************************/
void auto_walsh(values,num_steps,begin,end)
  real *values;
  int num_steps;
  real begin,end;
{
  int i;
  real tot_dist,step,cur_val;

  /* find the total distance */
  tot_dist = end-begin;

  /* figure out how big the step needs to be */
  step = tot_dist/(real)(num_steps-1);

  cur_val = begin;
  for(i=0;i<num_steps;i++){
    values[i] = cur_val;
    cur_val += step;
  }

}




/****************************************************************************
*
*                   Procedure walsh_update
*
* Arguments:  cell: pointer to cell_type
*          details: pointer to detail_type
*             step: integer
*         printing: char
*
* Returns: none
*
* Action:  Updates the geometry of the problem for step 'step in the walsh
*     diagram.
*
*  Note:  due to the potential of there being zeta optimization going on,
*    this routine makes a copy of the atoms the first time it is called and
*    resets all the zeta values to their initial values for each step.
*
*****************************************************************************/
void walsh_update(cell_type *cell,detail_type *details,int step,char printing)
{
  static atom_type *storage=0;
  static geom_frag_type *geom_frag_store=0;
  static xtal_defn_type xtal_storage;
  geom_frag_type *geom_frag,*new_frag,*prev_frag;
  int i;
  int num_atoms,num_vars,which_var,num_steps;
  real *values;
  atom_type *atom;

  num_atoms = cell->num_raw_atoms+cell->dim;
  num_vars = details->walsh_details.num_vars;
  num_steps = details->walsh_details.num_steps;

  /* some error checking */
  if( !num_vars ){
    NONFATAL_BUG("walsh_update called when there are no walsh variables.");
    return;
  }

  /* get a pointer for the variable values array */
  values = details->walsh_details.values;

  /* get memory for the atom store if we need it */
  if( !step && !storage ){
    storage = (atom_type *)calloc(num_atoms,sizeof(atom_type));
    if(!storage) fatal("Can't get space for the atomic storage array in walsh_update.");

    /* copy the atoms */
    bcopy((char *)cell->atoms,(char *)storage,num_atoms*sizeof(atom_type));

    /******

      if we are using crystallographic coordinates, then copy over
      the crystallograpic parms now (they can be varied in a Walsh
      diagram too.

    ********/
    if( cell->using_xtal_coords ){
      bcopy((char *)&(cell->xtal_defn),(char *)&xtal_storage,sizeof(xtal_defn_type));
    }

    /********

      if we're using geom_frags, copy them now

      the ordering of this is very important. The linked list is copied
       out forwards, so that it can be updated in order later.  the ordering
       of the atoms is crucial, because that is what will eventually determine
       the numbering of them in the cell->atoms array.  Changes here will
       completely munge projected DOSs etc.

    ********/
    geom_frag = cell->geom_frags;
    prev_frag = geom_frag_store = 0;
    while(geom_frag){
      new_frag = (geom_frag_type *)calloc(1,sizeof(geom_frag_type));
      if( !new_frag ) fatal("Can't get space for new_frag in walsh_update");
      /* copy the frag */
      bcopy(geom_frag,new_frag,sizeof(geom_frag_type));

      /* get its own atom list, etc. */
      new_frag->atoms = (atom_type *)calloc(new_frag->num_atoms,sizeof(atom_type));
      if( !new_frag->atoms) fatal("Can't get new_frag->atoms in walsh_update.");

      /* copy in the atoms */
      bcopy(geom_frag->atoms,new_frag->atoms,
            geom_frag->num_atoms*sizeof(atom_type));

      /*****

        put it in the list

      ******/
      if( !prev_frag ){
        /* the first frag */
        geom_frag_store = new_frag;
      } else{
        prev_frag->next = new_frag;
      }
      prev_frag = new_frag;

      /* move on to the next fragment */
      geom_frag = geom_frag->next;
    }
  }

  if( printing ){
    fprintf(output_file,"\n\n;-----------------------------------------------\n");

    /****

      print out some status information which will be used to generate the actual
      Walsh diagrams

      *****/
    fprintf(output_file,"#Walsh_step:  %d ",step);
    for(i=0;i<num_vars;i++){
      fprintf(output_file,"%lg ",values[i*num_steps + step]);
    }
    fprintf(output_file,"\n");
    fprintf(output_file,";-----------------------------------------------\n");
  }

  /******

    update the positions

  *******/

  /******
    first copy the atoms from the storage array to make sure that they all have
     the proper zeta values.
  ******/
  bcopy((char *)storage,(char *)cell->atoms,num_atoms*sizeof(atom_type));
  if( cell->using_xtal_coords ){
    bcopy((char *)&xtal_storage,(char *)&(cell->xtal_defn),sizeof(xtal_defn_type));
  }

  /* now loop over atoms */
  for(i=0;i<num_atoms;i++){
    atom = &(cell->atoms[i]);
    /*********
      check each coordinate to see if it is variable

        NOTE: this assumes that there can't be any initial values greater
         than 1000
    **********/
    if( cell->using_Zmat ){
      which_var = floor( fabs(atom->Zmat_loc.bond_length / 1000.0) );
      if( which_var ) atom->Zmat_loc.bond_length =
        values[(which_var-1)*num_steps + step];
      which_var = floor( fabs(atom->Zmat_loc.alpha / 1000.0) );
      if( which_var ) atom->Zmat_loc.alpha =
        values[(which_var-1)*num_steps + step];
      which_var = floor( fabs(atom->Zmat_loc.beta / 1000.0) );
      if( which_var ) atom->Zmat_loc.beta =
        values[(which_var-1)*num_steps + step];
    }
    else{
      which_var = floor( fabs(atom->loc.x / 1000.0) );
      if( which_var ) atom->loc.x = values[(which_var-1)*num_steps + step];
      which_var = floor( fabs(atom->loc.y / 1000.0) );
      if( which_var ) atom->loc.y = values[(which_var-1)*num_steps + step];
      which_var = floor( fabs(atom->loc.z / 1000.0) );
      if( which_var ) atom->loc.z = values[(which_var-1)*num_steps + step];
    }
  }


  /* check to see if the crystal definition needs to be updated */
  if( cell->using_xtal_coords ){
    for( i=0;i<3;i++ ){
      which_var = floor( fabs(cell->xtal_defn.axis_lengths[i] / 1000.0) );
      if( which_var )
        cell->xtal_defn.axis_lengths[i] = values[(which_var-1)*num_steps + step];
      which_var = floor( fabs(cell->xtal_defn.angles[i] / 1000.0) );
      if( which_var )
        cell->xtal_defn.angles[i] = values[(which_var-1)*num_steps + step];
    }
  }

  /* figure out where the atoms are */
  if( cell->using_Zmat ){
    eval_Zmat_locs(cell->atoms,cell->num_raw_atoms,cell->dim,printing);
  }

  if( cell->using_xtal_coords ){
    eval_xtal_coord_locs(unit_cell,printing);
  }

  /*****
    if we are doing a principle axis analysis, move the atoms to
    the principle axis frame.
  *****/
  if( details->find_princ_axes ){
    full_transform(cell->atoms,cell->COM,cell->princ_axes,num_atoms);
  }

  /******

    deal with the geom_frags

  ******/
  if( cell->geom_frags ){
    /*******

      first copy in each frag and check to see if it has variables that
      need to be updated.

    ********/
    geom_frag = cell->geom_frags;
    new_frag = geom_frag_store;
    while(geom_frag){
      if( !new_frag ) FATAL_BUG("Ran out of geom_frags in copy over.");
      bcopy(new_frag->atoms,geom_frag->atoms,
            geom_frag->num_atoms*sizeof(atom_type));
      for(i=0;i<geom_frag->num_atoms;i++){

        atom = &(geom_frag->atoms[i]);
        /*********
          check each coordinate to see if it is variable

          NOTE: this assumes that there can't be any initial values greater
          than 1000
          **********/
        if( geom_frag->using_Z_mat ){
          which_var = floor( fabs(atom->Zmat_loc.bond_length / 1000.0) );
          if( which_var ) atom->Zmat_loc.bond_length =
            values[(which_var-1)*num_steps + step];
          which_var = floor( fabs(atom->Zmat_loc.alpha / 1000.0) );
          if( which_var ) atom->Zmat_loc.alpha =
            values[(which_var-1)*num_steps + step];
          which_var = floor( fabs(atom->Zmat_loc.beta / 1000.0) );
          if( which_var ) atom->Zmat_loc.beta =
            values[(which_var-1)*num_steps + step];
        }
        else{
          which_var = floor( fabs(atom->loc.x / 1000.0) );
          if( which_var ) atom->loc.x = values[(which_var-1)*num_steps + step];
          which_var = floor( fabs(atom->loc.y / 1000.0) );
          if( which_var ) atom->loc.y = values[(which_var-1)*num_steps + step];
          which_var = floor( fabs(atom->loc.z / 1000.0) );
          if( which_var ) atom->loc.z = values[(which_var-1)*num_steps + step];
        }
      }
      geom_frag = geom_frag->next;
      new_frag = new_frag->next;
    }
    process_geom_frags(cell);
  }


  if( printing && !cell->using_Zmat){
    /************

      write out the new atomic positions

      *************/
    fprintf(output_file,"\n# ********* Atomic positions this Walsh Step:  *********\n");
    for(i=0;i<cell->num_atoms;i++){
      atom = &(cell->atoms[i]);
      fprintf(output_file,"%d %s %8.4f %8.4f %8.4g\n",
              i+1,atom->symb,atom->loc.x,atom->loc.y,atom->loc.z);
    }

    /* that's it! */
  }
}




/****************************************************************************
*
*                   Procedure walsh_output
*
* Arguments:  details: pointer to detail_type
*                cell: pointer to unit_cell_type
*            num_orbs: int
*            eigenset: pointer to eigenset_type
*             overlap: hermetian_matrix_type
*               hamil: hermetian_matrix_type
*          properties: prop_type
*orbital_lookup_table: pointer to int
*                step: int
*
* Returns: none
*
* Action:   This goes through the Walsh printing options in 'details and
*  puts the relevant information in the walsh output file.
*
*****************************************************************************/
void walsh_output(details,cell,num_orbs,eigenset,overlap,hamil,
                  properties,orbital_lookup_table,step)
  detail_type *details;
  cell_type *cell;
  int num_orbs;
  eigenset_type eigenset;
  hermetian_matrix_type overlap,hamil;
  prop_type properties;
  int *orbital_lookup_table;
  int step;
{
  static char first_call=1;
  int i;
  int things_so_far;
  int atom1,atom2;
  printing_info_type *p_info;


  /**********
    if this is the first time that this function was called in this particular run
    of the program then print out header information that tells what each of the
    columns in the file are.
  ***********/
  if( first_call ){
    fprintf(walsh_file,"# Walsh output for job: %s\n",details->title);
    fprintf(walsh_file,"# This is the key to the columns printed out below.\n");
    things_so_far = 1;

    for(i=0;i<details->walsh_details.num_vars;i++){
      fprintf(walsh_file,"# %d  Walsh Variable %d\n",things_so_far,i+1);
      things_so_far++;
    }

    /* step through the elements of the printing_info linked list */
    p_info = details->step_print_options;

    while(p_info){
      switch(p_info->which_to_print){
      case PRT_OP:
        switch(p_info->type){
        case P_DOS_ATOM:
          fprintf(status_file,"Atom Projection for Overlap Population ignored.\n");
          break;
        case P_DOS_ORB:
          fprintf(walsh_file,"# %d Overlap Population between orbitals %d and %d\n",
                  things_so_far,p_info->contrib1+1,p_info->contrib2+1);
          things_so_far++;
          break;
        }
        break;
      case PRT_ROP:
        switch(p_info->type){
        case P_DOS_ATOM:
          fprintf(walsh_file,"# %d Reduced Overlap Population between atoms %d and %d\n",
                  things_so_far,p_info->contrib1+1,p_info->contrib2+1);
          things_so_far++;
          break;
        case P_DOS_ORB:
          fprintf(status_file,"Orbital Projection for Reduced Overlap Population ignored.\n");
          break;
        }
        break;
      case PRT_OVERLAP:
        switch(p_info->type){
        case P_DOS_ATOM:
          fprintf(status_file,"Atom Projection for Overlap matrix ignored.\n");
          break;
        case P_DOS_ORB:
          fprintf(walsh_file,"# %d Overlap matrix element between orbitals %d and %d\n",
                  things_so_far,p_info->contrib1+1,p_info->contrib2+1);
          things_so_far++;
          break;
        }
        break;
      case PRT_HAMIL:
        switch(p_info->type){
        case P_DOS_ATOM:
          fprintf(status_file,"Atom Projection for Hamiltonian matrix ignored.\n");
          break;
        case P_DOS_ORB:
          fprintf(walsh_file,"# %d Hamiltonian matrix element between orbitals %d and %d\n",
                  things_so_far,p_info->contrib1+1,p_info->contrib2+1);
          things_so_far++;
          break;
        }
        break;
      case PRT_CHG_MAT:
        switch(p_info->type){
        case P_DOS_ATOM:
          fprintf(status_file,"Atom Projection for charge matrix ignored.\n");
          break;
        case P_DOS_ORB:
          fprintf(walsh_file,"# %d Charge matrix element between orbitals %d and %d\n",
                  things_so_far,p_info->contrib1+1,p_info->contrib2+1);
          things_so_far++;
          break;
        }
        break;
      case PRT_RCHG_MAT:
        switch(p_info->type){
        case P_DOS_ATOM:
          fprintf(walsh_file,"# %d Reduced charge matrix element between atoms %d and %d\n",
                  things_so_far,p_info->contrib1+1,p_info->contrib2+1);
          things_so_far++;
          break;
        case P_DOS_ORB:
          fprintf(status_file,"Orbital Projection for Reduced charge matrix ignored.\n");
          break;
        }
        break;
      case PRT_NET_CHG:
        switch(p_info->type){
        case P_DOS_ATOM:
          fprintf(walsh_file,"# %d Net charge on atom: %d \n",
                  things_so_far,p_info->contrib1+1);
          things_so_far++;
          break;
        case P_DOS_ORB:
          fprintf(status_file,"Orbital Projection for net charge\n");
          break;
        }
        break;
      case PRT_WAVE_FUNC:
        switch(p_info->type){
        case P_DOS_ATOM:
          fprintf(status_file,"Atom Projection for wave functions ignored.\n");
          break;
        case P_DOS_ORB:
          fprintf(walsh_file,"# %d Coefficient of AO %d in MO %d\n",
                  things_so_far,p_info->contrib1+1,p_info->contrib2+1);
          things_so_far++;
          break;
        }
        break;
      case PRT_DIST:
        switch(p_info->type){
        case P_DOS_ATOM:
          fprintf(walsh_file,"# %d Distance between atoms %d and %d\n",
                  things_so_far,p_info->contrib1+1,p_info->contrib2+1);
          things_so_far++;
          break;
        case P_DOS_ORB:
          fprintf(status_file,"Orbital Projection for distance matrix ignored.\n");
          break;
        }
        break;
      case PRT_ELECTROSTAT:
        fprintf(walsh_file,"# %d Electrostatic Repulsion Energy.\n",
                things_so_far);
        things_so_far++;
        break;
      case PRT_ENERGIES:
        fprintf(walsh_file,"# %d Total Energy.\n",
                things_so_far);
        things_so_far++;
        break;
      case PRT_ORB_ENERGY:
        fprintf(walsh_file,"# %d Energy of Orbital: %d.\n",
                things_so_far,p_info->contrib1+1);
        things_so_far++;
        break;
      case PRT_ORB_COEFF:
        fprintf(walsh_file,"# %d Coefficient of AO: %d in MO: %d.\n",
                things_so_far,p_info->contrib1+1,p_info->contrib2+1);
        things_so_far++;
        break;
      }
      p_info = p_info->next;
    }
    first_call = 0;
  }

  /*******************
    now actually print out the results that are desired
  ********************/

  /* print out the values of the Walsh variables */
  for(i=0;i<details->walsh_details.num_vars;i++){
    fprintf(walsh_file,"%8.6lf ",
            details->walsh_details.values[i*details->walsh_details.num_steps+step]);
  }

  /* Now step through the elements of the printing_info linked list */
  p_info = details->step_print_options;

  while(p_info){
    switch(p_info->which_to_print){
    case PRT_OP:
      switch(p_info->type){
      case P_DOS_ORB:
        fprintf(walsh_file,"%8.6lf ",properties.OP_mat[p_info->contrib1*num_orbs+
                                                        p_info->contrib2]);
        break;
      }
      break;
    case PRT_ROP:
      switch(p_info->type){
      case P_DOS_ATOM:
        if( p_info->contrib1 > p_info->contrib2 ){
          fprintf(walsh_file,"%8.6lf ",
                  properties.ROP_mat[p_info->contrib1*(p_info->contrib1+1)/2
                                     + p_info->contrib2]);
        } else{
          fprintf(walsh_file,"%8.6lf ",
                  properties.ROP_mat[p_info->contrib2*(p_info->contrib2+1)/2
                                     + p_info->contrib1]);
        }
        break;
      }
      break;
    case PRT_OVERLAP:
      switch(p_info->type){
      case P_DOS_ORB:
        fprintf(walsh_file,"%8.6lf ",
                HERMETIAN_R(overlap,(p_info->contrib1),(p_info->contrib2)));
                break;
      }
      break;
    case PRT_HAMIL:
      switch(p_info->type){
      case P_DOS_ORB:
        fprintf(walsh_file,"%8.6lf ",
                HERMETIAN_R(hamil,(p_info->contrib1),(p_info->contrib2)));
                break;
      }
      break;
    case PRT_CHG_MAT:
      switch(p_info->type){
      case P_DOS_ORB:
        fprintf(walsh_file,"%8.6lf ",properties.chg_mat[p_info->contrib1*num_orbs+
                                                         p_info->contrib2]);
        break;
      }
      break;
    case PRT_RCHG_MAT:
      switch(p_info->type){
      case P_DOS_ATOM:
        fprintf(walsh_file,"%8.6lf ",properties.Rchg_mat[p_info->contrib1*cell->num_atoms
                                                          +p_info->contrib2]);
          break;
      }
      break;
    case PRT_NET_CHG:
      switch(p_info->type){
      case P_DOS_ATOM:
        fprintf(walsh_file,"%8.6lf ",properties.net_chgs[p_info->contrib1]);
        break;
      }
      break;
    case PRT_WAVE_FUNC:
      switch(p_info->type){
      case P_DOS_ORB:
        fprintf(walsh_file,"%8.6lf ",EIGENVECT_R(eigenset,p_info->contrib2,
                                                 p_info->contrib1));
        break;
      }
      break;
    case PRT_DIST:
      switch(p_info->type){
      case P_DOS_ATOM:
        /********
          since the distance mat is stored as a symmetric matrix, we have to do
          a little extra work to get the element out
          ********/
        if( p_info->contrib1 < p_info->contrib2 ){
          atom1 = p_info->contrib1;
          atom2 = p_info->contrib2;
        } else{
          atom1 = p_info->contrib2;
          atom2 = p_info->contrib1;
        }
        fprintf(walsh_file,"%8.6lf ",
                cell->distance_mat[(atom2*(atom2+1))/2 + atom1]);
        break;
      }
      break;
    case PRT_ENERGIES:
      fprintf(walsh_file,"%8.6lf ",properties.total_E);
      break;
    case PRT_ORB_ENERGY:
      fprintf(walsh_file,"%8.6lf ",EIGENVAL(eigenset,p_info->contrib1));
      break;
    case PRT_ORB_COEFF:
      fprintf(walsh_file,"%8.6lf ",
              EIGENVECT_R(eigenset,p_info->contrib2,p_info->contrib1));
      break;
    case PRT_ELECTROSTAT:
      fprintf(walsh_file,"%8.6lf ",properties.electrostat_E);
      break;
    }
    p_info = p_info->next;
  }
  /* put in a carriage return to finish the line */
  fprintf(walsh_file,"\n");
}

