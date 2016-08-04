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
  Recent Edit History:
  23.06.98 gL:
    cleaned up stupid problems in grow_solid_with_surface
    which arose due to the changes in MO_center_list_type
  18.08.98 gL:
    turn off *all* selected atoms when either grow_solid
    or grow_solid_with_surface are called....
  26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)
  29.01.99 gL:
     changes to deal with new representation of lattice_vects
  21.05.1999 gL:
     YA core leak plugged
***/

/********

  This has everything that has to be kept separate for solids....
    this isn't a whole lot, because basically we can treat solids
    like molecules.

  written by greg Landrum,  September 1994

*********/


#include "viewkel.h"

/****************************************************************************
 *
 *                   Procedure grow_solid
 *
 * Arguments: num_args: int
 *             solid_p: array of pointers to char
 *
 * Returns: none
 *
 * Action: Prompts the user for the number of cells along each lattice
 *    vector to display, calculates the position of atoms along
 *    the lattice vectors, and adds them to the list of atoms.
 *
 ****************************************************************************/
void grow_solid(int num_args,char *solid_p[MAX_ARGS])
{
  int i,j,k,l;
  static int num_a=2, num_b=2, num_c=2;
  int frame;
  point_type dist_a,dist_b,dist_c;
  molec_type *solid;
  atom_type *temp_atoms;
  int old_num_atoms;
  int num_added;

  /* first get a pointer to the molecule that we are working with */
  solid = (molec_type *)solid_p[0];

  if( solid->num_dim < 1 ){
    display("Whoops...");
    error("grow_solid called with a molecule as an argument. This is a bug.");
    return;
  }


  display("Look in the xterm");
  printf("This crystal is %d dimensional.\n",solid->num_dim);
  readintparm("number along (a)",&num_a);
  if( num_a < 1 ){
    display("HEY!");
    error("Don't enter dumb values!");
    num_a = 1;
  }
  if( solid->num_dim > 1 ){
    readintparm("number along (b)",&num_b);
    if( num_b < 1 ){
      display("HEY!");
      error("Don't enter dumb values!");
      num_b = 1;
    }
  }
  else{
    num_b = num_c = 1;
  }
  if( solid->num_dim > 2 ){
    readintparm("number along (c)",&num_c);
    if( num_c < 1 ){
      display("HEY!");
      error("Don't enter dumb values!");
      num_c = 1;
    }
  }
  else{
    num_c = 1;
  }

  solid->num_along[0] = num_a;
  solid->num_along[1] = num_b;
  solid->num_along[2] = num_c;

  /* get the memory we'll need */
  if( solid->num_frames > 1){
    old_num_atoms = solid->num_atoms;
  } else {
    old_num_atoms = solid->num_atoms_in_cell;
  }
  solid->num_atoms = solid->num_atoms_in_cell*num_a*num_b*num_c;
  temp_atoms = (atom_type *)D_CALLOC(solid->num_atoms*solid->num_frames,
                                   sizeof(atom_type));
  if( !temp_atoms ) fatal("Memory allocation: can't get space for more atoms.");

  /* now grow the crystal */
  num_added = 0;
  for(frame = 0;frame < solid->num_frames; frame++){
    for(i=0;i<num_a;i++){
      dist_a.x = i*solid->orig_lattice[1].x;
      dist_a.y = i*solid->orig_lattice[1].y;
      dist_a.z = i*solid->orig_lattice[1].z;

      for(j=0;j<num_b;j++){
        dist_b.x = j*solid->orig_lattice[2].x;
        dist_b.y = j*solid->orig_lattice[2].y;
        dist_b.z = j*solid->orig_lattice[2].z;

        for(k=0;k<num_c;k++){
          dist_c.x = k*solid->orig_lattice[3].x;
          dist_c.y = k*solid->orig_lattice[3].y;
          dist_c.z = k*solid->orig_lattice[3].z;

          /* first copy in the old atom data */
          bcopy((char *)&(solid->atoms[frame*old_num_atoms]),
                (char *)&(temp_atoms[num_added]),
                solid->num_atoms_in_cell*sizeof(atom_type));

          /* now update the locations */
          for(l=0;l<solid->num_atoms_in_cell;l++){
            temp_atoms[num_added].loc.x =
              solid->atoms[frame*old_num_atoms+l].loc.x +
              dist_a.x + dist_b.x + dist_c.x;
            temp_atoms[num_added].loc.y =
              solid->atoms[frame*old_num_atoms+l].loc.y +
              dist_a.y + dist_b.y + dist_c.y;
            temp_atoms[num_added].loc.z =
              solid->atoms[frame*old_num_atoms+l].loc.z +
              dist_a.z + dist_b.z + dist_c.z;
            temp_atoms[num_added].num = num_added;

            /* we copied over some pointers too... this is bad */
            temp_atoms[num_added].linesto = 0;
#ifdef INCLUDE_ADF_PLOTS
            temp_atoms[num_added].displacements = 0;
#endif
            temp_atoms[num_added].p_surf = 0;
            /* make sure that copied atoms are not selected */
            temp_atoms[num_added].is_selected = 0;
            num_added++;
          }
        }
      }
    }
  }

  num_selected = 0;

  /* free up the atoms we currently have stored */
  for(i=0;i<old_num_atoms;i++){
    if( solid->atoms[i].p_surf ){
      D_FREE(solid->atoms[i].p_surf);
    }
    if( solid->atoms[i].linesto ){
      D_FREE(solid->atoms[i].linesto);
    }
  }
  D_FREE(solid->atoms);
  solid->atoms = temp_atoms;

  /* redetermine the bond locations */
  determine_connections(solid);

  /* update the box */
  if( num_a > 1 ) num_a--;
  if( num_b > 1 ) num_b--;
  if( num_c > 1 ) num_c--;
  solid->cell_box[0].x = solid->lattice_vect[0].x;
  solid->cell_box[0].y = solid->lattice_vect[0].y;
  solid->cell_box[0].z = solid->lattice_vect[0].z;
  solid->cell_box[1].x =
    num_a*solid->orig_lattice[1].x;
  solid->cell_box[1].y =
    num_a*solid->orig_lattice[1].y;
  solid->cell_box[1].z =
    num_a*solid->orig_lattice[1].z;
  solid->cell_box[3].x =
    num_b*solid->orig_lattice[2].x;
  solid->cell_box[3].y =
    num_b*solid->orig_lattice[2].y;
  solid->cell_box[3].z =
    num_b*solid->orig_lattice[2].z;
  solid->cell_box[4].x =
    num_b*solid->orig_lattice[3].x;
  solid->cell_box[4].y =
    num_b*solid->orig_lattice[3].y;
  solid->cell_box[4].z =
    num_b*solid->orig_lattice[3].z;
  V3Add(&solid->cell_box[1],&solid->cell_box[3],
        &solid->cell_box[2]);
  V3Add(&solid->cell_box[1],&solid->cell_box[4],
        &solid->cell_box[5]);
  V3Add(&solid->cell_box[2],&solid->cell_box[4],
        &solid->cell_box[6]);
  V3Add(&solid->cell_box[3],&solid->cell_box[4],
        &solid->cell_box[7]);
  for(i=1;i<8;i++){
    V3Add(&solid->cell_box[i],&solid->cell_box[0],
          &solid->cell_box[i]);
  }

  display("There ya go!");
}



/****************************************************************************
 *
 *                   Procedure grow_solid_with_surface
 *
 * Arguments: num_args: int
 *             solid_p: array of pointers to char
 *
 * Returns: none
 *
 * Action: This is the same as grow_solid (above), only the wavefunction
 *   coefficients are updates as well.
 *
 ****************************************************************************/
void grow_solid_with_surface(int num_args,char *solid_p[MAX_ARGS])
{
  int i,j,k,l,m,n,MO;
  static int num_a=2, num_b=2, num_c=2;
  point_type dist_a,dist_b,dist_c;
  MO_surface_type *MO_surf;
  MO_center_list_type *temp_centers;
  molec_type *solid;
  atom_type *temp_atoms;
  int num_added,old_num_atoms;
  float kdotr,cos_kdotr,sin_kdotr;

  /* first get a pointer to the molecule that we are working with */
  solid = (molec_type *)solid_p[0];
  MO_surf = (MO_surface_type *)solid_p[1];

  old_num_atoms = solid->num_atoms;
  if( solid->num_dim < 1 ){
    display("Whoops...");
    error("grow_solid_with_surface called with a molecule as an argument. This is a bug.");
    return;
  }

  display("Look in the xterm");
  printf("This crystal is %d dimensional.\n",solid->num_dim);
  printf(" Please enter the number of cells along each lattice direction on separate lines.\n");

  readintparm("number along (a)",&num_a);
  if( num_a < 1 ){
    display("HEY!");
    error("Don't enter dumb values!");
    num_a = 1;
  }
  if( solid->num_dim > 1 ){
    readintparm("number along (b)",&num_b);
    if( num_b < 1 ){
      display("HEY!");
      error("Don't enter dumb values!");
      num_b = 1;
    }
  }
  else{
    num_b = num_c = 1;
  }
  if( solid->num_dim > 2 ){
    readintparm("number along (c)",&num_c);
    if( num_c < 1 ){
      display("HEY!");
      error("Don't enter dumb values!");
      num_c = 1;
    }
  }
  else{
    num_c = 1;
  }

  solid->num_along[0] = num_a;
  solid->num_along[1] = num_b;
  solid->num_along[2] = num_c;

  /* clear out some old memory */
  for(i=0;i<MO_surf->num_centers;i++){
    D_FREE(MO_surf->MO_centers[i].AO_list);
  }

  D_FREE(MO_surf->MO_centers);

  /* get the memory we'll need */
  solid->num_atoms = solid->num_atoms_in_cell*num_a*num_b*num_c;
  temp_atoms = (atom_type *)D_CALLOC(solid->num_atoms,sizeof(atom_type));
  if( !temp_atoms )
    fatal("Memory allocation: can't get space for more atoms.");
  MO_surf->num_centers = MO_surf->num_centers_in_cell*num_a*num_b*num_c;
  temp_centers = (MO_center_list_type *)D_CALLOC(MO_surf->num_centers,
                                               sizeof(MO_center_list_type));
  if( !temp_centers ) fatal("Can't get space for more orbitals.");

  /* now grow the crystal */
  num_added = 0;
  for(i=0;i<num_a;i++){
    dist_a.x = i*solid->orig_lattice[1].x;
    dist_a.y = i*solid->orig_lattice[1].y;
    dist_a.z = i*solid->orig_lattice[1].z;

    for(j=0;j<num_b;j++){
      dist_b.x = j*solid->orig_lattice[2].x;
      dist_b.y = j*solid->orig_lattice[2].y;
      dist_b.z = j*solid->orig_lattice[2].z;

      for(k=0;k<num_c;k++){
        dist_c.x = k*solid->orig_lattice[3].x;
        dist_c.y = k*solid->orig_lattice[3].y;
        dist_c.z = k*solid->orig_lattice[3].z;

        /* first copy in the old atom data */
        bcopy((char *)solid->atoms,(char *)&(temp_atoms[num_added]),
              solid->num_atoms_in_cell*sizeof(atom_type));

        /* copy the old center data */
        bcopy((char *)MO_surf->raw_MO_centers,
              (char *)&(temp_centers[num_added]),
              MO_surf->num_centers_in_cell*
              sizeof(MO_center_list_type));

#ifdef DEBUG
        fprintf(stderr,"raw: C: %lf, Ci:% lf \n",
                MO_surf->raw_MO_centers[0].AO_list[0].coeff[0],
                MO_surf->raw_MO_centers[0].AO_list[0].coeffI[0]);
#endif


        /* now update the locations */
        for(l=0;l<solid->num_atoms_in_cell;l++){
          l = num_added % MO_surf->num_centers_in_cell;
          temp_centers[num_added].AO_list = (AO_list_type *)
            D_CALLOC(MO_surf->raw_MO_centers[l].num_AOs,
                     sizeof(AO_list_type));
          if( !temp_centers[num_added].AO_list )
            fatal("can't allocate new AO_list_type\n");

          bcopy((char *)MO_surf->raw_MO_centers[l].AO_list,
                (char *)temp_centers[num_added].AO_list,
                MO_surf->raw_MO_centers[l].num_AOs*sizeof(AO_list_type));


          /* we copied over some pointers too... this is bad */
          temp_atoms[num_added].linesto = 0;
#ifdef INCLUDE_ADF_PLOTS
          temp_atoms[num_added].displacements = 0;
#endif
          temp_atoms[num_added].loc.x =
            solid->atoms[l].loc.x + dist_a.x + dist_b.x + dist_c.x;
          temp_atoms[num_added].loc.y =
            solid->atoms[l].loc.y + dist_a.y + dist_b.y + dist_c.y;
          temp_atoms[num_added].loc.z =
            solid->atoms[l].loc.z + dist_a.z + dist_b.z + dist_c.z;
          temp_atoms[num_added].num = num_added;

          temp_centers[num_added].loc = &(temp_atoms[num_added].loc);

          for(MO=0;MO<MO_surf->num_MOs;MO++){
            /* figure out the phase */
            kdotr = MO_surf->kpoints[MO].x * i +
              MO_surf->kpoints[MO].y * j +
              MO_surf->kpoints[MO].z * k;
            cos_kdotr = cos(2.0*kdotr*PI);
            sin_kdotr = sin(2.0*kdotr*PI);

#ifdef DEBUG
            printf("(i,j,k): (%d %d %d), num: %d kdotr: %lf, cos: %lf, sin: %lf\n",
                   i,j,k,num_added,kdotr,cos_kdotr,sin_kdotr);
#endif
            /* update the phase of the wavefunction contributions */
            for(m=0;m<temp_centers[num_added].num_AOs;m++){
#ifdef DEBUG
              fprintf(stderr,"C: %lf, Ci:% lf \n",
                      temp_centers[num_added].AO_list[m].coeff[MO],
                      temp_centers[num_added].AO_list[m].coeffI[MO]);
#endif
              temp_centers[num_added].AO_list[m].coeff[MO]
                *= cos_kdotr;
              temp_centers[num_added].AO_list[m].coeffI[MO]
                *= sin_kdotr;
#ifdef DEBUG
              fprintf(stderr,"\tC: %lf, Ci:% lf \n",
                      temp_centers[num_added].AO_list[m].coeff[MO],
                      temp_centers[num_added].AO_list[m].coeffI[MO]);
#endif
            }
          }
#ifdef DEBUG
          for(n=0;j<temp_centers[n].num_AOs;n++){
            printf("\t% -6.4lf ",temp_centers[num_added].AO_list[n].coeff[MO_surf->active_MO]);
            if( !((n+1)%10) ) printf("\n");
          }
          printf("\n");
#endif

          num_added++;
        }
      }
    }
  }

#ifdef DEBUG
  printf("num_added: %d, num_centers: %d\n",num_added,MO_surf->num_centers);
  printf("\n-----\n");
  for(i=0;i<num_added;i++){
    for(j=0;j<temp_centers[i].num_AOs;j++){
      printf("% -6.4lf ",temp_centers[i].AO_list[j].coeff[MO_surf->active_MO]);
      if( !((j+1)%10) ) printf("\n");
    }
  }
  printf("\n");
#endif

  /* free up the atoms we currently have stored */
  for(i=0;i<old_num_atoms;i++){
    if( solid->atoms[i].p_surf ){
      D_FREE(solid->atoms[i].p_surf);
    }
    if( solid->atoms[i].linesto ){
      D_FREE(solid->atoms[i].linesto);
    }
  }
  D_FREE(solid->atoms);
  solid->atoms = temp_atoms;
  MO_surf->MO_centers = temp_centers;

#ifdef DEBUG
  printf("\n-----\n");
  for(i=0;i<MO_surf->num_centers;i++){
    for(j=0;j<MO_surf->MO_centers[i].num_AOs;j++){
      printf("% -6.4lf ",MO_surf->MO_centers[i].AO_list[j].coeff[MO_surf->active_MO]);
      if( !((j+1)%10) ) printf("\n");
    }
    printf("\n");
  }
#endif

  /* redetermine the bond locations */
  determine_connections(solid);

  /* redetermine the bounding box */
  determine_mol_bounds(MO_surf->molec,&(MO_surf->bmin),
                       &(MO_surf->bmax));

  /* update the box */
  if( num_a > 1 ) num_a--;
  if( num_b > 1 ) num_b--;
  if( num_c > 1 ) num_c--;
  solid->cell_box[0].x = solid->lattice_vect[0].x;
  solid->cell_box[0].y = solid->lattice_vect[0].y;
  solid->cell_box[0].z = solid->lattice_vect[0].z;
  solid->cell_box[1].x =
    num_a*solid->orig_lattice[1].x;
  solid->cell_box[1].y =
    num_a*solid->orig_lattice[1].y;
  solid->cell_box[1].z =
    num_a*solid->orig_lattice[1].z;
  solid->cell_box[3].x =
    num_b*solid->orig_lattice[2].x;
  solid->cell_box[3].y =
    num_b*solid->orig_lattice[2].y;
  solid->cell_box[3].z =
    num_b*solid->orig_lattice[2].z;
  solid->cell_box[4].x =
    num_b*solid->orig_lattice[3].x;
  solid->cell_box[4].y =
    num_b*solid->orig_lattice[3].y;
  solid->cell_box[4].z =
    num_b*solid->orig_lattice[3].z;
  V3Add(&solid->cell_box[1],&solid->cell_box[3],
        &solid->cell_box[2]);
  V3Add(&solid->cell_box[1],&solid->cell_box[4],
        &solid->cell_box[5]);
  V3Add(&solid->cell_box[2],&solid->cell_box[4],
        &solid->cell_box[6]);
  V3Add(&solid->cell_box[3],&solid->cell_box[4],
        &solid->cell_box[7]);
  for(i=1;i<8;i++){
    V3Add(&solid->cell_box[i],&solid->cell_box[0],
          &solid->cell_box[i]);
  }


  display("There ya go!");
}


