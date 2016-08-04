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
*     this file contains stuff for dealing with the distance matrix
*
*  created:  greg landrum  August 1993
*
*  Revision History:
*     27 January 2002: fixed dumb bug in the distance matrix printing code
*
*****************************************************************************/
#include "bind.h"

/* distances less than this trigger warnings */
#define TOO_SHORT 1.0

/****************************************************************************
*
*                   Procedure check_a_cell
*
* Arguments:  atoms: pointer to atom_type
*              vect: point_type
*         num_atoms: int
* closest_nn_contact: real
*        descriptor: pointer to char
*
* Returns: none
*
* Action:
*        Goes through 'atoms and checks the distances between
*   the unmoved atoms and the atoms in the cell translated by
*   'vect.  any distances less than 'closest_nn_contact are printed
*   out in the output file.
*
*****************************************************************************/
void check_a_cell(atoms,vect,num_atoms,closest_nn_contact,descriptor)
  atom_type *atoms;
  point_type vect;
  int num_atoms;
  real closest_nn_contact;
  char *descriptor;
{
  int i,j;
  real min_squared;
  real dist,temp;

  /* use the squared distance to avoid sqrts */
  min_squared = closest_nn_contact * closest_nn_contact;

  /* loop over all the atoms */
  for(i=0;i<num_atoms;i++){
    for(j=0;j<num_atoms;j++){
      temp = atoms[i].loc.x - (atoms[j].loc.x + vect.x);
      dist = temp*temp;
      temp = atoms[i].loc.y - (atoms[j].loc.y + vect.y);
      dist += temp*temp;
      temp = atoms[i].loc.z - (atoms[j].loc.z + vect.z);
      dist += temp*temp;

      /* is the distance within the tolerance? */
      if( dist <= min_squared ){
        fprintf(output_file,"Atom: %d  Atom: %d in cell %s  Distance %6.4lf A\n",i+1,j+1,
                descriptor,sqrt(dist));
      }

      /* check to see if this is close enough to trigger a warning */
      if( dist < TOO_SHORT && atoms[i].at_number > 0 && atoms[j].at_number > 0){
        fprintf(stderr,"!!! Warning !!! Distance between atoms %d and %d in cell %s\
(%6.4lf A) is suspicious.\n",i+1,j+1,descriptor,dist);
      }
    }
  }
}

/****************************************************************************
*
*                   Procedure check_nn_contacts
*
* Arguments:  details: pointer to detail_type
*                cell: pointer to cell_type
*
* Returns: none
*
* Action:
*     Goes through the unit cell and checks for contacts with nearest neighbors
*      that are less than 'details->close_nn_contact.
*      These contacts are printed to the output file.
*
*****************************************************************************/
void check_nn_contacts(cell,details)
  cell_type *cell;
  detail_type *details;
{
  int num_atoms;
  int i,j;
  int itab,jtab;
  point_type vect,cell_dim[3];
  char *symbols;
  int num_so_far;
  real dist,max_dist=-1;

  num_atoms = cell->num_atoms;

  /* primitive error checking */
  if( cell->dim <= 0 ){
    FATAL_BUG("check_nn_contacts called for non-extended system.");
  }

  /* find the dimensions of the unit cell */
  for(i=0;i<cell->dim;i++){
    itab = cell->tvects[i].begin;
    jtab = cell->tvects[i].end;
    cell_dim[i].x = cell->atoms[jtab].loc.x-cell->atoms[itab].loc.x;
    cell_dim[i].y = cell->atoms[jtab].loc.y-cell->atoms[itab].loc.y;
    cell_dim[i].z = cell->atoms[jtab].loc.z-cell->atoms[itab].loc.z;
  }

  /*******

    now do the nearest neighbor cells

  *******/

  /* 1 0 0 */
  vect.x = cell_dim[0].x; vect.y = cell_dim[0].y; vect.z = cell_dim[0].z;
  check_a_cell(cell->atoms,vect,num_atoms,details->close_nn_contact,"(1 0 0)");

  if( cell->dim != 1 ){
    /* 1 1 0 */
    vect.x = cell_dim[0].x + cell_dim[1].x;
    vect.y = cell_dim[0].y + cell_dim[1].y;
    vect.z = cell_dim[0].z + cell_dim[1].z;
    check_a_cell(cell->atoms,vect,num_atoms,details->close_nn_contact,"(1 1 0)");

    /* 0 1 0 */
    vect.x = cell_dim[1].x;
    vect.y = cell_dim[1].y;
    vect.z = cell_dim[1].z;
    check_a_cell(cell->atoms,vect,num_atoms,details->close_nn_contact,"(0 1 0)");

    /* -1 1 0 */
    vect.x = -cell_dim[0].x + cell_dim[1].x;
    vect.y = -cell_dim[0].y + cell_dim[1].y;
    vect.z = -cell_dim[0].z + cell_dim[1].z;
    check_a_cell(cell->atoms,vect,num_atoms,details->close_nn_contact,"(-1 1 0)");

    if( cell->dim != 2 ){
      /* 0 0 1 */
      vect.x = cell_dim[2].x;
      vect.y = cell_dim[2].y;
      vect.z = cell_dim[2].z;
      check_a_cell(cell->atoms,vect,num_atoms,details->close_nn_contact,"(0 0 1)");

      /* 1 1 1 */
      vect.x = cell_dim[0].x + cell_dim[1].x + cell_dim[2].x;
      vect.y = cell_dim[0].y + cell_dim[1].y + cell_dim[2].y;
      vect.z = cell_dim[0].z + cell_dim[1].z + cell_dim[2].z;
      check_a_cell(cell->atoms,vect,num_atoms,details->close_nn_contact,"(1 1 1)");

      /* 0 1 1 */
      vect.x = cell_dim[1].x + cell_dim[2].x;
      vect.y = cell_dim[1].y + cell_dim[2].y;
      vect.z = cell_dim[1].z + cell_dim[2].z;
      check_a_cell(cell->atoms,vect,num_atoms,details->close_nn_contact,"(0 1 1)");

      /* -1 1 1 */
      vect.x = -cell_dim[0].x + cell_dim[1].x + cell_dim[2].x;
      vect.y = -cell_dim[0].y + cell_dim[1].y + cell_dim[2].y;
      vect.z = -cell_dim[0].z + cell_dim[1].z + cell_dim[2].z;
      check_a_cell(cell->atoms,vect,num_atoms,details->close_nn_contact,"(-1 1 1)");

      /* 1 0 1 */
      vect.x = cell_dim[0].x + cell_dim[2].x;
      vect.y = cell_dim[0].y + cell_dim[2].y;
      vect.z = cell_dim[0].z + cell_dim[2].z;
      check_a_cell(cell->atoms,vect,num_atoms,details->close_nn_contact,"(1 0 1)");

      /* -1 0 1 */
      vect.x = -cell_dim[0].x + cell_dim[2].x;
      vect.y = -cell_dim[0].y + cell_dim[2].y;
      vect.z = -cell_dim[0].z + cell_dim[2].z;
      check_a_cell(cell->atoms,vect,num_atoms,details->close_nn_contact,"(-1 0 1)");

      /* 1 -1 1 */
      vect.x = cell_dim[0].x - cell_dim[1].x + cell_dim[2].x;
      vect.y = cell_dim[0].y - cell_dim[1].y + cell_dim[2].y;
      vect.z = cell_dim[0].z - cell_dim[1].z + cell_dim[2].z;
      check_a_cell(cell->atoms,vect,num_atoms,details->close_nn_contact,"(1 -1 1)");

      /* 0 -1 1 */
      vect.x = -cell_dim[1].x + cell_dim[2].x;
      vect.y = -cell_dim[1].y + cell_dim[2].y;
      vect.z = -cell_dim[1].z + cell_dim[2].z;
      check_a_cell(cell->atoms,vect,num_atoms,details->close_nn_contact,"(0 -1 1)");

      /* -1 -1 1 */
      vect.x = -cell_dim[0].x - cell_dim[1].x + cell_dim[2].x;
      vect.y = -cell_dim[0].y - cell_dim[1].y + cell_dim[2].y;
      vect.z = -cell_dim[0].z - cell_dim[1].z + cell_dim[2].z;
      check_a_cell(cell->atoms,vect,num_atoms,details->close_nn_contact,"(-1 -1 1)");
    }
  }
}



/****************************************************************************
*
*                   Procedure build_distance_matrix
*
* Arguments:  cell: pointer to cell_type
*          details: pointer to detail_type
*
* Returns: none
*
* Action:
*     generates the distance matrix for the unit cell.
*
*  see the file notes.outl for the representation of symmetric matrices.
*
*****************************************************************************/
void build_distance_matrix(cell,details)
  cell_type *cell;
  detail_type *details;
{
  int num_atoms;
  int i,j;
  point_type temp;
  char *symbols;
  int num_so_far;
  real dist,max_dist=-1;

  num_atoms = cell->num_atoms;

  /* space for the list of symbols */
  symbols = (char *)calloc(num_atoms*4,sizeof(char));
  if(!symbols) fatal("Can't get space for symbols");

  /********
    get space for the distance matrix
      (we only need half of the matrix, since it's symmetrical)
  ********/
  if( !cell->distance_mat ){
    cell->distance_mat = (real *)calloc((num_atoms*num_atoms)/2+num_atoms,
                                         sizeof(real));
    if(!cell->distance_mat) fatal("Can't allocate distance matrix.");
  }

  num_so_far = 0;
  for(i=0;i<num_atoms;i++){
    for(j=0;j<i;j++){
      temp.x = cell->atoms[i].loc.x - cell->atoms[j].loc.x;
      temp.y = cell->atoms[i].loc.y - cell->atoms[j].loc.y;
      temp.z = cell->atoms[i].loc.z - cell->atoms[j].loc.z;

      dist = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);

      cell->distance_mat[num_so_far++] = dist;

      /* check to see if the distance is too short */
      if( dist < TOO_SHORT && cell->atoms[i].at_number > 0 && cell->atoms[j].at_number > 0){
        fprintf(stderr,"!!! Warning !!! Distance between atoms %d and %d (%f A) \
is suspicious.\n",i+1,j+1,dist);
      }
      if( dist > max_dist ) max_dist = dist;

    }

    /* put in the diagonal element */
    cell->distance_mat[num_so_far++] = 0.0;

    /* copy the symbol into the list of symbols */
    symbols[i*4] = cell->atoms[i].symb[0];
    symbols[i*4+1] = cell->atoms[i].symb[1];
    symbols[i*4+2] = 0;
  }

  if( details->distance_mat_PRT ){
    print_sym_mat(cell->distance_mat,num_atoms,num_atoms,output_file,
                  "\n\n;****** DISTANCE MATRIX *********",symbols,details->line_width);

    /* check close contacts to nearest neighbors */
    if( cell->dim > 0 ){
      fprintf(output_file,"# Inter-Cell distances less than %6.4lf Angstroms:\n",
              details->close_nn_contact);
      check_nn_contacts(cell,details);

    }
  }
  free(symbols);

  if(details->dump_dist_mat) dump_distance_mats(cell,details);

}


/****************************************************************************
*
*                   Procedure display_lattice_parms
*
* Arguments:  cell: pointer to cell_type
*
* Returns: none
*
* Action:
*     dumps the lattice parameters to the output file and calculates the
*     reciprocal lattice vectors
*
*****************************************************************************/
void display_lattice_parms(cell)
  cell_type *cell;
{
  static char first_call=0;

  /*******

    if this is the first time that this was called, then
    subtract the dimension of the xtal from the number of atoms in the unit
    cell to account for the atoms included to determine translation vectors

  ********/
  if( !first_call ){
  }

  fprintf(output_file,"\n; ------  Lattice Parameters ------\n");
  fprintf(output_file,"#Dimensionality: %d\n",cell->dim);
  fprintf(output_file,"#Lattice Vectors\n");
  fprintf(output_file,"(a) %lf %lf %lf\n",
          cell->atoms[cell->tvects[0].end].loc.x-cell->atoms[cell->tvects[0].begin].loc.x,
          cell->atoms[cell->tvects[0].end].loc.y-cell->atoms[cell->tvects[0].begin].loc.y,
          cell->atoms[cell->tvects[0].end].loc.z-cell->atoms[cell->tvects[0].begin].loc.z);

  if( cell->dim > 1 ){
    fprintf(output_file,"(b) %lf %lf %lf\n",
            cell->atoms[cell->tvects[1].end].loc.x-
            cell->atoms[cell->tvects[1].begin].loc.x,
            cell->atoms[cell->tvects[1].end].loc.y-
            cell->atoms[cell->tvects[1].begin].loc.y,
            cell->atoms[cell->tvects[1].end].loc.z-
            cell->atoms[cell->tvects[1].begin].loc.z);

  }
  if( cell->dim > 2 ){
    fprintf(output_file,"(c) %lf %lf %lf\n",
            cell->atoms[cell->tvects[2].end].loc.x-
            cell->atoms[cell->tvects[2].begin].loc.x,
            cell->atoms[cell->tvects[2].end].loc.y-
            cell->atoms[cell->tvects[2].begin].loc.y,
            cell->atoms[cell->tvects[2].end].loc.z-
            cell->atoms[cell->tvects[2].begin].loc.z);
  }
  if( cell->dim >= 1 ){
    calc_reciprocal_lattice(cell);
  }

}






/****************************************************************************
*
*                   Procedure fill_distance_mat
*
* Arguments:  cell: pointer to cell_type
*        num_atoms: int
*         dist_mat: pointer to float
*             vect: pointer to point_type
*
*
* Returns: none
*
* Action:
*     generates the distance matrix for the cell at 'point
*
*
*****************************************************************************/
void fill_distance_matrix(cell_type *cell,int num_atoms,float *dist_mat,
                          point_type *vect)
{
  int num_so_far;
  point_type temp;
  real dist;
  int i,j;

  num_so_far = 0;
  for(i=0;i<num_atoms;i++){
    for(j=0;j<num_atoms;j++){
      temp.x = cell->atoms[i].loc.x -
        (cell->atoms[j].loc.x+vect->x);
      temp.y = cell->atoms[i].loc.y -
        (cell->atoms[j].loc.y+vect->y);
      temp.z = cell->atoms[i].loc.z -
        (cell->atoms[j].loc.z+vect->z);

      dist = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);

      dist_mat[num_so_far++] = (float)dist;
    }
  }
}

typedef struct{
  int x,y,z;
} cell_spec_type;


/****************************************************************************
*
*                   Procedure dump_distance_mats
*
* Arguments:  cell: pointer to cell_type
*          details: pointer to detail_type
*
* Returns: none
*
* Action:
*     generates the distance matrix for the unit cell and dumps
*      it to a binary output file.
*
*****************************************************************************/
void dump_distance_mats(cell,details)
  cell_type *cell;
  detail_type *details;
{
  int num_atoms;
  int i,j,itab,jtab;
  point_type cell_dim[3];
  point_type temp,vect;
  char tempfilename[240];
  int matfile;
  int num_so_far;
  real dist,max_dist=-1;
  float *dist_mat;
  cell_spec_type cell_vect;
  char *dummies;

  num_atoms = cell->num_atoms;

  for(i=0;i<cell->dim;i++){
    itab = cell->tvects[i].begin;
    jtab = cell->tvects[i].end;
    cell_dim[i].x = cell->atoms[jtab].loc.x-cell->atoms[itab].loc.x;
    cell_dim[i].y = cell->atoms[jtab].loc.y-cell->atoms[itab].loc.y;
    cell_dim[i].z = cell->atoms[jtab].loc.z-cell->atoms[itab].loc.z;
  }


  /* open the file */
  sprintf(tempfilename,"%s.DMAT",details->filename);
#ifndef USING_THE_MAC
  matfile = open(tempfilename,O_RDWR|O_APPEND|O_CREAT|O_TRUNC,S_IRUSR|S_IWUSR);
#else
  matfile = open(tempfilename,O_RDWR|O_TRUNC|O_APPEND|O_CREAT);
#endif
  if( matfile == -1 ){
    error("Can't open DMAT file for binary I/O");
    return;
  }

  /*****

    write the dimensionality of the system
    and number of atoms to the output file

  ******/
  write(matfile,(const char *)&(cell->dim),sizeof(int));
  write(matfile,(const char *)&(num_atoms),sizeof(int));


  /* get space for the dummies array */
  dummies = (char *)calloc(num_atoms,sizeof(char));
  if(!dummies) fatal("Can't allocate dummies array.");
  for(i=0;i<num_atoms;i++){
    if( cell->atoms[i].at_number == -1 ) dummies[i] = 1;
  }
  /* write the dummies array */
  write(matfile,(const char *)(dummies),num_atoms*sizeof(char));

  /********

    get space for the distance matrix

  ********/
  dist_mat = (float *)calloc((num_atoms*num_atoms),sizeof(float));
  if(!dist_mat) fatal("Can't allocate distance matrix.");

  num_so_far = 0;
  for(i=0;i<num_atoms;i++){
    for(j=0;j<num_atoms;j++){
      temp.x = cell->atoms[i].loc.x - cell->atoms[j].loc.x;
      temp.y = cell->atoms[i].loc.y - cell->atoms[j].loc.y;
      temp.z = cell->atoms[i].loc.z - cell->atoms[j].loc.z;

      dist = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);

      dist_mat[num_so_far++] = (float)dist;
    }
  }

  /* dump the distance matrix */
  write(matfile,(const char *)dist_mat,num_atoms*num_atoms*sizeof(float));

  /* do the nearest neighbors */
  if( cell->dim > 0 ){

    /* 1 0 0 */
    vect.x = cell_dim[0].x; vect.y = cell_dim[0].y; vect.z = cell_dim[0].z;
    fill_distance_matrix(cell,num_atoms,dist_mat,&vect);
    cell_vect.x=1;cell_vect.y=0;cell_vect.z=0;
    write(matfile,(const char *)(&cell_vect),sizeof(cell_spec_type));
    write(matfile,(const char *)dist_mat,
          num_atoms*num_atoms*sizeof(float));

    if( cell->dim != 1 ){
      /* 1 1 0 */
      vect.x = cell_dim[0].x + cell_dim[1].x;
      vect.y = cell_dim[0].y + cell_dim[1].y;
      vect.z = cell_dim[0].z + cell_dim[1].z;
      fill_distance_matrix(cell,num_atoms,dist_mat,&vect);
      cell_vect.x=1;cell_vect.y=1;cell_vect.z=0;
      write(matfile,(const char *)(&cell_vect),sizeof(cell_spec_type));
      write(matfile,(const char *)dist_mat,
            num_atoms*num_atoms*sizeof(float));

      /* 0 1 0 */
      vect.x = cell_dim[1].x;
      vect.y = cell_dim[1].y;
      vect.z = cell_dim[1].z;
      fill_distance_matrix(cell,num_atoms,dist_mat,&vect);
      cell_vect.x=0;cell_vect.y=1;cell_vect.z=0;
      write(matfile,(const char *)(&cell_vect),sizeof(cell_spec_type));
      write(matfile,(const char *)dist_mat,
            num_atoms*num_atoms*sizeof(float));

      /* -1 1 0 */
      vect.x = -cell_dim[0].x + cell_dim[1].x;
      vect.y = -cell_dim[0].y + cell_dim[1].y;
      vect.z = -cell_dim[0].z + cell_dim[1].z;
      fill_distance_matrix(cell,num_atoms,dist_mat,&vect);
      cell_vect.x=-1;cell_vect.y=1;cell_vect.z=0;
      write(matfile,(const char *)(&cell_vect),sizeof(cell_spec_type));
      write(matfile,(const char *)dist_mat,
            num_atoms*num_atoms*sizeof(float));

      if( cell->dim != 2 ){
        /* 0 0 1 */
        vect.x = cell_dim[2].x;
        vect.y = cell_dim[2].y;
        vect.z = cell_dim[2].z;
        fill_distance_matrix(cell,num_atoms,dist_mat,&vect);
        cell_vect.x=0;cell_vect.y=0;cell_vect.z=1;
        write(matfile,(const char *)(&cell_vect),sizeof(cell_spec_type));
        write(matfile,(const char *)dist_mat,
              num_atoms*num_atoms*sizeof(float));

        /* 1 1 1 */
        vect.x = cell_dim[0].x + cell_dim[1].x + cell_dim[2].x;
        vect.y = cell_dim[0].y + cell_dim[1].y + cell_dim[2].y;
        vect.z = cell_dim[0].z + cell_dim[1].z + cell_dim[2].z;
        fill_distance_matrix(cell,num_atoms,dist_mat,&vect);
        cell_vect.x=1;cell_vect.y=1;cell_vect.z=1;
        write(matfile,(const char *)(&cell_vect),sizeof(cell_spec_type));
        write(matfile,(const char *)dist_mat,
              num_atoms*num_atoms*sizeof(float));


        /* 0 1 1 */
        vect.x = cell_dim[1].x + cell_dim[2].x;
        vect.y = cell_dim[1].y + cell_dim[2].y;
        vect.z = cell_dim[1].z + cell_dim[2].z;
        fill_distance_matrix(cell,num_atoms,dist_mat,&vect);
        cell_vect.x=0;cell_vect.y=1;cell_vect.z=1;
        write(matfile,(const char *)(&cell_vect),sizeof(cell_spec_type));
        write(matfile,(const char *)dist_mat,
              num_atoms*num_atoms*sizeof(float));

        /* -1 1 0 */
        vect.x = -cell_dim[0].x + cell_dim[1].x + cell_dim[2].x;
        vect.y = -cell_dim[0].y + cell_dim[1].y + cell_dim[2].x;
        vect.z = -cell_dim[0].z + cell_dim[1].z + cell_dim[2].x;
        fill_distance_matrix(cell,num_atoms,dist_mat,&vect);
        cell_vect.x=-1;cell_vect.y=1;cell_vect.z=0;
        write(matfile,(const char *)(&cell_vect),sizeof(cell_spec_type));
        write(matfile,(const char *)dist_mat,
              num_atoms*num_atoms*sizeof(float));

        /* 1 0 1 */
        vect.x = cell_dim[0].x + cell_dim[2].x;
        vect.y = cell_dim[0].y + cell_dim[2].y;
        vect.z = cell_dim[0].z + cell_dim[2].z;
        fill_distance_matrix(cell,num_atoms,dist_mat,&vect);
        cell_vect.x=1;cell_vect.y=0;cell_vect.z=1;
        write(matfile,(const char *)(&cell_vect),sizeof(cell_spec_type));
        write(matfile,(const char *)dist_mat,
              num_atoms*num_atoms*sizeof(float));

        /* -1 0 1 */
        vect.x = -cell_dim[0].x + cell_dim[2].x;
        vect.y = -cell_dim[0].y + cell_dim[2].y;
        vect.z = -cell_dim[0].z + cell_dim[2].z;
        fill_distance_matrix(cell,num_atoms,dist_mat,&vect);
        cell_vect.x=-1;cell_vect.y=0;cell_vect.z=1;
        write(matfile,(const char *)(&cell_vect),sizeof(cell_spec_type));
        write(matfile,(const char *)dist_mat,
              num_atoms*num_atoms*sizeof(float));

        /* 1 -1 1 */
        vect.x = cell_dim[0].x - cell_dim[1].x + cell_dim[2].x;
        vect.y = cell_dim[0].y - cell_dim[1].y + cell_dim[2].y;
        vect.z = cell_dim[0].z - cell_dim[1].z + cell_dim[2].z;
        fill_distance_matrix(cell,num_atoms,dist_mat,&vect);
        cell_vect.x=1;cell_vect.y=-1;cell_vect.z=1;
        write(matfile,(const char *)(&cell_vect),sizeof(cell_spec_type));
        write(matfile,(const char *)dist_mat,
              num_atoms*num_atoms*sizeof(float));

        /* 0 -1 1 */
        vect.x = -cell_dim[1].x + cell_dim[2].x;
        vect.y = -cell_dim[1].y + cell_dim[2].y;
        vect.z = -cell_dim[1].z + cell_dim[2].z;
        fill_distance_matrix(cell,num_atoms,dist_mat,&vect);
        cell_vect.x=0;cell_vect.y=-1;cell_vect.z=1;
        write(matfile,(const char *)(&cell_vect),sizeof(cell_spec_type));
        write(matfile,(const char *)dist_mat,
              num_atoms*num_atoms*sizeof(float));


        /* -1 -1 0 */
        vect.x = -cell_dim[0].x - cell_dim[1].x + cell_dim[2].x;
        vect.y = -cell_dim[0].y - cell_dim[1].y + cell_dim[2].x;
        vect.z = -cell_dim[0].z - cell_dim[1].z + cell_dim[2].x;
        fill_distance_matrix(cell,num_atoms,dist_mat,&vect);
        cell_vect.x=-1;cell_vect.y=-1;cell_vect.z=0;
        write(matfile,(const char *)(&cell_vect),sizeof(cell_spec_type));
        write(matfile,(const char *)dist_mat,
              num_atoms*num_atoms*sizeof(float));

      }
    }
  }

  free(dist_mat);
  close(matfile);
}


