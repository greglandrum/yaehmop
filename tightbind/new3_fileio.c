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
*     this file contains stuff for dealing with file input from NEW3's twisted
*      input files.
*
*  created:  greg landrum  January 1993
*
*****************************************************************************/
#include "bind.h"


/****************************************************************************
*
*                   Procedure read_NEW3file
*
* Arguments:  cell: pointer to cell_type
*          details: pointer to detail_type
*             name: pointer to type char
*
* Returns: none
*
* Action: reads all the data out of the file 'infile
*
*****************************************************************************/
void read_NEW3file(cell,details,infile)
  cell_type *cell;
  detail_type *details;
  FILE *infile;
{
  char err_string[240];
  char instring[90];
  char foo_string[80];
  k_point_type *points;
  int num_k_points,max_k_points;

  int idle;
  int max_p_DOS;
  int use_gradients,avg_props;
  int which;
  real weight;
  real ecut,eerr;
  char found;
  int i,j;
  int temp;

  p_DOS_type *p_DOS,*temp_p_dos;
#ifdef USING_THE_MAC
  fprintf(stderr,"That doesn't look like a bind input file.  This version\n \
  of bind makes no attempt to read in new3 input files, so we're going to\n \
  error out now and let you try another input file\n");
  fatal("Bogus input file\n");

#else

  /* if we made it this far, the title of the file has already been read in.... */

  /* put some warnings in the status file */
  fprintf(status_file,"\n\n!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!\n");
  fprintf(status_file,"You are using a NEW3 style input file.\n");
  fprintf(status_file," Be sure that you put spaces between all numbers.\n");
  fprintf(status_file," Printing options are ignored, minimal printing is done.\n");
  fprintf(status_file," If you want more control over the calculation, please switch\n");
  fprintf(status_file,"  to using the newer input format.\n\n");

  skipcomments(infile,instring,FATAL);
  sscanf(instring,"%lf %d %d %d",&(details->rho),&(cell->num_atoms),
         &(details->just_geom),&(details->save_energies));

  /* allocate space for the atoms */
  cell->atoms = (atom_type *)calloc(cell->num_atoms,sizeof(atom_type));
  if(!cell->atoms){
    sprintf(err_string,"Can't allocate memory for: %d atoms.",cell->num_atoms);
    fatal("Can't allocate memory for the atoms.");
  }

  /* read in the individual atomic locations */
  for(i=0;i<cell->num_atoms;i++){
    skipcomments(infile,instring,FATAL);
    sscanf(instring,"%c%c  %lf %lf %lf",
           &(cell->atoms[i].symb[0]),&(cell->atoms[i].symb[1]),
           &(cell->atoms[i].loc.x),&(cell->atoms[i].loc.y),
           &(cell->atoms[i].loc.z));
  }
  cell->atoms[i].symb[2] = 0;

  fprintf(status_file,"Read: %d atoms\n",cell->num_atoms);

  /* more details */
  skipcomments(infile,instring,FATAL);
  sscanf(instring,"%d %d %d %d %d %d %d %d %d %d",&(cell->dim),
         &(cell->overlaps[0]),&(cell->overlaps[1]),&(cell->overlaps[2]),
         &(cell->tvects[0].begin),&(cell->tvects[0].end),
         &(cell->tvects[1].begin),&(cell->tvects[1].end),
             &(cell->tvects[2].begin),&(cell->tvects[2].end));

  /* check to see if this is specifying a molecular calculation */
  if(unit_cell->dim == 0){
    details->Execution_Mode = MOLECULAR;
  }

  /*********
    now check to see if the order of any of the end points
    needs to be changed
  *********/
  for(j=0;j<3 && j<cell->dim;j++){
    /* first decrement the tabs since we index arrays from 0 */
    cell->tvects[j].end--;
    cell->tvects[j].begin--;

    if(cell->tvects[j].end < cell->tvects[j].begin){
      temp = cell->tvects[j].end;
      cell->tvects[j].end = cell->tvects[j].begin;
      cell->tvects[j].begin = temp;
    }
    else if(cell->tvects[j].end == cell->tvects[j].begin &&
            cell->tvects[j].end ){
      /* both ends of the vector are the same */
      fatal("Begin and end of a translation vector are the same.");
    }
  }

  /* fill in the atomic parameters */
  fill_atomic_parms(cell->atoms,cell->num_atoms,infile);
  write_atom_coords(cell->atoms,cell->num_atoms,cell->using_Zmat,
                    cell->using_xtal_coords);
  write_atom_parms(details,cell->atoms,cell->num_atoms,1);



  fprintf(status_file,"Completed parameter acquisition.\n");

  /* this doesn't deal with actually reading in multiple occupations, etc. */
  skipcomments(infile,instring,FATAL);
/*
  sscanf(instring,"%d %d %lf %d %d %d %d %d %d %d %lf %lf",
         &(details->lower_level_PRT),&(details->upper_level_PRT),
         &(cell->num_electrons),&(details->num_bonds_OOP),
         &use_gradients,&(details->num_FMO),
         &(details->num_frags_FMO),&(details->num_occup_AVG),
         &avg_props,&idle,&ecut,&eerr);
*/
  details->gradients = (char)use_gradients;
  details->avg_props = (char)avg_props;

  if( details->avg_props ){
    /* read in the number of projected DOS curves */
    skipcomments(infile,instring,FATAL);

    sscanf(instring,"%d %d",&which,&(details->num_proj_DOS));
  }
  /*****
    add stuff for  overlap populations here
  *****/
  skipcomments(infile,instring,FATAL);

  /* get memory for the projected DOS data */
  details->proj_DOS = (p_DOS_type *)calloc(details->num_proj_DOS,sizeof(p_DOS_type));
  if( !(details->proj_DOS) ) fatal("Can't allocate memory for projected DOS structure.");

  /* now read in the individual DOS projections */
  for(i=0;i<details->num_proj_DOS;i++){
    skipcomments(infile,instring,FATAL);

    /* first get memory for the contributions to this projection */
    max_p_DOS = 2;

    p_DOS = &(details->proj_DOS[i]);

    /******
      insert code for reading in multiple contributions here
    *******/
    p_DOS->num_contributions = 1;
    sscanf(instring,"%d %d %lf",&(p_DOS->type),&which,&weight);

    p_DOS->weights = (real *)calloc(p_DOS->num_contributions,sizeof(real));
    if( !p_DOS->weights ) fatal("Can't allocate memory for projected DOS weightings.");
    p_DOS->weights[0] = weight;
    p_DOS->contributions = (int *)calloc(p_DOS->num_contributions,sizeof(int));
    if( !p_DOS->contributions )fatal("Can't allocate memory to hold projected DOS \
contributions.");

    /*******
      subtract one from the specified number, since arrays are indexed from
      zero in C
    *******/
    p_DOS->contributions[0] = which-1;
  }

  /*******

    subtract the dimension of the xtal from the number of atoms in the unit cell to
    account for the atoms included to determine translation vectors

  ********/
  cell->num_atoms -= cell->dim;

  /********
    skip over the energy window information
  *********/
  if( details->avg_props ){
    skipcomments(infile,instring,FATAL);
  }

  /**********

    do the k-points (if this is an extended calculation)

  **********/
  if( cell->dim > 0 ){
    /* first get some memory */
    max_k_points = 256;
    num_k_points = 0;
    points = (k_point_type *)calloc(max_k_points,sizeof(k_point_type));
    if(!points)fatal("Can't allocate memory for k point set.");

    while(skipcomments(infile,instring,IGNORE)>-1){
      /* read in the K point, ignoring the printing options */
      sscanf(instring,"%s %lf %lf %lf %d %d",
             foo_string,
             &(points[num_k_points].loc.x),&(points[num_k_points].loc.y),
             &(points[num_k_points].loc.z),
             &(points[num_k_points].weight),&(points[num_k_points].num_filled_bands));
      num_k_points++;

      /* check to see if there's still enough memory */
      if( num_k_points == max_k_points ){
        /* not enough space, reallocate the array */
        max_k_points += 256;
        points=(k_point_type *)realloc((char *)points,
                                       (unsigned)max_k_points*sizeof(k_point_type));
        if( !points )fatal("Can't reallocate k point array to get more space.");
      }
    }
    fprintf(status_file,"Read: %d k points\n",num_k_points);

    /* set the pointer in the details structure */
    details->K_POINTS = points;
    details->num_KPOINTS = num_k_points;
  }
  else{
    fprintf(status_file,"This will be a molecular calculation.\n");

    /* read out the printing options */
    skipcomments(infile,instring,IGNORE);
    points = (k_point_type *)calloc(1,sizeof(k_point_type));
    if( !points ) fatal("Can't allocate space for printing options.");

    num_k_points = 0;
    sscanf(instring,"%s %lf %lf %lf %d %d",
           foo_string,
           &(points[num_k_points].loc.x),&(points[num_k_points].loc.y),
           &(points[num_k_points].loc.z),
           &(points[num_k_points].weight),&(points[num_k_points].num_filled_bands));

    details->num_KPOINTS = 1;
    details->K_POINTS = points;
  }
#endif
}
