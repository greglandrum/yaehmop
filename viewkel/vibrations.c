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

#ifdef INCLUDE_ADF_PLOTS

/********

  this has got the stuff for dealing with vibrations

*********/
#include "viewkel.h"

/***
  Recent Edit History

   26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)

***/
/****************************************************************************
 *
 *                   Procedure read_vibration_data
 *
 * Arguments:   infile: pointer to type FILE
 *              molec: a pointer to molecule_type
 *
 * Returns: none
 *
 * Action: This routine reads in the information for a vibration of a molecule from
 *    'infile.  The data is stored in 'molec
 *
 *  The file format assumed is simple...
 *
 ****************************************************************************/
void read_vibration_data(FILE *infile,molec_type *molec)
{
  char instring[MAX_STR_LEN];

  int i,j;


  if(!infile)FATAL_BUG("No file passed to read_vibration_data!");
  if(!molec)FATAL_BUG("No molecule passed to read_vibration_data!");

  /******

    find and read out the number of atoms, checking for Walsh info along the way

  *******/
  skipcomments(infile,instring);
  upcase(instring);
  while(!strstr(instring,"VIBRATIONS")){
    if(skipcomments(infile,instring)<0){
      error("End of File hit while looking for vibration data");
      return;
    }
    upcase(instring);
  }
  skipcomments(infile,instring);
  sscanf(instring,"%d",&molec->num_vibrations);

  /* allocate space to store the vibrations now */
  for(i=0;i<molec->num_atoms;i++){
    molec->atoms[i].displacements = (point_type *)D_CALLOC(molec->num_vibrations,sizeof(point_type));
    if(!molec->atoms[i].displacements) fatal("Can't allocate space for displacements");
  }

  /* read in the vibrations */
  for(i=0;i<molec->num_vibrations;i++){
    for(j=0;j<molec->num_atoms;j++){
      if(skipcomments(infile,instring)<0) fatal("EOF hit reading vibration data");
      sscanf(instring,"%lf %lf %lf",&(molec->atoms[j].displacements[i].x),
             &(molec->atoms[j].displacements[i].y),&(molec->atoms[j].displacements[i].z));
    }
  }

  /* that's that */
}





/****************************************************************************
 *
 *                   Procedure new_vibration
 *
 * Arguments: filename: pointer to type char
 *
 * Returns: none
 *
 * Action: does everything to get space for and read in a new vibration
 *
 ****************************************************************************/
void new_vibration(char *filename)
{
  char file_name[80];
  char *theinline;
  FILE *infile;

  /* set up a new object to hold the vibration */
  makenewobject();
  whichobj = head->obj;

  /* now build the molecule primitive */
  whichobj->prim = (prim_type *)D_CALLOC(1,sizeof(prim_type));
  if( !whichobj->prim )fatal("Can't get space for vibration primitive.");
  whichobj->prim->which = MOLECULE;

  whichobj->prim->molec = (molec_type *)D_CALLOC(1,sizeof(molec_type));
  if( !whichobj->prim->molec )
    fatal("Can't get space for molecule.");
#ifndef USING_THE_MAC
  if( !filename ){
    display("Look in the xterm...");
#ifndef USE_READLINE
    printf("Enter the file name containing the vibration data: ");
    scanf("%s",file_name);
#else
    theinline= readline("Enter the file name containing the vibration data: ");
    add_history(theinline);
    if( theinline ){
      sscanf(theinline,"%s",file_name);
      free(theinline);
    } else {
      error("Bad file name");
      file_name[0] = 0;
    }
#endif
  }else{
    strcpy(file_name,filename);
  }

  /* open the file */
  infile = fopen(file_name,"r");
  if(!infile){
    printf("Problems opening file: %s\n",file_name);
    display("oooooops!");
    return;
  }
#else
        if(!filename){
                infile = choose_mac_file(file_name,MAC_FOPEN_OPEN_CD);
        } else{
                strcpy(file_name,filename);
                infile = fopen(file_name,"r");
        }
        if( !infile ){
             printf("Problems opening file: %s\n",file_name);
           display("oooooops!");
             return;
    }
#endif

        read_molecule_data(infile,whichobj->prim->molec);


  /* check to see if any atoms were actually read in.... */
  if(!whichobj->prim->molec->num_atoms){
    /* no... free the memory that we asked for */
    D_FREE(whichobj->prim->molec);
    D_FREE(whichobj->prim);
    D_FREE(whichobj);
    whichobj=0;
    head = head->next;
  }
  else{
    /* okay, read in the vibrational modes now */
    rewind(infile);
    read_vibration_data(infile,whichobj->prim->molec);

    if( !whichobj->prim->molec->num_vibrations ){
      error("Something went wrong whilst reading vibrations");
      return;
    }
    whichobj->prim->molec->active_vibn = 1;

    whichobj->scale.x=whichobj->scale.y=whichobj->scale.z=0.5;
    whichobj->trans.x=0;whichobj->trans.y=0;
    whichobj->trans.z=0;
    whichobj->cent.x = g_xmax/2;
    whichobj->cent.y = g_ymax/2;

#ifdef INTERACTIVE_USE
    /* create the button window */
    build_molec_button_window(&button_wins,whichobj->prim->molec);
    whichobj->prim->but_win = button_wins;
#endif
    /* set up some default parameters */
    whichobj->prim->molec->draw_connectors = 1;
    whichobj->prim->molec->fancy_lines = 1;
    whichobj->prim->molec->outlines_on = 1;
    whichobj->prim->molec->shading_on = 1;
    whichobj->prim->molec->hydrogens_on = 1;
    whichobj->prim->molec->dummies_on = 1;
    whichobj->prim->molec->line_width = 1;
    whichobj->prim->molec->rad_mult = 1.0;
    whichobj->prim->molec->tubes_on = 1;
    whichobj->prim->molec->bond_rad = 0.05;
    whichobj->prim->molec->vibration_scale = 10.0;
    strcpy(whichobj->prim->molec->filename,file_name);

  }

}

#endif
