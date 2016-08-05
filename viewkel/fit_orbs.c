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

#include "viewkel.h"
#include "orb_defs.h"

/***
  Recent Edit History:
   28.02.98 gL:
     modifications to reduce memory usage for MO storage.
   26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)
   14.02.2004 gL:
     support degree 7 radial functions.


***/


#ifdef USE_V5D
#include <v5d.h>
#else
#ifdef USE_HDF
#include <hdf.h>
#endif
#endif
/*******************

  this has got the stuff for doing fitting of orbitals

  created by greg Landrum  Feb. 1995

*******************/



/****************************************************************************
 *
 *                   Function choose_MO
 *
 * Arguments:  surf: pointer to MO_surface_type
 *
 * Returns: int
 *
 * Action:  allows the user to choose a new MO from the available list.
 *   returns the number of the MO selected, or -1 on failure.
 *
 ****************************************************************************/
int choose_MO(MO_surface_type *surf)
{
  int i;
  int which_MO;

  printf("There are %d MO's in the MO file. They are: \n",surf->num_MOs);
  for(i=0;i<surf->num_MOs;i++){
    if( !surf->adf_plot ){
      printf("\t%d: MO: %d at k point: (%6.4lf, %6.4lf, %6.4lf)\n",i+1,
             surf->MO_numbers[i],surf->kpoints[i].x,surf->kpoints[i].y,
             surf->kpoints[i].z);
    }else{
      printf("\t%d: MO: %d\n",i+1,i+1);
    }
  }
  which_MO = 0;
  while(!which_MO){
    printf("Which number MO would you like (a negative number to cancel)? ");

    which_MO = surf->active_MO+1;
    /*    scanf("%d\n",&which_MO); */
    readintparm("active MO",&which_MO);
    if( which_MO > surf->num_MOs ){
      printf("%d is an invalid choice, there are %d MOs in the file\n",
             which_MO,surf->num_MOs);
      which_MO = 0;
    }
  }
  which_MO--;

  return which_MO;
}


/****************************************************************************
 *
 *                   Function exclude_atoms
 *
 * Arguments:
 *
 * Returns: none
 *
 * Action:  allows the user to exclude particular atoms from the
 *  calculation of the MO surface
 *
 ****************************************************************************/
void exclude_atoms(int num_args,char **MO_surf_ptr)
{
  char instring[MAX_STR_LEN];
  int *atoms_to_exclude;
  int num_to_exclude;
  int i;

  MO_surface_type *MO_surf;

  display("Look in Xterm");
  MO_surf = (MO_surface_type *)MO_surf_ptr[0];

  printf("Please enter the numbers of the atoms you wish to exclude:\n");
  fgets(instring,MAX_STR_LEN,stdin);

  parse_integer_string(instring,&atoms_to_exclude,&num_to_exclude);

  /*******

    now set the exclusion...

    we set it in the raw_MO_centers array in case this is an extended
    system and the user has grown the crystal

  ********/
  for(i=0;i<num_to_exclude;i++){
    MO_surf->raw_MO_centers[atoms_to_exclude[i]-1].exclude = 1;
  }

  printf("%d Atoms excluded.\n",num_to_exclude);
  if( MO_surf->molec->num_dim >= 1 ){
    printf("If you have grown this crystal, you will need to grow it again\n");
    printf("  to make these changes take effect.  Sorry\n");
  }
}



/****************************************************************************
 *
 *                   Function include_atoms
 *
 * Arguments:
 *
 * Returns: none
 *
 * Action:  allows the user to include particular atoms for the
 *  calculation of the MO surface.  There's no reason to do this
 *  if they haven't alread excluded some.
 *
 ****************************************************************************/
void include_atoms(int num_args,char **MO_surf_ptr)
{
  char instring[MAX_STR_LEN];
  int *atoms_to_include;
  int num_to_include;
  int i;

  MO_surface_type *MO_surf;

  display("Look in Xterm");
  MO_surf = (MO_surface_type *)MO_surf_ptr[0];

  printf("Please enter the numbers of the atoms you wish to include:\n");
  fgets(instring,MAX_STR_LEN,stdin);

  parse_integer_string(instring,&atoms_to_include,&num_to_include);

  /*******

    now set the exclusion...

    we set it in the raw_MO_centers array in case this is an extended
    system and the user has grown the crystal

  ********/
  for(i=0;i<num_to_include;i++){
    MO_surf->raw_MO_centers[atoms_to_include[i]-1].exclude = 0;
  }

  printf("%d Atoms included.\n",num_to_include);
  if( MO_surf->molec->num_dim >= 1 ){
    printf("If you have grown this crystal, you will need to grow it again\n");
    printf("  to make these changes take effect.  Sorry\n");
  }
}



/****************************************************************************
 *
 *                   Function change_active_MO
 *
 * Arguments:
 *
 * Returns: none
 *
 * Action:  allows the user to change which MO surface they are
 *   looking at.
 *
 ****************************************************************************/
void change_active_MO(int num_args,char **MO_surf_ptr)
{
  MO_surface_type *MO_surf;
  int which_MO;

  display("Look in Xterm");
  MO_surf = (MO_surface_type *)MO_surf_ptr[0];
  which_MO = choose_MO(MO_surf);
  if( which_MO >= 0 ){
    MO_surf->active_MO = which_MO;
    display("Active MO changed!");
  }else{
    display("No change to active MO");
  }
}


/****************************************************************************
 *
 *                   Function read_MO_surface
 *
 * Arguments: infile: pointer to type FILE
 *              surf: pointer to MO_surface_type
 *
 * Returns: int
 *
 * Action:  Reads the information for an MO isosurface out of 'infile
 *     and puts it all into 'surf.
 *
 *   If there is a problem reading the file, then a nonzero number is
 *    returned to indicate failure.
 *
 ****************************************************************************/
int read_MO_surface(FILE *infile,MO_surface_type *surf)
{
  int *MO_numbers;
  point_type *MO_kpoints;
  char instring[MAX_STR_LEN],foostring[MAX_STR_LEN];
  char type[4];
  int i,j,k;
  int max_unique=5;
  int num_AOs,num_MOs;
  int AO_accum;
  int which_MO;

  int ns,np,nd,nf;
  float contraction;
  float exps,expp,expd1,expd2,expf1,expf2;
  float c_d1,c_d2,c_f1,c_f2;
  float Hs,Hp,Hd,Hf;

  int color;
  AO_list_type *AO;
  MO_center_list_type *center,*unique_center;


  /*******

    the beginning of the input file should have the parameters for all
    of the atoms.

  ********/
  skipcomments(infile,instring);
  upcase(instring);
  while(!strstr(instring,"#BEGIN_PARMS")){
    skipcomments(infile,instring);
    upcase(instring);
  }

  printf("Enter contraction coefficient [%lf]: ",DEFAULT_CONTRACT);

  fgets(foostring,80,stdin);
  if(foostring[0] != '\n' && foostring[0] != 0 ){
    contraction = (float)atof(foostring);
  } else{
    contraction = DEFAULT_CONTRACT;
  }

  surf->contraction = contraction;

  printf("\nEnter number of radial lookup table entries [%d]: ",
         DEFAULT_NUM_ENTRIES);

  fgets(foostring,80,stdin);

  if(foostring[0] != '\n' && foostring[0] != 0 ){
    surf->num_lookup_entries = atoi(foostring);
  } else{
    surf->num_lookup_entries = DEFAULT_NUM_ENTRIES;
  }


  printf("\nEnter minimum lookup table value [%lf]: ",
         DEFAULT_LOOKUP_MIN);

  fgets(foostring,80,stdin);
  if(foostring[0] != '\n' && foostring[0] != 0 ){
    surf->lookup_min = (float)atof(foostring);
  } else{
    surf->lookup_min = DEFAULT_LOOKUP_MIN;
  }

  printf("\nEnter maximum lookup table value [%lf]: ",
         DEFAULT_LOOKUP_MAX);

  fgets(foostring,80,stdin);
  if( foostring[0] != '\n' && foostring[0] != 0 ){
    surf->lookup_max = (float)atof(foostring);
  } else{
    surf->lookup_max = DEFAULT_LOOKUP_MAX;
  }

  /* loop until we hit the end of the parm section */
  AO_accum = 0;
  skipcomments(infile,instring);
  upcase(instring);
  color = 2;

  /* get space for the unique centers array */
  surf->unique_centers = (MO_center_list_type *)
    D_CALLOC(max_unique,sizeof(MO_center_list_type));
  if( !surf->unique_centers ) fatal("Can't get unique_center memory.");
  surf->num_unique = 0;

  /* read out the unique centers */
  while(!strstr(instring,"#END_PARMS")){
    sscanf(instring,"%s %d %lf %lf %d %lf %lf %d %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf",type,&ns,&exps,&Hs,&np,&expp,&Hp,&nd,&expd1,&Hd,&c_d1,&expd2,&c_d2,&nf,&expf1,&Hf,&c_f1,&expf2,&c_f2);

    color++;
    /* contract the orbitals */
    exps *= contraction;
    expp *= contraction;
    expd1 *= contraction;
    expd2 *= contraction;
    expf1 *= contraction;
    expf2 *= contraction;

    unique_center = &(surf->unique_centers[surf->num_unique]);
    unique_center->type = (char *)D_CALLOC(8,sizeof(char));
    if( !unique_center->type )
      fatal("Can't get space for type!");
    strcpy(unique_center->type,type);
    unique_center->atom_type = surf->num_unique;

    /********

      set up this unique atom

      here we just set the radial stuff so that we can use a lookup
      table for the radial values.  The angular operations
      are handled by the individual atoms.

    *********/
    unique_center->AO_list = (AO_list_type *)D_CALLOC(MAX_AOS,sizeof(AO_list_type));
    if(!unique_center->AO_list){
      fatal("can't allocate unique_center->AO_list");
    }

    /*******************

      s orbitals

    ********************/
    unique_center->num_AOs = 0;
    if( ns ){
      AO = &(unique_center->AO_list[unique_center->num_AOs]);
      switch(ns){
      case 1:
        AO->rad_func = eval_1_rad; break;
      case 2:
        AO->rad_func = eval_2_rad; break;
      case 3:
        AO->rad_func = eval_3_rad; break;
      case 4:
        AO->rad_func = eval_4_rad; break;
      case 5:
        AO->rad_func = eval_5_rad; break;
      case 6:
        AO->rad_func = eval_6_rad; break;
      case 7:
        AO->rad_func = eval_7_rad; break;
      }
      AO->zeta1 = exps;

      /* do the lookup_table */
      AO->rad_lookup_tbl = (lookup_table_type *)
        D_CALLOC(1,sizeof(lookup_table_type));
      if( !AO->rad_lookup_tbl ) fatal("Can't allocate lookup table\n");
      AO->rad_lookup_tbl->min_val = surf->lookup_min;
      AO->rad_lookup_tbl->max_val = surf->lookup_max;
      AO->rad_lookup_tbl->num_entries = surf->num_lookup_entries;
      AO->rad_lookup_tbl->step = (surf->lookup_max - surf->lookup_min) /
        (float)surf->num_lookup_entries;
      AO->rad_lookup_tbl->values = (float *)D_CALLOC(surf->num_lookup_entries,
                                                   sizeof(float));
      if( !AO->rad_lookup_tbl->values )
        fatal("Can't get lookup table entries memory");

      unique_center->num_AOs += 1;
    }

    /*******************

      p orbitals

    ********************/
    if( np ){
      for( j=unique_center->num_AOs; j<unique_center->num_AOs+3; j++ ){
        AO = &(unique_center->AO_list[j]);
        switch(np){
        case 1:
          AO->rad_func = eval_1_rad; break;
        case 2:
          AO->rad_func = eval_2_rad; break;
        case 3:
          AO->rad_func = eval_3_rad; break;
        case 4:
          AO->rad_func = eval_4_rad; break;
        case 5:
          AO->rad_func = eval_5_rad; break;
        case 6:
          AO->rad_func = eval_6_rad; break;
        case 7:
          AO->rad_func = eval_7_rad; break;
        }
        AO->zeta1 = expp;
        AO->zeta2 = 0;
        AO->C1 = 0;
        AO->C2 = 0;

        if( j == unique_center->num_AOs ){
          /* do the lookup_table */
          AO->rad_lookup_tbl = (lookup_table_type *)
            D_CALLOC(1,sizeof(lookup_table_type));
          if( !AO->rad_lookup_tbl ) fatal("Can't allocate lookup table\n");
          AO->rad_lookup_tbl->min_val = surf->lookup_min;
          AO->rad_lookup_tbl->max_val = surf->lookup_max;
          AO->rad_lookup_tbl->num_entries = surf->num_lookup_entries;
          AO->rad_lookup_tbl->step = (surf->lookup_max - surf->lookup_min) /
            (float)surf->num_lookup_entries;
          AO->rad_lookup_tbl->values =
            (float *)D_CALLOC(surf->num_lookup_entries,sizeof(float));
          if( !AO->rad_lookup_tbl->values )
            fatal("Can't get lookup table entries memory");
        }
        else{
          AO->rad_lookup_tbl =
            unique_center->AO_list[unique_center->num_AOs].rad_lookup_tbl;
        }


      }
      unique_center->num_AOs += 3;
    }

    /*******************

      d orbitals

      ********************/
    if( nd ){
      /* set the radial functions */
      for( j=unique_center->num_AOs; j<unique_center->num_AOs+5; j++ ){
        AO = &(unique_center->AO_list[j]);
        switch(nd){
        case 1:
          AO->rad_func = eval_1_rad; break;
        case 2:
          AO->rad_func = eval_2_rad; break;
        case 3:
          AO->rad_func = eval_3_rad; break;
        case 4:
          AO->rad_func = eval_4_rad; break;
        case 5:
          AO->rad_func = eval_5_rad; break;
        case 6:
          AO->rad_func = eval_6_rad; break;
        case 7:
          AO->rad_func = eval_7_rad; break;
        }
        AO->zeta1 = expd1;
        AO->zeta2 = expd2;
        AO->C1 = c_d1;
        AO->C2 = c_d2;

        if( j == unique_center->num_AOs ){
          /* do the lookup_table */
          AO->rad_lookup_tbl = (lookup_table_type *)
            D_CALLOC(1,sizeof(lookup_table_type));
          if( !AO->rad_lookup_tbl ) fatal("Can't allocate lookup table\n");
          AO->rad_lookup_tbl->min_val = surf->lookup_min;
          AO->rad_lookup_tbl->max_val = surf->lookup_max;
          AO->rad_lookup_tbl->num_entries = surf->num_lookup_entries;
          AO->rad_lookup_tbl->step = (surf->lookup_max - surf->lookup_min) /
            (float)surf->num_lookup_entries;
          AO->rad_lookup_tbl->values =
            (float *)D_CALLOC(surf->num_lookup_entries,sizeof(float));
          if( !AO->rad_lookup_tbl->values )
            fatal("Can't get lookup table entries memory");
        }
        else{
          AO->rad_lookup_tbl =
            unique_center->AO_list[unique_center->num_AOs].rad_lookup_tbl;
        }
      }
      unique_center->num_AOs += 5;
    }
    /********************

      f-orbitals

      ******************/

    if(nf){
      for(j=unique_center->num_AOs; j<unique_center->num_AOs+7; j++){
        AO = &(unique_center->AO_list[j]);
        switch(nf){
        case 1:
          AO->rad_func = eval_1_rad; break;
        case 2:
          AO->rad_func = eval_2_rad; break;
        case 3:
          AO->rad_func = eval_3_rad; break;
        case 4:
          AO->rad_func = eval_4_rad; break;
        case 5:
          AO->rad_func = eval_5_rad; break;
        case 6:
          AO->rad_func = eval_6_rad; break;
        case 7:
          AO->rad_func = eval_7_rad; break;
        }
        AO->zeta1 = expf1;
        AO->zeta2 = expf2;
        AO->C1 = c_f1;
        AO->C2 = c_f2;

        if(j==unique_center->num_AOs){
          AO->rad_lookup_tbl = (lookup_table_type *)D_CALLOC(1,sizeof(lookup_table_type));

          if(!AO->rad_lookup_tbl) fatal("Can't allocate lookup table\n");
          AO->rad_lookup_tbl->min_val = surf->lookup_min;
          AO->rad_lookup_tbl->max_val = surf->lookup_max;
          AO->rad_lookup_tbl->num_entries = surf->num_lookup_entries;
          AO->rad_lookup_tbl->step = (surf->lookup_max - surf->lookup_min)/(float)surf->num_lookup_entries;
          AO->rad_lookup_tbl->values = (float *)D_CALLOC(surf->num_lookup_entries,sizeof(float));
          if(!AO->rad_lookup_tbl->values) fatal("Can't get lookup table entries memory\n");
        }
        else{
          AO->rad_lookup_tbl = unique_center->AO_list[unique_center->num_AOs].rad_lookup_tbl;
        }
      }
      unique_center->num_AOs += 7;
    }


    surf->num_unique++;
    /* check to see if we need more memory */
    if(surf->num_unique == max_unique){
      max_unique++;
      surf->unique_centers = (MO_center_list_type *)
        D_REALLOC((char *)surf->unique_centers,
                (unsigned)(max_unique*sizeof(MO_center_list_type)));
      if( !surf->unique_centers )
        fatal("Can't D_REALLOC unique_center memory.");
    }

    /********

      now loop over all the centers and set up the
      angular parms of those that
      are this type.

      set the appropriate pointers to the lookup tables as well.

    *********/
    center = surf->MO_centers;
    for(i=0;i<surf->num_centers;i++){
#ifdef USING_THE_MAC
      strcpy(temptype,center[i].type);
      upcase(temptype);
      upcase(type);
      if(!strcmp(temptype,type)){
#else

      if( !strcasecmp(center[i].type,type) ){
#endif
        center[i].num_AOs = 0;
        if(center[i].AO_list) D_FREE(center[i].AO_list);
        if( unique_center->num_AOs ){
          center[i].AO_list = (AO_list_type *)
            D_CALLOC(unique_center->num_AOs,sizeof(AO_list_type));
          if(!center[i].AO_list) fatal("can't allocate center[i].AO_list");
        }
        /********

          this is the right type, add in the parms

        *********/

        center[i].color = color;
        /*******************

          s orbitals

        ********************/
        center[i].num_AOs = 0;
        if( ns ){
          AO = &(center[i].AO_list[center[i].num_AOs]);
          AO->ang_func = eval_s_ang;
          AO->rad_lookup_tbl =
            unique_center->AO_list[center[i].num_AOs].rad_lookup_tbl;
          center[i].num_AOs += 1;
        }

        /*******************

          p orbitals

          ********************/
        if( np ){
          for( j=center[i].num_AOs; j<center[i].num_AOs+3; j++ ){
            AO = &(center[i].AO_list[j]);
            AO->zeta1 = expp;
            AO->zeta2 = 0;
            AO->C1 = 0;
            AO->C2 = 0;
            AO->rad_lookup_tbl =
              unique_center->AO_list[j].rad_lookup_tbl;
          }
          j = center[i].num_AOs;
          center[i].AO_list[j].ang_func = eval_px_ang;
          center[i].AO_list[j+1].ang_func = eval_py_ang;
          center[i].AO_list[j+2].ang_func = eval_pz_ang;
          center[i].num_AOs += 3;
        }

        /*******************

          d orbitals

          ********************/
        if( nd ){
          for( j=center[i].num_AOs; j<center[i].num_AOs+5; j++ ){
            AO = &(center[i].AO_list[j]);
            AO->zeta1 = expd1;
            AO->zeta2 = expd2;
            AO->C1 = c_d1;
            AO->C2 = c_d2;
            AO->rad_lookup_tbl =
              unique_center->AO_list[j].rad_lookup_tbl;
          }
          j = center[i].num_AOs;
          center[i].AO_list[j].ang_func = eval_dx2y2_ang;
          center[i].AO_list[j+1].ang_func = eval_dz2_ang;
          center[i].AO_list[j+2].ang_func = eval_dxy_ang;
          center[i].AO_list[j+3].ang_func = eval_dxz_ang;
          center[i].AO_list[j+4].ang_func = eval_dyz_ang;
          center[i].num_AOs += 5;
        }
        AO_accum += center[i].num_AOs;
      }
    }
    skipcomments(infile,instring);
    upcase(instring);
  }


  /*******

    okay, that takes care of the parameters.. get the number of AO's now

  *******/
  while(!strstr(instring,"#NUM_AOS")){
    skipcomments(infile,instring);
    upcase(instring);
  }
  sscanf(instring,"%s %d",foostring,&num_AOs);

  /* error checking */
  if( num_AOs != AO_accum ){
    sprintf(foostring,
            "Num AO's in parm section: %d does not equal num AO's specified: %d\n",
            AO_accum,num_AOs);
    error(foostring);
  }

  /* read out the number of MO's in this file */
  while(!strstr(instring,"#NUM_MOS")){
    skipcomments(infile,instring);
    upcase(instring);
  }
  sscanf(instring,"%s %d",foostring,&num_MOs);

  /******

    read out the actual numbers and k points of the MO's
    so the user can be prompted intelligently

  *******/
  MO_numbers = (int *)D_CALLOC(num_MOs,sizeof(int));
  MO_kpoints = (point_type *)D_CALLOC(num_MOs,sizeof(point_type));
  if(!MO_numbers || !MO_kpoints )
    fatal("Can't get memory to store MO_numbers or MO_kpoints.");
  for(i=0;i<num_MOs;i++){
    skipcomments(infile,instring);
    sscanf(instring,"%d %lf %lf %lf",&(MO_numbers[i]),
           &(MO_kpoints[i].x),&(MO_kpoints[i].y),&(MO_kpoints[i].z));
  }

  /* set pointers to those arrays */
  surf->kpoints = MO_kpoints;
  surf->MO_numbers = MO_numbers;
  surf->num_MOs = num_MOs;

  /* get space to store the characters */
  surf->characters = (MO_character_type *)
    D_CALLOC(num_MOs,sizeof(MO_character_type));
  if( !surf->characters ) fatal("can't allocate MO characters.");

  /* ask the user which MO they wish to see... */
  which_MO = choose_MO(surf);
  if( which_MO < 0 ) return(-1);
  else{
    surf->active_MO = which_MO;
    surf->num_MOs = num_MOs;

#if 0
    /******

      we have to loop through all the centers and allocate
      space for the arrays of coefficients used to store the
      MO's

    ******/
    for( i=0; i<surf->num_centers; i++){
      AO = center[i].AO_list;
      for(j=0; j<center[i].num_AOs;j++){
        AO[j].coeff = (float *)D_CALLOC(num_MOs,sizeof(float));
        AO[j].coeffI = (float *)D_CALLOC(num_MOs,sizeof(float));
        if( !AO[j].coeffI  ) fatal("Can't get space for AO coeff list");
      }
    }

#endif
#if 0
    /* okay, skip forward until we find the MO the user wants */
    for(i=0;i<which_MO;i++){
      while(!strstr(instring,"#BEGIN_MO")){
        skipcomments(infile,instring);
        upcase(instring);
      }
      instring[0] = 0;
    }
#endif

    /* we're at the beginning, now read out the coefficients */
    center = surf->MO_centers;
    for(i=0;i<num_MOs;i++){
      while(!strstr(instring,"#BEGIN_MO")){
        skipcomments(infile,instring);
        upcase(instring);
      }
      sscanf(instring,"%s %d %d %d",foostring,&surf->characters[i].planes[X_AX],
             &surf->characters[i].planes[Y_AX],&surf->characters[i].planes[Z_AX]);
      for(j=0;j<surf->num_centers;j++){
        AO = center[j].AO_list;
        for(k=0;k<center[j].num_AOs;k++){
          skipcomments(infile,instring);
          upcase(instring);
          /* some error checking */

          if( strstr(instring,"#END_MO") ){
            return(-2);
          } else{
            sscanf(instring,"%lf %lf",&(AO[k].coeff[i]),&(AO[k].coeffI[i]));
          }
        }
      }
    }
  }

  /* if we made it to here, everything is fine... return now */
  return(0);

}



/****************************************************************************
 *
 *                   Procedure new_MO_surface
 *
 * Arguments: filename: pointer to type char
 *
 * Returns: none
 *
 * Action:  Does everything necessary to construct a new MO isosurface.
 *
 ****************************************************************************/
void new_MO_surface(char *filename)
{
  int i;
  char file_name[80],file_name_copy[80];
  char *theinline;
  char failed;
  FILE *infile;
  MO_surface_type *MO_surf;
  MO_center_list_type *center;

  failed = 0;

  /* set up a new object to hold the surface */
  makenewobject();
  whichobj = head->obj;

  /* now build the primitive */
  whichobj->prim = (prim_type *)D_CALLOC(1,sizeof(prim_type));
  if( !whichobj->prim )fatal("Can't get space for surface primitive.");
  whichobj->prim->which = MO_SURF;

  whichobj->prim->MO_surf = (MO_surface_type *)
    D_CALLOC(1,sizeof(MO_surface_type));
  if( !whichobj->prim->MO_surf)
    fatal("Can't get space for MO surface.");
  MO_surf = whichobj->prim->MO_surf;

#ifndef USING_THE_MAC
  if( !filename ){
    display("Look in the xterm...");
#ifndef USE_READLINE
    printf("Enter the input file name used to construct the MOs: ");
    scanf("%s\n",file_name);
#else
    theinline=
      readline("Enter the input file name used to construct the MOs: ");
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

#else
  if(!filename){
    infile = choose_mac_file(file_name,MAC_FOPEN_NOOPEN_CD);
    if(infile) fclose(infile);
  } else{
    strcpy(file_name,filename);
  }
#endif


  /*********

    open the output file to read out the atomic coordinates and set up
    the molecule structure associated with this surface.

  *********/
  strcpy(whichobj->prim->MO_surf->filename,file_name);
  strcpy(file_name_copy,file_name);
  strcat(file_name_copy,".out");
  infile = fopen(file_name_copy,"r");
  if(!infile){
    printf("Problems opening file: %s\n",file_name_copy);
    display("oooooops!");
    failed = 1;
  }

  /* read out the molecular data */
  if( !failed ){
    whichobj->prim->MO_surf->molec = (molec_type *)
      D_CALLOC(1,sizeof(molec_type));
    if( !whichobj->prim->MO_surf->molec )
      fatal("Can't get space to store the associated molecule.");

    read_molecule_data(infile,whichobj->prim->MO_surf->molec);


    /*********

      allocate space for the centers in the MO_surf structure,
      then loop through them

      *********/
    MO_surf->num_centers = MO_surf->molec->num_atoms;
    MO_surf->MO_centers = (MO_center_list_type *)
      D_CALLOC(MO_surf->num_centers,sizeof(MO_center_list_type));
    if( !MO_surf->MO_centers ) fatal("Can't get memory for MO_centers.");

    center = MO_surf->MO_centers;
    for(i=0;i<MO_surf->num_centers;i++){
      /* set the pointer to this center's location and symbol */
      center[i].loc = &(MO_surf->molec->atoms[i].loc);
      center[i].type = MO_surf->molec->atoms[i].type;
    }

    /*********

      open the MO file to read out the MO details

      *********/
    strcpy(file_name_copy,file_name);
    strcat(file_name_copy,".MO");
    infile = fopen(file_name_copy,"r");
    if(!infile){
      printf("Problems opening file: %s\n",file_name_copy);
      display("oooooops!");
      failed = 1;
    }

  }
  /* read out the MO details and make sure things went okay */
  if( !failed ){
    failed = read_MO_surface(infile,MO_surf);
  }

  if( failed ){
    /* free up the memory we asked for */
    if( MO_surf->molec ){
      D_FREE(MO_surf->molec);
    }
    if( MO_surf->MO_centers ){
      D_FREE(MO_surf->MO_centers);
    }

    D_FREE(MO_surf);
    D_FREE(whichobj->prim);
    D_FREE(whichobj);
    whichobj = 0;
    if( head ){
      head->obj = 0;
      head = head->next;
    }
  } else{
    whichobj->scale.x=whichobj->scale.y=whichobj->scale.z=0.5;
    whichobj->trans.x=0;whichobj->trans.y=0;
    whichobj->trans.z=0;
    whichobj->cent.x = g_xmax/2;
    whichobj->cent.y = g_ymax/2;
    MO_surf->adf_plot = 0;
    /* build the lookup table */
    build_radial_lookup_table(MO_surf);

#ifdef INTERACTIVE_USE
    /* create the button window */
    build_MO_surf_button_window(&button_wins,MO_surf);
    whichobj->prim->but_win = button_wins;
    build_molec_button_window(&whichobj->prim->but_win->child,MO_surf->molec);

#endif

    whichobj->prim->which = MO_SURF;
    MO_surf->surface_value = .08;
    MO_surf->surface_tolerance = .005;
    MO_surf->slop = 20.0;
    MO_surf->voxel_size = .2;
    MO_surf->search_radius = .1;
    MO_surf->do_shading = 1;
    MO_surf->display_molec = 1;

    /* copy the MO_list into a raw MO_list array */
    if(MO_surf->molec->num_dim > 0){
      MO_surf->raw_MO_centers =
        (MO_center_list_type *)D_CALLOC(MO_surf->num_centers,
                                      sizeof(MO_center_list_type));
      if( !MO_surf->raw_MO_centers )
        fatal("Can't get space for raw_MO_centers");

      bcopy(MO_surf->MO_centers,MO_surf->raw_MO_centers,
            MO_surf->num_centers*sizeof(MO_center_list_type));

      for(i = 0; i< MO_surf->num_centers; i++ ){
        MO_surf->raw_MO_centers[i].AO_list =
          (AO_list_type *)D_CALLOC(MO_surf->raw_MO_centers[i].num_AOs,
                                   sizeof(AO_list_type));
        if(!MO_surf->raw_MO_centers[i].AO_list)
          fatal("can't allocated MO_surf->raw_MO_centers[i].AO_list.");
        bcopy(MO_surf->MO_centers[i].AO_list,
              MO_surf->raw_MO_centers[i].AO_list,
              MO_surf->MO_centers[i].num_AOs*sizeof(AO_list_type));
      }
      MO_surf->num_centers_in_cell = MO_surf->num_centers;
    } else{
      MO_surf->raw_MO_centers = MO_surf->MO_centers;
    }

    /* set up some default parameters */
    MO_surf->molec->draw_connectors = 1;
    MO_surf->molec->fancy_lines = 1;
    MO_surf->molec->outlines_on = 1;
    MO_surf->molec->shading_on = 1;
    MO_surf->molec->hydrogens_on = 1;
    MO_surf->molec->dummies_on = 1;
    MO_surf->molec->line_width = 1;
    MO_surf->molec->rad_mult = 1.0;
    MO_surf->molec->tubes_on = 1;
    strcpy(MO_surf->plane_dirs,"XYZ");

    /* set up the gridded calculation */
    determine_mol_bounds(MO_surf->molec,&(MO_surf->bmin),
                         &(MO_surf->bmax));

  }
}

#ifdef INCLUDE_ADF_PLOTS
/****************************************************************************
 *
 *                   Function read_ADF_MO_data
 *
 * Arguments: infile: pointer to type FILE
 *           MO_surf: pointer to MO_surface_type
 *
 * Returns: int
 *
 * Action:  Reads out both atomic positions and the MO data for an ADF MO
 *       returns 0 on success, 1 otherwise
 *
 *
 ****************************************************************************/
int read_ADF_MO_data(FILE *infile,MO_surface_type *MO_surf)
{
  int i,j;
  char instring[MAX_STR_LEN],foostring[80];
  char tempstring[MAX_STR_LEN],*string_token;
  char *irrep_names;
  char *curr_irrep_name;
  molec_type *molec;
  int num_types;
  int *num_per_type;
  MO_center_list_type *curr_center,*center;
  AO_list_type *AO;
  int num_irreps;
  int atoms_so_far;
  int choice;
  int num_AOs_involved,*AOs_involved;
  int max_MOs,num_MOs,offset,which_MO;
  float *coeffs, *energies;

  /* quick error checking */
  if( !infile ) FATAL_BUG("Bogus file passed to read_ADF_MO_data");
  if( !MO_surf ) FATAL_BUG("Bogus MO_surf passed to read_ADF_MO_data");

  /* start out by getting space for the molecule */
  molec = (molec_type *)D_CALLOC(1,sizeof(molec_type));
  if(!molec)fatal("Can't get memory for the molecule.");
  MO_surf->molec = molec;

  MO_surf->adf_plot = 1;

  /* Get the radial lookup table parameters */
  printf("\nEnter number of radial lookup table entries [%d]: ",
         DEFAULT_NUM_ENTRIES);

  fgets(foostring,80,stdin);

  if(foostring[0] != '\n' && foostring[0] != 0 ){
    MO_surf->num_lookup_entries = atoi(foostring);
  } else{
    MO_surf->num_lookup_entries = DEFAULT_NUM_ENTRIES;
  }


  printf("\nEnter minimum lookup table value [%lf]: ",
         DEFAULT_LOOKUP_MIN);

  fgets(foostring,80,stdin);
  if(foostring[0] != '\n' && foostring[0] != 0 ){
    MO_surf->lookup_min = (float)atof(foostring);
  } else{
    MO_surf->lookup_min = DEFAULT_LOOKUP_MIN;
  }

  printf("\nEnter maximum lookup table value [%lf]: ",
         DEFAULT_LOOKUP_MAX);
  fgets(foostring,80,stdin);
  if( foostring[0] != '\n' && foostring[0] != 0 ){
    MO_surf->lookup_max = (float)atof(foostring);
  } else{
    MO_surf->lookup_max = DEFAULT_LOOKUP_MAX;
  }

  /* start reading out the file */
  skipcomments(infile,instring);
  upcase(instring);
  while(!strstr(instring,"#NUM_TYPES")){
    skipcomments(infile,instring);
    upcase(instring);
  }

  sscanf(instring,"%*s %d",&MO_surf->num_unique);
  num_types = MO_surf->num_unique;

  /* get space for the unique atoms */
  MO_surf->unique_centers = (MO_center_list_type *)D_CALLOC(MO_surf->num_unique,
                                                          sizeof(MO_center_list_type));
  if(!MO_surf->unique_centers)fatal("Can't allocate unique_centers");
  num_per_type = (int *)D_CALLOC(MO_surf->num_unique,sizeof(int));
  if(!num_per_type)fatal("Can't allocate num_per_type");

  for(i=0;i<MO_surf->num_unique;i++){
    curr_center = &(MO_surf->unique_centers[i]);
    curr_center->type = (char *)D_CALLOC(8,sizeof(char));
    if(!curr_center->type)fatal("can't allocated type");

    skipcomments(infile,instring);
    sscanf(instring,"%s %d",curr_center->type,&curr_center->num_AOs);

    curr_center->AO_list = (AO_list_type *)
      D_CALLOC(curr_center->num_AOs,sizeof(AO_list_type));
    if(!curr_center->AO_list) fatal("can't allocated curr_center->AO_list.");
#if 0
    if( curr_center->num_AOs >= MAX_AOS){
      FATAL_BUG("MAX_AOS too small for this system");
    }
#endif
    /* loop through the AOs and read the parameters for this type of atom */
    for(j=0;j<curr_center->num_AOs;j++){
      AO = &(curr_center->AO_list[j]);
      skipcomments(infile,instring);
      sscanf(instring,"%d %d %d %d %lf %lf",&AO->kx,&AO->ky,&AO->kz,&AO->kr,
             &AO->zeta1,&AO->norm_fact);

      /* get space for and set up the AO lookup tables */
      AO->rad_lookup_tbl = (lookup_table_type *)D_CALLOC(1,sizeof(lookup_table_type));
      if(!AO->rad_lookup_tbl) fatal("Can't allocate a lookup table");
      AO->rad_lookup_tbl->min_val = MO_surf->lookup_min;
      AO->rad_lookup_tbl->max_val = MO_surf->lookup_max;
      AO->rad_lookup_tbl->num_entries = MO_surf->num_lookup_entries;
      AO->rad_lookup_tbl->step = (MO_surf->lookup_max - MO_surf->lookup_min) /
        (float)MO_surf->num_lookup_entries;
      AO->rad_lookup_tbl->values = (float *)D_CALLOC(MO_surf->num_lookup_entries,
                                                   sizeof(float));
      if( !AO->rad_lookup_tbl->values )
        fatal("Can't get lookup table entries memory");
    }
    /* now read out the number of atoms of this type */
    skipcomments(infile,instring);
    sscanf(instring,"%d",&(num_per_type[i]));


  }

  /* read out the total number of atoms */
  while(!strstr(instring,"#NUM_ATOMS")){
    skipcomments(infile,instring);
    upcase(instring);
  }
  sscanf(instring,"%*s %d",&molec->num_atoms);

  molec->atoms = (atom_type *)D_CALLOC(molec->num_atoms,sizeof(atom_type));
  if(!molec->atoms)fatal("Can't allocate memory for the molec atoms");

  MO_surf->MO_centers = (MO_center_list_type *)D_CALLOC(molec->num_atoms,
                                                       sizeof(MO_center_list_type));
  if(!MO_surf->MO_centers)fatal("can't allocate centers");
  MO_surf->num_centers = molec->num_atoms;
  MO_surf->orig_num_centers = molec->num_atoms;
  MO_surf->raw_MO_centers = MO_surf->MO_centers;

  /* read in the atomic positions */
  for(i=0;i<molec->num_atoms;i++){
    skipcomments(infile,instring);
    sscanf(instring,"%lf %lf %lf",&molec->atoms[i].loc.x,
           &molec->atoms[i].loc.y,&molec->atoms[i].loc.z);
  }
  molec->num_frames = 1;
  molec->num_lines = (int *)D_CALLOC(molec->num_frames,sizeof(int));


  /*******

    connect the atoms with their lookup tables and fill in the other
    parameters which they will need

  ********/
  atoms_so_far = 0;
  for(i=0;i<num_types;i++){
    curr_center = &MO_surf->unique_centers[i];

    for(j=0;j<num_per_type[i];j++){
      center = &MO_surf->MO_centers[atoms_so_far];
      center->loc = &molec->atoms[atoms_so_far].loc;
      center->num_AOs = curr_center->num_AOs;

      center->AO_list = (AO_list_type *)
        D_CALLOC(center->num_AOs,sizeof(AO_list_type));
      if(!center->AO_list) fatal("can't allocated center->AO_list.");

      bcopy(curr_center->AO_list,center->AO_list,center->num_AOs*sizeof(AO_list_type));

      strcpy(molec->atoms[atoms_so_far].type,curr_center->type);
      molec->atoms[atoms_so_far].num = atoms_so_far;

      atoms_so_far++;
    }
  }

  /* fill in the atomic radii and colors */
  fill_atomic_info(molec->num_atoms,molec->atoms);
  molec->bond_tol = DEF_BOND_TOL;
  molec->old_bond_tol = DEF_BOND_TOL;
  determine_connections(molec);



  /* read the number of irreps */
  while(!strstr(instring,"#NUM_IRREPS")){
    skipcomments(infile,instring);
    upcase(instring);
  }
  sscanf(instring,"%*s %d",&num_irreps);

  irrep_names = (char *)D_CALLOC(num_irreps*80,sizeof(char));
  if(!irrep_names)fatal("Can't allocate irrep names");

  /*******

    scan through the file until we've found all the irreps,
    store the name of each so that the user can be intelligently
    prompted for which one they want

  ********/
  curr_irrep_name = irrep_names;
  for(i=0;i<num_irreps;i++){
    while(!strstr(instring,"#BEGIN_IRREP")){
      skipcomments(infile,instring);
      upcase(instring);
    }
    skipcomments(infile,instring);
    sscanf(instring,"%s",curr_irrep_name);
    curr_irrep_name+=80;
  }

  choice = 0;
  while(!choice){
    printf("there are %d irreps in the file\n",num_irreps);
    curr_irrep_name = irrep_names;
    for(i=0;i<num_irreps;i++){
      printf("%d: %s\n",i+1,curr_irrep_name);
      curr_irrep_name+=80;
    }
#if 0
    printf("which number irrep would you like to see (a negative number to cancel)? ");
    scanf("%d",&choice);
#endif
    strcpy(tempstring,"which irrep number you would like to see (a negative number to cancel)");
    readintparm(tempstring,&choice);
    if( choice > num_irreps || !choice ) printf("%d is an invalid choice\n",choice);
    if( choice < 0 ) return(1);
  }

  /* now rewind the file and find that irrep */
  rewind(infile);
  skipcomments(infile,instring);
  upcase(instring);
  for(i=0;i<choice;i++){
    while(!strstr(instring,"#BEGIN_IRREP")){
      skipcomments(infile,instring);
      upcase(instring);
    }
    if( i < choice-1 ){
      skipcomments(infile,instring);
      upcase(instring);
    }
  }
  /*******

    we're there, now skip past the name of the irrep then start reading
    in the info

  ********/
  skipcomments(infile,instring);
  skipcomments(infile,instring);
  sscanf(instring,"%d",&num_AOs_involved);
  AOs_involved = (int *)D_CALLOC(num_AOs_involved,sizeof(int));
  if(!AOs_involved)fatal("can't allocate AOs_involved");

  /* now read in the AOs involved */
  i=0;
  skipcomments(infile,instring);
  strcpy(tempstring,instring);
  string_token = strtok(tempstring," ");
  while(i<num_AOs_involved){
    if( string_token ){
      sscanf(string_token,"%d",&(AOs_involved[i]));
      AOs_involved[i]--;
      string_token = strtok(NULL," ");
      i++;
    }else{
      skipcomments(infile,instring);
      strcpy(tempstring,instring);
      string_token = strtok(tempstring," ");
    }
  }

  /* read in the number of MOs, then get each MO */
  skipcomments(infile,instring);
  sscanf(instring,"%d",&num_MOs);
  max_MOs = num_MOs;
  printf("There are probably more MOs (there are: %d) than you want to see in the file.\n",num_MOs);
  readintparm("num_MOs",&num_MOs);

  if( num_MOs < 1 || num_MOs > max_MOs ){
    fprintf(stderr,"Oh, you're a real joker.  num_MOs set to: %d\n",max_MOs);
    num_MOs = max_MOs;
    display("Nice try...");
  }

  MO_surf->num_MOs = num_MOs;
  if( num_MOs >= MAX_MOS ) FATAL_BUG("MAX_MOS is too small for this system");

  energies = (float *)D_CALLOC(num_MOs,sizeof(float));
  coeffs = (float *)D_CALLOC(num_MOs*num_AOs_involved,sizeof(float));
  if(!energies || !coeffs) fatal("Can't allocate either energies or coeffs");

  MO_surf->characters = (MO_character_type *)
    D_CALLOC(num_MOs,sizeof(MO_character_type));
  if( !MO_surf->characters ) fatal("can't allocate MO characters.");

  for(i=0;i<num_MOs;i++){
    skipcomments(infile,instring);
    sscanf(instring,"%lf",&energies[i]);

    skipcomments(infile,instring);
    strcpy(tempstring,instring);
    string_token = strtok(tempstring," ");
    j=0;
    while(j<num_AOs_involved){
      if( string_token ){
        sscanf(string_token,"%lf",&(coeffs[i*num_AOs_involved+j]));
        string_token = strtok(NULL," ");
        j++;
      }else{
        skipcomments(infile,instring);
        strcpy(tempstring,instring);
        string_token = strtok(tempstring," ");
      }
    }
  }

  /* a bit of error checking */
  skipcomments(infile,instring);
  upcase(instring);
  if(num_MOs == max_MOs && !strstr(instring,"#END_IRREP")){
    fprintf(stderr,"Found something unexpected at end of irrep section:\n");
    fprintf(stderr,"%s\n",instring);
    fprintf(stderr,"Continuing anyway... who knows what will happen\n");
    display("Danger Will Robinson!");
  }

  /********

    okay... now we have the coefficients, but they are in a screwy order
    (the same order as in the AOs_involved array).
    We need to copy them into the proper centers, do that now.

  ********/
  for(i=0;i<num_AOs_involved;i++){
    map_AO_number_to_center(AOs_involved[i],MO_surf,&center,&offset);
    for(j=0;j<num_MOs;j++){
      center->AO_list[offset].coeff[j] = coeffs[j*num_AOs_involved+i];
    }
  }

  which_MO = choose_MO(MO_surf);
  MO_surf->active_MO = which_MO;

  D_FREE(energies);
  D_FREE(coeffs);
  D_FREE(AOs_involved);
  D_FREE(num_per_type);
  D_FREE(irrep_names);
  return(0);

}

/****************************************************************************
 *
 *                   Procedure new_ADF_MO_surface
 *
 * Arguments: filename: pointer to type char
 *
 * Returns: none
 *
 * Action:  Does everything necessary to construct a new MO isosurface. for an ADF plot
 *
 ****************************************************************************/
void new_ADF_MO_surface(char *filename)
{
  char file_name[80];
  char *theinline;
  char instring[MAX_STR_LEN];
  char failed;
  FILE *infile;
  MO_surface_type *MO_surf;



  failed = 0;

  /* set up a new object to hold the surface */
  makenewobject();
  whichobj = head->obj;

  /* now build the primitive */
  whichobj->prim = (prim_type *)D_CALLOC(1,sizeof(prim_type));
  if( !whichobj->prim )fatal("Can't get space for surface primitive.");
  whichobj->prim->which = MO_SURF;

  whichobj->prim->MO_surf = (MO_surface_type *)
    D_CALLOC(1,sizeof(MO_surface_type));
  if( !whichobj->prim->MO_surf)
    fatal("Can't get space for MO surface.");
  MO_surf = whichobj->prim->MO_surf;

#ifndef USING_THE_MAC
  if( !filename ){
    display("Look in the xterm...");
#ifndef USE_READLINE
    printf("Enter the file name constructed from the Tape21: ");
    scanf("%s\n",file_name);
#else
    theinline=
      readline("Enter the file name constructed from the Tape21: ");
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

 #else
  if(!filename){
    infile = choose_mac_file(file_name,MAC_FOPEN_NOOPEN_CD);
    if(infile) fclose(infile);
  } else{
    strcpy(file_name,filename);
  }
#endif


  /*********

    open the file and make sure it's valid

  *********/
  strcpy(whichobj->prim->MO_surf->filename,file_name);
  infile = fopen(file_name,"r");
  if(!infile){
    printf("Problems opening file: %s\n",file_name);
    display("oooooops!");
    failed = 1;
  }
  if( !failed ){
    /* check for validity */
    skipcomments(infile,instring);
    upcase(instring);
    if( !strstr(instring,"#ADF_MO_DATA") ){
      printf("That doesn't look like an ADF MO data file... Cancelling Read\n");
      display("whoopsy!");
      failed = 1;
    }

    /* read out the data */
    rewind(infile);

    failed = read_ADF_MO_data(infile,MO_surf);
  }

  if( failed ){
    /* free up the memory we asked for */
    if( MO_surf->molec ){
      D_FREE(MO_surf->molec);
    }
    if( MO_surf->MO_centers ){
      D_FREE(MO_surf->MO_centers);
    }

    D_FREE(MO_surf);
    D_FREE(whichobj->prim);
    D_FREE(whichobj);
    whichobj = 0;
    if( head ){
      head->obj = 0;
      head = head->next;
    }
  } else{
    whichobj->scale.x=whichobj->scale.y=whichobj->scale.z=0.5;
    whichobj->trans.x=0;whichobj->trans.y=0;
    whichobj->trans.z=0;
    whichobj->cent.x = g_xmax/2;
    whichobj->cent.y = g_ymax/2;
    /* build the lookup table */
    build_radial_lookup_table(MO_surf);

#ifdef INTERACTIVE_USE
    /* create the button window */
    build_MO_surf_button_window(&button_wins,MO_surf);
    whichobj->prim->but_win = button_wins;
    build_molec_button_window(&whichobj->prim->but_win->child,MO_surf->molec);

#endif

    whichobj->prim->which = MO_SURF;
    MO_surf->surface_value = .08;
    MO_surf->surface_tolerance = .005;
    MO_surf->slop = 20.0;
    MO_surf->voxel_size = .2;
    MO_surf->search_radius = .1;
    MO_surf->do_shading = 1;
    MO_surf->display_molec = 1;

    MO_surf->raw_MO_centers = MO_surf->MO_centers;

    /* set up some default parameters */
    MO_surf->molec->draw_connectors = 1;
    MO_surf->molec->fancy_lines = 1;
    MO_surf->molec->outlines_on = 1;
    MO_surf->molec->shading_on = 1;
    MO_surf->molec->hydrogens_on = 1;
    MO_surf->molec->dummies_on = 1;
    MO_surf->molec->line_width = 1;
    MO_surf->molec->rad_mult = 1.0;
    MO_surf->molec->tubes_on = 1;
    MO_surf->adf_plot = 1;
    strcpy(MO_surf->plane_dirs,"XYZ");
    /* set up the gridded calculation */
    determine_mol_bounds(MO_surf->molec,&(MO_surf->bmin),
                         &(MO_surf->bmax));

  }
}
#endif
