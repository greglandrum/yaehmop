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
   31.03.98 gL:
     Added code to support sorting of energy levels so
     that they don't have to be ordered in energy in the
     .FMO file.  This is essential for dealing with ADF
     input files.
   26.04.98 gL:
     Added code to support FMO and energy level files which
     already contain the occupations of the MOs.  This is important
     for dealing with ADF files where things may be non-Aufbau.
   03.05.98 gL:
     bounding boxes adjusted
   26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)
***/


/********

  this has got the stuff for dealing with FMO (interaction) diagrams

*********/
#include "viewkel.h"



/****************************************************************************
 *
 *                   Procedure find_FMO_tic_sep
 *
 * Arguments: diagram: pointer to FMO_diagram_type
 *
 * Returns: none
 *
 * Action:  Determines a "reasonable" tic mark separation for the interaction
 *    diagram stored in 'diagram.
 *
 ****************************************************************************/
void find_FMO_tic_sep( FMO_diagram_type *diagram )
{
  float range,norm,tics,tic_sep,log_range;
  float ymax,ymin;

  range = fabs(diagram->min_y-diagram->max_y);

  log_range = log10(range);

  norm = pow(10.0,log_range-(float)((log_range >= 0.0 ) ? (int)log_range :
                                    ((int)log_range-1)));
  if (norm <= 2)
    tics = 0.2;
  else if (norm <= 5)
    tics = 0.5;
  else tics = 1.0;
  tic_sep = tics * pow(10.0,(float)((log_range >= 0.0 ) ? (int)log_range : ((int)log_range-1)));

  diagram->tic_sep_y = tic_sep;

  /* figure out the max and min values spanned by the tic marks */
  ymin = tic_sep * floor(diagram->min_y/tic_sep);
  ymax = tic_sep * ceil(diagram->max_y/tic_sep);
  diagram->num_tics_y = 1 + (ymax - ymin) / tic_sep;

  /* if the first tic mark is outside the range the user specified, skip it! */
  if( ymin < diagram->min_y ){
    diagram->tic_start_y = ymin + tic_sep;
    diagram->num_tics_y--;
  }
  else{
    diagram->tic_start_y = ymin;
  }
  /* same deal with the last tic mark */
  if( ymax > diagram->max_y ){
    diagram->num_tics_y--;
  }

/*
  diagram->max_y = ymax;
  diagram->min_y = ymin;
*/
}



/****************************************************************************
 *
 *                   Procedure fill_levels_with_electrons
 *
 * Arguments: levels: pointer to level_type
 *     num_electrons: float
 *            occups: point to float
 *
 * Returns: none
 *
 * Action: fills the energy levels in 'levels with electrons.  This assumes
 *  that the linked list 'levels is stored with the lowest energy
 *  levels first.
 *  if 'occups is nonzero, those will be used to fill the levels
 *
 ****************************************************************************/
void fill_levels_with_electrons(FMO_level_type *levels,float num_electrons,
                                float *occups)
{
  FMO_level_type *level;
  float electrons_left,electrons_required;
  int i,levels_past;

  /* step through the levels, filling as we go. */
  level = levels;
  levels_past = 0;
  electrons_left = num_electrons;
  while(level && electrons_left > 0.001){
    /* are there enough electrons left to completely fill this level? */
    if( !occups ){
      electrons_required = level->num_degen*2.0;
      if(electrons_left > electrons_required){
        level->num_electrons = electrons_required;
        electrons_left -= electrons_required;
        if( electrons_left < .001 ){
          electrons_left = 0.0;
          level->highest = 1;
        }
      }else{
        level->num_electrons = electrons_left;
        level->highest = 1;
        electrons_left = 0;
      }
    } else{
      level->num_electrons = 0;
      for(i=0; i<level->num_degen; i++){
        level->num_electrons += occups[levels_past];
        levels_past++;
      }
      if( occups[levels_past] == 0 ){
        level->highest = 1;
      }
    }
    level = level->next;
  }
}


/****************************************************************************
 *
 *                   Function map_orbital_to_level
 *
 * Arguments: levels: pointer to level_type
 *         which_orb: int
 *
 * Returns: pointer to FMO_level_type
 *
 * Action: finds which level in 'levels contains orbital 'which_orb. returns
 *   this level
 *
 ****************************************************************************/
FMO_level_type *map_orbital_to_level(FMO_level_type *levels,int which_orb)
{
  FMO_level_type *level;
  FMO_level_type *the_level;
  int num_orbs_traversed;
  int which_level;

  level = levels;
  which_level = -1;
  num_orbs_traversed = 0;

  /* loop until we either run out of levels or find the one we're looking for */
  while(which_level < 0 && level){
    num_orbs_traversed += level->num_degen;
    if( which_orb < num_orbs_traversed ){
      which_level = level->number;
      the_level = level;
    }
    else{
      level = level->next;
    }
  }

  /* check to see that we actually found a level */
  if( which_level < 0 ){
    fatal("map_orbital_to_level couldn't match an orbital.  This is probably a bug.");
  }

  return(the_level);
}

/* this function is used by qsort to sort energy levels */
int Ecompare(const void *E1p,const void *E2p)
{
  float v1,v2;
  v1 = *(float *)E1p;
  v2 = *(float *)E2p;
  return (int)(1000.0*(v1-v2));
}

/****************************************************************************
 *
 *                   Procedure preprocess_FMO_data
 *
 * Arguments: FMO_diagram: pointer to FMO_diagram_type
 *
 * Returns: none
 *
 * Action: Does the initial preprocessing of the FMO data. This includes
 *   sorting the energies into increasing order,
 *   converting the energies into levels (dealing with degeneracies), and
 *   filling in which levels should be connected.
 *
 ****************************************************************************/
void preprocess_FMO_data(FMO_diagram_type *FMO_diagram)
{
  int i,j,k,l;
  int itab;
  int num_levels;
  int which_main_level,which_frag_level;
  int which_frag,orbs_into_frag;
  FMO_fragment_type *FMO_frag;
  FMO_level_type *level,*next_level;
  FMO_level_type *main_level,*frag_level;
  FMO_connect_type *connector,*next_connector;
  float *chg_mat;
  float tot_contrib;
  int num_main_degen, num_frag_degen;

  /* !!! NOTE: this does not sort the charge matrix elements yet! */
  /* sort the energies */
  for(i=0;i<FMO_diagram->num_frags;i++){
    FMO_frag = &(FMO_diagram->frags[i]);
    qsort(FMO_frag->raw_energies,FMO_frag->num_orbs,
          sizeof(float),Ecompare);
  }
  qsort(FMO_diagram->raw_energies,FMO_diagram->tot_num_orbs,
        sizeof(float),Ecompare);

  /* now loop over the fragments and figure out which levels are degenerate */
  for(i=0;i<FMO_diagram->num_frags;i++){
    FMO_frag = &(FMO_diagram->frags[i]);
    num_levels = 0;

    /* if we've already allocated a list of levels, free it now */
    if( FMO_frag->levels ){
      level = FMO_frag->levels;
      while(level){
        next_level = level->next;
        D_FREE(level);
        level = next_level;
      }
    }

    level = 0;
    /*****
      go throught the levels backwards so that the ordering in the linked
      list is from lowest E to highest E. This makes filling the levels with
      electrons easier.
    *****/
    for(j=0;j<FMO_frag->num_orbs;j++){
      /* get space for the new level */
      if( level ){
        level->next = (FMO_level_type *)D_CALLOC(1,sizeof(FMO_level_type));
        if( !(level->next) )
          fatal("Can't get space for an FMO_level structure.");
        level = level->next;
      }
      else{
        /*****

          if this is the first level for this fragment, set a pointer to
          it in the fragment structure.

        *****/
        level = (FMO_level_type *)D_CALLOC(1,sizeof(FMO_level_type));
        if( !(level) ) fatal("Can't get space for an FMO_level structure.");
        FMO_frag->levels = level;
      }

      level->number = num_levels++;

      level->num_degen = 1;
      level->energy = FMO_frag->raw_energies[j];

      /* check for degeneracies */
      while( j+1 < FMO_frag->num_orbs &&
            fabs(FMO_frag->raw_energies[j+1] - level->energy) <
            FMO_diagram->degen_tol){
        level->num_degen++;
        j++;
      }
    }
  }

  /* repeat the same process for the main energy levels */
  if( FMO_diagram->levels ){
    level = FMO_diagram->levels;
    while(level){
      next_level = level->next;
      D_FREE(level);
      level = next_level;
    }
  }

  level = 0;
  num_levels = 0;
  for(j=0;j<FMO_diagram->tot_num_orbs;j++){
    /* get space for the new level */
   if( level ){
      level->next = (FMO_level_type *)D_CALLOC(1,sizeof(FMO_level_type));
      if( !(level->next) )
        fatal("Can't get space for an FMO_level structure.");
      level = level->next;
    }
    else{
      /* if this is the first level set a pointer to it */
      level = (FMO_level_type *)D_CALLOC(1,sizeof(FMO_level_type));
      if( !(level) ) fatal("Can't get space for an FMO_level structure.");
      FMO_diagram->levels = level;
    }

    level->number = num_levels++;
    level->num_degen = 1;
    level->energy = FMO_diagram->raw_energies[j];

    /* check for degeneracies */
    while( j < FMO_diagram->tot_num_orbs  &&
          fabs(FMO_diagram->raw_energies[j+1] - level->energy)
          < FMO_diagram->degen_tol){
      level->num_degen++;
      j++;
    }
  }

  /*****

    okay, now go through and set up the connections

  ******/
  if( FMO_diagram->num_frags ){
    /* if there currently is a connection list, free it now */
    if( FMO_diagram->connects ){
      connector = FMO_diagram->connects;
      while(connector){
        next_connector = connector->next;
        D_FREE(connector);
        connector = next_connector;
      }
    }

    chg_mat = FMO_diagram->chg_mat;
    connector = 0;
    FMO_diagram->connects = 0;

    for(i=0;i<FMO_diagram->tot_num_orbs;i++){

      /* map the orbital to a level */
      main_level = map_orbital_to_level(FMO_diagram->levels,i);
      which_main_level = main_level->number;
      num_main_degen = main_level->num_degen;

      which_frag = 0;
      FMO_frag = &(FMO_diagram->frags[which_frag]);
      orbs_into_frag = 0;

      for(j=0;j<FMO_diagram->tot_num_orbs;j++){
        /* check to see if we need to change to another fragment */
        if( orbs_into_frag == FMO_frag->num_orbs ){
          which_frag++;
          FMO_frag = &(FMO_diagram->frags[which_frag]);
          orbs_into_frag = 0;
        }


        /* map the orbital to a level */
        frag_level = map_orbital_to_level(FMO_frag->levels,orbs_into_frag);
        which_frag_level = frag_level->number;
        num_frag_degen = frag_level->num_degen;

        /******

          build up the total contribution between sets of degenerate
          levels

          *******/
        tot_contrib = 0.0;
        for(k=0;k<num_main_degen;k++){
          itab = (i+k)*FMO_diagram->tot_num_orbs;
          for(l=0;l<num_frag_degen;l++){
            tot_contrib += fabs(chg_mat[itab + j + l]);
          }
        }
        j += num_frag_degen-1;
        orbs_into_frag += num_frag_degen;

        /* don't make connectors that aren't going to be shown */
        if(tot_contrib >= FMO_diagram->min_contrib ){
          /****

            get memory for the new connector

            ****/
          next_connector =
            (FMO_connect_type *)D_CALLOC(1,sizeof(FMO_connect_type));
          if(!next_connector) fatal("Can't get memory for an FMO connector.");

          /* is this the first connector in the list? */
          if( !FMO_diagram->connects ){
            FMO_diagram->connects = next_connector;
            connector = next_connector;
          } else {
            connector->next = next_connector;
            connector = next_connector;
          }


          /* set up the connector */
          connector->which_frag = which_frag;
          connector->fragment_level = which_frag_level;
          connector->main_level = which_main_level;
          connector->contrib = tot_contrib;

          /* set the line-style */
          if( fabs(connector->contrib) < .25 ) connector->linestyle = 3;
          else if(fabs(connector->contrib) < .5) connector->linestyle = 2;
          else if(fabs(connector->contrib) < .75) connector->linestyle = 1;
          else connector->linestyle = 0;
        }
      }
      i += num_main_degen-1;
    }
  }
  /* now fill each fragment and the main set of levels with electrons */
  fill_levels_with_electrons(FMO_diagram->levels,
                             FMO_diagram->tot_num_electrons,
                             FMO_diagram->occups);
  for(i=0;i<FMO_diagram->num_frags;i++){
    fill_levels_with_electrons(FMO_diagram->frags[i].levels,
                               FMO_diagram->frags[i].num_electrons,
                               FMO_diagram->frags[i].occups);
  }
}



/****************************************************************************
 *
 *                   Procedure change_FMO_degen_tol
 *
 * Arguments:    num_args: int
 *                FMO_ptr: pointer to pointer to char
 *
 * Returns: none
 *
 * Action: changes the degeneracy tolerance in an FMO diagram
 *    this is intended to be called from a button window
 *
 ****************************************************************************/
void change_FMO_degen_tol(int num_args,char **FMO_ptr)
{
  FMO_diagram_type *diagram;

  diagram = (FMO_diagram_type *)FMO_ptr[0];

  /* get the new value */
  readfloatparm("FMO degeneracy tolerance",&diagram->degen_tol);

  /* update everything */
  preprocess_FMO_data(diagram);
}





/****************************************************************************
 *
 *                   Procedure read_FMO_data
 *
 * Arguments: infile: pointer to FILE
 *       FMO_diagram: pointer to FMO_diagram_type
 *
 * Returns: none
 *
 * Action: Actually reads the FMO information from the specified file.
 *
 ****************************************************************************/
void read_FMO_data(FILE *infile,FMO_diagram_type *FMO_diagram)
{
  FMO_fragment_type *FMO_frag;
  char instring[MAX_STR_LEN],tempstring[80];
  char read_occups;
  char *string_ptr;
  float *chg_mat;
  int i,j,itab;
  float min_E,max_E;

  min_E = 1e10;
  max_E = -1e10;

  if( !FMO_diagram ) FATAL_BUG("No FMO_diagram passed to read_FMO_data");

  /* make sure that this is really FMO data and that the file isn't empty */
  if( skipcomments(infile,instring) < 0 ){
    error("That file is empty!");
    display("Too Bad....");
    return;
  } else{
    upcase(instring);
    if(!strstr(instring,"FMO RESULTS")){
      error("That file doesn't contain FMO data!  Please give me the name of the .FMO file.");
      display("Too Bad....");
      return;
    }
  }

  /* check if we are gonna be reading occupations as well */
  if( strstr(instring,"OCCUP") ){
    read_occups = 1;
    printf("Reading occupation numbers from FMO file.\n");
  } else{
    read_occups = 0;
  }


  /********

    read out the information in the header

  *********/
  skipcomments(infile,instring);
  sscanf(instring,"%s %d",tempstring,&(FMO_diagram->tot_num_orbs));
  skipcomments(infile,instring);
  sscanf(instring,"%s %lf",tempstring,&(FMO_diagram->tot_num_electrons));
  skipcomments(infile,instring);
  sscanf(instring,"%s %d",tempstring,&(FMO_diagram->num_frags));

  if( FMO_diagram->num_frags < 1 ){
    printf("There are no fragments in the input file...\n");
    printf("\tI'll treat it like an energy level diagram then.\n");
    printf("\tI hope this is right.\n");
  }

  /* get space to store the fragments */
  if( FMO_diagram->num_frags){
    FMO_diagram->frags = (FMO_fragment_type *)D_CALLOC(FMO_diagram->num_frags,
                                                     sizeof(FMO_fragment_type));
    if( !FMO_diagram->frags ) fatal("Can't get memory for FMO_frags.");
  }

  /********

    loop and read out the number of orbitals and electrons in each fragment

  ********/
  for(i=0;i<FMO_diagram->num_frags;i++){
    skipcomments(infile,instring);
    sscanf(instring,"%d %lf",&(FMO_diagram->frags[i].num_orbs),
           &(FMO_diagram->frags[i].num_electrons));
  }

  /******

    loop over fragments, get space for, and read out orbital energies

  *******/
  for(i=0;i<FMO_diagram->num_frags;i++){
    FMO_frag = &(FMO_diagram->frags[i]);
    FMO_frag->raw_energies = (float *)D_CALLOC(FMO_frag->num_orbs,sizeof(float));
    if( !FMO_frag->raw_energies )
      fatal("Can't get memory for FMO fragment raw_energies.");

    if( read_occups ){
      FMO_frag->occups = (float *)D_CALLOC(FMO_frag->num_orbs,sizeof(float));
      if( !FMO_frag->occups )
        fatal("Can't get memory for FMO fragment occupations.");
    }

    /* read out the raw energies */
    for(j=0;j<FMO_frag->num_orbs;j++){
      skipcomments(infile,instring);
      if( FMO_frag->occups ){
        sscanf(instring,"%lf %lf",&(FMO_frag->raw_energies[j]),
               &(FMO_frag->occups[j]));
      } else{
        sscanf(instring,"%lf",&(FMO_frag->raw_energies[j]));
      }
      if( FMO_frag->raw_energies[j] > max_E )
        max_E = FMO_frag->raw_energies[j];
      if( FMO_frag->raw_energies[j] < min_E )
        min_E = FMO_frag->raw_energies[j];
    }
  }

  /*****

    now get space for and read out the molecular energies

  ******/
  FMO_diagram->raw_energies = (float *)D_CALLOC(FMO_diagram->tot_num_orbs,
                                              sizeof(float));
  if(!FMO_diagram->raw_energies)
    fatal("Can't get space for raw molecular energies.");

  if( read_occups ){
    FMO_diagram->occups = (float *)D_CALLOC(FMO_diagram->tot_num_orbs,
                                          sizeof(float));
    if( !FMO_diagram->occups )
      fatal("Can't allocate space for occupations.");
  }
  for(i=0;i<FMO_diagram->tot_num_orbs;i++){
    skipcomments(infile,instring);
    if( read_occups ){
      sscanf(instring,"%lf %lf",&(FMO_diagram->raw_energies[i]),
             &(FMO_diagram->occups[i]));
    }else{
      sscanf(instring,"%lf",&(FMO_diagram->raw_energies[i]));
    }

    if( FMO_diagram->raw_energies[i] > max_E )
      max_E = FMO_diagram->raw_energies[i];
    if( FMO_diagram->raw_energies[i] < min_E )
      min_E = FMO_diagram->raw_energies[i];

  }

  /*****

    now get space for and read out the elements of the charge matrix

  ******/
  if( FMO_diagram->num_frags ){
    FMO_diagram->chg_mat = (float *)D_CALLOC(FMO_diagram->tot_num_orbs*
                                           FMO_diagram->tot_num_orbs,
                                           sizeof(float));
    if( !FMO_diagram->chg_mat )
      fatal("Can't get memory for FMO charge matrix.");

    chg_mat = FMO_diagram->chg_mat;
    for(i=0;i<FMO_diagram->tot_num_orbs;i++){
      itab = i*FMO_diagram->tot_num_orbs;

      skipcomments(infile,instring);

      /* use strtok to read out the space delimited numbers */
      string_ptr = strtok(instring," ");
      if( string_ptr ){
        sscanf((const char *)string_ptr,"%lf",&(chg_mat[itab]));
      } else{
        fatal("EOF hit way too early in FMO file.");
      }

      for(j=1;j<FMO_diagram->tot_num_orbs;j++){
        string_ptr = strtok(0," ");
        if( string_ptr != 0 ){
          sscanf((const char *)string_ptr,"%lf",&(chg_mat[itab+j]));
        } else{
          fatal("ran out of entries when reading charge matrix elements.");
        }
      }
    }
  }
  /* set the min and max energies */
  FMO_diagram->min_y = min_E*1.1;
  FMO_diagram->max_y = max_E<0?max_E*.9:max_E*1.1;

  /* that's all we have to do */
}


/****************************************************************************
 *
 *                   Procedure new_FMO_diagram
 *
 * Arguments: none
 *
 * Returns: none
 *
 * Action: does everything to get space for and read in a new FMO interaction
 *   diagram
 *
 ****************************************************************************/
void new_FMO_diagram(char *filename)
{
  char file_name[80];
  char *theinline;
  char failed;
  FILE *infile;
  FMO_diagram_type *FMO_diagram;

  failed = 0;
  /* set up a new object to hold the graph */
  makenewobject();
  whichobj = head->obj;

  /* now build the graph primitive */
  whichobj->prim = (prim_type *)D_CALLOC(1,sizeof(prim_type));
  if( !whichobj->prim )fatal("Can't get space for graph primitive.");
  whichobj->prim->which = FMO_DIAGRAM;

  whichobj->prim->FMO_diagram = (FMO_diagram_type *)D_CALLOC(1,sizeof(FMO_diagram_type));
  if( !whichobj->prim->FMO_diagram )
    fatal("Can't get space for FMO interaction diagram.");

  FMO_diagram = whichobj->prim->FMO_diagram;
#ifndef USING_THE_MAC
  if( !filename ){
    display("Look in the xterm...");
#ifndef USE_READLINE
    printf("Enter the file name containing the FMO data: ");
    scanf("%s\n",file_name);
#else
    theinline= readline("Enter the file name containing the FMO data: ");
    add_history(theinline);
    if( theinline ){
      sscanf(theinline,"%s",file_name);
      free(theinline);
    } else {
      error("Bad file name");
      file_name[0] = 0;
    }
#endif

  } else{
    strcpy(file_name,filename);
  }

  /* open the file */
  infile = fopen(file_name,"r");
  if(!infile){
    printf("Problems opening file: %s\n",file_name);
    display("oooooops!");
    failed = 1;
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
    failed = 1;
  }
#endif

  strcpy(FMO_diagram->filename,file_name);

  if( !failed ){
    read_FMO_data(infile,FMO_diagram);
  }

  /* check to see if any data was actually read in.... */
  if(!FMO_diagram->raw_energies || failed ){
    /* no... free the memory that we asked for */
    D_FREE(FMO_diagram);
    D_FREE(whichobj->prim);
    D_FREE(whichobj);
    whichobj=0;
    head->obj = 0;
    head = head->next;
  }
  else{
    /******

      set some defaults and do the necessary preprocessing

    ******/
    FMO_diagram->min_contrib = FMO_MIN_CONTRIB;
    FMO_diagram->show_connects = 1;
    FMO_diagram->level_width = FMO_LEVEL_WIDTH;
    FMO_diagram->level_thickness = FMO_LEVEL_THICKNESS;
    FMO_diagram->electron_length = FMO_ELECTRON_LENGTH;
    FMO_diagram->degen_tol = FMO_DEGEN_TOL;
    if( FMO_diagram->num_frags >= 1 ){
      FMO_diagram->left_fragment = 1;
    } else FMO_diagram->left_fragment = -1;
    if( FMO_diagram->num_frags >= 2 ){
      FMO_diagram->right_fragment = 2;
    } else FMO_diagram->right_fragment = -1;
    FMO_diagram->do_y_tics = 1;
    FMO_diagram->show_box = 1;
    FMO_diagram->show_data = 1;

    preprocess_FMO_data(FMO_diagram);

    whichobj->scale.x=whichobj->scale.y=whichobj->scale.z=2.0*GRAPHICS_SCALE;
    whichobj->cent.x=-180.0;whichobj->cent.y=180.0;
    whichobj->cent.z=0*GRAPHICS_SCALE;
    whichobj->trans.x=0;whichobj->trans.y=0;
    whichobj->trans.z=0;

#ifdef INTERACTIVE_USE
    /* make the button window */
    build_FMO_button_window(&button_wins,FMO_diagram);
    whichobj->prim->but_win = button_wins;
#endif
  }
}

/****************************************************************************
 *
 *                   Procedure FMO_draw_electrons
 *
 * Arguments: diagram: pointer to FMO_diagram_type
 *              level: pointer to FMO_level_type
 *
 * Returns: none
 *
 * Action: draws in the electrons in level 'level.  This deals
 *   properly with partially occupied and degenerate levels.
 *
 ****************************************************************************/
void FMO_draw_electrons(FMO_diagram_type *diagram, FMO_level_type *level)
{
  int i;
  float level_center,xpos;
  float num_electrons;

/*  num_electrons = ceil(level->num_electrons);*/
  num_electrons = level->num_electrons;

  if( num_electrons/level->num_degen > 2 ){
    char errstring[80];
    sprintf(errstring,"Level at E= %lf, has a filling (%lf) which is larger than 2.0",
            level->energy,num_electrons);
    FATAL_BUG(errstring);
  }

  level_center = level->xmin + diagram->level_width/2.0;

  /* if the level is full, then this is easy */
  if( num_electrons == 2.0*level->num_degen ){
    for(i=0;i<level->num_degen;i++){
      xpos = level_center - diagram->level_width/4.0;
      g_line(xpos,level->yloc-diagram->electron_length/2.0,
             xpos,level->yloc+diagram->electron_length/2.0);
      xpos = level_center + diagram->level_width/4.0;
      g_line(xpos,level->yloc-diagram->electron_length/2.0,
             xpos,level->yloc+diagram->electron_length/2.0);

      level_center += (float)(diagram->level_width + FMO_DEGEN_LEVEL_SKIP);
    }
  } else {
    /********

      the level isn't fully occupied.... deal with it by looping
      through the sublevels twice, adding electrons as we go.

    ********/
    level_center = level->xmin + diagram->level_width/2.0;
    for(i=0;i<level->num_degen && num_electrons > 0.01;i++){
      xpos = level_center - diagram->level_width/4.0;
      g_line(xpos,level->yloc-diagram->electron_length/2.0,
             xpos,level->yloc+diagram->electron_length/2.0);

      level_center += (float)(diagram->level_width + FMO_DEGEN_LEVEL_SKIP);
      num_electrons -= 1.0;
    }

    /* loop the second time if we need to */
    if( num_electrons > 0.1 ){
      level_center = level->xmin + diagram->level_width/2.0;
      for(i=0;i<level->num_degen && num_electrons > 0.1;i++){
        xpos = level_center + diagram->level_width/4.0;
        g_line(xpos,level->yloc-diagram->electron_length/2.0,
               xpos,level->yloc+diagram->electron_length/2.0);

        level_center += (float)(diagram->level_width + FMO_DEGEN_LEVEL_SKIP);
        num_electrons -= 1.0;
      }
    }
  }
}



/****************************************************************************
 *
 *                   Procedure draw_FMO_diagram
 *
 * Arguments: prim: pointer to prim_type
 *             obj: pointer to object_type
 *
 * Returns: none
 *
 * Action: draws in an FMO interaction diagram
 *
 ****************************************************************************/
void draw_FMO_diagram(prim_type *prim,object_type *obj)
{
  int i;
  float xpos,ypos;
  float xcenter;
  point_type2D origin,dim,con_begin,con_end;
  char numstring[20];
  float xloc,yloc;
  float inv_yscale;
  float yscale;
  float yref;
  float xval,yval;
  int xoffset;

  FMO_level_type *level,*main_level,*fragment_level;
  FMO_connect_type *connector;
  FMO_diagram_type *diagram;

  /* a little error checking */
  if( !prim->FMO_diagram )
    fatal("draw_FMO_diagram called with invalid primitive.  This is a bug!");

  diagram = prim->FMO_diagram;

  /* check to see if we need to re-determine the location of tic marks */
  if( diagram->old_max_y != diagram->max_y || diagram->old_min_y != diagram->min_y){
    find_FMO_tic_sep(diagram);
    diagram->old_max_y = diagram->max_y;
    diagram->old_min_y = diagram->min_y;
  }

  /* determine the location of the origin on screen */
  origin.x = obj->cent.x + obj->trans.x + g_xmax / 2;
  origin.y = obj->cent.y - obj->trans.y + g_ymax / 2;
  dim.x = DEF_GRAPH_X * obj->scale.x;
  dim.y = DEF_GRAPH_Y * obj->scale.y;

  /* scaling terms */
  inv_yscale = (diagram->max_y - diagram->min_y) / DEF_GRAPH_Y;
  yscale = obj->scale.y/inv_yscale;

  /* find the point in data space which will appear at the origin */
  yref = diagram->min_y;

  /* determine the bounding box for this object */
  localmin.x = obj->bmin.x = origin.x;
  localmin.y = obj->bmin.y = origin.y - dim.y;
  localmax.x = obj->bmax.x = origin.x + dim.x;
  localmax.y = obj->bmax.y = origin.y;

  /* Y tics */
  if( diagram->do_y_tics ){
    for(i=0;i<floor(diagram->num_tics_y);i++){
      yloc = origin.y + yscale * (yref - diagram->tic_start_y -
                                        i * diagram->tic_sep_y);
      g_line(origin.x-obj->scale.x*TIC_DIM,yloc,origin.x,yloc);

      yval = (diagram->tic_start_y + i*diagram->tic_sep_y);
      if( fabs(yval) < 1e-12 ) yval = 0.0;
      sprintf(numstring,"%lg",yval);
      xloc = origin.x-obj->scale.x*TIC_DIM*1.1;

      g_right_text(xloc,yloc,numstring);
    }
  }

  if( diagram->ylegend[0] != 0 && diagram->do_y_tics ){
    yloc = origin.y - dim.y/2;
    xloc = origin.x-obj->scale.x*TIC_DIM*5.0;

    g_ylegend(xloc,yloc,diagram->ylegend);

  }


  /* Now do the title */
  if( diagram->title[0] != 0 && diagram->do_title){
    xloc = origin.x + dim.x/2;
    yloc = origin.y - dim.y;

    g_title(xloc,yloc,diagram->title);
  }

  /* do the levels of the molecule */
  xcenter = origin.x + dim.x/2;

  if( diagram->show_data ){
    g_change_linewidth((float)diagram->level_thickness);

    level = diagram->levels;

    while(level){
      /* figure out how far out from the center we have to move */
      xoffset = (level->num_degen*diagram->level_width +
                 (level->num_degen-1)*FMO_DEGEN_LEVEL_SKIP)/2;
      level->xmin = xcenter-xoffset;
      level->xmax = xcenter+xoffset;
      xpos = level->xmin;

      /* figure out the y position */
      level->yloc = ypos = origin.y + (yref - level->energy)*yscale;

      /* now draw the levels (if they're onscreen) */
      if( ypos > obj->bmin.y && ypos < obj->bmax.y ){
        for(i=0;i<level->num_degen;i++){
          g_line(xpos,ypos,xpos+diagram->level_width,ypos);
          xpos += diagram->level_width + FMO_DEGEN_LEVEL_SKIP;
        }
      }

      /* move onto the next level */
      level = level->next;
    }

    /* draw in the label */
    if( diagram->label ){
      yloc = origin.y+obj->scale.y*TIC_DIM*(1.1);
      g_center_text(xcenter,yloc,diagram->label);
    }


    /* move over to 1/5 the way from the left side and do the left fragment */
    xcenter = origin.x + dim.x/5;

    if( diagram->left_fragment-1 < diagram->num_frags && diagram->left_fragment > 0) {
      level = diagram->frags[diagram->left_fragment-1].levels;
    }else{
      level = 0;
    }

    while(level){
      /* figure out how far out from the center we have to move */
      xoffset = (level->num_degen*diagram->level_width +
                 (level->num_degen-1)*FMO_DEGEN_LEVEL_SKIP)/2;
      level->xmin = xcenter-xoffset;
      level->xmax = xcenter+xoffset;
      xpos = level->xmin;

      /* figure out the y position */
      level->yloc = ypos = origin.y + (yref - level->energy)*yscale;

      /* now draw the levels (if they're onscreen) */
      if( ypos > obj->bmin.y && ypos < obj->bmax.y ){
        for(i=0;i<level->num_degen;i++){
          g_line(xpos,ypos,xpos+diagram->level_width,ypos);
          xpos += diagram->level_width + FMO_DEGEN_LEVEL_SKIP;
        }
      }

      /* move onto the next level */
      level = level->next;
    }

    /* draw in the label */
    if( diagram->left_fragment-1 < diagram->num_frags && diagram->left_fragment > 0 &&
       diagram->frags[diagram->left_fragment-1].label ){
      yloc = origin.y+obj->scale.y*TIC_DIM*(1.1);
      g_center_text(xcenter,yloc,
                    diagram->frags[diagram->left_fragment-1].label);
    }


    /* move over to 4/5 the way from the left side and do the left fragment */
    xcenter = origin.x + 4*dim.x/5;

    if( diagram->right_fragment-1 < diagram->num_frags && diagram->right_fragment > 0){
      level = diagram->frags[diagram->right_fragment-1].levels;
    }else{
      level = 0;
    }

    while(level){
      /* figure out how far out from the center we have to move */
      xoffset = (level->num_degen*diagram->level_width +
                 (level->num_degen-1)*FMO_DEGEN_LEVEL_SKIP)/2;
      level->xmin = xcenter-xoffset;
      level->xmax = xcenter+xoffset;
      xpos = level->xmin;

      /* figure out the y position */
      level->yloc = ypos = origin.y + (yref - level->energy)*yscale;

      /* now draw the levels (if they're onscreen) */
      if( ypos > obj->bmin.y && ypos < obj->bmax.y ){
        for(i=0;i<level->num_degen;i++){
          g_line(xpos,ypos,xpos+diagram->level_width,ypos);
          xpos += diagram->level_width + FMO_DEGEN_LEVEL_SKIP;
        }
      }

      /* move onto the next level */
      level = level->next;
    }

    /* draw in the label */
    if( diagram->right_fragment-1 < diagram->num_frags && diagram->right_fragment > 0 &&
       diagram->frags[diagram->right_fragment-1].label ){
      yloc = origin.y+obj->scale.y*TIC_DIM*(1.1);
      g_center_text(xcenter,yloc,
                    diagram->frags[diagram->right_fragment-1].label);
    }


    /*************

      draw in the connectors now

    *************/
    if( diagram->show_connects ){
      g_change_linewidth(1);
      connector = diagram->connects;
      while(connector){
        /* find the fragment level to which this guy connects*/
        if( connector->which_frag == (diagram->left_fragment-1) ||
           connector->which_frag == (diagram->right_fragment-1) ){

          fragment_level = diagram->frags[connector->which_frag].levels;

          while(fragment_level && fragment_level->number != connector->fragment_level){
            fragment_level = fragment_level->next;
          }
          if( !fragment_level ) fatal("Can't find a fragment level for a connector.");

          /* find the main level to which this guy connects */
          main_level = diagram->levels;
          while(main_level && main_level->number != connector->main_level){
            main_level = main_level->next;
          }
          if( !main_level ) fatal("Can't find a main level for a connector.");

          /* determine the location of the connector end points */
          if( connector->which_frag == (diagram->left_fragment-1) ){
            con_begin.x = fragment_level->xmax;
            con_end.x = main_level->xmin;
          }else{
            con_begin.x = fragment_level->xmin;
            con_end.x = main_level->xmax;
          }

          /*******

            the y location of the end points is trickier, because we want
            the connectors to go to the edge of the bounding box if one of
            the levels to which they are connected isn't on screen

            ********/
          con_begin.y = fragment_level->yloc;
          con_end.y = main_level->yloc;

          /* this is the case where we don't have to draw the connector at all */
          if( (con_begin.y < obj->bmin.y || con_begin.y > obj->bmax.y) &&
             (con_end.y < obj->bmin.y || con_end.y > obj->bmax.y) ){
            con_begin.y = con_end.y = -1000;
          } else if( (con_begin.y > obj->bmin.y && con_begin.y < obj->bmax.y) &&
                    (con_end.y < obj->bmin.y || con_end.y > obj->bmax.y) ){
            /*****
              the beginning is onscreen, but the end point is off.
              leave the beginning alone and figure out where we should
              draw the end to.
              ******/
            if( con_end.y < obj->bmin.y ) yval = obj->bmin.y;
            else yval = obj->bmax.y;

            xval = con_begin.x + (yval - con_begin.y)*(con_end.x - con_begin.x)/
              (con_end.y-con_begin.y);

            con_end.x = xval;
            con_end.y = yval;
          } else if( (con_begin.y < obj->bmin.y || con_begin.y > obj->bmax.y) &&
                    (con_end.y > obj->bmin.y && con_end.y < obj->bmax.y) ){
            /*****
              the end is onscreen, but the beginning point is off.
              leave the end alone and figure out where we should
              draw the beginning to.
              ******/
            if( con_begin.y < obj->bmin.y ) yval = obj->bmin.y;
            else yval = obj->bmax.y;

            xval = con_end.x + (yval - con_end.y)*(con_end.x - con_begin.x)/
              (con_end.y-con_begin.y);

            con_begin.x = xval;
            con_begin.y = yval;
          }

          /* draw it in */
          g_change_linestyle(connector->linestyle);
          if(con_begin.y >= -100){
            g_line(con_begin.x,con_begin.y,con_end.x,con_end.y);
          }
        }
        connector = connector->next;
      }
    }
    g_change_linewidth(1);
  } /* end of if(diagram->show_data) */
  /*******

    if we need to, draw in the electron filling

  *******/
  if( diagram->electron_filling_mode != FMO_FILL_NONE ){
    g_change_linestyle(0);

    /* do the main energy levels */
    level = diagram->levels;

#if 0
    while( level && !(level->highest) ){
      if( diagram->electron_filling_mode == FMO_FILL_ALL ){
        if( level->yloc > obj->bmin.y && level->yloc < obj->bmax.y ){
          FMO_draw_electrons(diagram,level);
        }
      }

      level = level->next;
    }

    /* fill the highest occupied level */
    if( level ){
      if( level->yloc > obj->bmin.y && level->yloc < obj->bmax.y ){
        FMO_draw_electrons(diagram,level);
      }
    }

#endif
    while( level ){
      if( diagram->electron_filling_mode == FMO_FILL_ALL || level->highest){
        if( level->yloc > obj->bmin.y && level->yloc < obj->bmax.y ){
          FMO_draw_electrons(diagram,level);
        }
      }

      level = level->next;
    }

    /* do the energy levels of the fragments*/
    if( diagram->left_fragment-1 < diagram->num_frags && diagram->left_fragment > 0){
      level = diagram->frags[diagram->left_fragment-1].levels;
    }else{
      level = 0;
    }
#if 0
    while( level && !(level->highest) ){
      if( diagram->electron_filling_mode == FMO_FILL_ALL ){
        if( level->yloc > obj->bmin.y && level->yloc < obj->bmax.y ){
          FMO_draw_electrons(diagram,level);
        }
      }
      level = level->next;
    }
    if( level ){
      if( level->yloc > obj->bmin.y && level->yloc < obj->bmax.y ){
        FMO_draw_electrons(diagram,level);
      }
    }
#endif
    while( level ){
      if( diagram->electron_filling_mode == FMO_FILL_ALL || level->highest ){
        if( level->yloc > obj->bmin.y && level->yloc < obj->bmax.y ){
          FMO_draw_electrons(diagram,level);
        }
      }
      level = level->next;
    }

    if( diagram->right_fragment-1 < diagram->num_frags && diagram->right_fragment > 0){
      level = diagram->frags[diagram->right_fragment-1].levels;
    }else{
      level = 0;
    }
    while( level && !(level->highest) ){
      if( diagram->electron_filling_mode == FMO_FILL_ALL ){
        if( level->yloc > obj->bmin.y && level->yloc < obj->bmax.y ){
          FMO_draw_electrons(diagram,level);
        }
      }
      level = level->next;
    }
    if( level ){
      if( level->yloc > obj->bmin.y && level->yloc < obj->bmax.y ){
        FMO_draw_electrons(diagram,level);
      }
    }
  }

  /* draw a box around it */
  g_change_linestyle(0);
  if(diagram->show_box){
    g_rectangle(origin.x,origin.y,dim.x,dim.y);
  }
}
