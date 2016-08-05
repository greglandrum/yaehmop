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
   06.04.98 gL:
     Fixed stupid problem in show_selected_data... not
      being able to find all selected atoms is *not* a
      fatal_bug.  This happens all the time if molecules
      are deleted or crystals are grown... it's just not
      that big of a deal.  Changed the fatal_bug to an
      info message.
   09.07.98 gL:
      added code to generate rotation angles for bonds in VRML
   07.09.98 gL:
      dump_molecular_coords no longer dumps the coordinates of
      hidden atoms.
  09.09.98 gL:
    added invert_selected_atoms
  24.09.98 gL:
    initialize new fields of labels in show_selected_data
  26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)
  12.03.99 gL:
     fixed up some potential buffer overrun problems in the read
     functions.
  26.06.99 gL:
     hide_selected_atoms did not properly toggle the exclude variable
     in the MO_centers array.  This has been fixed.

***/


#include "viewkel.h"
#include "matrix_ops.h"

/* Procedure fatal_bug
 * prints an error message and terminates the program.
 *  this should be called when internal consistency error checking
 *   fails.
 */
void fatal_bug( char *errorstring, char *file, int line )
{
  fprintf( stderr, "FATAL ERROR: %s.\n",errorstring );
  fprintf( stderr, "The error occured at line: %d of file: %s\n",line,file );
  fprintf( stderr, "  This is a bug.  Please tell Greg that you saw this.\n");
  fprintf( stderr, "Execution Terminated.\n");

#ifndef USING_THE_MAC
  exit(-12);
#else
        scanf("/n/n/n");
#endif
}


/* Procedure fatal
 * prints an error message and terminates the program
 */
void fatal( char *errorstring )
{
  fprintf( stderr, "FATAL ERROR: %s.\nExecution Terminated.\n",
      errorstring );
#ifndef USING_THE_MAC
  exit(-1);
#else
        scanf("/n/n/n");
#endif
}

/* Procedure error
 * prints an error message
 */
void error( char *errorstring )
{
  fprintf( stderr, "ERROR: %s.\n", errorstring );
  return;
}

/****************************************************************************
*
*                   Procedure printmat
*
* Arguments: mat : matrix_type
* Returns: none
*
* Action: prints out matrix mat
*
*****************************************************************************/
void printmat( matrix_type *mat )
{
  int i,j;

  for( i=0; i<DIM; i++ ){
    for( j=0; j<DIM; j++ ) printf("\t%6.4g ",mat->matrix[i][j]);
    printf("\n");
  }
}



/****************************************************************************
 *
 *                   Procedure display
 *
 * Arguments:   string: pointer to type char
 * Returns:  none
 *
 * Action: Displays the string passed in in the bottom of the button window.
 *
 ****************************************************************************/
void display( char *string )
{
    int pos;

#ifdef X_GRAPHICS
    if( doing_X ){
      /* first clear out anything that was there */
      XFillRectangle(disp,gwin,blackgc,0,
                     g_ymax-BUTTONOFF-30,BUTWIDTH+30,30);

      /* now draw in the string (centered) */
      pos = BUTWIDTH / 2 - XTextWidth(big_font,string,strlen(string))/2;
      XDrawImageString(disp,gwin,bigtextgc,pos,g_ymax-BUTTONOFF,
                       string,strlen(string));

      /* flush the buffer to make sure that the string is displayed */
      XFlush(disp);
    }
#endif
#ifdef TEK_GRAPHICS
    if( doing_tek ){
      fprintf(stderr,"%s\n",string);
    }
#endif
    /* that's it */
 }



/****************************************************************************
 *
 *                   Procedure readfloatparm
 *
 * Arguments:   string: a pointer to char
 *                parm: a pointer to a float
 * Returns:  None
 *
 * Action: Reads in the floating point parameter 'parm from the user.
 *
 ****************************************************************************/
void readfloatparm( char *string, float *parm )
{
  char outstring[MAX_STR_LEN],instring[MAX_STR_LEN];


  /* display the prompt */
  display("Look in Xterm");
  strcpy(outstring,"\nEnter the new value of ");
  printf("%s%s [%lg] ",outstring,string,*parm);

  /* read in the value */
  fgets(instring,80,stdin);

  /* check to see if the current value should be kept */
  if(instring[0] != '\n' && instring[0] != 0){
    *parm = (float)atof(instring);
    strcpy(outstring,string);
    strcat(outstring," changed.");
  }
  else{
    strcpy(outstring,"No change to ");
    strcat(outstring,string);
  }
  display(outstring);
}
/****************************************************************************
 *
 *                   Procedure readintparm
 *
 * Arguments:   string: a pointer to char
 *                parm: a pointer to a int
 * Returns:  None
 *
 * Action: Reads in the integer parameter 'parm from the user.
 *
 ****************************************************************************/
void readintparm( char *string, int *parm )
{
  char outstring[MAX_STR_LEN],instring[MAX_STR_LEN];

  /* display the prompt */
  display("Look in Xterm");
  strcpy(outstring,"\nEnter the new value of ");
  printf("%s%s [%d] ",outstring,string,*parm);

  /* read in the value */
  fgets(instring,80,stdin);

  /* check to see if the current value should be kept */
  if(instring[0] != '\n' && instring[0] != 0){
    sscanf(instring,"%d",parm);
    strcpy(outstring,string);
    strcat(outstring," changed.");
  }
  else{
    strcpy(outstring,"No change to ");
    strcat(outstring,string);
  }
  display(outstring);
}


/****************************************************************************
 *
 *                   Procedure readcharparm
 *
 * Arguments:   string: a pointer to char
 *                parm: a pointer to a char
 * Returns:  None
 *
 * Action: Reads in the character parameter 'parm from the user.
 *
 ****************************************************************************/
void readcharparm( char *string, char *parm )
{
  char outstring[MAX_STR_LEN],instring[MAX_STR_LEN];

  /* display the prompt */
  display("Look in Xterm");
  strcpy(outstring,"\nEnter the new value of ");
  printf("%s%s [%c] ",outstring,string,*parm);

  /* read in the value */
  fgets(instring,80,stdin);

  /* check to see if the current value should be kept */
  if(instring[0] != '\n' && instring[0] != 0){
    sscanf(instring,"%c",parm);
    strcpy(outstring,string);
    strcat(outstring," changed.");
  }
  else{
    strcpy(outstring,"No change to ");
    strcat(outstring,string);
  }
  display(outstring);
}


/****************************************************************************
 *
 *                   Procedure readstringparm
 *
 * Arguments:   string: a pointer to char
 *          stringparm: array of pointers to char
 * Returns:  None
 *
 * Action: Reads in single line for 'stringparm from the user
 *
 ****************************************************************************/
void readstringparm( char *string, char **stringparm )
{
  char instring[240];

  /* display the prompt */
  display("Look in Xterm");
  printf("Enter the new value for %s [%s]: ",string,stringparm[0]);

  fgets(instring,240,stdin);

  /* check to see if we are done */
  if(instring[0] == '\n' || instring[0] == 0) return;
  else{
    strcpy(stringparm[0],instring);
    /* kill the carriage return at the end */
    stringparm[0][strlen(stringparm[0])-1] = 0;
  }

  display("Got it!");
}


/****************************************************************************
 *
 *                   Procedure readmultistringparm
 *
 * Arguments:   string: a pointer to char
 *           num_lines: int
 *          stringparm: array of pointers to char
 * Returns:  None
 *
 * Action: Reads in the 'num_lines long string 'stringparm from the user
 *
 ****************************************************************************/
void readmultistringparm( char *string, int num_lines, char **stringparm )
{
  int i,more;
  char instring[NORMAL_STR_LEN];

  /* display the prompt */
  display("Look in Xterm");
  printf("The %s consists of up to %d lines.\n",string,num_lines);
  printf("Enter each line. Hit <enter> to leave a line blank and skip the rest.\n");

  /* read in each of the strings */
  more = 1;
  for(i=0;i<num_lines;i++){
    /* read in the string */
    if( more ) fgets(instring,NORMAL_STR_LEN,stdin);

    /* check to see if we are done */
    if(instring[0] == '\n' || instring[0] == 0) more = 0;

    /* if we're not done, copy in the new string */
    if( more ){
      strcpy(stringparm[i],instring);

      /* kill the carriage return at the end */
      stringparm[i][strlen(stringparm[i])-1] = 0;
    }
    else {
      /* we are done, zero out the string */
      bzero(stringparm[i],NORMAL_STR_LEN*sizeof(char));
    }
  }

  printf("<done>\n");
  display("Got it!");

}

/****************************************************************************
*
*                   Function safe_strcpy
*
* Arguments:     str1, str2: pointers to type char
* Returns: pointer to char
*
* Action: This is a version of strcpy that sets the first
*  character of 'str1 to zero if 'str2 is NULL.  This is
*  to make the use of strtok a little easier.
*
*****************************************************************************/
char *safe_strcpy(char *str1,char *str2)
{
  if( !str1 ) fatal("safe_strcpy called with null str1");
  if( !str2 ){
    str1[0] = 0;
  } else{
    strcpy(str1,str2);
  }
  return(str1);
}


/****************************************************************************
 *
 *                   Procedure skipcomments
 *
 * Arguments: file : a pointer to file type
 *        instring : pointer to type char
 * Returns: an integer
 *
 * Action: Reads in lines from 'file' until one is hit that does not begin
 *     with a ; or a return. puts the first non-comment line into instring
 *     and then returns.  Returns negative if EOF is hit.
 *
 ****************************************************************************/
int skipcomments(FILE *file,char *string)
{
  char *error;

  /* use the first element of string to check for EOF */
  string[0] = 0;
  error = fgets(string,MAX_STR_LEN,file);
  while( error && (string[0] == '\n' || string[0] == ';')
        && string[0] != 0 ){
    string[0] = 0;
    error = fgets(string,MAX_STR_LEN,file);
  }

  if( error && string[0] != 0 ) return(0);
  else{
    return(-1);
  }

}

/*********
  converts a string to all uppercase
**********/
void upcase(char *string)
{
  int i,len,diff;

  diff = 'A' - 'a';

  len = strlen(string);
  for(i=0;i<len;i++){
    /* check to see if its a letter */
    if( string[i] >= 'a' && string[i] <= 'z' ){
      string[i] += diff;
    }
  }
}


void dump_molecular_coords(molec_type *molec)
{
  atom_type *atoms;
  int i;

  if(molec->num_frames > 1)
    atoms = &(molec->atoms[(molec->current_frame%molec->num_frames)*molec->num_atoms]);
  else atoms = molec->atoms;

  for(i=0;i<molec->num_atoms;i++){
    if(!atoms[i].exclude &&
       (molec->hydrogens_on || atoms[i].type[0] != 'H' ||
        atoms[i].type[1] != 0) &&
       (molec->dummies_on || atoms[i].type[0] != '&') ){
      printf("%d %s % -6.4lf % -6.4lf % -6.4lf\n",i+1,
             atoms[i].type,
             atoms[i].loc.x,atoms[i].loc.y,atoms[i].loc.z);
    }
  }
}

#ifdef NEED_DRAND
/*********************************

  function: my_drand

  argument: bound: float
  returns: float

  action: generates a random float between +/- 'bound

**********************************/
float my_drand(float bound)
{
  /* this can be used if drand48 exists and works */
  return(bound - 2.0*bound*drand48());
}
#endif
/****************************************************************************
 *
 *                   Procedure parse_integer_string
 *
 * Arguments: string: pointer to type char
 *            values: pointer to pointer to int
 *        num_values: pointer to int
 *
 *
 * Returns: none
 *
 * Action: reads a comma delimited list of integers out of 'string.
 *   returns these integers in 'values.  The special thing here is that
 *   this deals with entries of the form "3-7", parsing this into 3,4,5,6,7.
 *
 ****************************************************************************/
void parse_integer_string(char *string,int **values,int *num_values)
{
  int max_values;
  char local_string[MAX_STR_LEN],num_string[80];
  int i;
  int num1,num2;
  char foo_char;

  *num_values = 0;
  max_values = 10;

  /* get some initial memory */
  *values = (int *)D_CALLOC(max_values,sizeof(int));
  if( !(*values )) fatal("Can't allocate values in parse_integer_string");

  /* first copy the input string */
  safe_strcpy(local_string,string);

  /* now use strtok to chop it up */
  safe_strcpy(num_string,(char *)strtok(local_string,",\n"));

  while(num_string[0]){
    /* check to see if there's a - we need to deal with */
    if( strstr(num_string,"-") ){
      sscanf(num_string,"%d%c%d",&num1,&foo_char,&num2);
    }else{
      sscanf(num_string,"%d",&num1);
      num2 = num1;
    }
    for(i=num1;i<=num2;i++){
      (*values)[*num_values] = i;
      (*num_values)++;
      if(*num_values == max_values){
        max_values += 10;
        *values = (int *)D_REALLOC(*values,max_values*sizeof(int));
        if( !(*values )) fatal("Can't D_REALLOCate values in parse_integer_string");
      }
    }
    safe_strcpy(num_string,(char *)strtok(0,",\n"));
  }
}


/****************************************************************************
 *
 *                   Procedure hide_selected_atoms
 *
 * Arguments: num_selected: int
 *                     obj: pointer to object_type
 *
 *
 * Returns: none
 *
 * Action:  hides the atoms selected
 *
 ****************************************************************************/
void hide_selected_atoms(int num_selected,object_type *obj)
{
  atom_type **selected_atoms;
  molec_type *molec;
  int i,num_found,which;


  if(obj->prim->molec){
    molec = obj->prim->molec;
  }else if( obj->prim->MO_surf && obj->prim->MO_surf->molec ){
    molec = obj->prim->MO_surf->molec;
  } else{
    return;
  }

  selected_atoms = (atom_type **)D_CALLOC(num_selected,
                                        sizeof(atom_type *));
  if(!selected_atoms)fatal("Can't get memory for selected atoms array");

  num_found = 0;
  for(i=0;i<molec->num_atoms && num_found < num_selected;i++){
    if(molec->atoms[i].is_selected){
      /******

        insert a pointer to the atom into the selected_atoms
        array.  Put it in the slot corresponding to its selection
        order.

        ******/
      selected_atoms[molec->atoms[i].is_selected-1] = &(molec->atoms[i]);
      num_found++;
    }
  }
  if( num_found != num_selected )
    FATAL_BUG("Can't find all the selected atoms!");

  for(i=0;i<num_selected;i++){
    selected_atoms[i]->exclude = 1;
    selected_atoms[i]->is_selected = 0;
    if( obj->prim->MO_surf ){
      which = selected_atoms[i]->num;
      obj->prim->MO_surf->MO_centers[which].exclude = 1;
    }
  }

  D_FREE(selected_atoms);
}


/****************************************************************************
 *
 *                   Procedure show_all_atoms
 *
 * Arguments:   obj: pointer to object_type
 *
 *
 * Returns: none
 *
 * Action:  shows all the atoms
 *
 ****************************************************************************/
void show_all_atoms(object_type *obj)
{
  molec_type *molec;
  int i;

  if(obj->prim->molec){
    molec = obj->prim->molec;
  }else if( obj->prim->MO_surf && obj->prim->MO_surf->molec ){
    molec = obj->prim->MO_surf->molec;
  } else{
    return;
  }

  for(i=0;i<molec->num_atoms;i++){
    if(molec->atoms[i].exclude){
      num_selected++;
      molec->atoms[i].is_selected = num_selected;
    }
    molec->atoms[i].exclude = 0;
    if( obj->prim->MO_surf ) obj->prim->MO_surf->MO_centers[i].exclude = 0;
  }
}


/****************************************************************************
 *
 *                   Procedure show_selected_data
 *
 * Arguments: num_selected: int
 *                     obj: pointer to object_type
 *
 *
 * Returns: none
 *
 * Action:  shows relevant information about 'obj.
 *  At the moment this only deals with molecules, showing
 *   distances or angles, depending upon the number of
 *   atoms selected.
 *
 ****************************************************************************/
void show_selected_data(int num_selected,object_type *obj,int xloc,int yloc)
{
  atom_type **selected_atoms;
  molec_type *molec;
  label_type *label;
  atom_type *atoms;
  char label_string[80];
  int i,num_found;
  double val;
  float R0_val,valence;

  if(!num_selected) return;
  if(obj->prim->molec){
    molec = obj->prim->molec;
    if(molec->num_frames > 1)
      atoms = &(molec->atoms[molec->current_frame*molec->num_atoms]);
    else atoms = molec->atoms;
  }else if( obj->prim->MO_surf && obj->prim->MO_surf->molec ){
    molec = obj->prim->MO_surf->molec;
    atoms = molec->atoms;
  } else{
    return;
  }

  selected_atoms = (atom_type **)D_CALLOC(num_selected,
                                        sizeof(atom_type *));
  if(!selected_atoms)fatal("Can't get memory for selected atoms array");

  num_found = 0;
  for(i=0;i<molec->num_atoms && num_found < num_selected;i++){
    if(atoms[i].is_selected){
      /******

        insert a pointer to the atom into the selected_atoms
        array.  Put it in the slot corresponding to its selection
        order.

        ******/
      selected_atoms[atoms[i].is_selected-1] = &(atoms[i]);
      num_found++;
    }
  }
  if( num_found != num_selected ){
    num_selected = num_found;
    fprintf(stderr,"Info: number of selected atoms has changed.\n");
  }

  /* now process the selected atoms... */
  switch(num_selected){
  case 1:
    printf("Atom: %d %s (%lf %lf %lf)\n",selected_atoms[0]->num+1,
           selected_atoms[0]->type,selected_atoms[0]->loc.x,
           selected_atoms[0]->loc.y,selected_atoms[0]->loc.z);
    break;
  case 2:
    val = V3DistanceBetween2Points(&selected_atoms[0]->loc,
                                   &selected_atoms[1]->loc);
    printf("Dist between atoms %s(%d) and %s(%d): %lf\n",
           selected_atoms[0]->type,selected_atoms[0]->num+1,
           selected_atoms[1]->type,selected_atoms[1]->num+1,
           val);
#ifdef INCLUDE_BOND_VALENCE
    /* calculate the bond valence too, because we can */
    R0_val = 0.0;
    bond_length_to_bond_valence(selected_atoms[0],selected_atoms[1],
                                val,&R0_val,&valence);
    printf("\tBond valence: %6.3lf based on R0= %6.3lf\n",valence,R0_val);
#endif
    sprintf(label_string,"%6.3lf \\AA",val);
    new_label(label_string);
    label = head->obj->prim->label;
    head->obj->cent.x = xloc;
    head->obj->cent.y = yloc;
    label->num_atoms_labelled = 2;
    label->atoms_to_label = (atom_type **)D_CALLOC(2,sizeof(atom_type *));
    if(!label->atoms_to_label) fatal("can't get memory for atoms_to_label");
    label->atoms_to_label[0] = selected_atoms[0];
    label->atoms_to_label[1] = selected_atoms[1];
    label->show_lines = 1;
    break;
  case 3:
    val = V3AngleBetween3Points(&selected_atoms[0]->loc,
                                &selected_atoms[1]->loc,
                                &selected_atoms[2]->loc);
    val = val*180.0/M_PI;
    printf("Angle %s(%d)--%s(%d)--%s(%d): %lf\n",
           selected_atoms[0]->type,selected_atoms[0]->num+1,
           selected_atoms[1]->type,selected_atoms[1]->num+1,
           selected_atoms[2]->type,selected_atoms[2]->num+1,
           val);

    sprintf(label_string,"%6.1lf deg",val);
    new_label(label_string);
    label = head->obj->prim->label;
    head->obj->cent.x = xloc;
    head->obj->cent.y = yloc;
    label->num_atoms_labelled = 3;
    label->atoms_to_label = (atom_type **)D_CALLOC(3,sizeof(atom_type *));
    if(!label->atoms_to_label) fatal("can't get memory for atoms_to_label");
    label->atoms_to_label[0] = selected_atoms[0];
    label->atoms_to_label[1] = selected_atoms[1];
    label->atoms_to_label[2] = selected_atoms[2];
    label->show_lines = 1;

    break;
  case 4:
    val = V3DihedralAngle(&selected_atoms[0]->loc,
                          &selected_atoms[1]->loc,
                          &selected_atoms[2]->loc,
                          &selected_atoms[3]->loc);
    val = val*180.0/M_PI;
    printf("Dihedral %s(%d)--%s(%d)--%s(%d)--%s(%d): %lf\n",
           selected_atoms[0]->type,selected_atoms[0]->num+1,
           selected_atoms[1]->type,selected_atoms[1]->num+1,
           selected_atoms[2]->type,selected_atoms[2]->num+1,
           selected_atoms[3]->type,selected_atoms[3]->num+1,
           val);

    sprintf(label_string,"%6.1lf deg",val);
    new_label(label_string);
    label = head->obj->prim->label;
    head->obj->cent.x = xloc;
    head->obj->cent.y = yloc;
    label->num_atoms_labelled = 4;
    label->atoms_to_label = (atom_type **)D_CALLOC(4,sizeof(atom_type *));
    if(!label->atoms_to_label) fatal("can't get memory for atoms_to_label");
    label->atoms_to_label[0] = selected_atoms[0];
    label->atoms_to_label[1] = selected_atoms[1];
    label->atoms_to_label[2] = selected_atoms[2];
    label->atoms_to_label[3] = selected_atoms[3];
    label->show_lines = 1;
    break;

  default:
    printf("Too many atoms selected, I'm confused.\n");
  }
  if( selected_atoms ) D_FREE(selected_atoms);

}


/****************************************************************************
 *
 *                   Procedure unselect_all_atoms
 *
 * Arguments: num_selected: int
 *                     obj: pointer to object_type
 *
 *
 * Returns: none
 *
 * Action:  removes the selection information from
 *    all of the atoms.
 *
 ****************************************************************************/
void unselect_all_atoms(int num_selected,object_type *obj)
{
  int i;
  int num_atoms;
  molec_type *molec;
  atom_type *atoms;


  if( obj->prim->molec){
    molec = obj->prim->molec;
    if(molec->num_frames > 1)
      atoms = &(molec->atoms[molec->current_frame*molec->num_atoms]);
    else atoms = molec->atoms;
    num_atoms = molec->num_atoms;
  }
  else if( obj->prim->MO_surf && obj->prim->MO_surf->molec){
    atoms = obj->prim->MO_surf->molec->atoms;
    num_atoms = obj->prim->MO_surf->molec->num_atoms;
  } else{
    FATAL_BUG("unselect_all_atoms called without a molecule");
  }

  for(i=0;i<num_atoms;i++){
    atoms[i].is_selected = 0;
  }
}


/****************************************************************************
 *
 *                   Function invert_selected_atoms
 *
 * Arguments:  num_selected: pointer to int
 *           obj: pointer to object_type
 *
 *
 * Returns: none
 *
 * Action:   selected atoms are deselected and vice versa
 *
 *****************************************************************************/
void invert_selected_atoms(int *num_selected,object_type *obj)
{
  int i;
  int num,num_atoms;
  atom_type *atoms;
  molec_type *molec;

  if( obj->prim->molec){
    molec = obj->prim->molec;
    atoms = molec->atoms;
  }
  else if( obj->prim->MO_surf && obj->prim->MO_surf->molec){
    molec = obj->prim->MO_surf->molec;
    atoms = obj->prim->MO_surf->molec->atoms;
  } else{
    FATAL_BUG("unselect_all_atoms called without a molecule");
  }


  num = 1;

  if( molec->num_frames <= 1 ){
    num_atoms = molec->num_atoms;
  } else{
    num_atoms = molec->num_atoms*molec->num_frames;
  }

  /* step through the atoms */
  for(i=0;i<molec->num_atoms;i++){
    if( atoms[i].is_selected ){
      atoms[i].is_selected = 0;
    } else {
      atoms[i].is_selected = num;
      num++;
    }
  }
  num--;
  *num_selected = num;
}

#ifdef INCLUDE_ADF_PLOTS
void map_AO_number_to_center(int AO_num,MO_surface_type *MO_surf,
                             MO_center_list_type **center_p,int *offset)
{
  int i;
  int center_begin;
  int AOs_passed;
  MO_center_list_type *center;

  i=0;
  AOs_passed = 0;
  center = &MO_surf->MO_centers[i];
  center_begin = AOs_passed;
  AOs_passed += center->num_AOs;

  while(AO_num >= AOs_passed && i<MO_surf->num_centers){
    i++;
    center = &MO_surf->MO_centers[i];
    center_begin = AOs_passed;
    AOs_passed += center->num_AOs;
  }

  if(i==MO_surf->num_centers)
    FATAL_BUG("Can't map an AO_number to a center");

  *center_p = center;
  *offset = AO_num-center_begin;
}
#endif

void angles_from_bond_vect(point_type *p1,point_type *p2,
                           float *theta_y,float *theta_z,
                           float *len)
{
  point_type V;
  float dtot,dxz;
  float acos_y,acos_z;
  float t_y,t_z;

  V.x = p2->x-p1->x;
  V.y = p2->y-p1->y;
  V.z = p2->z-p1->z;

  dtot = V3Length(&V);
  dxz = sqrt(V.x*V.x+V.z*V.z);

  if( dxz > 0 ){
    acos_y = V.x / dxz;
    t_y = acos(acos_y);
  } else{
    t_y = 0;
  }

  acos_z = V.y / dtot;
  t_z = acos(acos_z);


  if( V.z < 0 ){
    t_y *= -1;
  }

  *len = dtot;
  *theta_y = t_y;
  *theta_z = t_z;
}
