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

  March '98: WG
  - f orbital number <-> name mapping [map_orb_num_to_name()]
  - f orbital labels for matrices printed in .out [print_labelled_mat()]
  - incl. f orbital parameters in .MO files [print_MOs()]
  - overlap_tab_from_vect() algorithm modified (Apr. 7th 1998)

  11.04.98, gL:
   modified fatal_bug and nonfatal_bug to take line number and file
     information.

  04.09.98, gL:
   added checks for invalid projected DOS curves and invalid
     COOP specifications to check_for_errors.

  26.09.98, gL:
   added some symmetry information to the MO output file.
   this makes dealing with it in viewkel a mite easier.

  26.10.98, gL:
   uh, I forgot something when I did that stuff above, so the MO
   symmetry printing wasn't working for walsh diagrams.  That's
   fixed now (the mirror plane toggles are static).  Yippee.

  24.02.99, gL:
    added a bogus etime function for linux installations which need it.

***/


/****************************************************************************
*
*     this file contains general purpose utility functions
*
*  created:  greg landrum  August 1993
*
*****************************************************************************/
#include "bind.h"


/* Procedure fatal_bug
 * prints an error message and terminates the program.
 *  this should be called when internal consistency error checking
 *   fails.
 */
void fatal_bug( char *errorstring, char *file, int line )
{
  fprintf( stderr, "FATAL ERROR: %s.\n",errorstring );
    fprintf( stderr, "The error occured at line: %d of file: %s\n",line,file );
  fprintf( stderr, "  This is a bug. Please report the error by e-mail to:\n");
  fprintf( stderr, "  \tyaehmop@xtended.chem.cornell.edu\n");
  fprintf( stderr, "\nExecution Terminated.\n");
  if( status_file ){
    fprintf( status_file, "FATAL ERROR: %s.\n",errorstring );
    fprintf( stderr, "The error occured at line: %d of file: %s\n",line,file );
    fprintf( status_file, "  This is a bug.  Please report the error by e-mail to:\n");
    fprintf( status_file, "  \tyaehmop@xtended.chem.cornell.edu\n");
    fprintf( status_file, "\nExecution Terminated.\n");
  }
  exit(-12);
}

/* Procedure nonfatal_bug
 * prints an error message
 *  this should be called when internal consistency error checking
 *   fails.
 */
void nonfatal_bug( char *errorstring, char *file, int line )
{
  fprintf( stderr, "ERROR: %s.\n",errorstring );
    fprintf( stderr, "The error occured at line: %d of file: %s\n",line,file );
  fprintf( stderr, "  This is a bug.  Please report the error by e-mail to:\n");
  fprintf( stderr, "  \tyaehmop@xtended.chem.cornell.edu\n");
  if(status_file){
    fprintf( status_file, "ERROR: %s.\n",errorstring );
    fprintf( stderr, "The error occured at line: %d of file: %s\n",line,file );
    fprintf( status_file, "  This is a bug.  Please report the error by e-mail to:\n");
    fprintf( status_file, "  \tyaehmop@xtended.chem.cornell.edu\n");
  }
}


/* Procedure fatal
 * prints an error message and terminates the program
 *
 * in case the job was queued, the error message is echoed to the status file
 */
void fatal( errorstring )
  char *errorstring;
{
#ifdef USING_THE_MAC
  char composite_string[255];
  const char *initstring = "FATAL ERROR: ";
  unsigned char len;
#endif
  fprintf( stderr, "FATAL ERROR: %s.\nExecution Terminated.\n",
      errorstring );
  if( status_file ){
    fprintf( status_file, "FATAL ERROR: %s.\nExecution Terminated.\n",
            errorstring );
    fflush(status_file);
  }
  if(output_file){
    fflush(output_file);
  }
#ifdef USING_THE_MAC
        bzero(composite_string,255*sizeof(char));
        strcpy((char *)composite_string,initstring);
        strcat((char *)composite_string,errorstring);
        ParamText(CtoPstr(composite_string),NULL,NULL,NULL);
        Alert(130, NULL);
#endif
  exit(-1);
}

/* Procedure error
 * prints an error message
 *
 * in case the job was queued, the error message is echoed to the status file
 */
void error( errorstring )
  char *errorstring;
{
#ifdef USING_THE_MAC
  char composite_string[255];
  const char *initstring = "ERROR: ";
  short disposition;
  unsigned char len;
#endif


  fprintf( stderr, "ERROR: %s.\n", errorstring );
  if( status_file ){
    fprintf( status_file, "ERROR: %s.\n", errorstring );
    fflush(status_file);
  }
  if( output_file) fflush(output_file);

#ifdef USING_THE_MAC
        bzero(composite_string,255*sizeof(char));
        strcpy((char *)composite_string,initstring);
        strcat((char *)composite_string,errorstring);

        ParamText(CtoPstr(composite_string),NULL,NULL,NULL);
        disposition = Alert(131, NULL);
        switch(disposition){
        case 2: fprintf(stderr,"Continuing\n");break;
        case 1: fprintf(stderr,"Terminated Execution\n");exit(555);
        }
#endif

  return;
}

/* used to deal with interrupts (^C's) */
void handle_sigint()
{
  fprintf(stderr,"Ack!  You've killed me!\n");
  fatal("User interrupt.");
}

/*********
  converts a string to all uppercase
**********/
void upcase(string)
  char *string;
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
char *safe_strcpy(str1,str2)
  char *str1, *str2;
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
*                   Procedure center_text_string
*
* Arguments:     src, dest: pointers to type char
*                 dest_len: int
* Returns: none
*
* Action: centers the string 'src in the string 'dest.
*     'dest should be at least 'dest_len+1 long and (obviously)
*     'dest_len should be longer than the string in 'src.
*
*****************************************************************************/
void center_text_string(char *src, char *dest, int dest_len)
{
  int src_len,shift,i;

  src_len = strlen(src);
  if( dest_len < src_len ) FATAL_BUG("center_text_string got a src longer than dest");
  shift = (int)floor((float)(dest_len-src_len)/2);

  for(i=0;i<shift;i++) dest[i] = ' ';
  strcpy(&(dest[shift]),src);
  for(i=shift+src_len;i<=dest_len;i++) dest[i] = ' ';
  dest[dest_len+1] = 0;
}


/****************************************************************************
*
*                   Procedure left_just_text_string
*
* Arguments:     src, dest: pointers to type char
*                 dest_len: int
* Returns: none
*
* Action: left justifies string 'src in the string 'dest.
*     'dest should be at least 'dest_len+1 long and (obviously)
*     'dest_len should be longer than the string in 'src.
*
*****************************************************************************/
void left_just_text_string(char *src, char *dest, int dest_len)
{
  int src_len,shift,i;

  src_len = strlen(src);
  if( dest_len < src_len ) FATAL_BUG("center_text_string got a src longer than dest");

  strcpy(dest,src);
  for(i=src_len;i<=dest_len;i++) dest[i] = ' ';
  dest[dest_len+1] = 0;
}



/****************************************************************************
*
*                   Procedure map_orb_num_to_name
*
* Arguments:     name: pointer to type char
*          orb_num: integer
*orbital_lookup_table: pointer to int
*                atoms: pointer to atom_type
*             num_atoms: integer
*
* Returns: none
*
* Action: Figures out the name of orbital 'orb_num
*
*****************************************************************************/
void map_orb_num_to_name(char *name,int orb_num,int *orbital_lookup_table,
                         int num_orbs,atom_type *atoms,int num_atoms)
{
  int i,j,orb_offset;
  atom_type *atom;

  if(orb_num >= num_orbs) FATAL_BUG("Bogus orb_num passed to map_orb_num_to_name.");

  for(i=0;i<num_atoms;i++){
    if(i == num_atoms-1 ){
      orb_offset = orb_num-orbital_lookup_table[i];
      atom = &(atoms[i]);
      break;
    }else{
      /* we need to find the next non-dummy atom */
      j=i+1;
      while(orbital_lookup_table[j]<0 && j<num_atoms) j++;
      if(j==num_atoms || orbital_lookup_table[j] > orb_num){
        orb_offset = orb_num-orbital_lookup_table[i];
        atom = &(atoms[i]);
        break;
      }
    }
  }



  if(atom->ns){
    if( orb_offset == 0){
      sprintf(name,"%2s(%3d) %ds",atom->symb,atom->which_atom+1,atom->ns);
      return;
    } else orb_offset--;
  }
  if( atom->np ){
    switch(orb_offset){
    case 0:
      sprintf(name,"%2s(%3d) %dpx",atom->symb,atom->which_atom+1,atom->np);
      return;
      break;
    case 1:
      sprintf(name,"%2s(%3d) %dpy",atom->symb,atom->which_atom+1,atom->np);
      return;
      break;
    case 2:
      sprintf(name,"%2s(%3d) %dpz",atom->symb,atom->which_atom+1,atom->np);
      return;
      break;
    default:
      orb_offset -= 3;
    }
  }
  if( atom->nd ){
    switch(orb_offset){
    case 0:
      sprintf(name,"%2s(%3d) %ddx2y2",atom->symb,atom->which_atom+1,atom->nd);
      return;
      break;
    case 1:
      sprintf(name,"%2s(%3d) %dz2",atom->symb,atom->which_atom+1,atom->nd);
      return;
      break;
    case 2:
      sprintf(name,"%2s(%3d) %ddxy",atom->symb,atom->which_atom+1,atom->nd);
      return;
      break;
    case 3:
      sprintf(name,"%2s(%3d) %ddxz",atom->symb,atom->which_atom+1,atom->nd);
      return;
      break;
    case 4:
      sprintf(name,"%2s(%3d) %ddyz",atom->symb,atom->which_atom+1,atom->nd);
      return;
      break;
    default:
      orb_offset -= 5;
    }
    if( atom->nf ){
      switch(orb_offset){
      case 0:
        sprintf(name,"%2s(%3d) %dz3",atom->symb,atom->which_atom+1,atom->nf);
        return;
        break;
      case 1:
        sprintf(name,"%2s(%3d) %dxz2",atom->symb,atom->which_atom+1,atom->nf);
        return;
        break;
      case 2:
        sprintf(name,"%2s(%3d) %dyz2",atom->symb,atom->which_atom+1,atom->nf);
        return;
        break;
      case 3:
        sprintf(name,"%2s(%3d) %dxyz",atom->symb,atom->which_atom+1,atom->nf);
        return;
        break;
      case 4:
        sprintf(name,"%2s(%3d) %dz(x2-y2)",atom->symb,atom->which_atom+1,atom->nf);
        return;
        break;
      case 5:
        sprintf(name,"%2s(%3d) %dx(x2-3y2)",atom->symb,atom->which_atom+1,atom->nf);
        return;
        break;
      case 6:
        sprintf(name,"%2s(%3d) %dy(3x2-y2)",atom->symb,atom->which_atom+1,atom->nf);
        return;
        break;
      default:
        orb_offset -= 7;
      }
    }
  }

  error("Can't map orbital number to atom.");
}





/****************************************************************************
*
*                   Procedure debugmat
*
* Arguments:     mat : pointer to type real
*     num_col,num_row: integers
*                 tol: a real
* Returns: none
*
* Action: prints out matrix 'mat to stdout 132 columns at a time
*   if a given number is less than 'tol then a 0 is printed instead
*
*****************************************************************************/
void debugmat( mat, num_row, num_col, tol )
  real *mat;
  int num_row,num_col;
  real tol;
{
  int i,j;
  int beg_col,end_col;
  int which_part,num_parts;
  real val;

  beg_col = 0;

  /* check to see if we're gonna have to break it into multiple parts */
  if( num_col < 8 ){
    end_col = num_col;
    num_parts = 1;
  }
  else{
    end_col = 8;
    num_parts = ceil((real)num_col/8.0);
  }

  for(which_part=0;which_part<num_parts;which_part++){
    for( i=-1; i<num_row; i++ ){

      /* start each row with a label */
      if( i >= 0 ) fprintf(stdout,"%4d",i+1);

      for( j=beg_col; j<end_col && j< num_col; j++ ){

        /* begin each column with a header */
        if( i == -1) fprintf(stdout,"  %-7d",j+1);
        else{
          val = mat[i*num_col+j];
          /* check to see if the value is less than the tolerance */
          if( fabs(val) <= tol ) val = 0.0;
          fprintf(stdout,"  %-6.4g",val);
        }
      }
      fprintf(stdout,"\n");
    }

    beg_col=end_col;
    end_col+=8;
  }
  fprintf(stdout,"\n");
  fprintf(stdout,"\n");
}


/****************************************************************************
*
*                   Procedure printmat
*
* Arguments:     mat : pointer to type real
*     num_col,num_row: integers
*             outfile: pointer to type FILE
*                 tol: a real
*             transpose: a char
*               width: int
* Returns: none
*
* Action: prints out matrix 'mat to file 'outfile 'width columns at a time
*   if a given number is less than 'tol then a 0 is printed instead
*   if 'transpose is nonzero, then the transpose of the matrix will
*    be printed instead.
*
*****************************************************************************/
void printmat( real *mat, int num_row, int num_col, FILE *outfile, real tol,
              char transpose,int width)
{
  int i,j;
  int beg_col,end_col;
  int beg_row,end_row;
  int which_part,num_parts;
  int cols_per_line;
  real val;

  /* figure out how many columns we get.... */
  cols_per_line = (int)floor((real)width/10);

  /* we need to subtract one from this to leave room for the row labels */
  cols_per_line--;

  beg_col = 0;
  beg_row = 0;

  /* check to see if we're gonna have to break it into multiple parts */
  if( !transpose ){
    if( num_col < cols_per_line ){
      end_col = num_col;
      num_parts = 1;
    }
    else{
      end_col = cols_per_line;
      num_parts = ceil((real)num_col/(real)cols_per_line);
    }
    for(which_part=0;which_part<num_parts;which_part++){
      for( i=-1; i<num_row; i++ ){

        /* start each row with a label */
        if( i >= 0 ){
          fprintf(outfile,"%4d",i+1);
        }
        else fprintf(outfile,"    ");

        for( j=beg_col; j<end_col && j< num_col; j++ ){

          /* begin each column with a header */
          if( i == -1) fprintf(outfile,"  %-7d",j+1);
          else{
            val = mat[i*num_col+j];
            /* check to see if the value is less than the tolerance */
            if( fabs(val) <= tol ) val = 0.0;
            fprintf(outfile,"  %-6.4lf",val);
          }
        }
        fprintf(outfile,"\n");
      }

      beg_col=end_col;
      end_col+=cols_per_line;
    }
  } else{
    if( num_row < cols_per_line ){
      end_row = num_row;
      num_parts = 1;
    }
    else{
      end_row = cols_per_line;
      num_parts = ceil((real)num_row/(real)cols_per_line);
    }
    for(which_part=0;which_part<num_parts;which_part++){
      for( i=-1; i<num_col; i++ ){

        /* start each row with a label */
        if( i >= 0 ){
          fprintf(outfile,"%4d",i+1);
        }
        else fprintf(outfile,"    ");

        for( j=beg_row; j<end_row && j< num_row; j++ ){

          /* begin each column with a header */
          if( i == -1) fprintf(outfile,"  %-7d",j+1);
          else{
            val = mat[j*num_col+i];
            /* check to see if the value is less than the tolerance */
            if( fabs(val) <= tol ) val = 0.0;
            fprintf(outfile,"  %-6.4lf",val);
          }
        }
        fprintf(outfile,"\n");
      }

      beg_row=end_row;
      end_row+=cols_per_line;
    }
  }
  fprintf(outfile,"\n");
  fprintf(outfile,"\n");
}


/****************************************************************************
*
*                   Procedure print_labelled_mat
*
* Arguments:     mat : pointer to type real
*     num_col,num_row: integers
*             outfile: pointer to type FILE
*                 tol: a real
*               atoms: pointer to atom_type
*           num_atoms: int
*orbital_lookup_table: pointer to int
*            num_orbs: int
*             transpose: a char
*           label_which: a char
*               width: int
* Returns: none
*
* Action: prints out matrix 'mat to file 'outfile 'width columns at a time
*   if a given number is less than 'tol then a 0 is printed instead
*   if 'transpose is nonzero, then the transpose of the matrix will
*    be printed instead.
*
*****************************************************************************/
void print_labelled_mat( real *mat, int num_row, int num_col, FILE *outfile, real tol,
                        atom_type *atoms,int num_atoms, int *orbital_lookup_table,
                        int num_orbs,char transpose,char label_which,int width)
{
  char name_string[80],just_string[80];
  int i,j;
  int beg_col,end_col;
  int beg_row,end_row;
  int which_part,num_parts;
  int cols_per_line;
  real val;

  /* figure out how many columns we get.... */
  cols_per_line = (int)floor((real)width/18);

  /* we need to subtract one from this to leave room for the row labels */
  cols_per_line--;


  beg_col = 0;
  beg_row = 0;

  /* check to see if we're gonna have to break it into multiple parts */
  if( !transpose ){
    if( num_col < cols_per_line ){
      end_col = num_col;
      num_parts = 1;
    }
    else{
      end_col = cols_per_line;
      num_parts = ceil((real)num_col/(real)cols_per_line);
    }
    for(which_part=0;which_part<num_parts;which_part++){
      for( i=-1; i<num_row; i++ ){

        /* start each row with a label */
        if( i >= 0 ){
          if( label_which == LABEL_ROWS || label_which == LABEL_BOTH ){
            map_orb_num_to_name(name_string,i,orbital_lookup_table,
                                num_orbs,atoms,num_atoms);
            left_just_text_string(name_string,just_string,18);
            fprintf(outfile,"%s",just_string);
          }else{
            fprintf(outfile,"%18d",i+1);
          }
        }
        else fprintf(outfile,"               ");

        for( j=beg_col; j<end_col && j< num_col; j++ ){

          /* begin each column with a header */
          if( i == -1){
            if( label_which == LABEL_COLS || label_which == LABEL_BOTH ){
              map_orb_num_to_name(name_string,j,orbital_lookup_table,
                                  num_orbs,atoms,num_atoms);
              center_text_string(name_string,just_string,18);
              fprintf(outfile,"%s",just_string);
            }
            else{
              sprintf(name_string,"%d",j+1);
              center_text_string(name_string,just_string,18);
              fprintf(outfile,"%s",just_string);
            }
          }
          else{
            val = mat[i*num_col+j];
            /* check to see if the value is less than the tolerance */
            if( fabs(val) <= tol ) val = 0.0;
            fprintf(outfile," %-18.4lf",val);
          }
        }
        fprintf(outfile,"\n");
      }

      beg_col=end_col;
      end_col+=cols_per_line;
    }
  } else{
    if( num_row < cols_per_line ){
      end_row = num_row;
      num_parts = 1;
    }
    else{
      end_row = cols_per_line;
      num_parts = ceil((real)num_row/(real)cols_per_line);
    }
    for(which_part=0;which_part<num_parts;which_part++){
      for( i=-1; i<num_col; i++ ){

        /* start each row with a label */

        if( i >= 0 ){
          if( label_which == LABEL_COLS || label_which == LABEL_BOTH ){
            map_orb_num_to_name(name_string,i,orbital_lookup_table,
                                num_orbs,atoms,num_atoms);
            left_just_text_string(name_string,just_string,18);
            fprintf(outfile,"%s",just_string);
          }else{
            fprintf(outfile,"%18d",i+1);
          }
        }
        else fprintf(outfile,"               ");

        for( j=beg_row; j<end_row && j< num_row; j++ ){

          /* begin each column with a header */
          if( i == -1){
            if( label_which == LABEL_ROWS || label_which == LABEL_BOTH ){
              map_orb_num_to_name(name_string,j,orbital_lookup_table,
                                  num_orbs,atoms,num_atoms);
              center_text_string(name_string,just_string,18);
              fprintf(outfile,"%s",just_string);
            }else{
              sprintf(name_string,"%d",j+1);
              center_text_string(name_string,just_string,18);
              fprintf(outfile,"%s",just_string);
            }
          }else{
            val = mat[j*num_col+i];
            /* check to see if the value is less than the tolerance */
            if( fabs(val) <= tol ) val = 0.0;
            fprintf(outfile," %-18.4lf",val);
          }
        }
        fprintf(outfile,"\n");
      }

      beg_row=end_row;
      end_row+=cols_per_line;
    }
  }
  fprintf(outfile,"\n");
  fprintf(outfile,"\n");
}


/****************************************************************************
*
*                   Procedure print_sym_mat
*
* Arguments:     mat : pointer to type real
*     num_col,num_row: integers
*             outfile: pointer to type FILE
*              legend: (optional) pointer to char
*              titles: (optional) pointer to char
*               width: int
* Returns: none
*
* Action: prints out symmetric matrix 'mat to file 'outfile
*         'width columns at a time
*
*
*   if 'legend is non-null it is printed above the array
*
*   if 'titles is non-null, then the elements of that array will be
*    used to label the columns.  NOTE that each label should be no
*    more than 4 characters long.  otherwise the columns will be
*    labelled with numbers;
*
*  see the file notes.outl for the representation of symmetric matrices.
*
*****************************************************************************/
void print_sym_mat( real *mat, int num_row, int num_col, FILE *outfile, char *legend,
                   char *titles,int width)
{
  int i,j;
  int beg_col,end_col;
  int which_part,num_parts;
  int num_so_far;
  int cols_per_line;

  /* figure out how many columns we get.... */
  cols_per_line = (int)floor((real)width/10);

  /* we need to subtract one from this to leave room for the row labels */
  cols_per_line--;


  /* check to see if we're gonna have to break it into multiple parts */
  beg_col = 0;
  if( num_col < cols_per_line ){
    end_col = num_col;
    num_parts = 1;
  }
  else{
    end_col = cols_per_line;
    num_parts = ceil((real)num_col/(real)cols_per_line);
  }

  /* print the legend */
  if( legend ) fprintf(outfile,"%s\n",legend);

  num_so_far = 0;
  for(which_part=0;which_part<num_parts;which_part++){
    fprintf(outfile,"\n");
    for( i=beg_col-1; i<num_row; i++ ){

      /* start each row with a label */
      if( i >= beg_col ){
        if( titles ){
          if( i<9 ){
            fprintf(outfile,"%4s(%4d)  ",&(titles[4*i]),i+1);
          } else if(i<99){
            fprintf(outfile,"%4s(%4d)  ",&(titles[4*i]),i+1);
          } else {
            fprintf(outfile,"%4s(%4d) ",&(titles[4*i]),i+1);
          }
        }
        else fprintf(outfile,"  %-7d",i+1);
      }
      else fprintf(outfile,"            ");

      /* begin each column with a header */
      if( i == beg_col-1){
        for(j=beg_col; j<end_col; j++ ){
          if( titles ){
            if( j<9 ){
              fprintf(outfile,"%4s(%4d)  ",&(titles[4*j]),j+1);
            } else if(j<99){
              fprintf(outfile,"%4s(%4d)  ",&(titles[4*j]),j+1);
            } else {
              fprintf(outfile,"%4s(%4d)  ",&(titles[4*j]),j+1);
            }
          }
          else fprintf(outfile,"  %-7d",j+1);
        }
      }
      else{
        /* go ahead and print out the values for this row */
        for( j=beg_col; j<end_col && j<=i; j++ ){
          fprintf(outfile,"  %-10.4f",mat[i*(i+1)/2 +j]);
        }
      }
      fprintf(outfile,"\n");
    }


    beg_col=end_col;
    end_col+=cols_per_line;
    if( end_col > num_col ) end_col = num_col;

  }
}

/****************************************************************************
 *
 *                   Procedure skipcomments
 *
 * Arguments: file : a pointer to file type
 *        instring : pointer to type char
 *          toggle : a char
 * Returns: an integer
 *
 * Action: Reads in lines from 'file' until one is hit that does not begin
 *     with a ; or a return. puts the first non-comment line into instring
 *     and then returns.
 *
 *    if 'toggle is set to FATAL then hitting EOF will result in
 *      program termination with a call to fatal.
 *    if 'toggle is set to ERROR then EOF results in a call to error then
 *        the function returns.
 *    if 'toggle is set to IGNORE then EOF is ignored.
 *
 *   in any case, if the function returns and EOF has been hit the return
 *    value is -1.
 *
 ****************************************************************************/
int skipcomments(FILE *file,char *string,char toggle)
{
  int i;

  /* use the first element of string to check for EOF */
  string[0] = 0;
  fgets(string,MAX_STR_LEN,file);

  /*******
    deal with the fact that the string may contain only spaces, which
    we will want to skip
  ********/
  i = 0;
  while(string[i] == ' ') i++;
  while( string[i] == '\n' || string[i] == ';'
        && string[i] != 0 ){
    string[0] = 0;
    fgets(string,MAX_STR_LEN,file);
    i = 0;
    while(string[i] == ' ') i++;
  }

  if( string[0] != 0 ) return(0);
  else{
    switch(toggle){
    case FATAL:
      fatal("End of File (EOF) hit in skipcomments.");
      break;
    case ERROR:
      error("End of File (EOF) hit in skipcomments, execution continuing.");
      break;
    case IGNORE:
      break;
    }
    return(-1);
  }
}


/****************************************************************************
*
*                   Procedure find_atoms_orbs
*
* Arguments: num_orbs: int
*           num_atoms: int
*                atom: int
*orbital_lookup_table: pointer to int
*           begin,end: pointers to int
*
* Returns: none
*
* Action: Finds the beginning and end of the orbitals of atom 'atom in
*  'orbital_lookup_table.  This does check for dummy atoms.
*
*****************************************************************************/
void find_atoms_orbs(num_orbs,num_atoms,atom,orbital_lookup_table,begin,end)
  int num_orbs,num_atoms;
  int atom;
  int *orbital_lookup_table;
  int *begin,*end;
{
  char done;
  int next_real_atom;

  /* first check to see if this is a dummy atom */
  if( orbital_lookup_table[atom] < 0 ){
    /* it is, set begin and end, then return */
    *begin = -1;
    *end = -1;
  }else{
    /* okay, the beginning is easy */
    *begin = orbital_lookup_table[atom];

    /**********
      finding the end is a little harder, there are a couple of possible cases
      to deal with.

      if atom is the last one, the end is just num_orbs
      if the next atom is a dummy then we need to look at the one after that.
    ***********/
    next_real_atom = atom + 1;
    done = 0;
    while(next_real_atom < num_atoms && !done ){
      if( orbital_lookup_table[next_real_atom] >= 0 ){
        /* this is the real atom, set end and finish */
        *end = orbital_lookup_table[next_real_atom];
        done = 1;
      }else{
        /* it's a dummy atom, increment and try the next */
        next_real_atom++;
      }
    }
    /*******
      if done isn't set, it's because we ran off the end of the list, so just
      set end to num_orbs.
    *******/
    if( !done ){
      *end = num_orbs;
    }
  }
  /* that's all there is to do! */
}


/****************************************************************************
 *
 *                   Function overlap_tab_from_vect
 *
 * Arguments:  vect: pointer to point_type
 *             cell: pointer to cell_type
 *
 * Returns: int
 *
 * Action:  Determines which R space overlap matrix should be used to
 *   to evaluate the overlap between atoms in the unit cell and the
 *   cell at the end of the lattice vector 'vect.
 *
 *****************************************************************************/

int overlap_tab_from_vect(vect,cell)
  point_type *vect;
  cell_type *cell;
{
  int x,y,z;
  int L,M,N;
  int val;

  /* get the number of overlaps in each direction */
  x = cell->overlaps[0];
  y = cell->overlaps[1];
  z = cell->overlaps[2];

  /* get the components of the lattice vector */
  L = (int)vect->x;
  M = (int)vect->y;
  N = (int)vect->z;

  if(N == 0){
    if( M == 0 ){
      /* case 1 */
      val = ABS(L);
    }
    else{
      /* case 2 */
      if( M < 0 && L != 0 ){
        val = (x + L) * y + ABS(M) + x;
      }
      else{
        val = (x - L) * y + ABS(M) + x;
      }
    }
  }
  else{
    /* case 3 */
    if( N > 0 ){
      val = (y * (2*x + 1) + x + 1) +
        (abs(N) - 1)*(2*y + 1)*(2*x + 1) +
          (x - L)*(2*y + 1) + (y - M);
    }
    if( N < 0 ){
      val = (y * (2*x + 1) + x + 1) +
        (abs(N) - 1)*(2*y + 1)*(2*x + 1) +
          (x + L)*(2*y + 1) + (y + M);
    }
  }
  return(val);
}


/****************************************************************************
*
*                   Procedure check_for_errors
*
* Arguments:   cell: pointer to cell_type
*           details: pointer to detail_type
*          num_orbs: integer
* Returns: none
*
* Action:  This does any internal consistancy checking that needs
*  to be taken care of after the input file has been read in.
*
*  The point of this is to prevent trap as many as possible things
*   which can give rise to core dumps.
*
*  This is not comprehensive, there will always be something that the
*   user can do that will screw things up.
*
*****************************************************************************/
void check_for_errors(cell,details,num_orbs)
  cell_type *cell;
  detail_type *details;
  int num_orbs;
{
  int i,j;
  COOP_type *COOP_ptr,*COOP_ptr2;

  /* did they specify the number of electrons? */
  if( cell->num_electrons < 0.0 ){
    fatal("Negative number of electrons specified.");
  }

  /*******

    if they are just doing a geometry calculation or a band structure then the
    number of electrons is not needed.

  ********/
  if( cell->num_electrons == 0.0 && !details->just_geom ){
    if( details->avg_props || details->Execution_Mode == MOLECULAR ){
       if( cell->num_electrons == 0.0 && details->avg_props ){
         fatal("You forgot to specify the number of electrons.");
       }
     }
  }

  /* mabye there are too many electrons? */
  if( cell->num_electrons > (real)(2.0*num_orbs)){
    fatal("num_electrons is greater than twice the num_orbs.\n\tThis is probably wrong");
  }

  /* what about K points for average properties calculations? */
  if( details->avg_props && details->Execution_Mode != MOLECULAR ){
    if( details->num_KPOINTS == 0 ){
      fatal("No K points given for an extended average properties calculation.");
    }
  }

  /* did they specify a lattice? */
  if( details->Execution_Mode != MOLECULAR ){
    if( cell->dim <= 0 ){
      fatal("You either forgot to specify the lattice parameters or specified them wrong.");
    }
  }

  for(i=0;i<details->num_proj_DOS;i++){
    switch(details->proj_DOS[i].type){
    case P_DOS_ORB:
      for(j=0;j<details->proj_DOS[i].num_contributions; j++){
        if( details->proj_DOS[i].contributions[j] > num_orbs)
          fatal("Orbital projected DOS specified which is larger than num_orbs.");
      }
      break;
    case P_DOS_ATOM:
      for(j=0;j<details->proj_DOS[i].num_contributions; j++){
        if( details->proj_DOS[i].contributions[j] > cell->num_atoms)
          fatal("Atom projected DOS specified which is larger than num_atoms.");
      }
      break;
    case P_DOS_FMO:
      if( !details->num_FMO_frags )
        fatal("FMO projected DOS specified without an FMO specification.");
      for(j=0;j<details->proj_DOS[i].num_contributions; j++){
        if( details->proj_DOS[i].contributions[j] > num_orbs)
          fatal("FMO projected DOS specified which is larger than num_orbs.");
      }
      break;
    }
  }

  COOP_ptr = details->the_COOPS;
  while(COOP_ptr){
    COOP_ptr2 = COOP_ptr;
    while(COOP_ptr2){
      switch(COOP_ptr2->type){
      case P_DOS_ORB:
        if( COOP_ptr2->contrib1 > num_orbs ||
            COOP_ptr2->contrib2 > num_orbs )
          fatal("Orbital COOP contribution larger than num_orbs.");
        break;
      case P_DOS_ATOM:
        if( COOP_ptr2->contrib1 > cell->num_atoms ||
            COOP_ptr2->contrib2 > cell->num_atoms )
          fatal("Atom COOP contribution larger than num_atoms.");
        break;
      case P_DOS_FMO:
        if( !details->num_FMO_frags )
          fatal("FMO COOP specified without an FMO specification.");
        if( COOP_ptr2->contrib1 > num_orbs ||
            COOP_ptr2->contrib2 > num_orbs )
          fatal("FMO COOP contribution larger than num_orbs.");
        break;
      }
      COOP_ptr2 = COOP_ptr2->next_to_avg;
    }
    COOP_ptr = COOP_ptr->next_type;
  }

  /* did they combine goofy stuff with the just_avgE option? */
  if( details->just_avgE ){
    if( details->num_FMO_frags )
      fatal("Can't do FMO analysis with Just Average Energy specified.");
    if( details->num_proj_DOS )
      fatal("Can't do projected DOS with Just Average Energy specified.");
    if( details->the_COOPS )
      fatal("Can't do COOPs with Just Average Energy specified.");
    if( details->chg_mat_PRT || details->OP_mat_PRT ||
       details->ROP_mat_PRT || details->Rchg_mat_PRT ||
       details->wave_fn_PRT || details->net_chg_PRT ||
       details->avg_OP_mat_PRT || details->avg_ROP_mat_PRT ){
      fatal("There's a bad printing option combined with Just Average Energy keyword.");
    }
  }
}



/****************************************************************************
*
*                   Procedure build_orbital_lookup_table
*
* Arguments:   cell: pointer to cell_type
*          num_orbs: pointer to int
*        orbital_lookup_table: pointer to pointer to int
* Returns: none
*
* Action:  Counts the number of orbitals and builds the orbital lookup table.
*
*****************************************************************************/
void build_orbital_lookup_table(cell,num_orbs,orbital_lookup_table)
  cell_type *cell;
  int *num_orbs;
  int **orbital_lookup_table;
{
  int i;
  int num_so_far,atoms_so_far;
  int num_atoms;
  geom_frag_type *geom_frag;

  /*******

    figure out how many atoms there really are (including
    the atoms in geom_frags.

  *******/
  num_atoms = cell->num_raw_atoms;
  geom_frag = cell->geom_frags;
  while(geom_frag){
    num_atoms += geom_frag->num_atoms;
    geom_frag = geom_frag->next;
  }

/*  cell->num_atoms = num_atoms; */

  /* get space for the orbital lookup table */
  *orbital_lookup_table = (int *)calloc(cell->num_atoms,sizeof(int));
  if( !(*orbital_lookup_table) ) fatal("Can't allocate orbital lookup table.");

  /********

    now get the number of orbitals in the unit cell and fill the orbital
    lookup table

  ********/
  num_so_far = 0;
  for(i=0;i<cell->num_raw_atoms;i++){
    /* trap dummy atoms */
    if( cell->atoms[i].at_number < 0 ){
      (*orbital_lookup_table)[i] = -1;
    }
    else{
      (*orbital_lookup_table)[i] = num_so_far;
      if( cell->atoms[i].ns != 0 ) num_so_far++;
      if( cell->atoms[i].np != 0 ) num_so_far += 3;
      if( cell->atoms[i].nd != 0 ) num_so_far += 5;
      if( cell->atoms[i].nf != 0 ) num_so_far += 7;
    }
  }
  /* add the geom_frags to the lookup table */
  geom_frag = cell->geom_frags;
  atoms_so_far = cell->num_raw_atoms;
  while(geom_frag){
    for(i=0;i<geom_frag->num_atoms;i++){
      /* trap dummy atoms */
      if( geom_frag->atoms[i].at_number < 0 ){
        (*orbital_lookup_table)[atoms_so_far] = -1;
      }
      else{
        (*orbital_lookup_table)[atoms_so_far] = num_so_far;
        if( geom_frag->atoms[i].ns != 0 ) num_so_far++;
        if( geom_frag->atoms[i].np != 0 ) num_so_far += 3;
        if( geom_frag->atoms[i].nd != 0 ) num_so_far += 5;
        if( geom_frag->atoms[i].nf != 0 ) num_so_far += 7;
      }
      atoms_so_far++;
    }
    geom_frag = geom_frag->next;
  }


  *num_orbs = num_so_far;
  fprintf(status_file,"There are %d orbitals in the unit cell.\n",num_so_far);
}


/****************************************************************************
*
*                   Procedure print_MOs
*
* Arguments:  details: pointer to detail_type
*          num_orbs: pointer to int
*          eigenset: eigenset_type
*          kpoint: int
*    unique_atoms: pointer to atom_type
* num_unique_atoms: int
*        num_atoms: int
* orbital_lookup_table: pointer to int
*
* Returns: none
*
* Action:  Writes the MO output file.
*
*****************************************************************************/
void print_MOs(details,num_orbs,eigenset,kpoint,unique_atoms,num_unique_atoms,
               num_atoms,orbital_lookup_table)
  detail_type *details;
  int num_orbs;
  eigenset_type eigenset;
  int kpoint;
  atom_type *unique_atoms;
  int num_unique_atoms;
  int num_atoms;
  int *orbital_lookup_table;
{
  static char first_call=1;
  static int x_mirror_present,y_mirror_present,z_mirror_present;
  char filename[240];

  int i,j,k;
  int curr_atom,begin_atom,end_atom;
  int num_walsh_steps;
  sym_op_type *symm_op;
  int ops_passed;
  real this_character;



  /*******

    check to see if this is the first call... if so we need to open
    the output file and write out the header.

  ********/
  if( first_call ){
    /* open the file */
    strcpy(filename,details->filename);
    strcat(filename,".MO");
    MO_file = fopen(filename,"w+");
    if(!MO_file) fatal("Can't open MO file.");

    /* write the atomic parms */
    fprintf(MO_file,"#begin_parms\n");
    for(i=0;i<num_unique_atoms;i++){
      fprintf(MO_file,"%s %d %lf %lf %d %lf %lf %d %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf\n",
              unique_atoms[i].symb,
              unique_atoms[i].ns,unique_atoms[i].exp_s,unique_atoms[i].coul_s,
              unique_atoms[i].np,unique_atoms[i].exp_p,unique_atoms[i].coul_p,
              unique_atoms[i].nd,unique_atoms[i].exp_d,unique_atoms[i].coul_d,
              unique_atoms[i].coeff_d1,unique_atoms[i].exp_d2,unique_atoms[i].coeff_d2,
              unique_atoms[i].nf,unique_atoms[i].exp_f,unique_atoms[i].coul_f,
              unique_atoms[i].coeff_f1,unique_atoms[i].exp_f2,unique_atoms[i].coeff_f2);
    }
    fprintf(MO_file,"#end_parms\n");

    fprintf(MO_file,"#Num_AOs %d\n",num_orbs);


    /* now write out which MO's are being printed */
    if( details->walsh_details.num_steps == 0 ){
      num_walsh_steps = 1;
    } else{
      num_walsh_steps = details->walsh_details.num_steps;
    }

    if( details->Execution_Mode != MOLECULAR ){
      fprintf(MO_file,"#Num_MOs %d\n",
              num_walsh_steps*details->num_MOs_to_print*details->num_KPOINTS);
      for(i=0;i<num_walsh_steps;i++){
        for(j=0;j<details->num_KPOINTS;j++){
          for(k=0;k<details->num_MOs_to_print;k++){
            fprintf(MO_file,"%d %lf %lf %lf %d\n",details->MOs_to_print[k]+1,
                    details->K_POINTS[j].loc.x,details->K_POINTS[j].loc.y,
                    details->K_POINTS[j].loc.z,i+1);
          }
        }
      }
    } else{
      fprintf(MO_file,"#Num_MOs %d\n",details->num_MOs_to_print*num_walsh_steps);
      for(i=0;i<num_walsh_steps;i++){
        for(j=0;j<details->num_MOs_to_print;j++){
          fprintf(MO_file,"%d 0 0 0 %d\n",details->MOs_to_print[j]+1,i+1);
        }
      }
    }

    first_call = 0;
    /***
      if we did symmetry analysis, dump some info about that
      into the MO file.
    ***/
    x_mirror_present = -1;
    y_mirror_present = -1;
    z_mirror_present = -1;
    if( details->use_symmetry && details->Execution_Mode==MOLECULAR){
      symm_op = sym_ops_present;
      ops_passed = 0;
      while(symm_op){
        if(symm_op->type == Mirror ){
          if( symm_op->axis.x == 1.0 ) x_mirror_present = ops_passed;
          else if (symm_op->axis.y == 1.0 ) y_mirror_present = ops_passed;
          else if (symm_op->axis.z == 1.0 ) z_mirror_present = ops_passed;
        }
        symm_op = symm_op->next;
        ops_passed++;
      }
    }
  }


/* round a to nearest int */
#define ROUND(a)        (int)floor((a)+0.5)

  /* okay... we're set, write the MO's that we need to */
  for(i=0;i<details->num_MOs_to_print;i++){
    fprintf(MO_file,"#begin_mo");
    if( x_mirror_present > -1 ){
      this_character = details->characters[x_mirror_present*num_orbs+
                                          details->MOs_to_print[i]];

      if( fabs(1-fabs(this_character)) > 0.001 ){
        this_character = 0;
      }
    } else{
      this_character = 0;
    }
    fprintf(MO_file," %d",ROUND(this_character));
    if( y_mirror_present > -1 ){
      this_character = details->characters[y_mirror_present*num_orbs+
                                          details->MOs_to_print[i]];
      if( fabs(1-fabs(this_character)) > 0.001 ){
        this_character = 0;
      }
    } else{
      this_character = 0;
    }
    fprintf(MO_file," %d",ROUND(this_character));
    if( z_mirror_present > -1 ){
      this_character = details->characters[z_mirror_present*num_orbs+
                                          details->MOs_to_print[i]];
      if( fabs(1-fabs(this_character)) > 0.001 ){
        this_character = 0;
      }
    } else{
      this_character = 0;
    }
    fprintf(MO_file," %d\n",ROUND(this_character));

    for( curr_atom=0; curr_atom<num_atoms; curr_atom++ ){
      find_atoms_orbs(num_orbs,num_atoms,curr_atom,orbital_lookup_table,
                      &(begin_atom),&(end_atom));
      if( begin_atom >= 0 ){
        fprintf(MO_file,"; %d\n",curr_atom+1);
        for(j=begin_atom; j<end_atom; j++){
          fprintf(MO_file,"%6.4lf",EIGENVECT_R(eigenset,details->MOs_to_print[i],j));
          if( details->Execution_Mode != MOLECULAR ){
            fprintf(MO_file," %6.4lf",EIGENVECT_I(eigenset,details->MOs_to_print[i],j));
          }
          fprintf(MO_file,"\n");
        }
      }
    }
    fprintf(MO_file,"#end_mo\n");
  }
}

#ifdef NEED_DSIGN

/****************************************************************************
*
*                   Procedure d_sign
*
* Arguments:  a,b: doubles
*
* Returns: double
*
* Action:  This mimics the fortran 77 function sign.
*
*
*        b\a|  <0  |   >=0  |
*        ---|------|--------|
*         <0|  a   |  -a    |
*        >=0| -a   |   a    |
*
*
*****************************************************************************/
double d_sign(a,b)
  double a,b;
{
  double x;

  if( a >= 0 ) x = a;
  else x = -a;

  if( b >= 0 ) return x;
  else return -x;
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
void parse_integer_string(string,values,num_values)
  char *string;
  int **values,*num_values;
{
  int max_values;
  char local_string[400],num_string[80];
  int i;
  int num1,num2;
  char foo_char;

  *num_values = 0;
  max_values = 10;

  /* if space in *values has already been allocated, blow it out now */
  if( *values ) free(*values);

  /* get some initial memory */
  *values = (int *)calloc(max_values,sizeof(int));
  if( !(*values )) fatal("Can't allocate values in parse_integer_string");

  /* first copy the input string */
  safe_strcpy(local_string,string);

  /* now use strtok to chop it up */
  safe_strcpy(num_string,strtok(local_string,",\n"));

  while(num_string[0]){
    /* check to see if there's a - we need to deal with */
    if( strstr(num_string,"-") ){
      sscanf(num_string,"%d%c%d",&num1,&foo_char,&num2);
      /* deal with negative numbers */
      if(num1<0){
        num2 = num1;
      }
    }else{
      sscanf(num_string,"%d",&num1);
      num2 = num1;
    }
    for(i=num1;i<=num2;i++){
      (*values)[*num_values] = i;
      (*num_values)++;
      if(*num_values == max_values){
        max_values += 10;
        *values = (int *)my_realloc(*values,max_values*sizeof(int));
        if( !(*values )) fatal("Can't reallocate values in parse_integer_string");
      }
    }
    safe_strcpy(num_string,strtok(0,",\n"));
  }
}

/****************************************************************************
 *
 *                   Procedure dump_hermetian_mat
 *
 * Arguments: file: integer
 *             mat: pointer to real
 *        num_orbs: integer
 *
 *
 *
 * Returns: none
 *
 * Action: Does a binary dump of the contents of 'mat into 'file
 *
 ****************************************************************************/
void dump_hermetian_mat(file,mat,num_orbs)
  int file;
  real *mat;
  int num_orbs;
{
  write(file,(const char *)mat,num_orbs*num_orbs*sizeof(real));
}

/****************************************************************************
 *
 *                   Procedure dump_sparse_mat
 *
 * Arguments: file: pointer to FILE
 *             mat: pointer to real
 *        num_orbs: integer
 *         cut_off: real
 *
 *
 *
 * Returns: none
 *
 * Action: Dumps the contents of 'mat (which is assumed to be hermetian)
 *     to 'file in the sparse matrix
 *     format Heinrich Roder likes.
 *
 ****************************************************************************/
void dump_sparse_mat(file,mat,num_orbs,cut_off)
  FILE *file;
  real *mat;
  int num_orbs;
  real cut_off;
{
  int i,j,itab,jtab;
  int num_non_zero,num_written;
  int *nonzero_elements_per_row;
  int *position_in_row;

  nonzero_elements_per_row = (int *)calloc(num_orbs,sizeof(int));
  if( !nonzero_elements_per_row )
    fatal("Can't get memory for nonzero_elements_per_row");
  position_in_row = (int *)malloc(num_orbs*num_orbs*sizeof(int));
  if( !position_in_row )
    fatal("Can't get memory for position_in_row");

  /**********

    count up the nonzero elements (total and number in each row)

  ***********/
  num_non_zero = 0;
  for(i=0;i<num_orbs;i++){
    itab = i * num_orbs;
    for( j=i+1;j<num_orbs;j++){
      if( fabs(mat[itab+j]) > cut_off || fabs(mat[j*num_orbs+i]) > cut_off ){
        num_non_zero+=2;
        nonzero_elements_per_row[i]++;
        nonzero_elements_per_row[j]++;
      }
    }
    if( fabs(mat[itab+i]) > cut_off ){
      num_non_zero++;
      nonzero_elements_per_row[i]++;
    }
  }
  fprintf(stderr,"%d of %d hermetian matrix elements were found to be nonzero\n",
          num_non_zero,num_orbs*num_orbs);

  fprintf(file,"%d\n",num_orbs);
  fprintf(file,"%d\n",num_non_zero);
  /* write the number of non zero elements per row */
  for(i=0;i<num_orbs;i++){
    fprintf(file,"%d\n",nonzero_elements_per_row[i]);
  }

  /* write the nonzero elements */
  num_written = 0;
  for(i=0;i<num_orbs;i++){
    itab = i*num_orbs;
    for(j=0;j<num_orbs;j++){
      jtab = j*num_orbs;
      if( fabs(mat[itab+j]) > cut_off || fabs(mat[jtab+i]) > cut_off){
        if( j>i ){
          fprintf(file,"%lf %lf\n",mat[itab+j],mat[jtab+i]);

        } else{
          if( j != i ){
            fprintf(file,"%lf %lf\n",mat[jtab+i],mat[itab+j]);
          } else{
            fprintf(file,"%lf 0.0\n",mat[jtab+i]);
          }
        }
        position_in_row[num_written++] = j;
      }
    }
  }
  if( num_written != num_non_zero ){
    fprintf(stderr,"num_written (%d) doesn't match num_non_zero (%d)\n",
            num_written,num_non_zero);
  }

  /* now write out the positions of the nonzero elements */
  for(i=0;i<num_non_zero;i++){
    fprintf(file,"%d\n",position_in_row[i]);
  }

  free(position_in_row);
  free(nonzero_elements_per_row);
}





/****************************************************************************
 *
 *                   Procedure charge_to_num_electrons
 *
 * Arguments: cell: pointer to cell_type
 *
 * Returns: none
 *
 * Action: Figures out how many electrons are in the unit cell based
 *   upon the specified charge.
 *
 ****************************************************************************/
void charge_to_num_electrons(cell)
  cell_type *cell;
{
  int i;
  real accum;

  if( cell->geom_frags )
    fatal("Geom Frags and Charge keywords are incompatible (for now).");

  accum = 0.0;
  for(i=0;i<cell->num_atoms;i++){
    accum += (real)(cell->atoms[i].num_valence);
  }
  accum -= cell->charge;
  cell->num_electrons = accum;
}


#ifdef NEED_ETIME
#ifdef UNDERSCORE_FORTRAN
int etime_()
#else
int etime()
#endif
{
  return 0;
}
#endif
