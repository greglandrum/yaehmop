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
*     this file contains general purpose utility functions
*
*  created:  greg landrum  August 1993
*
*****************************************************************************/
#include <stdio.h>

#define FATAL 0
#define ERROR 1
#define IGNORE 2
#define MAX_STR_LEN 2048

/* Procedure fatal
 * prints an error message and terminates the program
 *
 * in case the job was queued, the error message is echoed to the status file
 */
void fatal( errorstring )
  char *errorstring;
{
  fprintf( stderr, "FATAL ERROR: %s.\nExecution Terminated.\n",
      errorstring );
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
  fprintf( stderr, "ERROR: %s.\n", errorstring );
  return;
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
int skipcomments(file,string,toggle)
  FILE *file;
  char *string;
  char toggle;
{

  /* use the first element of string to check for EOF */
  string[0] = 0;
  fgets(string,MAX_STR_LEN,file);
  while( string[0] == '\n' || string[0] == ';'
        && string[0] != 0 ){
    string[0] = 0;
    fgets(string,MAX_STR_LEN,file);
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
