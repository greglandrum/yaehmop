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
*     this is the main procedure for the program bind
*
*  read the file README.bind to find out more
*
*  created:  greg landrum  August 1993
*
*****************************************************************************/
/***
  Edit History:

  March '98: WG
    - print FMO-FMO COOP's to .out

  07.05.98 gL:
    added netCDF support.
   14.02.2004 gL:
     added a "10th anniversary" message :-)

***/
#include "bind.h"

// For booleans
#include "stdbool.h"

FILE *COHP_file;
const char greetings[]="Welcome to the 10th Anniversary edition of YAeHMOP!\n";

void main(argc, argv)
  int argc;
  char **argv;
{
  FILE *temp_file;
  char file_name[500];
  FILE *the_file=0;
  bool use_stdin_stdout = false;

  if( argc == 2 && strcmp(argv[1], "-v") == 0){
    fprintf(stdout, "version: %s\n", VERSION_STRING);
    exit(0);
  }

  /* make sure the program was called with the right arguments */
  if( argc < 2){
    fprintf(stderr,"Usage: bind <inputfile> [paramfile]\n");
    exit(-1);
  }

  // If the argument is --use_stdin_stdout, use stdin and stdout! Makes sense...
  if (strcmp(argv[1], "--use_stdin_stdout") == 0)
    use_stdin_stdout = true;
  if(!use_stdin_stdout)
    strncpy(file_name,argv[1],500);
  else
    strcpy(file_name,"stdin");

  // If there is a third argument, it is the parameter file
  char* parm_file_name = NULL;
  if (argc > 2)
    parm_file_name = argv[2];

  /* install the sig_int handler */
  signal(SIGINT,handle_sigint);

  fprintf(stderr,greetings);

  run_bind(file_name, use_stdin_stdout, parm_file_name);

  exit(0);
}
