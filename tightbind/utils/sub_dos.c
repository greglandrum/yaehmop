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
/************************************************************************
  This is a program for subtracting 2 DOS curves from each other and
   writing the results to an output file.

   This currently does very little error checking...

   Created by greg Landrum March 1994
************************************************************************/
#include "fit_props.h"

#ifndef USING_THE_MAC
void main(argc, argv)
  int argc;
  char **argv;
#else
void main()
#endif
{
  int i,j,k;
  FILE *infile1,*infile2,*outfile,*the_file;
  char instring1[MAX_STR_LEN],com_string[80],foostring1[80],foostring2[80];
  char instring2[MAX_STR_LEN];
  real E_min,E_max,E_step,curr_E,E_diff;
  real E_min2,E_max2,E_step2;
  real broad1,broad2;
  int num_E_steps;
  int curve1,curve2;
  real val1, val2;
  char window_by_hand;

#ifdef USING_THE_MAC
  int argc;
  char argv[4][80];

        /* set up some stuff for Sioux */
        //SIOUXSettings.standalone = FALSE;
        SIOUXSettings.asktosaveonclose = FALSE;
        SIOUXSettings.autocloseonquit = FALSE;
        printf("Starting bind.\n");

  the_file = choose_mac_file(argv[1],MAC_FOPEN_OPEN_CD);
  if( !the_file ) {
          fatal("User cancelled intial file open");
  }
  the_file = choose_mac_file(argv[2],MAC_FOPEN_OPEN_CD);
  if( !the_file ) {
          fatal("User cancelled second file open");
  }

  printf("Enter the name of the output file: ");
  scanf("%s",argv[3]);

  argc = 4;
  /* get the command line arguments */
//  argc = ccommand(&argv);

#endif

  /* open the files */
  if( argc < 4 ){
    fprintf(stderr,"Usage: sub_dos input_file1 input_file2 output_file\n");
    fatal("Invalid Usage.");
  }

  /*******

    open the input and output files.

  *******/
  strcpy(instring1,argv[1]);
  infile1 = fopen(instring1,"r");
  if(!infile1){
    fprintf(stderr,"Can't open input file: %s\n", instring1);
    fatal("Can't open file.");
  }

  strcpy(instring1,argv[2]);
  infile2 = fopen(instring1,"r");
  if(!infile2){
    fprintf(stderr,"Can't open input file: %s\n", instring1);
    fatal("Can't open file.");
  }

  strcpy(instring1,argv[3]);
  outfile = fopen(instring1,"w+");
  if(!outfile){
    fprintf(stderr,"Can't open output file: %s\n", instring1);
    fatal("Can't open file.");
  }

  /* get the numbers of the curves to subtract */
  printf("NOTE:  curve 0 is the total DOS, curve 1 is the first projection.\n");
  printf("Enter the number of the curve to use from file 1: ");
  scanf("%d",&curve1);
  printf("Enter the number of the curve to use from file 2: ");
  scanf("%d",&curve2);

  /******

    read out the energy windows

    we're assuming that the energy window starts on the
    second line of each file.

  ******/
  printf("Reading out the energy window.\n");
  skipcomments(infile1,instring1);
  skipcomments(infile1,instring1);
  skipcomments(infile2,instring2);
  skipcomments(infile2,instring2);

  sscanf(instring1,"%s %lf",foostring1,&E_min);
  sscanf(instring2,"%s %lf",foostring1,&E_min2);
  skipcomments(infile1,instring1);
  skipcomments(infile2,instring2);
  sscanf(instring1,"%s %lf",foostring1,&E_max);
  sscanf(instring2,"%s %lf",foostring1,&E_max2);
  skipcomments(infile1,instring1);
  skipcomments(infile2,instring2);
  sscanf(instring1,"%s %lf",foostring1,&E_step);
  sscanf(instring2,"%s %lf",foostring1,&E_step2);
  skipcomments(infile1,instring1);
  skipcomments(infile2,instring2);
  sscanf(instring1,"%s %lf",foostring1,&broad1);
  sscanf(instring2,"%s %lf",foostring1,&broad2);

  /* make sure that the energy windows are the same */
  window_by_hand = 0;
  if( E_min != E_min2 ){
    fprintf(stderr,"Energy min's (%lf and %lf) aren't equal.\n",E_min,E_min2);
    window_by_hand = 1;
  }
  if( E_min != E_min2 ){
    fprintf(stderr,"Energy max's (%lf and %lf) aren't equal.\n",E_max,E_max2);
    window_by_hand = 1;
  }
  if( E_step != E_step2 ){
    fprintf(stderr,"Energy steps (%lf and %lf) aren't equal.\n",E_step,
            E_step2);
    window_by_hand = 1;
  }
  if( broad1 != broad2 ){
    fprintf(stderr,"The broadening values (%lf and %lf) aren't equal.\n",
            broad1,broad2);
    fprintf(stderr,"It's okay to continue, but this is a questionably valid procedure.\n");
  }



  if( window_by_hand ){
    fprintf(stderr,
            "Enter the window by hand.  This is a dubious thing to do.\n");
    printf("Enter E min: ");
    scanf("%lf",&E_min);
    printf("Enter E max: ");
    scanf("%lf",&E_max);
    printf("Enter Energy Step: ");
    scanf("%lf",&E_step);
  }

  printf("E min = %lf\n",E_min);
  printf("E max = %lf\n",E_max);
  printf("E step = %lf\n",E_step);

  num_E_steps = ceil(fabs(E_max-E_min)/E_step)+1;

  /* now read in the file */
  skipcomments(infile1,instring1);
  while(instring1[0] == '#') skipcomments(infile1,instring1);
  skipcomments(infile2,instring2);
  while(instring2[0] == '#') skipcomments(infile2,instring2);


  fprintf(outfile,"# Subtracted Density of States data\n");
  fprintf(outfile,"#E_min: %lf\n",E_min);
  fprintf(outfile,"#E_max: %lf\n",E_max);
  fprintf(outfile,"#E_step: %lf\n",E_step);
  fprintf(outfile,"#Broadening: %lf\n",broad1);


  /* read in lines until we hit the end of the dos data */
  curr_E = E_min;
  while(instring1[0] != '#' ){
    /* find the appropriate numbers */
    strcpy(foostring1,strtok(instring1," "));
    for(j=0;j<curve1;j++)
      strcpy(foostring1,strtok(0," "));
    sscanf(foostring1,"%lf",&val1);

    strcpy(foostring2,strtok(instring2," "));
    for(j=0;j<curve2;j++)
      strcpy(foostring2,strtok(0," "));
    sscanf(foostring2,"%lf",&val2);
/*
printf("%lf %lf %lf\n",val1,val2,curr_E);
*/
    fprintf(outfile,"%lf %lf\n",val1-val2,curr_E);
    skipcomments(infile1,instring1);
    skipcomments(infile2,instring2);
    curr_E += E_step;
  }

  /* now find the integration data */
  while(!strstr(instring1,"INTEGRATION")) skipcomments(infile1,instring1);
  while(!strstr(instring2,"INTEGRATION")) skipcomments(infile2,instring2);

  /* read in lines until we hit the end of the dos data */
  fprintf(outfile,"# END OF DOS\n");
  fprintf(outfile,"# BEGIN INTEGRATION\n");
  curr_E = E_min;
  skipcomments(infile1,instring1);
  skipcomments(infile2,instring2);
  while(instring1[0] != '#' && instring1[0] != '\n'){

    strcpy(foostring1,strtok(instring1," "));
    for(j=0;j<curve1;j++)
      strcpy(foostring1,strtok(0," "));
    sscanf(foostring1,"%lf",&val1);

    strcpy(foostring2,strtok(instring2," "));
    for(j=0;j<curve2;j++)
      strcpy(foostring2,strtok(0," "));
    sscanf(foostring2,"%lf",&val2);

    sscanf(foostring1,"%lf",&val1);
    sscanf(foostring2,"%lf",&val2);

    fprintf(outfile,"%lf %lf\n",val1-val2,curr_E);
    skipcomments(infile1,instring1);
    skipcomments(infile2,instring2);
    curr_E += E_step;
  }
}




