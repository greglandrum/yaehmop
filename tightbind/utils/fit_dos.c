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
  This is program for reading DOS information out of output files
   and generating a broadened set of points suitable for printing.

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
  FILE *infile,*outfile,*the_file;
  char instring[MAX_STR_LEN],com_string[80],type[80],foostring[80];
  char file_name[240];
  point_type *points,*tpoints;
  real *integration;
  int *num_states,*t_states;

  int max_p,num_p,num_so_far;
  real E_min,E_max,E_step,curr_E,E_diff;
  real DOS_here;
  real broadening;
  int num_E_steps;
  int done,num,num_curves,max_curves;
  int points_per_DOS;
  int num_DOS,max_DOS;
  int which_int;
  real fermi;
  int eof_hit;
  real norm_fact;
  int num_walsh_steps;
  int which_walsh_step;
  int temp;

#ifdef USING_THE_MAC
  int argc;
  char argv[4][80];

        /* set up some stuff for Sioux */
        //SIOUXSettings.standalone = FALSE;
        SIOUXSettings.asktosaveonclose = FALSE;
        SIOUXSettings.autocloseonquit = FALSE;
        printf("Starting fit_dos.\n");

  the_file = choose_mac_file(argv[1],MAC_FOPEN_OPEN_CD);
  if( !the_file ) {
          fatal("User cancelled intial file open");
  } else{
          argc = 2;
  }

  /* get the command line arguments */
//  argc = ccommand(&argv);

#endif


  /* open the files */
  if( argc < 2 ){
    fprintf(stderr,"Usage: fit_dos <input_file>\n");
    fatal("Invalid Usage.");
  }

  /*******

    figure out the names of the input and output files

  *******/
  strcpy(instring,argv[1]);
  strcat(instring,".out");
  infile = fopen(instring,"r");
  if(!infile){
    fprintf(stderr,"Can't open input file: %s\n", instring);
    fatal("Can't open file.");
  }


  E_step = ENERGY_STEP;
  broadening = BROADENING;

  /* prompt for the energy window */
  printf("Enter E min: ");
  scanf("%lf",&E_min);
  printf("Enter E max: ");
  scanf("%lf",&E_max);
  printf("Enter broadening: ");
  scanf("%lf",&broadening);
  printf("Enter Energy Step: ");
  scanf("%lf",&E_step);

  /* determine the gaussian normalization factor */
  norm_fact = sqrt(broadening/(real)PI);

  /* get some memory */
  max_p = 10000;
  num_p = 0;
  num_so_far = 0;
  points = (point_type *)calloc(max_p,sizeof(point_type));
  if( !points ) fatal("Can't get space to store the points.");

  max_curves = 10;
  num_curves = 0;
  num_states = (int *)calloc(max_curves,sizeof(int));

  if( !num_states ) fatal("Can't get space for num_states.");


  /*******

    if the file included a walsh diagram, the information for the
     walsh stuff should be at the beginning of the file.

    Look for this information and
    loop until we find a string that begins with # and contains the
    string DENSITY
  ********/
  eof_hit = skipcomments(infile,instring,IGNORE);
  upcase(instring);
  done = 0;
  num_walsh_steps = 0;
  while( eof_hit >= 0 && !done ){
    if( instring[0] == '#' ){
      if(strstr(instring,"WALSH INFO")&&!strstr(instring,"JOB_TITLE")){
        /***
          read out the walsh information
        ****/

        /* the first line has the number of variables */
        skipcomments(infile,instring,FATAL);
        /***
          the next line has the number of steps, which is what we
          are looking for.
        ****/
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%d",&num_walsh_steps);
        done = 1;
      } else if( strstr(instring,"DENSITY") &&!strstr(instring,"JOB_TITLE")){
        done = 1;
      } else{
        eof_hit = skipcomments(infile,instring,IGNORE);
        upcase(instring);
      }
    }else{
      eof_hit = skipcomments(infile,instring,IGNORE);
      upcase(instring);
    }
  }
  if( eof_hit < 0 ){
    fatal("End of file hit before DOS data was found.\n");
  }

  /*****

    if there's walsh data in the file, then prompt the user
    for which step they want to fit.

  ******/
  if( num_walsh_steps ){
    which_walsh_step = 0;
    while(!which_walsh_step){
      printf("There are %d walsh steps in the file.\n",num_walsh_steps);
      printf(" Which would you like to fit? ");
      scanf("%d",&which_walsh_step);
      if( which_walsh_step > num_walsh_steps || which_walsh_step <= 0){
        fprintf(stderr,"%d is a bogus value.\n");
        which_walsh_step = 0;
      }
    }

    /* now skip to the DOS data for this walsh step */
    i = 0;
    while(i != which_walsh_step){
      while(instring[0] != '#' || !strstr(instring,"WALSH_STEP")&&!strstr(instring,"JOB_TITLE")){
        eof_hit = skipcomments(infile,instring,IGNORE);
        if( eof_hit < 0 ){
          fatal("End of file hit before the walsh step was found.");
        }
        upcase(instring);
      }
      eof_hit = skipcomments(infile,instring,FATAL);
      upcase(instring);
      i++;
    }
    /* we're at the right step... move on to the DOS data */
    while(instring[0] != '#' || !strstr(instring,"DENSITY") &&!strstr(instring,"JOB_TITLE")){
      eof_hit = skipcomments(infile,instring,IGNORE);
      if( eof_hit < 0 ){
        fatal("End of file hit before the DOS data was found.");
      }
      upcase(instring);
    }
  }

  /* open the output file */
  if( num_walsh_steps ){
    sprintf(instring,"%s.step%d.DOS",argv[1],which_walsh_step);
  }else{
    strcpy(instring,argv[1]);
    strcat(instring,".DOS");
  }
  outfile = fopen(instring,"w+");
  if(!outfile){
    fprintf(stderr,"Can't open output file: %s\n", instring);
    fatal("Can't open file.");
  }


  num_DOS = 1;

  fprintf(outfile,"#DOS DATA\n");

  /* write out the energy window, step size, and broadening */
  fprintf(outfile,"#E_min: %lf\n",E_min);
  fprintf(outfile,"#E_max: %lf\n",E_max);
  fprintf(outfile,"#E_step: %lf\n",E_step);
  fprintf(outfile,"#Broadening: %lf\n",broadening);

  /* read out the number of states. */
  skipcomments(infile,instring,FATAL);
  sscanf(instring,"%d",&(num_states[0]));
  num_curves++;
  skipcomments(infile,instring,FATAL);

  /******

    okay, we're at the beginning of the DOS data, now read it all in.
     Stop when we hit a line beginning with a #

  *******/
  skipcomments(infile,instring,FATAL);
  while(instring[0] != '#'){
    sscanf(instring,"%lf %lf",&(points[num_p].height),&(points[num_p].energy));
    num_p++;
    num_so_far++;

    /* check to see if we need more memory */
    if( num_p == max_p ){
      max_p += 100;
      tpoints = (point_type *)calloc(max_p,sizeof(point_type));
      if( !tpoints ) fatal("Can't get additional space to store the points.");
      /* free up the memory that we are currently using */
      bcopy(points,tpoints,num_p*sizeof(point_type));
      free(points);
      points = tpoints;
    }
    skipcomments(infile,instring,FATAL);
  }
  points_per_DOS = num_p;


  /*********

    now look for the projected DOS curves

  ***********/
  done = 0;
  while( !done ){
    /* loop until we find a string that begins with # */
    skipcomments(infile,instring,IGNORE);
    while( instring[0] != '#' )  skipcomments(infile,instring,IGNORE);
    /* check to see if this is a projected DOS. */
    if( !strstr(instring,"DENSITY") || !strstr(instring,"PROJECTED") &&!strstr(instring,"JOB_TITLE")){
      done = 1;
    }
    if( !done ){
      num_DOS++;

      skipcomments(infile,instring,FATAL);
      /* read out the number of states. */
      sscanf(instring,"%d",&(num_states[num_curves]));
      skipcomments(infile,instring,FATAL);
      num_curves++;
      /* check to see if we need to get more space for the num_states array */
      if( num_curves == max_curves ){
        t_states = num_states;
        max_curves += 10;
        num_states = (int *)calloc(max_curves,sizeof(int));
        if( !num_states ) fatal("Can't reallocate num_states.");

        /* copy over the old data */
        bcopy((char *)t_states,(char *)num_states,num_curves*sizeof(int));

        /* free up the old memory */
        free(t_states);
      }


      /******
        deal with the possibility of multiple contributions to the projected
        DOS.
      ******/
      while(instring[0] != '#'){
        sscanf(instring,"%s %d",type,&num);
        fprintf(outfile,"# Projection %d  Type: %s Number: %d.\n",
                num_DOS-1,type,num);
        skipcomments(infile,instring,IGNORE);
      }

      /******
        okay, we're at the beginning of the DOS data, now read it all in.
        Stop when we hit a line beginning with a #
      *******/
      skipcomments(infile,instring,FATAL);
      num_so_far = 0;
      while(instring[0] != '#'){
        sscanf(instring,"%lf %lf",&(points[num_p].height),&(points[num_p].energy));
        num_p++;
        num_so_far++;

        /* check to see if we need more memory */
        if( num_p == max_p ){
          max_p += 100;
          tpoints = (point_type *)calloc(max_p,sizeof(point_type));
          if( !tpoints ) fatal("Can't get additional space to store the points.");
          /* free up the memory that we are currently using */
          bcopy(points,tpoints,num_p*sizeof(point_type));
          free(points);
          points = tpoints;
        }
        skipcomments(infile,instring,FATAL);
      }
    }
  }
  /******

    now that we have all the points, generate the curve by smoothing
    the data we have with Gaussians.

  *******/
  num_E_steps = (int)ceil(fabs(E_max-E_min)/E_step)+1;


  /* get space to store the integrations of the curves */
  integration = (real *)calloc(num_DOS*num_E_steps,sizeof(real));
  if( !integration ) fatal("Can't get space for integration data.");

  /* loop over energy values */
  for(i=0, curr_E=E_min;i<num_E_steps;i++, curr_E += E_step ){
    /* loop over DOS points at each energy value */
    num_p = 0;
    for( j=0; j<num_DOS; j++ ){
      DOS_here = 0;

      /* add the value of the integration at the previous energy step */
      if( i != 0 ){
        integration[i*num_DOS+j] = integration[(i-1)*num_DOS+j];
      }

      /* now loop over all the points in this DOS curve */
      for(k=0;k<points_per_DOS;k++,num_p++){
        E_diff = curr_E - points[num_p].energy;
        E_diff *= E_diff;
        if( E_diff <= 25.0 ){
          DOS_here += points[num_p].height*exp(-broadening*E_diff);
        }
      }
      DOS_here *= norm_fact;
      fprintf(outfile,"%lf ",DOS_here);

      /* add on the contribution to the DOS here */
#if 0
      integration[i*num_DOS+j] += DOS_here * E_step / (real)num_states[j];
#endif
      integration[i*num_DOS+j] += 2.0*DOS_here * E_step / (real)num_states[j];
    }
    fprintf(outfile,"%lf\n",curr_E);
  }

  /******

    now print out the integration data

  ******/
  fprintf(outfile,"# END OF DOS\n");
  fprintf(outfile,"# BEGIN INTEGRATION\n");
  for(i=0, curr_E=E_min;i<num_E_steps;i++, curr_E += E_step ){
    for( j=0; j<num_DOS; j++ ){
      fprintf(outfile,"%lf ",integration[i*num_DOS+j]);
    }
    fprintf(outfile,"%lf\n",curr_E);
  }

  /****
    find and read out the fermi energy
  *****/
  eof_hit = skipcomments(infile,instring,IGNORE);
  upcase(instring);

  while(eof_hit >= 0 && !strstr(instring,"FERMI_ENERGY")&&!strstr(instring,"JOB_TITLE")){
    eof_hit = skipcomments(infile,instring,IGNORE);
    upcase(instring);
  }

  if( eof_hit >= 0 ){
    sscanf(instring,"%s %lf",foostring,&fermi);
  }
  else{
    fermi = 0;
  }
  fprintf(outfile,"\n#FERMI_ENERGY: %lf\n",fermi);

  printf("Done!\n\n");
}




