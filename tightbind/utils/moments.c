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
   and generating a series of moments

   This currently does very little error checking...

   Created by greg Landrum March 1996
************************************************************************/
#include "fit_props.h"

void main(argc,argv)
  int argc;
  char **argv;
{
  int i,j,k;
  FILE *infile,*outfile;
  char instring[MAX_STR_LEN],com_string[80],type[80],foostring[80];
  point_type *points,*tpoints;
  real *moments;
  int *num_states,*t_states;

  int max_p,num_p,num_so_far;
  int jtab;
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
  int num_moments;
  int which_dos;
  real E_min,E_max,E_scale,scaled_E,E_accum;
  real a_val,b_val;

  /* open the files */
  if( argc < 3 ){
    fprintf(stderr,"Usage: moments <input_file> num_moments\n");
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

  sscanf(argv[2],"%d",&num_moments);


  /* get some memory */
  max_p = 100;
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
      if(strstr(instring,"WALSH INFO")){
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
      } else if( strstr(instring,"DENSITY") ){
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
      while(instring[0] != '#' || !strstr(instring,"WALSH_STEP")){
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
    while(instring[0] != '#' || !strstr(instring,"DENSITY") ){
      eof_hit = skipcomments(infile,instring,IGNORE);
      if( eof_hit < 0 ){
        fatal("End of file hit before the DOS data was found.");
      }
      upcase(instring);
    }
  }

  /* open the output file */
  if( num_walsh_steps ){
    sprintf(instring,"%s.step%d.MOM",argv[1],which_walsh_step);
  }else{
    strcpy(instring,argv[1]);
    strcat(instring,".MOM");
  }
  outfile = fopen(instring,"w+");
  if(!outfile){
    fprintf(stderr,"Can't open output file: %s\n", instring);
    fatal("Can't open file.");
  }


  num_DOS = 1;

  E_min = 1e8;
  E_max = -1e8;

  fprintf(outfile,"#MOMENTS DATA\n");
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

    if( points[num_p].energy < E_min ) E_min = points[num_p].energy;
    if( points[num_p].energy > E_max ) E_max = points[num_p].energy;

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


  /* figure out the energy scaling */
  E_scale = 2.0/(E_max-E_min);

  fprintf(outfile,"#NUM_MOMENTS: %d\n",num_moments);
  /* write out the energy window, and scaling term */
  fprintf(outfile,"#E_min: %lf\n",E_min);
  fprintf(outfile,"#E_max: %lf\n",E_max);
  fprintf(outfile,"#E_scale: %lf\n",E_scale);

  /*********

    now look for the projected DOS curves

  ***********/
  done = 0;
  while( !done ){
    /* loop until we find a string that begins with # */
    skipcomments(infile,instring,IGNORE);
    while( instring[0] != '#' )  skipcomments(infile,instring,IGNORE);
    /* check to see if this is a projected DOS. */
    if( !strstr(instring,"DENSITY") || !strstr(instring,"PROJECTED") ){
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


  fprintf(outfile,"#NUM_DOSs %d\n",num_DOS);

  /* get space to store the individual moments */
  moments = (real *)calloc(num_DOS,sizeof(real));
  if( !moments ) fatal("Can't get space for moments data.");

  a_val = (E_max - E_min)/2.0;
  b_val = (E_max + E_min)/2.0;


  /* loop over energy values */
  for(i=0; i<num_moments;i++){
    fprintf(outfile,"%d ",i);
    printf("%d\n",i);
    bzero(moments,num_DOS*sizeof(real));
    num_p = 0;
    /* loop over DOS points at each energy value */
    for(k=0;k<points_per_DOS;k++){

      /* figure out the scaled version of E^n */
      scaled_E = (points[k].energy - b_val)/a_val;

/*fprintf(stderr,"%d %d %lf\n",i,k,scaled_E);*/

      E_accum = 1;
      for(j=1;j<=i;j++) E_accum *= scaled_E;
      scaled_E = E_accum;



      /* loop over the DOS values */
      for( j=0; j<num_DOS; j++ ){
        jtab = j*points_per_DOS;
        moments[j] += scaled_E*points[jtab+k].height;
      }
    }

    /* now write the values of this moment */
    for( j=0;j<num_DOS;j++){
      fprintf(outfile,"%lf ",moments[j]);
    }
    fprintf(outfile,"\n");
  }

  fprintf(outfile,"#END MOMENTS\n");
  fclose(outfile);

  if( argc > 3 ){
    for(k=0;k<points_per_DOS;k++){

      /* figure out the scaled version of E^n */
      scaled_E = (points[k].energy - b_val)/a_val;
      fprintf(stderr,"%lf\n",scaled_E);
    }
  }


}





