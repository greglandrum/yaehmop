/************************************************************************
  This is program for reading DOS information out of output files
   and generating a projected energy...

   This currently does very little error checking...

   Created by greg Landrum February, 1996
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
  int *num_states,*t_states;
  point_type *projection,*total;
  int max_p,num_p,num_so_far;
  int done,num,num_curves,max_curves;
  int points_per_DOS;
  int num_DOS,max_DOS;
  real num_electrons,electrons_left;
  real fermi,perch_accum, avgE_accum;
  double height_accum;
  int eof_hit;
  real norm_fact;
  int num_walsh_steps;
  int which_walsh_step;
  int temp;

  /* open the files */
  if( argc < 2 ){
    fprintf(stderr,"Usage: perch <input_file>\n");
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

  num_DOS = 1;

  /* read out the number of states. */
  skipcomments(infile,instring,FATAL);
  sscanf(instring,"%d",&(num_states[0]));
  num_curves++;
  skipcomments(infile,instring,FATAL);

  /******

    okay, we're at the beginning of the total DOS data, now read
    it all in.
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

  /****
    find and read out the fermi energy
  *****/
  eof_hit = skipcomments(infile,instring,IGNORE);
  upcase(instring);

  while(eof_hit >= 0 && !strstr(instring,"FERMI_ENERGY")){
    eof_hit = skipcomments(infile,instring,IGNORE);
    upcase(instring);
  }

  if( eof_hit >= 0 ){
    sscanf(instring,"%s %lf",foostring,&fermi);
  }
  else{
    fprintf(stderr,"Can't find the Fermi E in input file.\n");
    fprintf(stderr,"Please enter a value: ");
    scanf("%lf",&fermi);
  }

  fprintf(stderr,"The Fermi Energy is: %lf\n",fermi);



  /* some error checking */
  for(i=0;i<points_per_DOS;i++){
    height_accum = 0.0;
    for(j=1;j<num_curves;j++){
      height_accum += points[j*points_per_DOS+i].height;
    }
    if( height_accum != points[i].height ){
      fprintf(stderr,"Point: %d (E=%lf) doesn't add up (diff = %lf)\n",
              i,points[i].energy,points[i].height-height_accum);
    }
  }




  /**********

    we've got all the projections, loop over each one
    up to the Fermi Energy and add up the
    energy contributed by that projection.

  **********/



  total = points;
  for(i=1;i<num_curves;i++){
    j = 0;
    projection = &(points[i*points_per_DOS]);
    perch_accum = 0.0;
    avgE_accum = 0.0;
    while(total[j].energy <= fermi){
      perch_accum += 2.0*projection[j].height*total[j].energy;
      avgE_accum += 2.0*total[j].height*total[j].energy;

      j++;
    }

    printf("Curve: %d \t PERCH: %lf  \tAvg_E: %lf\n",i,perch_accum,
           avgE_accum);
  }
}





