/************************************************************************
  This is program for reading COOP information out of output files
   and generating a broadened set of points suitable for printing.

   This currently does very little error checking...

   Created by greg Landrum June 1994
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
  point_type *points,*tpoints;
  real *integration;
  int max_p,num_p,num_so_far;
  real E_min,E_max,E_step,curr_E,E_diff;
  real COOP_here;
  real broadening;
  int num_E_steps;
  int done,num;
  int points_per_COOP;
  int num_COOP,max_COOP;
  int which_int;
  int eof_hit;
  real fermi;
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
        printf("Starting fit_coop.\n");

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
    fprintf(stderr,"Usage: fit_coop <input_file>\n");
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
  max_p = 100;
  num_p = 0;
  num_so_far = 0;
  points = (point_type *)calloc(max_p,sizeof(point_type));
  if( !points ) fatal("Can't get space to store the points.");



  /*******

    if the file included a walsh diagram, the information for the
     walsh stuff should be at the beginning of the file.

    Look for this information and
    loop until we find a string that begins with # and contains the
    string COOP
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
      } else if( strstr(instring,"COOP")&&!strstr(instring,"JOB_TITLE")){
        done = 1;
      }
      else{
        eof_hit = skipcomments(infile,instring,IGNORE);
        upcase(instring);
      }

    }
    else{
      eof_hit = skipcomments(infile,instring,IGNORE);
      upcase(instring);
    }
  }
  if( eof_hit < 0 ){
    fatal("End of file hit before initial COOP data was found.\n");
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

    /* now skip to the COOP data for this walsh step */
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
    /* we're at the right step... move on to the COOP data */
    while(instring[0] != '#' || !strstr(instring,"COOP") &&!strstr(instring,"JOB_TITLE")){
      eof_hit = skipcomments(infile,instring,IGNORE);
      if( eof_hit < 0 ){
        fatal("End of file hit before the COOP data was found.");
      }
      upcase(instring);
    }
  }

  /* open the output file */
  if( num_walsh_steps ){
    sprintf(instring,"%s.step%d.COOP",argv[1],which_walsh_step);
  }else{
    strcpy(instring,argv[1]);
    strcat(instring,".COOP");
  }
  outfile = fopen(instring,"w+");
  if(!outfile){
    fprintf(stderr,"Can't open output file: %s\n", instring);
    fatal("Can't open file.");
  }

  fprintf(outfile,"#COOP DATA\n");

  /* write out the energy window, step size, and broadening */
  fprintf(outfile,"#E_min: %lf\n",E_min);
  fprintf(outfile,"#E_max: %lf\n",E_max);
  fprintf(outfile,"#E_step: %lf\n",E_step);
  fprintf(outfile,"#Broadening: %lf\n",broadening);


  /* read out the number of curves */
  skipcomments(infile,instring,FATAL);
  sscanf(instring,"%d",&num_COOP);



  /******

    okay, we're at the beginning of the projected COOP data,
    read in each of the curves

  *******/
  skipcomments(infile,instring,FATAL);
  num_p = 0;
  points_per_COOP = 0;
  for(i=0;i<num_COOP;i++){

    /* move to the beginning of the data */
    while(!(instring[0] == '#' && strstr(instring,"BEGIN"))&&!strstr(instring,"JOB_TITLE")){
      skipcomments(infile,instring,FATAL);
    }

    /* now read in the data until we a line containing END */
    skipcomments(infile,instring,FATAL);

    while(!(instring[0] == '#' && strstr(instring,"END"))&&!strstr(instring,"JOB_TITLE")){
      sscanf(instring,"%lf %lf",&(points[num_p].height),&(points[num_p].energy));
      if( !i )points_per_COOP++;
      num_p++;

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

  /******

    now that we have all the points, generate the curve by smoothing
    the data we have with Gaussians.

  *******/
  num_E_steps = ceil(fabs(E_max-E_min)/E_step) + 1;


  /* get space to store the integrations of the curves */
  integration = (real *)calloc(num_COOP*num_E_steps,sizeof(real));
  if( !integration ) fatal("Can't get space for integration data.");

  /* loop over energy values */
  for(i=0, curr_E=E_min;i<num_E_steps;i++, curr_E += E_step ){
    /* loop over COOP points at each energy value */
    num_p = 0;
    for( j=0; j<num_COOP; j++ ){
      COOP_here = 0;

      /* add the value of the integration at the previous energy step */
      if( i != 0 ){
        integration[i*num_COOP+j] = integration[(i-1)*num_COOP+j];
      }

      /* now loop over all the points in this COOP curve */
      for(k=0;k<points_per_COOP;k++,num_p++){
        E_diff = curr_E - points[num_p].energy;
        E_diff *= E_diff;
        if( E_diff <= 25.0 ){
          COOP_here += points[num_p].height*exp(-broadening*E_diff);
        }
      }
      COOP_here *= norm_fact;
      fprintf(outfile,"% -6.4lf ",COOP_here);

      /* add on the contribution to the COOP here */
      integration[i*num_COOP+j] += COOP_here * E_step;
    }
    fprintf(outfile,"% -6.4lf\n",curr_E);
  }

  /******

    now print out the integration data

  ******/
  fprintf(outfile,"# END OF COOP\n");
  fprintf(outfile,"# BEGIN INTEGRATION\n");
  for(i=0, curr_E=E_min;i<num_E_steps;i++, curr_E += E_step ){
    for( j=0; j<num_COOP; j++ ){
      fprintf(outfile,"% -6.4lf ",integration[i*num_COOP+j]);
    }
    fprintf(outfile,"% -6.4lf\n",curr_E);
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




