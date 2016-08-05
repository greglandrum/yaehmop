/************************************************************************
  This is program for reading molecular orbital energy levels and
    symmetries out of output files and generating an output file
    which can be used to generate a walsh diagram.

   This currently does very little error checking...

   Created by greg Landrum June 1994
************************************************************************/
#include "fit_props.h"
#include "fit_walsh.h"



/****************************************************************************
*
*                   Function compare_characters
*
* Arguments: symm1,symm2: pointers to reals
*               num_symm: int
*
* Returns: real
*
* Action: This is used when characters cannot be matched exactly between
*  steps, it just compares the signs of all the characters.
*
*****************************************************************************/
real compare_characters(symm1,symm2,num_symm)
  real *symm1,*symm2;
  int num_symm;
{
  int i;
  real dev;

  dev = 0.0;
  for(i=0;i<num_symm; i++){
    /* assign a large penalty to characters which change sign */
    if( fabs(symm1[i]) > 1e-3 && fabs(symm2[i]) > 1e-3 ){
      if( (symm1[i] > 0 && symm2[i] < 0) ||
          (symm1[i] < 0 && symm2[i] > 0) ){
        dev += 100;
      }
    }
  }
  return(dev);
}


/****************************************************************************
*
*                   Function compare_symms
*
* Arguments: symm1,symm2: pointers to reals
*               num_symm: int
*
* Returns: real
*
* Action: Returns the absolute value of the
*   deviation between the two arrays 'symm1 and 'symm2
*
*****************************************************************************/
real compare_symms(symm1,symm2,num_symm)
  real *symm1,*symm2;
  int num_symm;
{
  int i;
  real dev;

  dev = 0.0;

  for(i=0;i<num_symm; i++){
    dev += (symm1[i]-symm2[i])*(symm1[i]-symm2[i]);
#if 0
    /* assign a large penalty to characters which change sign */
    if( fabs(symm1[i]) > 1e-3 && fabs(symm2[i]) > 1e-3 ){
      if( (symm1[i] > 0 && symm2[i] < 0) ||
          (symm1[i] < 0 && symm2[i] > 0) ){
        dev += 100;
      }
    }
#endif
  }
  return(sqrt(dev));
}


/****************************************************************************
*
*                   Procedure construct_lines
*
* Arguments: points:  pointer to walsh_point_type
* num_orbs,num_symm,num_steps: integers
*
* Returns: none
*
* Action: This matches up energy levels with the same symmetry from
*  one step to the next.   As long as no degeneracies are present, this
*  will always work.  However, in the case of degeneracies, things won't
*  be so clear-cut.  The way this is dealt with is by connecting degenerate
*  points with whichever point is closest.
*
*****************************************************************************/
void construct_lines(points,num_orbs,num_symm,num_steps)
  walsh_point_type *points;
  int num_orbs,num_symm,num_steps;
{
  int i,j,k;
  int *matched,*found;
  char foundit;
  real closest,dist,ediff;
  int which_is_closest;

  /* get space to use to keep track of which orbitals have been matched. */
  matched = (int *)calloc(num_orbs,sizeof(int));
  if(!matched) fatal("Can't get space for matched array.");
  found = (int *)calloc(num_orbs,sizeof(int));
  if(!found) fatal("Can't get space for matched array.");

  for(i=0;i<num_steps-1;i++){

    /* zero out the matched array */
    bzero(matched,num_orbs*sizeof(int));
    bzero(found,num_orbs*sizeof(int));

#if 1
    /* now loop over the orbitals here and match them */
    for(j=0;j<num_orbs;j++){
      foundit = 0;
      for(k=0;k<num_orbs && !foundit;k++){
        /********
          if this orbital has already been matched, we don't have to
          do anything
          *********/
        if( !matched[k] ){
          /* it hasn't been, see if it matches the orbital we are trying to match */
          if( compare_symms(points[i*num_orbs+j].symmetries,
                            points[(i+1)*num_orbs+k].symmetries,num_symm) < 1e-6 ){
            foundit = 1;
            points[i*num_orbs+j].next = &(points[(i+1)*num_orbs+k]);
            found[j] = 1;
fprintf(stderr,"%d: %d->%d\n",i,j,k);
            matched[k] = 1;
          }
        }
      }
    }
#endif
    /*******

      now loop through again and match up all the orbitals that didn't have
      exact matches in the previous loop.

      ********/
    for( j=0;j<num_orbs;j++){
      /******
        if we didn't find it, then point this orbital at the closest
        non-matched orbital at the next step
        ******/
      if( !found[j] ){
        printf("No exact: %d, step: %d\n",j,i);
        closest = 1e6;
        which_is_closest = 1000;
        for(k=0;k<num_orbs;k++){
          if( !matched[k] ){
            dist =  compare_symms(points[i*num_orbs+j].symmetries,
                                  points[(i+1)*num_orbs+k].symmetries,num_symm);
            /* weight the distance in symmetry space by the energy difference */
            ediff = fabs(points[i*num_orbs+j].energy - points[(i+1)*num_orbs+k].energy);

            dist = dist*(1+ediff)*(1+ediff);
            if( dist < closest ){
fprintf(stderr,"%d: %d->%d\n",i,j,k);
              closest = dist;
              which_is_closest = k;
            }
          }
        }
        if( which_is_closest >= num_orbs ){
          fprintf(stderr,"Can't match orbital %d in step %d\n",j,i);
          fatal("Can't match an orbital.");
        }
        points[i*num_orbs+j].next = &(points[(i+1)*num_orbs+which_is_closest]);
        found[j] = 1;
        matched[which_is_closest] = 1;
      }
    }
  }
}


/****************************************************************************
*
*                   Procedure extract_walsh_data
*
* Arguments:  infile: pointer to type FILE
*           p_points:  pointer to pointer to walsh_point_type
*    p_xvals,p_tot_E:  pointers to pointers to real
* p_num_orbs,p_num_symm,p_num_steps: pointers to int
*
* Returns: none
*
* Action: reads all the needed information from 'infile
*
*****************************************************************************/
void extract_walsh_data(infile,p_points,p_xvals,p_tot_E,p_num_orbs,p_num_symm,
                        p_num_steps)
  FILE *infile;
  walsh_point_type **p_points;
  real **p_xvals,**p_tot_E;
  int *p_num_orbs,*p_num_symm,*p_num_steps;
{
  char instring[240],com_string[80];
  int i,j,k;
  int which_var,num_vars,num_p;
  int num_orbs,num_symm,num_steps;
  real *xvals,*tot_E;
  walsh_point_type *points;

  /*****
    first read until we find out the number of orbitals in the
    unit cell.
  ******/
  skipcomments(infile,instring,FATAL);
  upcase(instring);
  while(instring[0] != '#' || !strstr(instring,"NUM_ORBITALS")){
    skipcomments(infile,instring,FATAL);
    upcase(instring);
  }
  sscanf(instring,"%s %d",com_string,&num_orbs);

  rewind(infile);

  /*****

    okay, now find out how many Walsh variables and steps there are.

  ******/
  skipcomments(infile,instring,FATAL);
  upcase(instring);
  while(instring[0] != '#' || !strstr(instring,"WALSH") || strstr(instring,"JOB_TITLE")){
    skipcomments(infile,instring,FATAL);
    upcase(instring);
  }
  skipcomments(infile,instring,FATAL);
  sscanf(instring,"%d",&num_vars);
  skipcomments(infile,instring,FATAL);
  sscanf(instring,"%d",&num_steps);

  /* now prompt to see which of the variables to use to index the plot. */
  printf("There are %d Walsh variables in this file.\n",num_vars);
  which_var = -1;
  while(which_var < 1 || which_var > num_vars){
    printf("Which should be used as the x axis in the plot? ");
    scanf("%d",&which_var);
    if( which_var < 1 || which_var > num_vars )
      printf("That was a bogus value, try again.\n");
  }
  which_var -= 1;

  /* find out how many symmetry elements need to be dealt with */
  skipcomments(infile,instring,FATAL);
  upcase(instring);
  while(instring[0] != '#' || !(strstr(instring,"WALSH") &&
                                strstr(instring,"SYMMETRY")) ){
    skipcomments(infile,instring,FATAL);
    upcase(instring);
  }
  skipcomments(infile,instring,FATAL);
  sscanf(instring,"%d",&num_symm);


  /* get some memory */
  num_p = num_steps * num_orbs;
  points = (walsh_point_type *)calloc(num_p,sizeof(walsh_point_type));
  if( !points ) fatal("Can't get space to store the points.");

  xvals = (real *)calloc(num_steps,sizeof(real));
  if( !xvals )fatal("Can't get space to store xvals.");

  tot_E = (real *)calloc(num_steps,sizeof(real));
  if( !tot_E )fatal("Can't get space to store total energies.");

  /******
    now go through and get space to store the characters
    with respect to the sym ops.
  *******/
  for(i=0;i<num_p;i++){
    points[i].symmetries = (real *)calloc(num_symm,sizeof(real));
    if( !points[i].symmetries )
      fatal("Can't get memory for characters.");
  }

  /********
    loop through the steps contained in the output file.
  ********/
  for(i=0;i<num_steps;i++){

    /* first find the information at the beginning of the step */
    while(instring[0] != '#' || !strstr(instring,"WALSH_STEP")){
      skipcomments(infile,instring,FATAL);
      upcase(instring);
    }

    /* skip over the initial string and the number of the step */
    strtok(instring," ");
    strtok(0," ");

    /* now read out the values of the variable being used as the ordinate */
    j=0;
    while(j != which_var ){
      strtok(0," ");
      j++;
    }
    sscanf(strtok(0," "),"%lf",&(xvals[i]));

    /* okay, read ahead until we hit the energies */
    while(instring[0] != '#' || !strstr(instring,"ENERGIES") ||
          strstr(instring,"FRAGMENT")){
      skipcomments(infile,instring,FATAL);
      upcase(instring);
    }

    /* read'em out */
    for(j=0;j<num_orbs;j++){
      skipcomments(infile,instring,FATAL);
      sscanf(instring,"%s %lf",com_string,&(points[i*num_orbs+j].energy));
    }


    /* get the total energy */
    skipcomments(infile,instring,FATAL);
    sscanf(instring,"%s %lf",com_string,&(tot_E[i]));

    /* find the characters */
    while(instring[0] != '#' || !strstr(instring,"CHARAC")){
      skipcomments(infile,instring,FATAL);
      upcase(instring);
    }

    /* now read them out too (using strtok again) */
    for(j=0;j<num_orbs;j++){
      skipcomments(infile,instring,FATAL);

      /* skip over the initial string */
      strtok(instring," ");

      for(k=0;k<num_symm;k++){
        sscanf(strtok(0," "),"%lf",&(points[i*num_orbs+j].symmetries[k]));
      }
    }
  }

  /* set all the pointers, then return */
  *p_points = points;
  *p_xvals = xvals;
  *p_tot_E = tot_E;
  *p_num_orbs = num_orbs;
  *p_num_symm = num_symm;
  *p_num_steps = num_steps;
}




void main(argc,argv)
  int argc;
  char **argv;
{
  int i,j,k;
  FILE *infile,*outfile;
  char instring[240],com_string[80],type[80];
  walsh_point_type *points,*next;
  real *xvals,*tot_E;
  real *lines;
  int max_p,num_p,num_so_far;
  real broadening;
  int done,num;
  int num_orbs;
  int num_steps,num_symm,step;

  /* open the files */
  if( argc < 2 ){
    fprintf(stderr,"Usage: fit_walsh <input_file>\n");
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

  strcpy(instring,argv[1]);
  strcat(instring,".WALSH");
  outfile = fopen(instring,"w+");
  if(!outfile){
    fprintf(stderr,"Can't open output file: %s\n", instring);
    fatal("Can't open file.");
  }


  /* parse the input file */
  extract_walsh_data(infile,&points,&xvals,&tot_E,&num_orbs,&num_symm,&num_steps);

  /* form the lines */
  construct_lines(points,num_orbs,num_symm,num_steps);

  /********

    put the lines into a form which is easier to print

  ********/
  lines = (real *)calloc(num_orbs*num_steps,sizeof(real));
  if(!lines) fatal("Can't get space for lines.");

  /* loop through the points array and construct the lines */
  for(i=0;i<num_orbs;i++){
    lines[i] = points[i].energy;
    next = points[i].next;
    for(step=1; step < num_steps; step++){
      if( !next )fatal("Something screwy happened with line construction. This is a bug.");
      lines[step*num_orbs+i] = next->energy;
      next = next->next;
    }
  }

  /* put some status information in the output file. */
  fprintf(outfile,"#WALSH_DATA\n");
  fprintf(outfile,"; There are:\n");
  fprintf(outfile,"%d orbitals.\n",num_orbs);
  fprintf(outfile,"%d steps were taken.\n",num_steps);

  /* print out the lines */
  for( i=0;i<num_steps;i++){
    for(j=0;j<num_orbs;j++){
      fprintf(outfile,"%lf ",lines[i*num_orbs+j]);
    }
    fprintf(outfile,"%lf ",tot_E[i]);
    fprintf(outfile,"%lf\n",xvals[i]);
  }
}
