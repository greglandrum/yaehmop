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

  This program reads FCO data out of a binary file, broadens it,
   and writes the results to a text file.

   Created by greg Landrum June 1996
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fcntl.h>
#include "fit_props.h"


#define E_CUTOFF 3
#define NUM_TABLE_ENTRIES 10000
#define ZERO_TOL 1e-5

typedef struct{
  int num_entries;
  real min_val,max_val;
  real step;
  real *values;
} lookup_table_type;

#define READ_FROM_LOOKUP_TBL(__tbl__,__x__) \
(( (__x__) < __tbl__.max_val && (__x__) > __tbl__.min_val ) ? \
 __tbl__.values[(int)floor(((__x__)-__tbl__.min_val)/__tbl__.step + 0.5)] : 0.0)
#define INSERT_INTO_LOOKUP_TBL(__tbl__,__x__,__val__) \
 ((__x__<__tbl__.max_val && __x__>__tbl__.min_val) ?  \
  __tbl__.values[(int)floor((__x__-__tbl__.min_val)/__tbl__.step + 0.5)]=__val__ : 0)

lookup_table_type gauss_lookup_table;


typedef struct{
  real value,total_E,frag_E;
  int which_frag;
} FCO_point_type;


void build_lookup_tbl(lookup_table_type *table)
{
  int i;
  real step, val,loc;

  fprintf(stderr,"Building lookup table.");
  table->min_val = 0;
  /* use the zero tolerance to set the limits on the table */
  table->max_val = -log((real)ZERO_TOL);
  table->num_entries = NUM_TABLE_ENTRIES;
  step = (table->max_val - table->min_val) / table->num_entries;
  table->step = step;

  table->values = (real *)calloc(table->num_entries,sizeof(real));
  if(!table->values)fatal("Can't allocate lookup table entries");

  loc = table->min_val;
  for(i=0;i<table->num_entries;i++){
    table->values[i] = exp(-loc);
    loc += step;
  }
}



int sort_FCO_helper(const void *p1,const void *p2)
{
  FCO_point_type *FCO1,*FCO2;
  FCO1 = (FCO_point_type *)p1;
  FCO2 = (FCO_point_type *)p2;

  /* sort by both total E and fragment E */
  if( FCO1->total_E > FCO2->total_E ) return(1);
  else if( FCO1->total_E < FCO2->total_E ) return(-1);
  else if( FCO1->frag_E > FCO2->frag_E ) return(1);
  else if( FCO1->frag_E < FCO2->frag_E ) return(-1);
  else return(0);
}



#ifndef USING_THE_MAC
void main(argc, argv)
  int argc;
  char **argv;
#else
void main()
#endif
{
  int FCO_file;
  int test_int;
  FILE *outfile,*the_file;
  char filename[240];
  int num_Kpoints,num_orbs,num_FCO_frags;
  real *total_Es,*frag_Es,*charge_mat;
  real *results;
  FCO_point_type *FCO_point_array,curr_FCO_p,*this_point;
  int FCO_points_so_far;
  int *FCO_num_orbs;
  real tot_K_weight;
  real E_min,E_max,E_step,broadening,norm_fact,E_diff;
  real curr_tot_E, curr_frag_E;
  real frag_E_diff,tot_E_diff,value;
  real *frag_vals;
  real total_DOS,tempval;
  int i,j,k,itab,jtab,ktab,points_per_side;
  int which_frag,orbs_this_frag;

long int calls_to_exp = 0;
 #ifdef USING_THE_MAC
  int argc;
  char argv[4][80];

        /* set up some stuff for Sioux */
        //SIOUXSettings.standalone = FALSE;
        SIOUXSettings.asktosaveonclose = FALSE;
        SIOUXSettings.autocloseonquit = FALSE;
        printf("Starting fit_FCO.\n");

  the_file = choose_mac_file(argv[1],MAC_FOPEN_OPEN_CD);
  if( !the_file ) {
          fatal("User cancelled intial file open");
  } else{
          argc = 2;
  }

  /* get the command line arguments */
//  argc = ccommand(&argv);

#endif

  if(argc < 2){
    fatal("Usage: fit_FCO <infile>");
  }


  /* open the input file */
  sprintf(filename,"%s.FCO",argv[1]);
#ifndef USING_THE_MAC
  FCO_file = open(filename,O_RDONLY,"r");
#else
  FCO_file = open(filename,O_RDONLY);
#endif
  if( FCO_file == -1 ){
    error("Can't open file for binary I/O.");
    return;
  }

  sprintf(filename,"%s.GRID",argv[1]);
  outfile = fopen(filename,"w+");
  if( !outfile ) fatal("Can't open output file.");

  /* read out the header information */


  /* first read the test int to make sure things are hunky dorey */
  read(FCO_file,(const char *)&test_int,sizeof(int));
  if( test_int != 4231 ){
    fprintf(stderr,"fit_FCO just failed while reading in the test int value.\n");
    fprintf(stderr,"  the point of this is to make sure the binary file you give\n");
    fprintf(stderr,"  is okay.  This does not seem to be the case.\n\n");
    fprintf(stderr,"Please make sure that you typed the file name properly\n");
    fprintf(stderr,"  and that you are on the same kind of system that you\n");
    fprintf(stderr,"  used to write the file originally.\n");
    fatal("Can't deal with file\n");
  }
  read(FCO_file,(const char *)&num_Kpoints,sizeof(int));
  read(FCO_file,(const char *)&tot_K_weight,sizeof(real));
  read(FCO_file,(const char *)&num_orbs,sizeof(int));
  read(FCO_file,(const char *)&num_FCO_frags,sizeof(int));

  /* get memory for the FCO_num_orbs now... */
  frag_vals = (real *)calloc(num_FCO_frags,sizeof(real));
  FCO_num_orbs = (int *)calloc(num_FCO_frags,sizeof(int));
  if( !FCO_num_orbs || !frag_vals) fatal("Can't get memory for FCO_num_orbs");

  read(FCO_file,(const char *)FCO_num_orbs,num_FCO_frags*sizeof(int));


  /* get memory to store everything we're gonna need */
  total_Es = (real *)calloc(num_orbs,sizeof(real));
  frag_Es = (real *)calloc(num_orbs,sizeof(real));
  charge_mat = (real *)calloc(num_orbs*num_orbs,sizeof(real));
  if( !charge_mat || !frag_Es ) fatal("Can't get K point memory");

  FCO_point_array = (FCO_point_type *)calloc(num_orbs*num_orbs*num_Kpoints,
                                             sizeof(FCO_point_type));
  if(!FCO_point_array)fatal("Can't get memory for FCO_point_array");


  /* prompt for the energy window */
  printf("\nEnter E min: ");
  scanf("%lf",&E_min);
  printf("Enter E max: ");
  scanf("%lf",&E_max);
  printf("Enter broadening: ");
  scanf("%lf",&broadening);
  printf("Enter Energy Step: ");
  scanf("%lf",&E_step);

  fprintf(stderr,"Reading...\n");

  /*******

     okay, go ahead, read out all the data and store it in the FCO_point_array...
     woo hoo!

  ********/
  FCO_points_so_far=0;
  for(k=0;k<num_Kpoints;k++){
    read(FCO_file,(const char *)total_Es,num_orbs*sizeof(real));
    read(FCO_file,(const char *)frag_Es,num_orbs*sizeof(real));
    read(FCO_file,(const char *)charge_mat,num_orbs*num_orbs*sizeof(real));

    /* loop over total orbs */
    for(i=0;i<num_orbs;i++){
      itab = i*num_orbs;
      /* loop over fragment orbs */
      which_frag = 0;
      orbs_this_frag = 0;
      for(j=0;j<num_orbs;j++){
        if( orbs_this_frag >= FCO_num_orbs[which_frag] ){
          which_frag++;
          orbs_this_frag = 0;
        }

        /* make sure the point falls in our energy window */
        if( total_Es[i] <= E_max && total_Es[i] >= E_min &&
           frag_Es[j] <= E_max && frag_Es[j] >= E_min ){
          FCO_point_array[FCO_points_so_far].which_frag = which_frag;
          FCO_point_array[FCO_points_so_far].total_E =
            total_Es[i];
          FCO_point_array[FCO_points_so_far].frag_E =
            frag_Es[j];
          FCO_point_array[FCO_points_so_far].value =
            charge_mat[itab+j];
          FCO_points_so_far++;
        }
        orbs_this_frag++;
      }
    }
  }
#if 0
  fprintf(stderr,"Sorting...\n");
  /* okay, we've got all the data, now sort it */
  qsort((void *)FCO_point_array,FCO_points_so_far,
        sizeof(FCO_point_type),sort_FCO_helper);

  /**************

    the data is now sorted in order of increasing total_E, with
    each total_E row containing columns of frag_Es in increasing
    order.... now we need to go through and broaden the data and
    write it to the output file.

    ****************/
#endif
  /* write the header now */
  fprintf(outfile,"#FCO_DATA\n");
  fprintf(outfile,"#MIN_X %lf\n",E_min);
  fprintf(outfile,"#MAX_X %lf\n",E_max);
  fprintf(outfile,"#STEP_X %lf\n",E_step);
  fprintf(outfile,"#MIN_Y %lf\n",E_min);
  fprintf(outfile,"#MAX_Y %lf\n",E_max);
  fprintf(outfile,"#STEP_Y %lf\n",E_step);
  fprintf(outfile,"#NUM_CURVES %d\n",num_FCO_frags);

  points_per_side = 0;
  curr_tot_E = E_min;
  while(curr_tot_E <= E_max){
    curr_tot_E += E_step;
    points_per_side++;
  }
  fprintf(outfile,"#NUM_X %d\n",points_per_side);
  fprintf(outfile,"#NUM_Y %d\n",points_per_side);

#if 1
  build_lookup_tbl(&gauss_lookup_table);
#endif
  fprintf(stderr,"\nBroadening...\n");
  /* determine the gaussian normalization factor */
  norm_fact = sqrt(broadening/(real)M_PI)/(2.0*tot_K_weight);


  fprintf(outfile,"\n#BEGIN_DATA\n");

  /********

    get space to store the results.

    it's important that this be done with a calloc, we're
    gonna rely on the fact that uninitialized values are
    zero to make this mess more efficient.

    This is one of the rare cases where I think efficiency
    outweighs readability concerns.  That's because the
    naive implementation of this takes *way* too long.

  ********/
  results = (real *)calloc(points_per_side*points_per_side*num_FCO_frags,
                           sizeof(real));
  if( !results ) fatal("Can't allocate results matrix.");


  for(i=0;i<FCO_points_so_far;i++){
    this_point = &(FCO_point_array[i]);
    itab = this_point->which_frag*points_per_side*points_per_side;

    /* loop over total_E until we're in the energy window */
    curr_tot_E = E_min;
    for(k=0;k<points_per_side;k++,curr_tot_E += E_step){
      tot_E_diff = fabs(this_point->total_E - curr_tot_E);
      if( tot_E_diff < E_CUTOFF ){
        ktab = itab + k*points_per_side;
        /* loop over frag_E until we're within the E cutoff window */
        curr_frag_E = E_min;
        for(j=0;j<points_per_side;j++,curr_frag_E += E_step){
          E_diff = tot_E_diff + fabs(this_point->frag_E - curr_frag_E);

          if( E_diff <= E_CUTOFF ){
            tempval = broadening*E_diff*E_diff;
            results[ktab + j] +=
              this_point->value * norm_fact *
#if 1
                READ_FROM_LOOKUP_TBL(gauss_lookup_table,tempval);
#else
            exp(-tempval);
#endif
calls_to_exp++;
          }
        }
      }  else if( curr_tot_E > this_point->total_E) break;
    }
  }

fprintf(stderr,"num_calls_to_exp: %ld\n",calls_to_exp);
  /*******

    now write the values

  *******/
  fprintf(stderr,"Writing...\n");
  for(i=0;i<points_per_side;i++){
    for(j=0;j<points_per_side;j++){
      for(k=0;k<num_FCO_frags;k++){
        value = results[(k*points_per_side + j)*points_per_side + i];
        if( value >= ZERO_TOL ){
          fprintf(outfile,"%lf ",value);
        }else{
          fprintf(outfile,"-.001 ");
        }
      }
      fprintf(outfile,"\n");
    }
    fprintf(outfile,"\n");
  }


  /* write the broadened total and fragment DOSs now */
  fprintf(outfile,"\n#BEGIN_DOS\n");
  curr_tot_E = E_min;
  while(curr_tot_E <= E_max){
    total_DOS = 0;
    for(i=0;i<num_FCO_frags;i++) frag_vals[i] = 0;

    for(i=0;i<FCO_points_so_far;i++){
      E_diff = fabs(FCO_point_array[i].total_E - curr_tot_E);
      if( E_diff < 3.0 ){
        total_DOS += norm_fact * exp(-broadening*E_diff*E_diff);
      }
      E_diff = fabs(FCO_point_array[i].frag_E - curr_tot_E);
      if( E_diff < 3.0 ){
        frag_vals[FCO_point_array[i].which_frag] +=
            norm_fact * exp(-broadening*E_diff*E_diff);
      }
    }
    if( total_DOS > 0.0001 )
      fprintf(outfile,"%6.3lg ",total_DOS);
    else
      fprintf(outfile,"0 ");

    for(i=0;i<num_FCO_frags;i++){
      if( fabs(frag_vals[i]) > 0.0001 )
        fprintf(outfile,"%6.3lg ",frag_vals[i]);
      else
        fprintf(outfile,"0 ");

    }
    fprintf(outfile,"\n");

    curr_tot_E += E_step;
  }
  fprintf(outfile,"\n");
  close(FCO_file);
  fclose(outfile);
}



