/************************************************************************

  This program reads the values for a hermetian matrix out of a
   binary file, and then generates a .ps file which shows in a
   block way.

   Created by greg Landrum March 1996
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fcntl.h>

#include "fit_props.h"

typedef struct{
  int x,y,z;
} cell_spec_type;

typedef struct{
  int atom1,atom2;
  cell_spec_type cell;
  real dist;
} dist_type;


int sort_distances_helper(dist1,dist2)
  void *dist1, *dist2;
{
  real diff;

  diff = (real)((dist_type *)dist1)->dist -
    (real)((dist_type *)dist2)->dist;
  if( diff > 0 ) return(1);
  else if( diff < 0 ) return(-1);
  else return( 0 );
}

void main(argc, argv)
  int argc;
  char **argv;
{
  FILE *the_file;
  int infile;
  int read_err;
  real min_dist,max_dist,dist_tol;
  float *dist_mat;
  real curr_dist;
  dist_type *dist_array;
  cell_spec_type cell_spec;
  int num_in_array,max_in_array;
  int i,j;
  int num_dim,num_atoms;
  int num_so_far;
  int curr_type;
  char *dummies;

  if(argc < 2){
    fatal("Usage: cooperate infile [min_dist max_dist tolerance]");
  }
  if( argc >= 3){
    sscanf(argv[2],"%lf",&min_dist);
  }else{
    min_dist = .01;
  }
  if( argc >= 4){
    sscanf(argv[3],"%lf",&max_dist);
  }else{
    max_dist = 4.0;
  }
  if( argc >= 5 ){
    sscanf(argv[4],"%lf",&dist_tol);
  } else{
    dist_tol = .01;
  }

  /* open the file */
  infile = open(argv[1],O_RDONLY,"r");
  if( infile == -1 ){
    fatal("Can't open file for binary I/O.");
  }

  /* read out the header */
  read(infile,(const char *)&num_dim,sizeof(int));
  read(infile,(const char *)&num_atoms,sizeof(int));
  fprintf(stderr,"The crystal is %d dimensional and has %d atoms\n",
          num_dim,num_atoms);

  /* get space for the dummies array */
  dummies = (char *)calloc(num_atoms,sizeof(char));
  if(!dummies) fatal("Can't get space for dummies.");

  /* read out the dummies array */
  read_err = read(infile,(const char *)dummies,num_atoms*sizeof(char));
  if( read_err <= 0 ) fatal("Can't read dummy array.");

  /* get space for the distance matrix */
  dist_mat = (float *)malloc(num_atoms*num_atoms*sizeof(float));
  if( !dist_mat ){
    fatal("Can't allocate space for dist_mat");
  }
  /* get some initial space for the dist_array */
  max_in_array = num_atoms*num_atoms;
  dist_array = (dist_type *)calloc(max_in_array,sizeof(dist_type));
  if( !dist_array ) fatal("Can't get space for dist_array");


  /* read out the first distance matrix*/
  read_err = read(infile,(const char *)dist_mat,num_atoms*num_atoms*sizeof(float));
  if( read_err <= 0 ) fatal("Can't read first distance matrix");

  /********

    do the home unit cell, which is slightly different from
    the others in that i < j and we are guaranteed to have
    enough space in the dist_array.

  ********/
  num_in_array = 0;
  for(i=0;i<num_atoms;i++){
    for(j=0;j<i;j++){
      num_so_far = i*num_atoms + j;
      if( dist_mat[num_so_far] >= min_dist &&
         dist_mat[num_so_far] <= max_dist ){
        dist_array[num_in_array].atom1 = i;
        dist_array[num_in_array].atom2 = j;
        dist_array[num_in_array].dist = (real)dist_mat[num_so_far];
        num_in_array++;
      }
    }
  }

  /* read out the location of the next distance matrix*/
  read_err = read(infile,(const char *)&cell_spec,sizeof(cell_spec_type));
  while(read_err > 0 ){
    /* get the next matrix if there is one  */
    read_err = read(infile,(const char *)dist_mat,num_atoms*num_atoms*sizeof(float));
    if(read_err > 0){

      for(i=0;i<num_atoms;i++){
        for(j=0;j<=i;j++){
          num_so_far = i*num_atoms + j;
          if( dist_mat[num_so_far] >= min_dist &&
             dist_mat[num_so_far] <= max_dist ){
            dist_array[num_in_array].atom1 = i;
            dist_array[num_in_array].atom2 = j;
            dist_array[num_in_array].dist = (real)dist_mat[num_so_far];
            dist_array[num_in_array].cell.x = cell_spec.x;
            dist_array[num_in_array].cell.y = cell_spec.y;
            dist_array[num_in_array].cell.z = cell_spec.z;
            num_in_array++;
            if( num_in_array == max_in_array ){
              max_in_array += num_atoms;
              dist_array =
                (dist_type *)realloc((char *)dist_array,
                                     max_in_array*sizeof(dist_type));
              if( !dist_array ) fatal("Can't realloc dist_array");
            }
          }
        }
      }
      read_err = read(infile,(const char *)&cell_spec,sizeof(cell_spec_type));
    }
  }

  close(infile);

  /*********

    hokay, we've read in and stored all the entries we need,
    now sort the array and then print it out

  *********/
  qsort((void *)dist_array,num_in_array,sizeof(dist_type),sort_distances_helper);


  curr_type = 1;
  curr_dist = dist_array[0].dist;
  for(i=0;i<num_in_array;i++){

    if( !(dummies[dist_array[i].atom2]) && !(dummies[dist_array[i].atom1]) ){
      /**********

        check to see if we need to change "types"  i.e. if the distance
        difference is greater than dist_tol

        ***********/
      if( fabs(dist_array[i].dist - curr_dist) > dist_tol ){
        curr_type++;
        curr_dist = dist_array[i].dist;
      }
      printf("Atom XXX%d   %d %d    %d %d %d \t; %f\n",curr_type,dist_array[i].atom1+1,
             dist_array[i].atom2+1,dist_array[i].cell.x,dist_array[i].cell.y,
             dist_array[i].cell.z,dist_array[i].dist);
    }

  }



}




