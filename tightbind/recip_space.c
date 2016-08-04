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
*     this file contains stuff for dealing with reciprocal space,
*      i.e. calculating reciprocal lattice vectors and automagically
*      determining k-point sets.
*
*  created:  greg landrum  December 1997
*
*****************************************************************************/
#include "bind.h"
#include "symmetry.h"


/****************************************************************************
*
*                   Procedure calc_reciprocal_lattice
*
* Arguments:  cell: pointer to cell_type
*
* Returns: none
*
* Action:
*     Uses the direct lattice vectors in 'cell to calculate the
*     corresponding reciprocal lattice.
*
*   The relationships here are, of course:
*    b0 = (a1 x a2) / V
*    b1 = (a2 x a0) / V
*    b2 = (a0 x a1) / V
*  where V = a0 . (a1 x a2) is the cell volume, a[0..2] is the direct latticce
*   and b[0..2] is the reciprocal lattice
*
*****************************************************************************/
void calc_reciprocal_lattice(cell_type *cell)
{
  int i,j;
  point_type tvects[3],rvects[3];
  point_type tempvect;
  real V;

  /* for molecules we shouldn't even be here */
  if( cell->dim < 1 ) return;

  /* generate a copy of the direct lattice */
  for(i=0;i<cell->dim;i++){
    tvects[i].x = cell->atoms[cell->tvects[i].end].loc.x -
      cell->atoms[cell->tvects[i].begin].loc.x;
    tvects[i].y = cell->atoms[cell->tvects[i].end].loc.y -
      cell->atoms[cell->tvects[i].begin].loc.y;
    tvects[i].z = cell->atoms[cell->tvects[i].end].loc.z -
      cell->atoms[cell->tvects[i].begin].loc.z;
  }

  /* for one dimensional systems, there is nothing to do */
  if(cell->dim == 1 ){
    V = sqrt(dot_prod(&tvects[0],&tvects[0]));
    cell->recip_vects[0].x = tvects[0].x / V;
    cell->recip_vects[0].y = tvects[0].y / V;
    cell->recip_vects[0].z = tvects[0].z / V;
    return;
  }

  /*********

    for two dimensional systems, this is easiest if we
    generate a bogus third direct vector which is normalized and
    perpendicular to both a0 and a1:
     a2 = a0 x a1 / |a0 x a1|

  ***********/
  if( cell->dim == 2 ){
    cross_prod(&tvects[0],&tvects[1],&tempvect);
    normalize_vector(&tempvect,&tvects[2]);
  }

  /* calculate the cell volume */
  cross_prod(&tvects[1],&tvects[2],&tempvect);
  V = dot_prod(&tvects[0],&tempvect);

  /* now generate the reciprocal lattice vectors */
  cross_prod(&tvects[1],&tvects[2],&rvects[0]);
  scale_vector(&rvects[0],1.0/V);
  cross_prod(&tvects[2],&tvects[0],&rvects[1]);
  scale_vector(&rvects[1],1.0/V);
  if( cell->dim == 3 ){
    cross_prod(&tvects[0],&tvects[1],&rvects[2]);
    scale_vector(&rvects[2],1.0/V);
  }

  for(i=0;i<cell->dim;i++){
    cell->recip_vects[i].x = rvects[i].x;
    cell->recip_vects[i].y = rvects[i].y;
    cell->recip_vects[i].z = rvects[i].z;
  }
#ifdef DEBUG_SYMM
  printf("Reciprocal Lattice:\n");
  for(i=0;i<cell->dim;i++){
    printf("b[%d] = (%lf %lf %lf)\n",i,rvects[i].x,rvects[i].y,
           rvects[i].z);
  }

  /* check orthonormality */
  for(i=0;i<cell->dim;i++){
    for(j=0;j<cell->dim;j++){
      V = dot_prod(&tvects[i],&rvects[j]);
      printf("a[%d] . b[%d] = %lf\n",i,j,V);
    }
  }
#endif
}

/****************************************************************************
*
*                   Procedure gen_k_point_mesh
*
* Arguments:  points: pointer to pointer to point_type
*            num_per_vect: int[3]
*            num_generated: pointer to int
*            include_high_symm_p: char
*             dim: int
*            offset: real
*
* Returns: none
*
* Action:
*    generates a 'dim dimensional mesh of k points with 'num_per_vect along
*    each direction.
*    if 'include_high_symm_p is nonzero, the zone edges and center
*     will be included, otherwise the points are offset from the
*     edges and center by 'offset.
*    the k points are generated in fractional coordinates in reciprocal
*      space.
*
*****************************************************************************/
void gen_k_point_mesh(point_type **points,int num_per_vect[3],
                      int *num_generated, char include_high_symm_p,int dim,
                      real offset)
{
  int i,j,k;
  int num_p;
  real step[3];
  point_type loc,origin;
  point_type *p_array;
  int num_so_far;

  if(dim < 1) return;

  /* figure out the number of points to be generated and allocate space */
  num_p = 2 * num_per_vect[0];
  if( dim > 1 ){
    num_p *= 2*num_per_vect[1];
  }else{
    num_per_vect[1] = num_per_vect[2] = 0;
  }
  if( dim > 2 ){
    num_p *= 2*num_per_vect[2];
  }else{
    num_per_vect[2] = 0;
  }
  p_array = (point_type *)my_calloc(num_p,sizeof(point_type));
  if(!p_array) fatal("can't allocate p_array");

  /* figure out the spacing between points */
  for(i=0;i<dim;i++){
    if( include_high_symm_p ) step[i] = 0.5 / (real)(num_per_vect[i]-1);
    else step[i] = (0.5 - 2.0*offset) /  (real)(num_per_vect[i]-1);
  }

  if( include_high_symm_p ) offset = 0;

  /*********

    now generate the points themselves

    because of the simple-minded generation algorithm, there are
    going to be duplicate points generated if we are using the high
    symmetry points.  Since duplicates will be removed later, however,
    this doesn't cause problems.

    in order to get multiplicities right later, it's going to
    be important to generate *all* the points, even though we
    are absolutely guaranteed not to need k and -k.

    because we set num_per_vect to 0 for directions beyond the
    dimensionality of the crystal above, we can just loop like this
    is a three dimensional crystal here.

  **********/
  num_so_far = 0;
  for(i=0;i<num_per_vect[0];i++){
    for(j=0;j<num_per_vect[1];j++){
      for(k=0;k<num_per_vect[2];k++){
        p_array[num_so_far].x = offset + (real)i*step[0];
        p_array[num_so_far].y = offset + (real)j*step[1];
        p_array[num_so_far].z = offset + (real)k*step[2];
        num_so_far++;
        p_array[num_so_far].x = (offset + (real)i*step[0]);
        p_array[num_so_far].y = (offset + (real)j*step[1]);
        p_array[num_so_far].z = -(offset + (real)k*step[2]);
        num_so_far++;
        p_array[num_so_far].x = (offset + (real)i*step[0]);
        p_array[num_so_far].y = -(offset + (real)j*step[1]);
        p_array[num_so_far].z = -(offset + (real)k*step[2]);
        num_so_far++;
        p_array[num_so_far].x = (offset + (real)i*step[0]);
        p_array[num_so_far].y = -(offset + (real)j*step[1]);
        p_array[num_so_far].z = offset + (real)k*step[2];
        num_so_far++;

        p_array[num_so_far].x = -(offset + (real)i*step[0]);
        p_array[num_so_far].y = -(offset + (real)j*step[1]);
        p_array[num_so_far].z = -(offset + (real)k*step[2]);
        num_so_far++;
        p_array[num_so_far].x = -(offset + (real)i*step[0]);
        p_array[num_so_far].y = -(offset + (real)j*step[1]);
        p_array[num_so_far].z = (offset + (real)k*step[2]);
        num_so_far++;
        p_array[num_so_far].x = -(offset + (real)i*step[0]);
        p_array[num_so_far].y = (offset + (real)j*step[1]);
        p_array[num_so_far].z = (offset + (real)k*step[2]);
        num_so_far++;
        p_array[num_so_far].x = -(offset + (real)i*step[0]);
        p_array[num_so_far].y = (offset + (real)j*step[1]);
        p_array[num_so_far].z = -(offset + (real)k*step[2]);
        num_so_far++;

      }
                        if(dim == 2){
      p_array[num_so_far].x = (offset + (real)i*step[0]);
      p_array[num_so_far].y = (offset + (real)j*step[1]);
      num_so_far++;
      p_array[num_so_far].x = (offset + (real)i*step[0]);
      p_array[num_so_far].y = -(offset + (real)j*step[1]);
      num_so_far++;

      p_array[num_so_far].x = -(offset + (real)i*step[0]);
      p_array[num_so_far].y = -(offset + (real)j*step[1]);
      num_so_far++;
      p_array[num_so_far].x = -(offset + (real)i*step[0]);
      p_array[num_so_far].y = (offset + (real)j*step[1]);
      num_so_far++;
      }
    }
          if( dim == 1 ){
    p_array[num_so_far].x = (offset + (real)i*step[0]);
    num_so_far++;

    p_array[num_so_far].x = -(offset + (real)i*step[0]);
    num_so_far++;
    }
  }

        *points = p_array;
  *num_generated = num_so_far;
}


/****************************************************************************
*
*                   Procedure automagic_k_points
*
* Arguments:  details: pointer to detail_type
*            cell: pointer to cell_type
*
* Returns: none
*
* Action:
*   Automagically generates a k-points set for 'cell
*
*****************************************************************************/
void automagic_k_points(detail_type *details,cell_type *cell)
{
  int i,j;
  point_type *raw_points;
  real *degen_p;
  int num_raw_points,num_points;
  real tot_weight;
  int nonorthogonal_lattice;
        int num_duplicates;

  /* start by generating a reducible mesh */
  gen_k_point_mesh(&raw_points,details->points_per_axis,&num_raw_points,
                   details->use_high_symm_p,cell->dim,details->k_offset);

  /* okay, now we need to generate the reciprocal lattice */
  calc_reciprocal_lattice(cell);

  /********

    that's in hand, now we can reduce the k_point set
      start with the easy ones: removing identical points,
      then move on to using symmetry to reduce the set.

  *********/
  degen_p = (real *)my_calloc(num_raw_points,sizeof(real));
  if(!degen_p)fatal("can't allocate degen_p array.");
          num_duplicates = 0;
  for(i=0;i<num_raw_points;i++){
    if( degen_p[i] == 0 ){
      for(j=i+1;j<num_raw_points;j++){
        if( POINTS_ARE_THE_SAME(&raw_points[i],&raw_points[j],
                                details->symm_tol) ){
          degen_p[j] = -1;
          num_duplicates++;
        }
      }
    }
  }

#ifdef DEBUG_SYMM
  printf("Here's the reducible k-point mesh (%d duplicates removed) for (%d %d %d) points per axis,\n",
         num_duplicates,details->points_per_axis[0],details->points_per_axis[1],
         details->points_per_axis[2]);
  printf("and an offset of %lf\n",details->k_offset);
  for(i=0;i<num_raw_points;i++){
    if( degen_p[i] == 0 ){
      printf("% 6.4lf % 6.4lf % 6.4lf\n",raw_points[i].x,raw_points[i].y,
             raw_points[i].z);
    }
  }
#endif

  /*******

    if the reciprocal lattice is not orthogonal, removing symmetry
    equivalent k-points requires a lot more smarts that the program
    currently has, so check for this, warn the user, and
    go with half the Brillouin zone.

  *******/
  nonorthogonal_lattice = check_for_orthogonal_basis(cell->recip_vects,cell->dim,
                                                details->symm_tol);
  if( nonorthogonal_lattice ){
    printf("WARNING: The lattice is non-orthogonal.\n");
    printf("\tI am not currently smart enough to generate a\n");
    printf("\tminimal k-point set for non-orthogonal lattices.\n");
    printf("\tI am going to make an attempt to reduce the zone somewhat,\n");
    printf("\tbut you should check the resulting k-point set *carefully*.\n");
    printf("\tFor a more efficent sampling, you will have to use one of the\n");
    printf("\tdistributed k-point programs.\n");
    fprintf(status_file,"WARNING: The lattice is non-orthogonal.\n");
    fprintf(status_file,"\tI am not currently smart enough to generate a\n");
    fprintf(status_file,"\tminimal k-point set for non-orthogonal lattices.\n");
    fprintf(status_file,
            "\tI am going to make an attempt to reduce the zone somewhat,\n");
    fprintf(status_file,
            "\tbut you should check the resulting k-point set *carefully*.\n");
    fprintf(status_file,
            "\tfor a more efficent sampling, you will have to use one of the\n");
    fprintf(status_file,"\tdistributed k-point programs.\n");

    fprintf(output_file,";WARNING: The lattice is non-orthogonal.\n");
    fprintf(output_file,";\tI am not currently smart enough to generate a\n");
    fprintf(output_file,";\tminimal k-point set for non-orthogonal lattices.\n");
    fprintf(output_file,
            ";\tI am going to make an attempt to reduce the zone somewhat,\n");
    fprintf(output_file,
            ";\tbut you should check the resulting k-point set *carefully*.\n");
    fprintf(output_file,
            ";\tfor a more efficent sampling, you will have to use one of the\n");
    fprintf(output_file,";\tdistributed k-point programs.\n");
    reduce_kpoints(details,cell,raw_points,num_raw_points,degen_p,
                              &num_points,nonorthogonal_lattice);
  }else{
#ifdef DEBUG_SYMM
    printf("Good, the lattice is orthogonal.\n");
#endif
    reduce_kpoints(details,cell,raw_points,num_raw_points,degen_p,
                              &num_points,nonorthogonal_lattice);
  }
#ifdef DEBUG_SYMM
  printf("Of the %d unique points originally present, %d survived:\n",
                                  num_raw_points-num_duplicates,num_points);
  tot_weight = 0;
  for(i=0;i<num_raw_points;i++){
    if( degen_p[i] > 0 ){
      printf("(%6.4lf %6.4lf %6.4lf) %4.2lf\n",raw_points[i].x,raw_points[i].y,
             raw_points[i].z,degen_p[i]);
      tot_weight += degen_p[i];
    }
  }
  printf("The sum of the k-weights is: %6.3lf\n",tot_weight);
#endif


}
