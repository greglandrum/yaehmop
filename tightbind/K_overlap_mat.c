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
*     this file contains stuff for the overlap matrices in K space
*
*  created:  greg landrum  September 1993
*
*****************************************************************************/
#include "bind.h"


/******
  the sizes of the transformation matrices are increased a little bit
  to make the code more readable
*******/
#define P_SIZE 9+BEGIN_P
#define D_SIZE 25+BEGIN_D
#define F_SIZE 49+BEGIN_F




/****************************************************************************
*
*                   Procedure build_K_overlap_FAT
*
* Arguments:    cell: pointer to cell type
*             kpoint: pointer to k_point type
*  overlapR,overlapK: hermetian_matrix_type
*           num_orbs: int
*
* Returns: none
*
* Action: This just performs the weighted sum of the R-overlaps for the given
*   k-point.
*
****************************************************************************/
void build_k_overlap_FAT(cell,kpoint,overlapR,overlapK,num_orbs)
  cell_type *cell;
  k_point_type *kpoint;
  hermetian_matrix_type overlapR,overlapK;
  int num_orbs;
{
  int i,j,k,l,m;
  int itab,jtab,ktab;
  int ltab,mtab;

  point_type kpointloc;

  real kdotR,temp;
  real cos_term,sin_term;
  real *which_overlap;


  kpointloc.x = TWOPI*kpoint->loc.x;
  kpointloc.y = TWOPI*kpoint->loc.y;
  kpointloc.z = TWOPI*kpoint->loc.z;


  /*****
    we'll use which_overlap to keep track of our location within the
    R space overlap matrix.
  ******/
  which_overlap = overlapR.mat;

  /* copy the unit cell overlap values into the k space matrix */
  for(l=0;l<num_orbs;l++){
    ltab = l*num_orbs;
    for(m=0;m<=l;m++){
      mtab = m*num_orbs;
      overlapK.mat[ltab+m] = 0.0;
      overlapK.mat[mtab+l] = overlapR.mat[mtab+l];
    }
  }

  /* sum up the individual overlaps, just like when overlapR was built */
  for(i=1;i<=cell->overlaps[0];i++){
    which_overlap += num_orbs*num_orbs;

    kdotR = kpointloc.x*(real)i;
    cos_term = cos(kdotR);
    sin_term = sin(kdotR);
    for(l=0;l<num_orbs-1;l++){
      ltab = l*num_orbs;
      for(m=l+1;m<num_orbs;m++){
        mtab = m*num_orbs;
        /* real part */
        overlapK.mat[ltab+m]+= cos_term*which_overlap[ltab+m];
        /* imaginary part */
        overlapK.mat[mtab+l] -= sin_term*which_overlap[mtab+l];
      }
      /* diagonal element (real) */
      overlapK.mat[ltab+l] += cos_term*which_overlap[ltab+l];
    }
    /* the last diagonal element */
    overlapK.mat[num_orbs*num_orbs-1] += cos_term*which_overlap[num_orbs*num_orbs-1];
  }

  if( cell->dim != 1 ){

    /******
      take care of the remainder of the layer containing the unit cell
      for each overlap the procedure is the same as above
    *******/
    for(i=0;i<=2*cell->overlaps[0];i++){
      itab = cell->overlaps[0]-i;
      for(j=1;j<=cell->overlaps[1];j++){
        which_overlap += num_orbs*num_orbs;
        kdotR = kpointloc.x*(real)itab+kpointloc.y*(real)j;
        cos_term = cos(kdotR);
        sin_term = sin(kdotR);
        for(l=0;l<num_orbs-1;l++){
          ltab = l*num_orbs;
          for(m=l+1;m<num_orbs;m++){
            mtab = m*num_orbs;
            overlapK.mat[ltab+m]+= cos_term*which_overlap[ltab+m];
            overlapK.mat[mtab+l] -= sin_term*which_overlap[mtab+l];
          }
          overlapK.mat[ltab+l] += cos_term*which_overlap[ltab+l];
        }
        overlapK.mat[num_orbs*num_orbs-1] += cos_term*which_overlap[num_orbs*num_orbs-1];
      }
    }
  }
  if( cell->dim == 3 ){
    /* do the layers above the one containing the unit cell */
    for(i=1;i<=cell->overlaps[2];i++){
      itab = i;
      for(j=0;j<=2*cell->overlaps[0];j++){
        jtab = cell->overlaps[0] - j;
        for(k=0;k<=2*cell->overlaps[1];k++){
          ktab = cell->overlaps[1] - k;
          which_overlap += num_orbs*num_orbs;
          kdotR = kpointloc.x*(real)jtab+kpointloc.y*(real)ktab+
            kpointloc.z*(real)i;
          cos_term = cos(kdotR);
          sin_term = sin(kdotR);
          for(l=0;l<num_orbs-1;l++){
            ltab = l*num_orbs;
            for(m=l+1;m<num_orbs;m++){
              mtab = m*num_orbs;
              overlapK.mat[ltab+m]+= cos_term*which_overlap[ltab+m];
              overlapK.mat[mtab+l] -= sin_term*which_overlap[mtab+l];
            }
            overlapK.mat[ltab+l] += cos_term*which_overlap[ltab+l];
          }
          overlapK.mat[num_orbs*num_orbs-1] += cos_term*which_overlap[num_orbs*num_orbs-1];
        }
      }
    }
  }


#ifdef PRINTMAT
fprintf(output_file,"---------- S(k) ------\n");
printmat(overlapK.mat,num_orbs,num_orbs,output_file,1e-6,details->line_width);
#endif

  /* that's it! */
}



/****************************************************************************
*
*                   Procedure build_K_overlap_THIN
*
* Arguments:    cell: pointer to cell type
*            details: pointer to detail_type
*             kpoint: pointer to k_point type
*  overlapR,overlapK: pointer to real
*           num_orbs: int
*
* Returns: none
*
* Action: This generates the overlap matrix for a given k-point
*
****************************************************************************/
void build_k_overlap_THIN(cell,details,kpoint,overlapR,overlapK,num_orbs)
  cell_type *cell;
  detail_type *details;
  k_point_type *kpoint;
  hermetian_matrix_type overlapR,overlapK;
  int num_orbs;
{
  int i,j,k,l,m;
  int itab,jtab,ktab;
  int ltab,mtab;

  point_type kpointloc;

  real kdotR;
  real cos_term,sin_term;
  real *which_overlap;

  static char first_call=1;
  static point_type cell_dim[3];
  static point_type distances;
  real temp,min=100.0;
  int min_dir;

  /* this is some stuff that only needs to be done once */
  if( first_call ){
    first_call = 0;
    /* find the dimensions of the unit cell */
    for(i=0;i<cell->dim;i++){
      itab = cell->tvects[i].begin;
      jtab = cell->tvects[i].end;
      cell_dim[i].x = cell->atoms[jtab].loc.x-cell->atoms[itab].loc.x;
      cell_dim[i].y = cell->atoms[jtab].loc.y-cell->atoms[itab].loc.y;
      cell_dim[i].z = cell->atoms[jtab].loc.z-cell->atoms[itab].loc.z;

      /* we use this to keep track of the shortest distance */
      temp = sqrt(cell_dim[i].x*cell_dim[i].x+
                  cell_dim[i].y*cell_dim[i].y+
                  cell_dim[i].z*cell_dim[i].z);
      if( temp < min ){
        min = temp;
        min_dir = i;
      }
    }

    /* a quick value for rho */
    if( fabs(details->rho) <= 1e-3 )
      details->rho = min*(cell->overlaps[min_dir]+1)+.01;
    fprintf(output_file,"\n\n; RHO = %lf\n",details->rho);
  }

  kpointloc.x = TWOPI*kpoint->loc.x;
  kpointloc.y = TWOPI*kpoint->loc.y;
  kpointloc.z = TWOPI*kpoint->loc.z;


  /******

    first do the unit cell

  ******/
  distances.x=distances.y=distances.z=0.0;
  calc_R_overlap(overlapR.mat,cell,details,
                 num_orbs,distances,TRUE,orbital_lookup_table);

  /*********
    copy the unit cell overlap values into the k space matrix
    and put 1's on the diagonal and copy the elements across the diagonal
  **********/
  for(j=0;j<num_orbs;j++){
    jtab = j*num_orbs;
    overlapK.mat[jtab+j] = 1.0;
    for(k=0;k<j;k++){
      ktab = k*num_orbs;
      overlapK.mat[jtab+k] = 0.0;
      overlapK.mat[ktab+j]=overlapR.mat[jtab+k];
    }
  }

  /*******

    now move into the other cells.

  ********/
  for(i=1;i<=cell->overlaps[0];i++){
    distances.x = i*cell_dim[0].x;
    distances.y = i*cell_dim[0].y;
    distances.z = i*cell_dim[0].z;
    calc_R_overlap(overlapR.mat,cell,details,
                   num_orbs,distances,FALSE,orbital_lookup_table);

    kdotR = kpointloc.x*(real)i;
    cos_term = cos(kdotR);
    sin_term = sin(kdotR);
    for(l=0;l<num_orbs-1;l++){
      ltab = l*num_orbs;
      for(m=l+1;m<num_orbs;m++){
        mtab = m*num_orbs;
        /* real part */
        overlapK.mat[ltab+m]+= cos_term*overlapR.mat[ltab+m];
        /* imaginary part */
        overlapK.mat[mtab+l] -= sin_term*overlapR.mat[mtab+l];
      }
      /* diagonal element (real) */
      overlapK.mat[ltab+l] += cos_term*overlapR.mat[ltab+l];
    }
    /* the last diagonal element */
    overlapK.mat[num_orbs*num_orbs-1] += cos_term*overlapR.mat[num_orbs*num_orbs-1];
  }

  if( cell->dim != 1 ){

    /******
      take care of the remainder of the layer containing the unit cell
      for each overlap the procedure is the same as above
    *******/
    for(i=0;i<=2*cell->overlaps[0];i++){
      itab = cell->overlaps[0] - i;
      for(j=1;j<=cell->overlaps[1];j++){
        distances.x = itab*cell_dim[0].x + j*cell_dim[1].x;
        distances.y = itab*cell_dim[0].y + j*cell_dim[1].y;
        distances.z = itab*cell_dim[0].z + j*cell_dim[1].z;

        calc_R_overlap(overlapR.mat,cell,details,
                       num_orbs,distances,FALSE,orbital_lookup_table);

        kdotR = kpointloc.x*(real)itab+kpointloc.y*(real)j;
        cos_term = cos(kdotR);
        sin_term = sin(kdotR);
        for(l=0;l<num_orbs-1;l++){
          ltab = l*num_orbs;
          for(m=l+1;m<num_orbs;m++){
            mtab = m*num_orbs;
            overlapK.mat[ltab+m]+= cos_term*overlapR.mat[ltab+m];
            overlapK.mat[mtab+l] -= sin_term*overlapR.mat[mtab+l];
          }
          overlapK.mat[ltab+l] += cos_term*overlapR.mat[ltab+l];
        }
        overlapK.mat[num_orbs*num_orbs-1] += cos_term*overlapR.mat[num_orbs*num_orbs-1];
      }
    }
  }
  if( cell->dim == 3 ){
    /* do the layers above the one containing the unit cell */
    for(i=1;i<=cell->overlaps[2];i++){
      itab = i;
      for(j=0;j<=2*cell->overlaps[0];j++){
        jtab = cell->overlaps[0] - j;
        for(k=0;k<=2*cell->overlaps[1];k++){
          ktab = cell->overlaps[1] - k;

          distances.x = itab*cell_dim[2].x + jtab*cell_dim[0].x +
            ktab*cell_dim[1].x;
          distances.y = itab*cell_dim[2].y + jtab*cell_dim[0].y +
            ktab*cell_dim[1].y;
          distances.z = itab*cell_dim[2].z + jtab*cell_dim[0].z +
            ktab*cell_dim[1].z;

          calc_R_overlap(overlapR.mat,cell,details,
                         num_orbs,distances,FALSE,orbital_lookup_table);

          kdotR = kpointloc.x*(real)jtab+kpointloc.y*(real)ktab+
            kpointloc.z*(real)i;
          cos_term = cos(kdotR);
          sin_term = sin(kdotR);
          for(l=0;l<num_orbs-1;l++){
            ltab = l*num_orbs;
            for(m=l+1;m<num_orbs;m++){
              mtab = m*num_orbs;
              overlapK.mat[ltab+m]+= cos_term*overlapR.mat[ltab+m];
              overlapK.mat[mtab+l] -= sin_term*overlapR.mat[mtab+l];
            }
            overlapK.mat[ltab+l] += cos_term*overlapR.mat[ltab+l];
          }
          overlapK.mat[num_orbs*num_orbs-1] += cos_term*overlapR.mat[num_orbs*num_orbs-1];
        }
      }
    }
  }

#ifdef PRINTMAT
fprintf(output_file,"---------- S(k) ------\n");
printmat(overlapK.mat,num_orbs,num_orbs,output_file,1e-6,details->line_width);
#endif
}



/****************************************************************************
 *
 *                   Procedure build_all_K_overlaps
 *
 * Arguments:    cell: pointer to cell type
 *            details: pointer to detail_type
 *  overlapR,overlapK: hermetian_matrix_type
 *           num_orbs: int
 *
 * Returns: none
 *
 * Action: This performs the weighted sum of the R-overlaps for all
 *   k-points.   The R overlap matrices are generated on the fly to conserve
 *   memory
 *
 ****************************************************************************/
void build_all_K_overlaps(cell,details,overlapR,overlapK,num_orbs,
                          tot_overlaps,orbital_lookup_table)
  cell_type *cell;
  detail_type *details;
  hermetian_matrix_type overlapR,overlapK;
  int num_orbs,tot_overlaps;
  int *orbital_lookup_table;
{
  k_point_type *kpoint;
  int which_k;
  int i,j,k,l,m;
  int itab,jtab,ktab;
  int ltab,mtab;

  point_type kpointloc;

  real kdotR,temp;
  real cos_term,sin_term;
  real *which_Koverlap;
  int overlaps_so_far;


  overlaps_so_far = 0;
  R_space_overlap_matrix(cell,details,overlapR,num_orbs,tot_overlaps,
                         orbital_lookup_table,overlaps_so_far);
  overlaps_so_far++;
  /* copy the unit cell overlap values into the k space matrix */
  for(l=0;l<num_orbs;l++){
    ltab = l*num_orbs;
    for(m=0;m<=l;m++){
      mtab = m*num_orbs;
      for(which_k=0;which_k<details->num_KPOINTS;which_k++){
        kpoint = &(details->K_POINTS[which_k]);
        which_Koverlap = &(overlapK.mat[which_k*(num_orbs*num_orbs)]);
        kpointloc.x = TWOPI*kpoint->loc.x;
        kpointloc.y = TWOPI*kpoint->loc.y;
        kpointloc.z = TWOPI*kpoint->loc.z;


        which_Koverlap[ltab+m] = 0.0;
        which_Koverlap[mtab+l] = overlapR.mat[mtab+l];
      }
    }
  }

  /* sum up the individual overlaps, just like when overlapR was built */
  for(i=1;i<=cell->overlaps[0];i++){
    R_space_overlap_matrix(cell,details,overlapR,num_orbs,tot_overlaps,
                           orbital_lookup_table,overlaps_so_far);
    overlaps_so_far++;

    for(which_k=0;which_k<details->num_KPOINTS;which_k++){
      kpoint = &(details->K_POINTS[which_k]);
      which_Koverlap = &(overlapK.mat[which_k*(num_orbs*num_orbs)]);
      kpointloc.x = TWOPI*kpoint->loc.x;
      kpointloc.y = TWOPI*kpoint->loc.y;
      kpointloc.z = TWOPI*kpoint->loc.z;
      kdotR = kpointloc.x*(real)i;
      cos_term = cos(kdotR);
      sin_term = sin(kdotR);

      for(l=0;l<num_orbs-1;l++){
        ltab = l*num_orbs;
        for(m=l+1;m<num_orbs;m++){
          mtab = m*num_orbs;

          /* real part */
          which_Koverlap[ltab+m]+= cos_term*overlapR.mat[ltab+m];
          /* imaginary part */
          which_Koverlap[mtab+l] -= sin_term*overlapR.mat[mtab+l];
        }
        /* diagonal element (real) */
        which_Koverlap[ltab+l] += cos_term*overlapR.mat[ltab+l];
      }
      /* the last diagonal element */
      which_Koverlap[num_orbs*num_orbs-1] +=
        cos_term*overlapR.mat[num_orbs*num_orbs-1];
    }
  }
  if( cell->dim != 1 ){

    /******
      take care of the remainder of the layer containing the unit cell
      for each overlap the procedure is the same as above
      *******/
    for(i=0;i<=2*cell->overlaps[0];i++){
      itab = cell->overlaps[0]-i;
      for(j=1;j<=cell->overlaps[1];j++){
        R_space_overlap_matrix(cell,details,overlapR,num_orbs,tot_overlaps,
                               orbital_lookup_table,overlaps_so_far);
        overlaps_so_far++;

        for(which_k=0;which_k<details->num_KPOINTS;which_k++){
          kpoint = &(details->K_POINTS[which_k]);
          which_Koverlap = &(overlapK.mat[which_k*(num_orbs*num_orbs)]);
          kpointloc.x = TWOPI*kpoint->loc.x;
          kpointloc.y = TWOPI*kpoint->loc.y;
          kpointloc.z = TWOPI*kpoint->loc.z;
          kdotR = kpointloc.x*(real)itab+kpointloc.y*(real)j;
          cos_term = cos(kdotR);
          sin_term = sin(kdotR);

          for(l=0;l<num_orbs-1;l++){
            ltab = l*num_orbs;
            for(m=l+1;m<num_orbs;m++){
              mtab = m*num_orbs;

              which_Koverlap[ltab+m]+= cos_term*overlapR.mat[ltab+m];
              which_Koverlap[mtab+l] -= sin_term*overlapR.mat[mtab+l];
            }
            which_Koverlap[ltab+l] += cos_term*overlapR.mat[ltab+l];
          }
          which_Koverlap[num_orbs*num_orbs-1] +=
            cos_term*overlapR.mat[num_orbs*num_orbs-1];
        }
      }
    }
  }
  if( cell->dim == 3 ){
    /* do the layers above the one containing the unit cell */
    for(i=1;i<=cell->overlaps[2];i++){
      itab = i;
      for(j=0;j<=2*cell->overlaps[0];j++){
        jtab = cell->overlaps[0] - j;
        for(k=0;k<=2*cell->overlaps[1];k++){
          ktab = cell->overlaps[1] - k;

          R_space_overlap_matrix(cell,details,overlapR,num_orbs,tot_overlaps,
                                 orbital_lookup_table,overlaps_so_far);
          overlaps_so_far++;

          for(which_k=0;which_k<details->num_KPOINTS;which_k++){
            kpoint = &(details->K_POINTS[which_k]);
            which_Koverlap = &(overlapK.mat[which_k*(num_orbs*num_orbs)]);
            kpointloc.x = TWOPI*kpoint->loc.x;
            kpointloc.y = TWOPI*kpoint->loc.y;
            kpointloc.z = TWOPI*kpoint->loc.z;
            kdotR = kpointloc.x*(real)jtab+kpointloc.y*(real)ktab+
              kpointloc.z*(real)i;

            cos_term = cos(kdotR);
            sin_term = sin(kdotR);

            for(l=0;l<num_orbs-1;l++){
              ltab = l*num_orbs;
              for(m=l+1;m<num_orbs;m++){
                mtab = m*num_orbs;
                which_Koverlap[ltab+m]+= cos_term*overlapR.mat[ltab+m];
                which_Koverlap[mtab+l] -= sin_term*overlapR.mat[mtab+l];
              }
              which_Koverlap[ltab+l] += cos_term*overlapR.mat[ltab+l];
            }
            which_Koverlap[num_orbs*num_orbs-1] +=
              cos_term*overlapR.mat[num_orbs*num_orbs-1];
          }
        }
      }
    }
  }

#ifdef PRINTMAT
  fprintf(output_file,"---------- S(k) ------\n");
  printmat(overlapK.mat,num_orbs,num_orbs,output_file,1e-6,details->line_width);
#endif

  /* that's it! */
}
