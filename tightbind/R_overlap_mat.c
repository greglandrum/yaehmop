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
*     this file contains stuff for the overlap matrices in R space
*
*  created:  greg landrum  August 1993
*
*****************************************************************************/
/***
  Edit History:

  March '98: WG
  - f orbitals added [calc_R_overlap()]
***/

#include "bind.h"


/* some calculation savers! */
#define SQRT3 1.732050807568877
#define SQRT6 2.449489742783178
#define SQRT10 3.16227766016838
#define SQRT15 3.872983346207417
#define AUI 1.889644746



/******
  the sizes of the transformation matrices are increased a little bit
  to make the code more readable
*******/
#define P_SIZE 9+BEGIN_P
#define D_SIZE 25+BEGIN_D
#define F_SIZE 49+BEGIN_F


/****************************************************************************
*
*                   Procedure zero_overlaps
*
* Arguments: overlap: pointer to type real
*            details: pointer to detail type
*           num_orbs: int
*           num_atoms: int
*     doing_unit_cell: char
* orbital_lookup_table: pointer to int.
*
* Returns: none
*
* Action:  zeroes out the selected members of the overlap matrix 'overlap.
*
****************************************************************************/
void zero_overlaps(real *overlap,detail_type *details,int num_orbs,int num_atoms,
                   char doing_unit_cell,int *orbital_lookup_table)
{
  int i,j,k;
  int begin1,end1,begin2,end2;
  int contrib1,contrib2;

  if( !details->overlaps_off )
    FATAL_BUG("Bad overlaps_off passed to zero_overlaps.");

  for(i=0;i<details->num_overlaps_off;i++){
    if( (details->overlaps_off[i].inter_cell && !doing_unit_cell) ||
       (!details->overlaps_off[i].inter_cell && doing_unit_cell) ){
      contrib1 = details->overlaps_off[i].which1;
      contrib2 = details->overlaps_off[i].which2;
      if(details->overlaps_off[i].type == P_DOS_ORB){
        overlap[contrib1*num_orbs+contrib2] = 0.0;
        overlap[contrib2*num_orbs+contrib1] = 0.0;
      } else if(details->overlaps_off[i].type == P_DOS_ATOM){
        find_atoms_orbs(num_orbs,num_atoms,contrib1,orbital_lookup_table,
                        &begin1,&end1);
        find_atoms_orbs(num_orbs,num_atoms,contrib2,orbital_lookup_table,
                        &begin2,&end2);

        /* ignore dummy atoms */
        if(begin1 >=0 && begin2 >= 0 ){
          for(j=begin1;j<end1;j++){
            for(k=begin2;k<end2;k++){
              overlap[j*num_orbs+k] = 0.0;
              overlap[k*num_orbs+j] = 0.0;
            }
          }
        }
      }
    }
  }
}




/****************************************************************************
*
*                   Procedure calc_R_overlap
*
* Arguments: overlap: pointer to type real
*               cell: pointer to cell type
*            details: pointer to detail type
*           num_orbs: int
*          distances: point_type
*    doing_unit_cell: BOOLEAN
* orbital_lookup_table: pointer to int.
*
* Returns: none
*
* Action: Some day I'll figure out exactly what's going on here and
*   rewrite this comment, until then it's gonna have to go largely
*   uncommented.  Most of the code is snarfed directly from the original
*   (FORTRAN) version, except that I don't go out of my way to make it
*   efficient at the price of readability since this is not where the
*   program is gonna spend its time.
*
*
*  Hugh and I figured out some of what happens here. (but we didn't make the
*    code any more palatable).
*
*   Between each pair of atoms mov is called.  Mov evaluates
*     the overlap matrix elements between those atoms in a basis
*     of sigma, pi, and delta (with phi if those icky f orbitals
*     are being used).  The ugly transformation matrix crap that
*     is at the top of the double loop over all atoms constructs
*     the matrix to transform from this simple (spherical polar)
*     basis back into the cartesian (sensible, but hard to do the
*     math) basis.  I still have no idea how mov works and God
*     forbid that I should *ever* have to look at lovlap.
*
****************************************************************************/
void calc_R_overlap(real *overlap,cell_type *cell,detail_type *details,
                    int num_orbs,point_type distances,char doing_unit_cell,
                    int *orbital_lookup_table)
{
  point_type dist_vect;
  real tot_dist,xy_dist,cosA,sinA,cosB,sinB;
  real s2B,c2B,c2A,sA3,c3B,s3B;
  real p_trans_mat[P_SIZE],d_trans_mat[D_SIZE],f_trans_mat[F_SIZE];
  real sigma,pi,delta,phi;
  real temp;

  int i,j,k,j_end;
  int i_tab,j_tab;
  int i_orb,j_orb;
  int la,lb;
  int q_num1,q_num2;

  /* zero out the transformation matrices just to make sure */
  bzero(p_trans_mat,(P_SIZE)*sizeof(real));
  bzero(d_trans_mat,(D_SIZE)*sizeof(real));
  bzero(f_trans_mat,(F_SIZE)*sizeof(real));
  bzero(overlap,num_orbs*num_orbs*sizeof(real));

  /*
    printf("add: %6.4lf %6.4lf %6.4lf\n",distances.x,distances.y,distances.z);
    printf("MOVLAP\n");
    */
  j_end = cell->num_atoms;
  for(i=0;i<cell->num_atoms;i++){

    /******
      check to see if we are evaluating the overlap matrix
      within the unit cell. If so, only some of the elements need
      to be evaluated.
      ******/
    if(doing_unit_cell) j_end = i;

    /*******
      since i refers to the other atoms in the unit cell, a tab
      needs to be kept in order to actually find where the
      orbitals for the i'th atom begin
      ********/

    i_tab = orbital_lookup_table[i];

    /* trap dummy atoms */
    if(i_tab >= 0 ){
      for(j=0;j<j_end;j++){
        j_tab = orbital_lookup_table[j];
        if(j_tab >= 0){
          /* this is the distance vector between the two atoms */
          dist_vect.x = cell->atoms[i].loc.x - cell->atoms[j].loc.x + distances.x;
          dist_vect.y = cell->atoms[i].loc.y - cell->atoms[j].loc.y + distances.y;
          dist_vect.z = cell->atoms[i].loc.z - cell->atoms[j].loc.z + distances.z;
          /*
            fprintf(stderr,"dist: %6.4lf %6.4lf %6.4lf\n",dist_vect.x,dist_vect.y,dist_vect.z);
            */
          temp = dist_vect.x*dist_vect.x + dist_vect.y*dist_vect.y;

          /* distance in the xy plane */
          xy_dist = sqrt(temp);
          /* total distance */
          tot_dist = sqrt(temp + dist_vect.z*dist_vect.z);


          /************

            if the total distance is equal to zero, then set the overlap
            term to zero.  This allows "fake" atoms to be placed on top
            of real atoms to expand the basis set without causing problems

                    OR

            if the total distance is greater than the cut off, then
            we can skip evaluating this puppy...


            throughout this we make the implicit assumption that the
            overlap matrix has been initially zeroed.

          *********/
          if( tot_dist >= 1e-6 && tot_dist <= details->rho ){

            /********

              set up the cosines and sines of the angles A & B

              A is the angle between the distance vector and the z axis
              B is the angle between the xy projection of the distance
              vector and the x axis

              *********/
            if( xy_dist < 1e-5 ){
              cosB = 1.0;
              sinB = 0.0;
              sinA = 0.0;
            }
            else{
              cosB = dist_vect.x/xy_dist;
              sinB = dist_vect.y/xy_dist;
              sinA = xy_dist/tot_dist;
            }
            cosA = dist_vect.z/tot_dist;

            /*******
              build the p projection matrix
              (order:  x,y,z)
              *******/
            p_trans_mat[0+BEGIN_P] = sinA*cosB;
            p_trans_mat[1+BEGIN_P] =sinA*sinB;
            p_trans_mat[2+BEGIN_P] =cosA;
            p_trans_mat[3+BEGIN_P] =cosA*cosB;
            p_trans_mat[4+BEGIN_P] =cosA*sinB;
            p_trans_mat[5+BEGIN_P] =-sinA;
            p_trans_mat[6+BEGIN_P] =-sinB;
            p_trans_mat[7+BEGIN_P] =cosB;
            p_trans_mat[8+BEGIN_P] =0.0;

            /*******
              build the d projection matrix
              (order:  x2-y2,z2,xy,xz,yz)
              *******/
            if( cell->atoms[i].nd || cell->atoms[j].nd ){

              /* some useful definitions */
              c2A = (cosA*cosA)-(sinA*sinA);
              c2B = (cosB*cosB)-(sinB*sinB);
              s2B = 2*sinB*cosB;


              d_trans_mat[0+BEGIN_D] =SQRT3*.5*sinA*sinA*c2B;
              d_trans_mat[1+BEGIN_D] =1.0-1.5*sinA*sinA;
              d_trans_mat[2+BEGIN_D] =SQRT3*cosB*sinB*sinA*sinA;
              d_trans_mat[3+BEGIN_D] =SQRT3*cosA*sinA*cosB;
              d_trans_mat[4+BEGIN_D] =SQRT3*cosA*sinA*sinB;
              d_trans_mat[5+BEGIN_D] =cosA*sinA*c2B;
              d_trans_mat[6+BEGIN_D] =-SQRT3*cosA*sinA;
              d_trans_mat[7+BEGIN_D] =cosA*sinA*s2B;
              d_trans_mat[8+BEGIN_D] =cosB*c2A;
              d_trans_mat[9+BEGIN_D] =sinB*c2A;
              d_trans_mat[10+BEGIN_D] =-sinA*s2B;
              d_trans_mat[11+BEGIN_D] =0.0;
              d_trans_mat[12+BEGIN_D] =sinA*c2B;
              d_trans_mat[13+BEGIN_D] =-p_trans_mat[4+BEGIN_P];
              d_trans_mat[14+BEGIN_D] =p_trans_mat[3+BEGIN_P];
           }

            /* only do these if both atoms have d or f orbitals */

            if((cell->atoms[i].nd || cell->atoms[i].nf) && (cell->atoms[j].nd || cell->atoms[j].nf)){
              d_trans_mat[15+BEGIN_D] =.50*(1.0+cosA*cosA)*c2B;
              d_trans_mat[16+BEGIN_D] =.50*SQRT3*sinA*sinA;
              d_trans_mat[17+BEGIN_D] =cosB*sinB*(1.0+cosA*cosA);
              d_trans_mat[18+BEGIN_D] =-cosA*sinA*cosB;
              d_trans_mat[19+BEGIN_D] =-cosA*sinA*sinB;
              d_trans_mat[20+BEGIN_D] =-cosA*s2B;
              d_trans_mat[21+BEGIN_D] =0.0;
              d_trans_mat[22+BEGIN_D] =cosA*c2B;
              d_trans_mat[23+BEGIN_D] =p_trans_mat[1+BEGIN_P];
              d_trans_mat[24+BEGIN_D] =-p_trans_mat[0+BEGIN_P];
            }
            /*******
              build the f projection matrix
              *******/

            if( cell->atoms[i].nf || cell->atoms[j].nf ){

              /* some useful definitions */
              sA3 = sinA*sinA*sinA;
              c3B = (c2B*cosB)-(s2B*sinB);
              s3B = (c2B*sinB)+(s2B*cosB);

              f_trans_mat[0+BEGIN_F] =0.5*cosA*(5*cosA*cosA-3);
              f_trans_mat[1+BEGIN_F] =SQRT6*0.25*cosB*sinA*(5*cosA*cosA-1);
              f_trans_mat[2+BEGIN_F] =SQRT6*0.25*sinB*sinA*(5*cosA*cosA-1);
              f_trans_mat[3+BEGIN_F] =SQRT15*sinB*cosB*cosA*sinA*sinA;
              f_trans_mat[4+BEGIN_F] =SQRT15*0.5*c2B*cosA*sinA*sinA;
              f_trans_mat[5+BEGIN_F] =SQRT10*0.25*c3B*sA3;
              f_trans_mat[6+BEGIN_F] =SQRT10*0.25*s3B*sA3;
              f_trans_mat[7+BEGIN_F] =(-SQRT6)*0.25*sinA*(5*cosA*cosA-1);
              f_trans_mat[8+BEGIN_F] =0.25*cosB*cosA*(15*cosA*cosA-11);
              f_trans_mat[9+BEGIN_F] =0.25*sinB*cosA*(15*cosA*cosA-11);
              f_trans_mat[10+BEGIN_F] =SQRT10*0.25*s2B*sinA*(3*cosA*cosA-1);
              f_trans_mat[11+BEGIN_F] =SQRT10*0.25*c2B*sinA*(3*cosA*cosA-1);
              f_trans_mat[12+BEGIN_F] =SQRT15*0.25*c3B*cosA*sinA*sinA;
              f_trans_mat[13+BEGIN_F] =SQRT15*0.25*s3B*cosA*sinA*sinA;
              f_trans_mat[14+BEGIN_F] =0.0;
              f_trans_mat[15+BEGIN_F] =(-0.25)*sinB*(5*cosA*cosA-1);
              f_trans_mat[16+BEGIN_F] =0.25*cosB*(5*cosA*cosA-1);
              f_trans_mat[17+BEGIN_F] =SQRT10*0.5*c2B*sinA*cosA;
              f_trans_mat[18+BEGIN_F] =(-SQRT10)*sinA*cosA*sinB*cosB;
              f_trans_mat[19+BEGIN_F] =(-SQRT15)*0.25*s3B*sinA*sinA;
              f_trans_mat[20+BEGIN_F] =SQRT15*0.25*sinA*sinA*cosB*(4*cosB*cosB-3);
            }

            /*   only do these if both atoms have d or f orbitals */

            if( (cell->atoms[i].nd || cell->atoms[i].nf) && (cell->atoms[j].nd || cell->atoms[j].nf)){
              f_trans_mat[21+BEGIN_F] =0.0;
              f_trans_mat[22+BEGIN_F] =SQRT10*0.5*sinB*cosA*sinA;
              f_trans_mat[23+BEGIN_F] =(-SQRT10)*0.5*cosB*cosA*sinA;
              f_trans_mat[24+BEGIN_F] =c2B*(cosA*cosA-sinA*sinA);
              f_trans_mat[25+BEGIN_F] =(-s2B)*(cosA*cosA-sinA*sinA);
              f_trans_mat[26+BEGIN_F] =(-SQRT6)*0.5*s3B*sinA*cosA;
              f_trans_mat[27+BEGIN_F] =SQRT6*0.5*cosA*sinA*cosB*(4*cosB*cosB-3);
              f_trans_mat[28+BEGIN_F] =SQRT15*0.5*cosA*sinA*sinA;
              f_trans_mat[29+BEGIN_F] =SQRT10*0.25*cosB*sinA*(1-3*cosA*cosA);
              f_trans_mat[30+BEGIN_F] =SQRT10*0.25*sinA*sinB*(1-3*cosA*cosA);
              f_trans_mat[31+BEGIN_F] =sinB*cosB*cosA*(3*cosA*cosA-1);
              f_trans_mat[32+BEGIN_F] =0.5*c2B*cosA*(3*cosA*cosA-1);
              f_trans_mat[33+BEGIN_F] =SQRT6*0.25*c3B*sinA*(1+cosA*cosA);
              f_trans_mat[34+BEGIN_F] =SQRT6*0.25*s3B*sinA*(1+cosA*cosA);
            }

            /*   only do these if both atoms have f orbitals */

            if( cell->atoms[i].nf && cell->atoms[j].nf){
              f_trans_mat[35+BEGIN_F] =SQRT10*0.25*sinA*(cosA*cosA-1);
              f_trans_mat[36+BEGIN_F] =SQRT15*0.25*cosB*cosA*sinA*sinA;
              f_trans_mat[37+BEGIN_F] =SQRT15*0.25*sinB*cosA*sinA*sinA;
              f_trans_mat[38+BEGIN_F] =(-SQRT6)*0.5*sinB*cosB*sinA*(1+cosA*cosA);
              f_trans_mat[39+BEGIN_F] =(-SQRT6)*0.25*c2B*sinA*(1+cosA*cosA);
              f_trans_mat[40+BEGIN_F] =0.25*c3B*cosA*(3+cosA*cosA);
              f_trans_mat[41+BEGIN_F] =0.25*s3B*cosA*(3+cosA*cosA);
              f_trans_mat[42+BEGIN_F] =0.0;
              f_trans_mat[43+BEGIN_F] =(-SQRT15)*0.25*sinB*sinA*sinA;
              f_trans_mat[44+BEGIN_F] =SQRT15*0.25*cosB*sinA*sinA;
              f_trans_mat[45+BEGIN_F] =(-SQRT6)*0.5*c2B*sinA*cosA;
              f_trans_mat[46+BEGIN_F] =SQRT6*sinB*cosB*sinA*cosA;
              f_trans_mat[47+BEGIN_F] =(-0.25)*s3B*(1+3*cosA*cosA);
              f_trans_mat[48+BEGIN_F] =0.25*c3B*(1+3*cosA*cosA);
            }
            /* AUI is 1/BOHR where BOHR is the Bohr radius */
            tot_dist *= AUI;

            /*************

              Now actually evaluate the overlaps

              **************/

            /*-----------
              < S(i) | S(j) >
              -------------*/
            q_num1 = cell->atoms[i].ns;
            q_num2 = cell->atoms[j].ns;
            if( q_num1 && q_num2 ){
              la = 0;
              lb = 0;

              mov(&sigma,&pi,&delta,&phi,i,j,tot_dist,q_num1,q_num2,
                  la,lb,cell->atoms);

              /* fprintf(stderr,"S-S: (sigma) %.12f\n",sigma);  */

              overlap[i_tab*num_orbs+j_tab] = sigma;
            }

            /*-----------
              < P(i) | S(j) >
              -------------*/
            q_num1 = cell->atoms[i].np;
            if( q_num1 && q_num2 ){
              la = 1;
              lb = 0;

              mov(&sigma,&pi,&delta,&phi,i,j,tot_dist,q_num1,q_num2,
                  la,lb,cell->atoms);

              sigma *= -1;

              /* fprintf(stderr,"P-S: (sigma) %.12f\n",sigma);  */

              for(j_orb=BEGIN_P;j_orb<=END_P;j_orb++){
                overlap[(i_tab+j_orb)*num_orbs+j_tab] = p_trans_mat[j_orb]*sigma;
              }
            }

            /*-----------
              < P(i) | P(j) >
              -------------*/
            q_num2 = cell->atoms[j].np;
            if( q_num1 && q_num2 ){
              la = 1;
              lb = 1;

              mov(&sigma,&pi,&delta,&phi,i,j,tot_dist,q_num1,q_num2,
                  la,lb,cell->atoms);

              sigma *= -1;

              /* fprintf(stderr,"P-P: (sigma,pi) %.12f %.12f\n",sigma,pi);  */

              for(j_orb=BEGIN_P;j_orb<=END_P;j_orb++){
                for(i_orb=j_orb;i_orb<=END_P;i_orb++){
                  overlap[(i_tab+i_orb)*num_orbs+j_tab+j_orb] =
                    p_trans_mat[j_orb]*p_trans_mat[i_orb]*sigma +
                      (p_trans_mat[i_orb+3]*p_trans_mat[j_orb+3] +
                       p_trans_mat[i_orb+6]*p_trans_mat[j_orb+6]) * pi;
                  overlap[(i_tab+j_orb)*num_orbs+j_tab+i_orb] =
                    overlap[(i_tab+i_orb)*num_orbs+j_tab+j_orb];

/*                  fprintf(stderr,"P-P: orbital %d and %d\t%f \n",j_orb+1,i_orb+1,overlap[(i_tab+j_orb)*num_orbs + j_tab + i_orb]);     */
                }
              }
            }

            /*-----------
              < S(i) | P(j) >
              -------------*/
            q_num1 = cell->atoms[i].ns;
            if( q_num1 && q_num2 ){
              la = 0;
              lb = 1;

              mov(&sigma,&pi,&delta,&phi,i,j,tot_dist,q_num1,q_num2,
                  la,lb,cell->atoms);

              /* fprintf(stderr,"S-P: (sigma) %.12f\n",sigma);  */

              for(j_orb=BEGIN_P;j_orb<=END_P;j_orb++){
                overlap[i_tab*num_orbs+j_tab+j_orb] = p_trans_mat[j_orb]*sigma;
              }
            }

            /*-----------
              < S(i) | D(j) >
              -------------*/
            q_num2 = cell->atoms[j].nd;
            if( q_num1 && q_num2 ){
              la = 0;
              lb = 2;

              mov(&sigma,&pi,&delta,&phi,i,j,tot_dist,q_num1,q_num2,
                  la,lb,cell->atoms);

              /* fprintf(stderr,"S-D: (sigma) %.12f\n",sigma);   */

              for(j_orb=BEGIN_D;j_orb<=END_D;j_orb++){
                overlap[i_tab*num_orbs+j_tab+j_orb] = d_trans_mat[j_orb]*sigma;
              }
            }

            /*-----------
              < P(i) | D(j) >
              -------------*/
            q_num1 = cell->atoms[i].np;
            if( q_num1 && q_num2 ){
              la = 1;
              lb = 2;

              mov(&sigma,&pi,&delta,&phi,i,j,tot_dist,q_num1,q_num2,
                  la,lb,cell->atoms);

              sigma *= -1;

              /* fprintf(stderr,"P-D: (sigma,pi) %.12f %.12f\n",sigma,pi);   */

              for(i_orb=BEGIN_P;i_orb<=END_P;i_orb++){
                for(j_orb=BEGIN_D;j_orb<=END_D;j_orb++){
                  overlap[(i_tab+i_orb)*num_orbs+j_tab+j_orb] =
                    p_trans_mat[i_orb]*d_trans_mat[j_orb]*sigma +
                      (d_trans_mat[j_orb+5]*p_trans_mat[i_orb+3] +
                       d_trans_mat[j_orb+10]*p_trans_mat[i_orb+6]) * pi;
                }
              }
            }

            /*-----------
              < D(i) | S(j) >
              -------------*/
            q_num1 = cell->atoms[i].nd;
            q_num2 = cell->atoms[j].ns;
            if( q_num1 && q_num2 ){
              la = 2;
              lb = 0;

              mov(&sigma,&pi,&delta,&phi,i,j,tot_dist,q_num1,q_num2,
                  la,lb,cell->atoms);

              /* fprintf(stderr,"D-S: (sigma) %.12f\n",sigma);   */

              for(j_orb=BEGIN_D;j_orb<=END_D;j_orb++){
                overlap[(i_tab+j_orb)*num_orbs+j_tab] = d_trans_mat[j_orb]*sigma;
              }
            }

            /*-----------
              < D(i) | P(j) >
              -------------*/
            q_num2 = cell->atoms[j].np;
            if( q_num1 && q_num2 ){
              la = 2;
              lb = 1;

              mov(&sigma,&pi,&delta,&phi,i,j,tot_dist,q_num1,q_num2,
                  la,lb,cell->atoms);

              pi *= -1;

              /* fprintf(stderr,"D-P: (sigma,pi) %.12f %.12f\n",sigma,pi);   */

              for(i_orb=BEGIN_P;i_orb<=END_P;i_orb++){
                for(j_orb=BEGIN_D;j_orb<=END_D;j_orb++){
                  overlap[(i_tab+j_orb)*num_orbs+j_tab+i_orb] =
                    p_trans_mat[i_orb]*d_trans_mat[j_orb]*sigma +
                      (p_trans_mat[i_orb+3]*d_trans_mat[j_orb+5] +
                       d_trans_mat[j_orb+10]*p_trans_mat[i_orb+6]) * pi;
                }
              }
            }

            /*-----------
              < D(i) | D(j) >
              -------------*/
            q_num2 = cell->atoms[j].nd;
            if( q_num1 && q_num2 ){
              la = 2;
              lb = 2;

              mov(&sigma,&pi,&delta,&phi,i,j,tot_dist,q_num1,q_num2,
                  la,lb,cell->atoms);

              pi *= -1;

              /* fprintf(stderr,"D-D: (sigma,pi,delta) %.12f %.12f %.12f\n",
                      sigma,pi,delta);  */

              for(i_orb=BEGIN_D;i_orb<=END_D;i_orb++){
                for(j_orb=BEGIN_D;j_orb<=END_D;j_orb++){
                  overlap[(i_tab+j_orb)*num_orbs+(j_tab+i_orb)] =
                    d_trans_mat[i_orb]*d_trans_mat[j_orb]*sigma +
                      (d_trans_mat[i_orb+5]*d_trans_mat[j_orb+5] +
                       d_trans_mat[j_orb+10]*d_trans_mat[i_orb+10])*pi+
                         (d_trans_mat[i_orb+15]*d_trans_mat[j_orb+15] +
                          d_trans_mat[i_orb+20]*d_trans_mat[j_orb+20])*delta;

                  /*
                    fprintf(stderr,"D-D element(%d,%d)= %f\n",
                    i_orb,j_orb,overlap[(i_tab+j_orb)*num_orbs+j_tab+i_orb]);
                    */

                  overlap[(i_tab+i_orb)*num_orbs+j_tab+j_orb] =
                    overlap[(i_tab+j_orb)*num_orbs+j_tab+i_orb];
                }
              }
            }

            /*------------
              < S(i) | F(j) >
              -------------*/
            q_num1 = cell->atoms[i].ns;
            q_num2 = cell->atoms[j].nf;
            if(q_num1 && q_num2){
              la = 0;
              lb = 3;

              mov(&sigma,&pi,&delta,&phi,i,j,tot_dist,
                  q_num1,q_num2,la,lb,cell->atoms);

              /* fprintf(stderr, "S-F:  (sigma) %f \n",sigma);  */

              for(j_orb=BEGIN_F;j_orb<=END_F;j_orb++){
                overlap[(i_tab*num_orbs)+j_tab+j_orb]=
                  f_trans_mat[j_orb]*sigma;

                /*
                  fprintf(stderr,"S-F element(%d)= %f\n",
                  j_orb,overlap[(i_tab*num_orbs)+j_tab+j_orb]);
                  */

              }
            }

            /*------------
              < P(i) | F(j) >
              -------------*/
            q_num1 = cell->atoms[i].np;
            q_num2 = cell->atoms[j].nf;
            if(q_num1 && q_num2){
              la = 1;
              lb = 3;

              mov(&sigma,&pi,&delta,&phi,i,j,tot_dist,
                  q_num1,q_num2,la,lb,cell->atoms);

              sigma *= -1;

              /* fprintf(stderr, "P-F:  (sigma,pi) %f %f \n",
                      sigma, pi);  */

              for(i_orb=BEGIN_P; i_orb<=END_P; i_orb++){
                for(j_orb=BEGIN_F; j_orb<=END_F; j_orb++){
                  overlap[(i_tab+i_orb)*num_orbs + j_tab + j_orb]=
                    p_trans_mat[i_orb]*f_trans_mat[j_orb]*sigma +
                      (p_trans_mat[i_orb+3]*f_trans_mat[j_orb+7] +
                       p_trans_mat[i_orb+6]*f_trans_mat[j_orb+14])*pi;

                  /*
                    fprintf(stderr,"P-F element(%d,%d)= %f\n",
                    i_orb,j_orb,overlap[(i_tab+i_orb)*num_orbs+j_tab+j_orb]);
                    */

                }
              }
            }

            /*------------
              < D(i) | F(j) >
              ------------*/
            q_num1 = cell->atoms[i].nd;
            q_num2 = cell->atoms[j].nf;
            if(q_num1 && q_num2){
              la = 2;
              lb = 3;

              mov(&sigma,&pi,&delta,&phi,i,j,tot_dist,
                  q_num1,q_num2,la,lb,cell->atoms);

              pi *= -1;

              /* fprintf(stderr, "D-F:  (sigma,pi,delta) %f %f %f \n",
                      sigma,pi,delta);   */

              for(i_orb=BEGIN_D;i_orb<=END_D;i_orb++){
                for(j_orb=BEGIN_F;j_orb<=END_F;j_orb++){
                  overlap[(i_tab+i_orb)*num_orbs + j_tab + j_orb]=
                    d_trans_mat[i_orb]*f_trans_mat[j_orb]*sigma +
                      (d_trans_mat[i_orb+5]*f_trans_mat[j_orb+7] +
                       d_trans_mat[i_orb+10]*f_trans_mat[j_orb+14])*pi +
                         (d_trans_mat[i_orb+15]*f_trans_mat[j_orb+28] +
                          d_trans_mat[i_orb+20]*f_trans_mat[j_orb+21])*delta;

                  /*
                    fprintf(stderr,"D-F element(%d,%d)= %f\n",
                    i_orb,j_orb,overlap[(i_tab+i_orb)*num_orbs+j_tab+j_orb]);
                    */

                }
              }
            }

            /*------------
              < F(i) | F(j) >
              ------------*/
            q_num1 = cell->atoms[i].nf;
            q_num2 = cell->atoms[j].nf;
            if(q_num1 && q_num2){
              la = 3;
              lb = 3;

              mov(&sigma,&pi,&delta,&phi,i,j,tot_dist,
                  q_num1,q_num2,la,lb,cell->atoms);

              sigma *= -1;

              /* fprintf(stderr, "F-F:  (sigma,pi,delta,phi) %f %f %f %f \n",
                      sigma,pi,delta,phi);   */

              for(i_orb=BEGIN_F;i_orb<=END_F; i_orb++){
                for(j_orb=BEGIN_F;j_orb<=END_F; j_orb++){
                  overlap[(i_tab+j_orb)*num_orbs + j_tab + i_orb]=
                    f_trans_mat[i_orb]*f_trans_mat[j_orb]*sigma +
                      (f_trans_mat[i_orb+7]*f_trans_mat[j_orb+7] +
                       f_trans_mat[i_orb+14]*f_trans_mat[j_orb+14])*pi +
                         (f_trans_mat[i_orb+21]*f_trans_mat[j_orb+21] +
                          f_trans_mat[i_orb+28]*f_trans_mat[j_orb+28])*delta +
                            (f_trans_mat[i_orb+35]*f_trans_mat[j_orb+35]+
                             f_trans_mat[i_orb+42]*f_trans_mat[j_orb+42])*phi;

                  /*
                    fprintf(stderr,"F-F element(%d,%d)= %f\n",
                    i_orb,j_orb,overlap[(i_tab+j_orb)*num_orbs
                    +j_tab+i_orb]);
                    */
                  overlap[(i_tab+i_orb)*num_orbs + j_tab + j_orb]=
                    overlap[(i_tab+j_orb)*num_orbs + j_tab + i_orb];
                }
              }
            }

            /*------------
              < F(i) | S(j) >
              -----------*/
            q_num1 = cell->atoms[i].nf;
            q_num2 = cell->atoms[j].ns;
            if(q_num1 && q_num2){
              la = 3;
              lb = 0;

              mov(&sigma,&pi,&delta,&phi,i,j,tot_dist,
                  q_num1,q_num2,la,lb,cell->atoms);

              sigma *= -1;

              /* fprintf(stderr, "F-S: (sigma) %f \n",sigma);  */

              for(i_orb=BEGIN_F;i_orb<=END_F;i_orb++){
                overlap[(i_tab+i_orb)*num_orbs + j_tab]=
                  f_trans_mat[i_orb]*sigma;
              }
            }

            /*------------
              < F(i) | P(j) >
              ------------*/
            q_num1 = cell->atoms[i].nf;
            q_num2 = cell->atoms[j].np;
            if(q_num1 && q_num2){
              la = 3;
              lb = 1;

              mov(&sigma,&pi,&delta,&phi,i,j,tot_dist,
                  q_num1,q_num2,la,lb,cell->atoms);

              sigma *= -1;

              /* fprintf(stderr,"F-P: (sigma,pi) %f %f \n",sigma,pi);  */

              for(i_orb=BEGIN_F;i_orb<=END_F;i_orb++){
                for(j_orb=BEGIN_P;j_orb<=END_P;j_orb++){
                  overlap[(i_tab+i_orb)*num_orbs + j_tab + j_orb]=
                    f_trans_mat[i_orb]*p_trans_mat[j_orb]*sigma +
                      (f_trans_mat[i_orb+7]*p_trans_mat[j_orb+3] +
                       f_trans_mat[i_orb+14]*p_trans_mat[j_orb+6])*pi;
                }
              }
            }

            /*-------------
              < F(i) | D(j) >
              ------------*/
            q_num1 = cell->atoms[i].nf;
            q_num2 = cell->atoms[j].nd;
            if(q_num1 && q_num2){
              la = 3;
              lb = 2;

              mov(&sigma,&pi,&delta,&phi,i,j,tot_dist,
                  q_num1,q_num2,la,lb,cell->atoms);

              sigma *= -1; delta *= -1;

              /* fprintf(stderr, "F-D: (sigma,pi,delta) %f %f %f \n",sigma,pi,delta); */

              for(i_orb=BEGIN_F;i_orb<=END_F;i_orb++){
                for(j_orb=BEGIN_D;j_orb<=END_D;j_orb++){
                  overlap[(i_tab+i_orb)*num_orbs + j_tab + j_orb]=
                    f_trans_mat[i_orb]*d_trans_mat[j_orb]*sigma +
                      (f_trans_mat[i_orb+7]*d_trans_mat[j_orb+5] +
                       f_trans_mat[i_orb+14]*d_trans_mat[j_orb+10])*pi +
                         (f_trans_mat[i_orb+21]*d_trans_mat[j_orb+20] +
                          f_trans_mat[i_orb+28]*d_trans_mat[j_orb+15])*delta;
                }
              }
            }
            /* END OF MATRIX ELEMENTS */
          }
        }
      }
    }
  }
#ifdef PRINTMAT
fprintf(status_file,"---- S(R) (%lf %lf %lf) rho: %lf----\n",
        distances.x,distances.y,distances.z,details->rho);
printmat(overlap,num_orbs,num_orbs,status_file,1e-10,0,details->line_width);
#endif

  /* do the symmetrization thing */
  if( !doing_unit_cell ){
    for(i=0;i<num_orbs;i++){
      if(i<num_orbs-1){
        for(j=i+1;j<num_orbs;j++){
          overlap[j*num_orbs+i] -= overlap[i*num_orbs+j];
          overlap[i*num_orbs+j] = overlap[j*num_orbs+i]+
            2*overlap[i*num_orbs+j];
        }
      }
      overlap[i*num_orbs+i]*=2;
    }
  }

  /* if we are zero'ing out any overlaps, do so now */
  if( details->num_overlaps_off && details->overlaps_off ){
    zero_overlaps(overlap,details,num_orbs,cell->num_atoms,doing_unit_cell,
                  orbital_lookup_table);
  }


//fprintf(stdout,".\n");


#if 0
#ifdef PRINTMAT
fprintf(status_file,"---- S(R) (%lf %lf %lf) rho: %lf----\n",
        distances.x,distances.y,distances.z,details->rho);
printmat(overlap,num_orbs,num_orbs,status_file,1e-6,0,details->line_width);
#endif
#endif
}


/****************************************************************************
 *
 *                   Procedure R_space_overlap_matrix
 *
 * Arguments:  cell: pointer to cell type
 *          details: pointer to detail type
 *          overlap: hermetian_matrix_type
 *         num_orbs: int
 *     tot_overlaps: int
 * orbital_lookup_table: pointer to int.
 *
 * Returns: none
 *
 * Action:
 *     generates the total overlap matrix in real space
 *
 *  Below is an example of what happens for one overlap in every direction
 *   *'ed cells are visited.  The unit cell contains an X (and is visited)
 *
 *
 *    /-----------\       /-----------\       /-----------\
 *    |   |   |   |       | * | * | * |       | * | * | * |
 *    |---|---|---|       |---|---|---|       |---|---|---|
 *    |   |   |   |       |   | X | * |       | * | * | * |
 *    |---|---|---|       |---|---|---|       |---|---|---|    ^ a2
 *    |   |   |   |       |   |   |   |       | * | * | * |    |
 *    \-----------/       \-----------/       \-----------/    |
 *                                                                            -->a1
 *              Layer -1                   Layer 0                     Layer 1
 *
 *
 *****************************************************************************/
void R_space_overlap_matrix(cell,details,overlap,num_orbs,tot_overlaps,
                            orbital_lookup_table,which_one)
  cell_type *cell;
  detail_type *details;
  hermetian_matrix_type overlap;
  int num_orbs,tot_overlaps;
  int *orbital_lookup_table;
  int which_one;
{
  char err_string[240];
  int overlaps_so_far;
  int i,j,k;
  int itab,jtab,ktab;
  int overlap_tab;

  char found = 0;
  point_type cell_dim[3];
  point_type distances;
  real temp,min=100.0;
  int min_dir;

  /* if this is the first call for this cycle, find the dimensions */
  if( cell->dim > 0 && which_one==0){
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
  if( fabs(details->rho) <= 1e-3 )
    details->rho = 10.0;

  /* initialize the overlap matrix to zeroes (just in case) */
  if( details->store_R_overlaps ){
    for(i=0;i<tot_overlaps;i++){
      itab = i*num_orbs*num_orbs;
      for(j=0;j<num_orbs;j++){
        jtab = j*num_orbs;
        for(k=0;k<num_orbs;k++){
          overlap.mat[itab+jtab+k] = 0.0;
        }
      }
    }
  }
  overlaps_so_far = 0;
  overlap_tab = 0;

  /******

    first do the unit cell

  ******/
  distances.x=distances.y=distances.z=0.0;
  if( details->store_R_overlaps || overlaps_so_far == which_one ){
    calc_R_overlap(&(overlap.mat[overlap_tab]),cell,details,
                   num_orbs,distances,TRUE,orbital_lookup_table);
    found = 1;
    /* put 1's on the diagonal and copy the elements across the diagonal */
    for(j=0;j<num_orbs;j++){
      jtab = j*num_orbs;
      overlap.mat[overlap_tab+jtab+j] = 1.0;
      for(k=0;k<j;k++){
        ktab = k*num_orbs;
        overlap.mat[overlap_tab+ktab+j]=overlap.mat[overlap_tab+jtab+k];
        overlap.mat[overlap_tab+jtab+k] =0.0;
      }
    }
  }
  overlaps_so_far++;
  if( details->store_R_overlaps ) overlap_tab += num_orbs*num_orbs;

  /*******

    now move into the other cells.

  ********/
  if( cell-> dim > 0){
    for(i=1;i<=cell->overlaps[0];i++){
      if( details->store_R_overlaps || overlaps_so_far == which_one ){
        distances.x = i*cell_dim[0].x;
        distances.y = i*cell_dim[0].y;
        distances.z = i*cell_dim[0].z;
        calc_R_overlap(&(overlap.mat[overlap_tab]),cell,details,
                       num_orbs,distances,FALSE,orbital_lookup_table);
        found = 1;
      }
      overlaps_so_far++;
      if( details->store_R_overlaps ) overlap_tab += num_orbs*num_orbs;
    }
  }

  /* if the crystal is 1-D then we're done. */
  if( cell->dim > 1 ){

    /* take care of the remainder of the layer containing the unit cell */
    for(i=0;i<=2*cell->overlaps[0];i++){
      itab = cell->overlaps[0] - i;
      for(j=1;j<=cell->overlaps[1];j++){
        if( details->store_R_overlaps || overlaps_so_far == which_one ){
          distances.x = itab*cell_dim[0].x + j*cell_dim[1].x;
          distances.y = itab*cell_dim[0].y + j*cell_dim[1].y;
          distances.z = itab*cell_dim[0].z + j*cell_dim[1].z;

          calc_R_overlap(&(overlap.mat[overlap_tab]),cell,details,
                         num_orbs,distances,FALSE,orbital_lookup_table);
          found = 1;
        }
        overlaps_so_far++;
        if( details->store_R_overlaps ) overlap_tab += num_orbs*num_orbs;
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
          if( details->store_R_overlaps || overlaps_so_far == which_one ){
            distances.x = itab*cell_dim[2].x + jtab*cell_dim[0].x +
              ktab*cell_dim[1].x;
            distances.y = itab*cell_dim[2].y + jtab*cell_dim[0].y +
              ktab*cell_dim[1].y;
            distances.z = itab*cell_dim[2].z + jtab*cell_dim[0].z +
              ktab*cell_dim[1].z;

            calc_R_overlap(&(overlap.mat[overlap_tab]),cell,details,
                           num_orbs,distances,FALSE,orbital_lookup_table);
            found = 1;
          }
          overlaps_so_far++;
          if( details->store_R_overlaps ) overlap_tab += num_orbs*num_orbs;
        }
      }
    }
  }

  /* check to see if enough overlaps were done */
  if( details->store_R_overlaps && overlaps_so_far != tot_overlaps ){
    sprintf(err_string,"*** Warning *** Number of overlaps done: %d does \
not match calculated number: %d.  This is a bug.\n",overlaps_so_far,tot_overlaps);
    FATAL_BUG(err_string);
  }
  if( !details->store_R_overlaps && !found ){
    sprintf(err_string,"Something ain't right! so_far: %d tot: %d which:%d\n",
            overlaps_so_far,tot_overlaps,which_one);
    NONFATAL_BUG(err_string);
  }

  if( details->store_R_overlaps )
    fprintf(status_file,"\n\nDone evaluating overlap integrals.\n");
}







