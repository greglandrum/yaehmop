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

#include "viewkel.h"


#define AUI 1.889644746

/* calculation savers */
float sqrt3=1.7320508075688772;
float sqrt5=2.2360679774997898;
float sqrt7=2.6457513110645907;
float sqrt15=3.8729833462074170;
float sqrt42=6.4807406984078604;
float sqrt70=8.3666002653407556;
float sqrt105=10.2469507659595980;

float fourpi=4.0*M_PI;

extern double drand48();

/************************************

  This has got the routines for evaluating the values
  and derivatives of wavefunctions.

  it's broken into 2 parts, radial functions and angular
   functions

************************************/

/***
  Recent Edit History:

   26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)

   14.02.2004 gL:
     Updated radial function calculators.

***/

#ifdef INCLUDE_ADF_PLOTS
/*************

  ADF plots only only need a single function for all
  possible angular and radial functions... here they are:

*************/
void eval_ADF_ang(point_type *loc,int kx,int ky,int kz,float *ang_val,
                  point_type *ang_grad)
{
  int i;
  float accum;
  float x,y,z;

  x = loc->x*AUI; y = loc->y*AUI; z = loc->z*AUI;

  if( kx > 0 ){
    accum = x;
    for(i=1;i<kx;i++)accum *= x;
  } else accum = 1.0;
  *ang_val = accum;
  if( ky > 0 ){
    accum = y;
    for(i=1;i<ky;i++)accum *= y;
    *ang_val *= accum;
  }
  if( kz > 0 ){
    accum = z;
    for(i=1;i<kz;i++)accum *= z;
    *ang_val *= accum;
  }
}

void eval_ADF_rad(point_type *loc,float r,int kr,float zeta,float *rad_val,
                  point_type *rad_grad)
{
  int i;
  float accum;

  if( kr > 0 ){
    accum = r;
    for(i=1;i<kr;i++) accum *= r;
  }else accum = 1.0;

  *rad_val = accum*exp(-zeta*r);
}
#endif


/************************************

       These are the angular functions.

  as arguments all of these functions take:

             loc: pointer to point_type
            r,r2: floats
             val: pointer to float
        gradient: pointer to point_type

   they return the value of the angular part of the wavefunction in 'val
   at location 'loc (with 'r = spherical coordiate r value
      and 'r2 = 'r^2)
   and the value of the gradient of the wavefunc in 'gradient

************************************/



void eval_s_ang(point_type *loc,float r,float r2,float *val,
                point_type *gradient)
{

  /* this one is easy.... */
  *val = 1.0;

#if 0
  gradient->x = 0.0;
  gradient->y = 0.0;
  gradient->z = 0.0;
#endif
}

void eval_px_ang(point_type *loc,float r,float r2,float *val,
                point_type *gradient)
{

  *val = sqrt3*loc->x/r;

#if 0
  gradient->x = sqrt3/r - loc->x*(*val)/r2;
  gradient->y = -loc->y*(*val)/r2;
  gradient->z = -loc->z*(*val)/r2;
#endif
}

void eval_py_ang(point_type *loc,float r,float r2,float *val,
                point_type *gradient)
{
  *val = sqrt3*loc->y/r;

#if 0
  gradient->x = -loc->x*(*val)/r2;
  gradient->y = sqrt3/r - loc->y*(*val)/r2;
  gradient->z = -loc->z*(*val)/r2;
#endif

}

void eval_pz_ang(point_type *loc,float r,float r2,float *val,
                point_type *gradient)
{

  *val = sqrt3*loc->z/r;

#if 0
  gradient->x = -loc->x*(*val)/r2;
  gradient->y = -loc->y*(*val)/r2;
  gradient->z = sqrt3/r - loc->z*(*val)/r2;
#endif
}

void eval_dx2y2_ang(point_type *loc,float r,float r2,float *val,
                    point_type *gradient)
{
  *val = 0.5*sqrt15*(loc->x*loc->x-loc->y*loc->y)/r2;

#if 0
  gradient->x = loc->x/r2*(sqrt15-2.0*(*val));
  gradient->y = -loc->y/r2*(sqrt15+2.0*(*val));
  gradient->x = -2.0*loc->z/r2*(*val);
#endif

}

void eval_dz2_ang(point_type *loc,float r,float r2,float *val,
                  point_type *gradient)
{
  *val = sqrt5*(loc->z*loc->z - 0.5*(loc->x*loc->x+loc->y*loc->y))/r2;

#if 0
  gradient->x = -loc->x/r2*(sqrt5 + 2.0*(*val));
  gradient->y = -loc->y/r2*(sqrt5 + 2.0*(*val));
  gradient->z = 2.0*loc->z/r2*(sqrt5 - *val);
#endif
}

void eval_dxy_ang(point_type *loc,float r,float r2,float *val,
                point_type *gradient)
{

  *val = sqrt15*(loc->x*loc->y)/r2;

#if 0
  gradient->x = (sqrt15*loc->y - 2*loc->x*(*val))/r2;
  gradient->y = (sqrt15*loc->x - 2*loc->y*(*val))/r2;
  gradient->z = -2.0*loc->z/r2*(*val);
#endif
}

void eval_dxz_ang(point_type *loc,float r,float r2,float *val,
                  point_type *gradient)
{

  *val = sqrt15*(loc->x*loc->z)/r2;

#if 0
  gradient->x = (sqrt15*loc->z - 2*loc->x*(*val))/r2;
  gradient->y = -2.0*loc->y/r2*(*val);
  gradient->z = (sqrt15*loc->x - 2*loc->z*(*val))/r2;
#endif
}

void eval_dyz_ang(point_type *loc,float r,float r2,float *val,
                point_type *gradient)
{

  *val = sqrt15*(loc->y*loc->z)/r2;

#if 0
  gradient->x = -2.0*loc->x/r2*(*val);
  gradient->y = (sqrt15*loc->z - 2*loc->y*(*val))/r2;
  gradient->z = (sqrt15*loc->y - 2*loc->z*(*val))/r2;
#endif
}

void eval_fz3_ang(point_type *loc, float r, float r2, float *val,
                point_type *gradient)
{

  *val = sqrt7*loc->z*(loc->z*loc->z - 1.5*(loc->x*loc->x + loc->y*loc->y))/(r*r2);

}

void eval_fxz2_ang(point_type *loc, float r, float r2, float *val,
                point_type *gradient)
{

  *val = sqrt42*loc->x*(loc->z*loc->z - 0.25*(loc->x*loc->x + loc->y*loc->y))/(r*r2);

}

void eval_fyz2_ang(point_type *loc, float r, float r2, float *val,
                point_type *gradient)
{

  *val = sqrt42*loc->y*(loc->z*loc->z - 0.25*(loc->x*loc->x + loc->y*loc->y))/(r*r2);

}

void eval_fxyz_ang(point_type *loc, float r, float r2, float *val,
                point_type *gradient)
{

  *val = sqrt105*loc->x*loc->y*loc->z/(2*r*r2);

}

void eval_fzx2_zy2_ang(point_type *loc, float r, float r2, float *val,
                point_type *gradient)
{

  *val = sqrt105*loc->z*(loc->x*loc->x - loc->y*loc->y)/(2*r*r2);

}

void eval_fx3_ang(point_type *loc, float r, float r2, float *val,
                point_type *gradient)
{

  *val = sqrt70*loc->x*(loc->x*loc->x - 3*loc->y*loc->y)/(4*r*r2);

}

void eval_fy3_ang(point_type *loc, float r, float r2, float *val,
                point_type *gradient)
{

  *val = sqrt70*loc->y*(3*loc->x*loc->x - loc->y*loc->y)/(4*r*r2);

}


/************************************

       These are the radial functions.

  as arguments all of these functions take:

             loc: pointer to point_type
            r,r2: floats
            zeta1,zeta2: floats
            C1,C2: floats
             val: pointer to float
        gradient: pointer to point_type

   they return the value of the radial part of the wavefunction in 'val
   at distance 'r,
   and the value of the gradient of the wavefunc in 'gradient

   they are sorted in order of increasing princ. quantum number


************************************/

/************************************************************************************
 *
 * Rnl(r) = ( sqrt((2*zeta) ** (2n+1)) / sqrt(fact(2*n)) )/4*pi * r**(n-1) * exp(-zeta*r)
 *
 * Thanks to Bob Kematick (kmatick@binghamton.edu) for providing this code
 *
************************************************************************************/
void eval_n_rad(int n,point_type *loc,float r,float r2,
                           float zeta1,float zeta2,float C1,
                        float C2,float *val,point_type *gradient)
{
  float norm_factor,v1,v2;
  /*****  a factorial table  *****/
  double f[21]={ 1.0, 1.0, 2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,
                 362880.0,3628800.0,39916800.0,479001600.0,6227020800.0 ,
                 87178291200.0 ,1307674368000.0 ,20922789888000.0 ,355687428096000.0 ,
                 6402373705728000.0 ,121645100408832000.0 ,2432902008176640000.0 };



  if( C1==0.0 || C2==0.0 ){
    norm_factor = sqrt( pow(2.0*zeta1,2*n+1) ) / sqrt(f[2*n]) / fourpi;
    *val=norm_factor * pow(r,n-1) * exp(-zeta1*r);
  } else{
    norm_factor = sqrt( pow(2.0*zeta1,2*n+1) ) / sqrt(f[2*n]) / fourpi;
    v1=norm_factor * pow(r,n-1) * exp(-zeta1*r);
    norm_factor = sqrt( pow(2.0*zeta2,2*n+1) ) / sqrt(f[2*n]) / fourpi;
    v2=norm_factor * pow(r,n-1) * exp(-zeta2*r);
    *val=(C1*v1+C2*v2);
  }

#if 0
  gradient->x = (*val)*loc->x/r*(2.0/r - zeta1);
  gradient->y = (*val)*loc->y/r*(2.0/r - zeta1);
  gradient->z = (*val)*loc->z/r*(2.0/r - zeta1);
#endif
}

#ifndef ORIG_NORMALIZATION
void eval_1_rad(point_type *loc,float r,float r2,
                float zeta1,float zeta2,float C1,
                float C2,float *val,point_type *gradient)
{
  eval_n_rad(1,loc,r, r2, zeta1,zeta2, C1, C2,val,gradient);

}

void eval_2_rad(point_type *loc,float r,float r2,
                float zeta1,float zeta2,float C1,
                float C2,float *val,point_type *gradient)
{
  eval_n_rad(2,loc,r, r2, zeta1,zeta2, C1, C2,val,gradient);

}

void eval_3_rad(point_type *loc,float r,float r2,
                float zeta1,float zeta2,float C1,
                float C2,float *val,point_type *gradient)
{
  eval_n_rad(3,loc,r, r2, zeta1,zeta2, C1, C2,val,gradient);

}

void eval_4_rad(point_type *loc,float r,float r2,
                float zeta1,float zeta2,float C1,
                float C2,float *val,point_type *gradient)
{
  eval_n_rad(4,loc,r, r2, zeta1,zeta2, C1, C2,val,gradient);

}

void eval_5_rad(point_type *loc,float r,float r2,
                float zeta1,float zeta2,float C1,
                float C2,float *val,point_type *gradient)
{
  eval_n_rad(5,loc,r, r2, zeta1,zeta2, C1, C2,val,gradient);

}
void eval_6_rad(point_type *loc,float r,float r2,
                float zeta1,float zeta2,float C1,
                float C2,float *val,point_type *gradient)
{
  eval_n_rad(6,loc,r, r2, zeta1,zeta2, C1, C2,val,gradient);

}

void eval_7_rad(point_type *loc,float r,float r2,
                float zeta1,float zeta2,float C1,
                float C2,float *val,point_type *gradient)
{
  eval_n_rad(7,loc,r, r2, zeta1,zeta2, C1, C2,val,gradient);
}
#else


void eval_1_rad(point_type *loc,float r,float r2,
                float zeta1,float zeta2,float C1,
                float C2,float *val,point_type *gradient)
{
  float norm_factor = 2.0*sqrt(zeta1*zeta1*zeta1);

  *val = norm_factor*exp(-zeta1*r);

#if 0
  gradient->x = -zeta1*loc->x/r*(*val);
  gradient->y = -zeta1*loc->y/r*(*val);
  gradient->z = -zeta1*loc->z/r*(*val);
#endif
}

void eval_2_rad(point_type *loc,float r,float r2,
                float zeta1,float zeta2,float C1,
                float C2,float *val,point_type *gradient)
{
  float norm_factor = 2.0*sqrt(pow(zeta1,5.0)/3.0)/fourpi;
  *val = norm_factor*r*exp(-zeta1*r);

#if 0
  gradient->x = (*val)*loc->x/r*(1.0/r - zeta1);
  gradient->y = (*val)*loc->y/r*(1.0/r - zeta1);
  gradient->z = (*val)*loc->z/r*(1.0/r - zeta1);
#endif
}


void eval_3_rad(point_type *loc,float r,float r2,
                float zeta1,float zeta2,float C1,
                float C2,float *val,point_type *gradient)
{
  float norm_factor;

  if( C1 == 0.0 || C2 == 0.0 ){
    norm_factor = 2.0*sqrt(pow(zeta1,7.0)*10.0)/(15.0*fourpi);
    *val = norm_factor*r*r*exp(-zeta1*r);
  } else{
    *val = (sqrt(8.0/45.0)/fourpi)*r*r*(C1*pow(zeta1,3.5)*exp(-zeta1*r) + C2*pow(zeta2,3.5)*exp(-zeta2*r));
  }


#if 0
  gradient->x = (*val)*loc->x/r*(2.0/r - zeta1);
  gradient->y = (*val)*loc->y/r*(2.0/r - zeta1);
  gradient->z = (*val)*loc->z/r*(2.0/r - zeta1);
#endif
}

void eval_4_rad(point_type *loc,float r,float r2,
                float zeta1,float zeta2,float C1,
                float C2,float *val,point_type *gradient)
{
  float norm_factor;

  if( C1 == 0.0 || C2 == 0.0 ){
    norm_factor = 2.0*sqrt(pow(zeta1,9.0)*35.0)/(105.0*fourpi);
    *val = norm_factor*pow(r,3.0)*exp(-zeta1*r);
  } else{
    *val = 2.0*sqrt(35.0)*pow(r,3.0)*(C1*pow(zeta1,4.5)*exp(-zeta1*r) +
                                      C2*pow(zeta2,4.5)*exp(-zeta2*r))/(105.0*fourpi);

  }
#if 0
  gradient->x = (*val)*loc->x/r*(2.0/r - zeta1);
  gradient->y = (*val)*loc->y/r*(2.0/r - zeta1);
  gradient->z = (*val)*loc->z/r*(2.0/r - zeta1);
#endif
}

void eval_5_rad(point_type *loc,float r,float r2,
                float zeta1,float zeta2,float C1,
                float C2,float *val,point_type *gradient)
{
  float norm_factor;

  if( C1==0.0 || C2==0.0 ){
    norm_factor = 2.0*sqrt(pow(zeta1,11.0)*14.0)/(315*fourpi);
    *val = norm_factor*pow(r,4.0)*exp(-zeta1*r);
  } else{
    *val = 2.0*sqrt(14.0)*pow(r,4.0)*(C1*pow(zeta1,5.5)*exp(-zeta1*r) +
                                      C2*pow(zeta2,5.5)*exp(-zeta2*r))/(315.0*fourpi);

  }
#if 0
  gradient->x = (*val)*loc->x/r*(2.0/r - zeta1);
  gradient->y = (*val)*loc->y/r*(2.0/r - zeta1);
  gradient->z = (*val)*loc->z/r*(2.0/r - zeta1);
#endif
}

void eval_6_rad(point_type *loc,float r,float r2,
                float zeta1,float zeta2,float C1,
                float C2,float *val,point_type *gradient)
{
  float norm_factor;

  if( C1==0.0 || C2==0.0 ){
    norm_factor = 2.0*sqrt(pow(zeta1,13.0)*462.0)/(10395.0*fourpi);
    *val = norm_factor*pow(r,5.0)*exp(-zeta1*r);
  }else{
    *val = 2.0*sqrt(462.0)*pow(r,5.0)*(C1*pow(zeta1,6.5)*exp(-zeta1*r) +
                                       C2*pow(zeta2,6.5)*exp(-zeta2*r))/(10395.0*fourpi);

  }
#if 0
  gradient->x = (*val)*loc->x/r*(2.0/r - zeta1);
  gradient->y = (*val)*loc->y/r*(2.0/r - zeta1);
  gradient->z = (*val)*loc->z/r*(2.0/r - zeta1);
#endif
}

void eval_7_rad(point_type *loc,float r,float r2,
                float zeta1,float zeta2,float C1,
                float C2,float *val,point_type *gradient)
{
  eval_n_rad(7,loc,r, r2, zeta1,zeta2, C1, C2,val,gradient);
}
#endif



/****************************************************************************
 *
 *                   Procedure build_radial_lookup_table
 *
 * Arguments: MO_surf: pointer to MO_surface_type
 *
 * Returns: none
 *
 * Action:  loops through all the AOs in 'MO_surf and builds
 *    their radial lookup tables
 *
 ****************************************************************************/
void build_radial_lookup_table(MO_surface_type *MO_surf)
{
  int i,j,step;
  float r_gap,r;
  float rad_val;
  point_type rad_grad,loc;

  AO_list_type *AO;
  MO_center_list_type *center;
  lookup_table_type *lookup_tbl;

  /* do the calculation savers */
  sqrt3 = sqrt(3.0);
  sqrt5 = sqrt(5.0);
  sqrt15 = sqrt(15.0);

  /* we just need to loop over the unique centers */
  fprintf(stderr,"Building radial lookup table.\n");
  for(i=0;i<MO_surf->num_unique;i++){
    center = &(MO_surf->unique_centers[i]);
    /* loop over the AOs */
    for(j=0;j<center->num_AOs;j++){
      AO = &(center->AO_list[j]);
      /* now fill the lookup table */
      lookup_tbl = AO->rad_lookup_tbl;

      /* first make sure that it hasn't been filled yet */
      if( !(lookup_tbl->filled) ){
fprintf(stderr,".");
        r_gap = (lookup_tbl->max_val - lookup_tbl->min_val) /
          (float)lookup_tbl->num_entries;
        for(step=0; step<lookup_tbl->num_entries; step++){
          r  = lookup_tbl->min_val + (float)step*r_gap;

          /*********

            we need to convert r into atomic units

            AUI is 1/BOHR where BOHR is the Bohr radius

          *********/
          r *= AUI;

#ifndef INCLUDE_ADF_PLOTS
          AO->rad_func(&loc,r,r*r,AO->zeta1,AO->zeta2,
                       AO->C1,AO->C2,&rad_val,&rad_grad);
#else
          if(!MO_surf->adf_plot){
            AO->rad_func(&loc,r,r*r,AO->zeta1,AO->zeta2,
                         AO->C1,AO->C2,&rad_val,&rad_grad);
          }else{
            eval_ADF_rad(&loc,r,AO->kr,AO->zeta1,&rad_val,&rad_grad);
          }
#endif

          lookup_tbl->values[step] = rad_val;
        }
        lookup_tbl->filled = 1;
      }
    }
fprintf(stderr,"\n");
  }
}



/****************************************************************************
 *
 *                   Procedure calc_MO_value
 *
 * Arguments:  which_MO: int
 *            MO_info: pointer to MO_info_type
 *            centers: pointer to MO_center_list_type
 *        num_centers: int
 *
 * Returns: none
 *
 * Action:  Calculates the value and gradient of the MO defined by
 *     'centers at 'MO_info->loc.
 *
 ****************************************************************************/
void calc_MO_value(int which_MO,MO_info_type *MO_info,
                   MO_center_list_type *centers,int num_centers,
                   char ADF_plot)
{
  static int first_call = 1;
  int i,j;
  AO_list_type *AO;
  point_type loc;
  point_type ang_grad;
  float rad_val, ang_val;
  float r,r2;
  float val,valI;

  /*****

    check to see if we need to evaluate the calculation savers

  ******/
  if( first_call ){
    sqrt3 = sqrt(3.0);
    sqrt5 = sqrt(5.0);
    sqrt15 = sqrt(15.0);
    first_call = 0;
  }

  MO_info->val = 0.0;
  val = 0.0;valI = 0.0;
  MO_info->grad.x = 0.0;MO_info->grad.y = 0.0;MO_info->grad.z = 0.0;

  /* loop over the centers */
  for(i=0;i<num_centers;i++){
    if( !centers[i].exclude ){
      /*****

        figure out the location of the particle relative to this
        center.

        ******/
      loc.x = MO_info->loc.x - centers[i].loc->x;
      loc.y = MO_info->loc.y - centers[i].loc->y;
      loc.z = MO_info->loc.z - centers[i].loc->z;

      r2 = loc.x*loc.x + loc.y*loc.y + loc.z*loc.z;
      r = sqrt(r2);

      /* loop over AO's on this center */
      if( r != 0.0 ){
        for(j=0;j<centers[i].num_AOs;j++){
          AO = &(centers[i].AO_list[j]);
          /*****

            evaluate the contributions of this AO to the MO and the
            gradient of the MO.

            ******/
          if( fabs(AO->coeff[which_MO]) > .001 ){
#ifndef INCLUDE_ADF_PLOTS
            AO->ang_func(&loc,r,r2,&ang_val,&ang_grad);
#else
            if( !ADF_plot ) AO->ang_func(&loc,r,r2,&ang_val,&ang_grad);
            else{
              /* evaluate the angular part */
              eval_ADF_ang(&loc,AO->kx,AO->ky,AO->kz,&ang_val,&ang_grad);
              /* multiply it by the normalization term */
              ang_val *= AO->norm_fact;
            }
#endif

            /* pull the radial value out of the lookup table */
            rad_val =
              READ_FROM_LOOKUP_TBL(AO->rad_lookup_tbl,r);
            val += AO->coeff[which_MO]*ang_val*rad_val;
            if(!ADF_plot) valI += AO->coeffI[which_MO]*ang_val*rad_val;
#if 0
            MO_info->grad.x += AO->coeff*(ang_val*rad_grad.x +
                                          rad_val*ang_grad.x);
            MO_info->grad.y += AO->coeff*(ang_val*rad_grad.y +
                                          rad_val*ang_grad.y);
            MO_info->grad.z += AO->coeff*(ang_val*rad_grad.z +
                                          rad_val*ang_grad.z);

#endif
          }
        }
      }
    }
  }
  MO_info->val = val;
}


