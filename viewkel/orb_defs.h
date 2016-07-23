/*******************************************************
*      Copyright (C) 1995, 1998, 1999 Greg Landrum
*
*  This file is part of yaehmop.
*
*   This is free software.
* 
*  Permission is granted to modify, or otherwise fold, spindle, and mutilate this
*    code provided all copyright notices are left intact.
*
*  This code may be distributed to your heart's content, in whatever form,
*    provided no fee is charged for the distribution, all copyright notices are
*    left intact, and the source is distributed (without fee) along with any
*    binaries to anyone who requests it.
*
*  There are, of course, no warranties at all on this program.
*
********************************************************************/


/*********

  This file contains the prototypes for orbital functions

**********/

/***
  Edit History :

  Wingfield Glassey 11th April 1998

  - prototypes for eval_{fz3,fxz2,fyz2,fxyz,fzx2-zy2,fx3,fy3}_ang() defined

   14.02.2004 gL:
   support degree 7 radial functions.
    
***/


#ifndef PROTO
# if defined(_NO_PROTO) || defined(_alpha) || defined(MIPSEL)
#  define PROTO(x) ()
# else /* defined(_NO_PROTO) || defined(__alpha) || defined(MIPSEL) */
#  define PROTO(x) x
# endif /* defined(_NO_PROTO) || defined(__alpha) || defined(MIPSEL) */
#endif /* PROTO */

#ifdef INCLUDE_ADF_PLOTS
extern void eval_ADF_ang PROTO((point_type *loc,int kx,int ky,int kz,float *ang_val,
				point_type *ang_grad));
extern void eval_ADF_rad PROTO((point_type *loc,float r,int kr,float zeta,float *rad_val,
				point_type *rad_grad));
#endif
extern void eval_s_ang PROTO((point_type *loc,float r,float r2,
			      float *val,point_type *gradient));
extern void eval_px_ang PROTO((point_type *loc,float r,float r2,
			       float *val,point_type *gradient));
extern void eval_py_ang PROTO((point_type *loc,float r,float r2,
			       float *val,point_type *gradient));
extern void eval_pz_ang PROTO((point_type *loc,float r,float r2,
			       float *val,point_type *gradient));
extern void eval_dx2y2_ang PROTO((point_type *loc,float r,float r2,
				  float *val,point_type *gradient));
extern void eval_dz2_ang PROTO((point_type *loc,float r,float r2,
				float *val,point_type *gradient));
extern void eval_dxy_ang PROTO((point_type *loc,float r,float r2,
				float *val,point_type *gradient));
extern void eval_dxz_ang PROTO((point_type *loc,float r,float r2,
				float *val,point_type *gradient));
extern void eval_dyz_ang PROTO((point_type *loc,float r,float r2,
				float *val,point_type *gradient));
extern void eval_fz3_ang PROTO((point_type *loc,float r,float r2,
				float *val,point_type *gradient));
extern void eval_fxz2_ang PROTO((point_type *loc,float r,float r2,
				float *val,point_type *gradient));
extern void eval_fyz2_ang PROTO((point_type *loc,float r,float r2,
				float *val,point_type *gradient));
extern void eval_fxyz_ang PROTO((point_type *loc,float r,float r2,
				float *val,point_type *gradient));
extern void eval_fzx2_zy2_ang PROTO((point_type *loc,float r,float r2,
				float *val,point_type *gradient));
extern void eval_fx3_ang PROTO((point_type *loc,float r,float r2,
				float *val,point_type *gradient));
extern void eval_fy3_ang PROTO((point_type *loc,float r,float r2,
				float *val,point_type *gradient));
extern void eval_1_rad PROTO((point_type *loc,float r,float r2,
			      float zeta1,float zeta2,float C1,
			      float C2,float *val,point_type *gradient));
extern void eval_2_rad PROTO((point_type *loc,float r,float r2,
			      float zeta1,float zeta2,float C1,
			      float C2,float *val,point_type *gradient));
extern void eval_3_rad PROTO((point_type *loc,float r,float r2,
			      float zeta1,float zeta2,float C1,
			      float C2,float *val,point_type *gradient));
extern void eval_4_rad PROTO((point_type *loc,float r,float r2,
			      float zeta1,float zeta2,float C1,
			      float C2,float *val,point_type *gradient));
extern void eval_5_rad PROTO((point_type *loc,float r,float r2,
			      float zeta1,float zeta2,float C1,
			      float C2,float *val,point_type *gradient));
extern void eval_6_rad PROTO((point_type *loc,float r,float r2,
			      float zeta1,float zeta2,float C1,
			      float C2,float *val,point_type *gradient));
extern void eval_7_rad PROTO((point_type *loc,float r,float r2,
			      float zeta1,float zeta2,float C1,
			      float C2,float *val,point_type *gradient));

