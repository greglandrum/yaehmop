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
#if defined(_NO_PROTO) || defined(_alpha) || defined(MIPSEL)
#define PROTO(x) ()
#else /* defined(_NO_PROTO) || defined(__alpha) || defined(MIPSEL) */
#define PROTO(x) x
#endif /* defined(_NO_PROTO) || defined(__alpha) || defined(MIPSEL) */
#endif /* PROTO */

#ifdef INCLUDE_ADF_PLOTS
extern void eval_ADF_ang PROTO((point_type * loc, int kx, int ky, int kz,
                                float *ang_val, point_type *ang_grad));
extern void eval_ADF_rad PROTO((point_type * loc, float r, int kr, float zeta,
                                float *rad_val, point_type *rad_grad));
#endif
extern void eval_s_ang
PROTO((point_type * loc, float r, float r2, float *val, point_type *gradient));
extern void eval_px_ang
PROTO((point_type * loc, float r, float r2, float *val, point_type *gradient));
extern void eval_py_ang
PROTO((point_type * loc, float r, float r2, float *val, point_type *gradient));
extern void eval_pz_ang
PROTO((point_type * loc, float r, float r2, float *val, point_type *gradient));
extern void eval_dx2y2_ang
PROTO((point_type * loc, float r, float r2, float *val, point_type *gradient));
extern void eval_dz2_ang
PROTO((point_type * loc, float r, float r2, float *val, point_type *gradient));
extern void eval_dxy_ang
PROTO((point_type * loc, float r, float r2, float *val, point_type *gradient));
extern void eval_dxz_ang
PROTO((point_type * loc, float r, float r2, float *val, point_type *gradient));
extern void eval_dyz_ang
PROTO((point_type * loc, float r, float r2, float *val, point_type *gradient));
extern void eval_fz3_ang
PROTO((point_type * loc, float r, float r2, float *val, point_type *gradient));
extern void eval_fxz2_ang
PROTO((point_type * loc, float r, float r2, float *val, point_type *gradient));
extern void eval_fyz2_ang
PROTO((point_type * loc, float r, float r2, float *val, point_type *gradient));
extern void eval_fxyz_ang
PROTO((point_type * loc, float r, float r2, float *val, point_type *gradient));
extern void eval_fzx2_zy2_ang
PROTO((point_type * loc, float r, float r2, float *val, point_type *gradient));
extern void eval_fx3_ang
PROTO((point_type * loc, float r, float r2, float *val, point_type *gradient));
extern void eval_fy3_ang
PROTO((point_type * loc, float r, float r2, float *val, point_type *gradient));
extern void eval_1_rad
PROTO((point_type * loc, float r, float r2, float zeta1, float zeta2, float C1,
       float C2, float *val, point_type *gradient));
extern void eval_2_rad
PROTO((point_type * loc, float r, float r2, float zeta1, float zeta2, float C1,
       float C2, float *val, point_type *gradient));
extern void eval_3_rad
PROTO((point_type * loc, float r, float r2, float zeta1, float zeta2, float C1,
       float C2, float *val, point_type *gradient));
extern void eval_4_rad
PROTO((point_type * loc, float r, float r2, float zeta1, float zeta2, float C1,
       float C2, float *val, point_type *gradient));
extern void eval_5_rad
PROTO((point_type * loc, float r, float r2, float zeta1, float zeta2, float C1,
       float C2, float *val, point_type *gradient));
extern void eval_6_rad
PROTO((point_type * loc, float r, float r2, float zeta1, float zeta2, float C1,
       float C2, float *val, point_type *gradient));
extern void eval_7_rad
PROTO((point_type * loc, float r, float r2, float zeta1, float zeta2, float C1,
       float C2, float *val, point_type *gradient));
