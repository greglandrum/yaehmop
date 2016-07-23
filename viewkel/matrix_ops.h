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
/*******************************************************
*
*  This file is part of yaehmop.
*
*  There are, of course, no warranties at all on this program.
*
********************************************************************/
/***
  Recent Edit History:

  04.09.98 gL:
    Conditional definitions of ABS

***/

/****************************

  These are operations on matrices and vectors adapted from:
     2d and 3d Vector C Library
     by Andrew Glassner
     from "Graphics Gems", Academic Press, 1990

*****************************/

#ifndef _MATRIX_OPS_
#define _MATRIX_OPS_

#ifndef _MY_DEFINES_
#include "defines.h"
#endif

/*********************/
/* 3d geometry types */
/*********************/

typedef point_type Point3;
typedef point_type Vector3;

typedef struct IntPoint3Struct {/* 3d integer point */
  int x, y, z;
} IntPoint3;

/***********************/
/* one-argument macros */
/***********************/

/* absolute value of a */
#ifndef ABS
#define ABS(a) (((a) < 0) ? -(a) : (a))
#endif

/* round a to nearest int */
#define ROUND(a) floor((a) + 0.5)

/* take sign of a, either -1, 0, or 1 */
#define ZSGN(a) (((a) < 0) ? -1 : (a) > 0 ? 1 : 0)

/* take binary sign of a, either -1, or 1 if >= 0 */
#define SGN(a) (((a) < 0) ? -1 : 1)

/* shout if something that should be true isn't */
#define ASSERT(x)                                                              \
  if (!(x))                                                                    \
    fprintf(stderr, " Assert failed: x\n");

/* square a */
#define SQR(a) ((a) * (a))

/***********************/
/* two-argument macros */
/***********************/

/* find minimum of a and b */
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

/* find maximum of a and b */
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

/* swap a and b (see Gem by Wyvill) */
#define SWAP(a, b)                                                             \
  {                                                                            \
    a ^= b;                                                                    \
    b ^= a;                                                                    \
    a ^= b;                                                                    \
  }

/* linear interpolation from l (when a=0) to h (when a=1)*/
/* (equal to (a*h)+((1-a)*l) */
#define LERP(a, l, h) ((l) + (((h) - (l)) * (a)))

/* clamp the input to the specified range */
#define CLAMP(v, l, h) ((v) < (l) ? (l) : (v) > (h) ? (h) : v)

/****************************/
/* memory allocation macros */
/****************************/

/* create a new instance of a structure (see Gem by Hultquist) */
#define NEWSTRUCT(x) (struct x *)(D_MALLOC((unsigned)sizeof(struct x)))

/* create a new instance of a type */
#define NEWTYPE(x) (x *)(D_MALLOC((unsigned)sizeof(x)))

/********************/
/* useful constants */
/********************/

#define PITIMES2 6.283185 /* 2 * pi */
#define PIOVER2 1.570796  /* pi / 2 */
#define EVAL 2.718282     /* the venerable e */
#define SQRT2 1.414214    /* sqrt(2) */
#define SQRT3 1.732051    /* sqrt(3) */
#define GOLDEN 1.618034   /* the golden ratio */
#define DTOR 0.017453     /* convert degrees to radians */
#define RTOD 57.29578     /* convert radians to degrees */

/************/
/* booleans */
/************/

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif
#define ON 1
#define OFF 0
typedef int boolean;  /* boolean data type */
typedef boolean flag; /* flag data type */

#endif
