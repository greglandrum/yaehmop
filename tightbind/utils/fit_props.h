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
  This is the include file for the property fitting programs

   Created by greg Landrum March 1994
************************************************************************/
#include <stdio.h>
#include <math.h>

/******
  These are used by the fileio routines.
*******/
#define FATAL 0
#define ERROR 1
#define IGNORE 2

/********
  The program will compile to use doubles instead of floats to represent
everything.
  If for some reason you want to use floats, then you must put -DUSE_FLOATS in
the
  CFLAGS field of the makefile
*********/
#ifndef USE_FLOATS
typedef double real;
#else
typedef float real;
#endif

#ifndef USE_BZERO
#define bzero(a, b) (memset((void *)(a), 0, (b)))
#define bcopy(a, b, c) (memcpy((void *)(b), (const void *)(a), (c)))
#endif

/* this is the step size used to generate the smoothed data (in eV) */
#define ENERGY_STEP .01

/* this is the broadening factor */
#define BROADENING .01

#define PI 3.141592653589793
#ifndef M_PI
#define M_PI 3.141592653589793
#endif

#define MAX_STR_LEN 2048

/******
  this is the data type that will be used to store the properties data
  within the program.
******/
typedef struct {
  real height;
  real energy;
} point_type;
