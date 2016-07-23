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
/******************************

  contour.h

  This file contains the definitions required to draw contour plots

  Nabbed from gnuplot 3.5 by gL, June 1996

*******************************/

/***
  Recent modification history

  26.09.98 gL:
    added a prev pointer to iso_curve_type to make
      that a doubly linked list (to enable use of
      symmetry without sorting in the evaluation of
      MO planes).
***/
#ifndef _CONTOUR_
#define _CONTOUR_

#ifndef _MY_DEFINES_
#include "defines.h"
#endif

#ifndef _BASIC_OBJECTS_
#include "basic_objects.h"
#endif

struct curve_points {
  struct curve_points *next_cp; /* pointer to next plot in linked list */
  int p_max;                    /* how many points are allocated */
  int p_count;                  /* count of points in points */
  point_type *points;
};

struct gnuplot_contours {
  struct gnuplot_contours *next;
  point_type *coords;
  char isNewLevel;
  char label[12];
  int num_pts;
};

typedef struct gnuplot_contours gnuplot_contour_type;

typedef struct iso_curve {
  struct iso_curve *next, *prev;
  int p_max;   /* how many points are allocated */
  int p_count; /* count of points in points */
  point_type *points;
} iso_curve_type;

#endif
