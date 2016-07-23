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

  basic_objects.h

  This file contains the basic objects (points, lines, etc.)
   needed by viewkel.

  Created by gL June 1995

*******************************/

/*******
  Recent modification history

********/

#ifndef _BASIC_OBJECTS_

#define _BASIC_OBJECTS_

#ifndef _MY_DEFINES_
#include "defines.h"
#endif

/********

  this is the node for the matrix stack

********/
typedef struct matrix_type_def {
  float matrix[DIM][DIM];
  struct matrix_type_def *next;
} matrix_type;

/**********

  a 3D point

**********/
typedef struct { float x, y, z; } point_type;

/***********

  a 2D point

*************/
typedef struct { float x, y; } point_type2D;

/************

  used to draw lines between atoms

************/
typedef struct {
  long int end1, end2;
  char custom, drawn, tube, breaking, arrow, dashed;
  float length;
  int thickness;
  int type;
} line_type;

/************

  a triangle with vertices as tabs (used to draw polygonal surfaces)


*************/
typedef struct {
  int vertices[3];
  point_type center, normal;
  char perp;
  int color;
} triangle_type;

/***********

  triangle vertex, with position and normal

************/
typedef struct vertex { point_type position, normal; } vertex_type;

/************

  a triangle with real vertices specified (used to draw polygonal surfaces)

*************/
typedef struct {
  point_type vertices[3];
  point_type center, normal;
  int color;
} explicit_triangle_type;

/* a particle */
typedef struct particle_type_def {
  char on_surf;
  char phase;
  int num_neighbors;
  int num;
  float rad;
  float val;
  point_type loc;
  struct particle_type_def *next;
  struct particle_type_def *neighbors[6];
} particle_type;

#endif
