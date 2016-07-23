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

/********

  this has got the stuff for dealing with coordination
  polyhedra

  Created by greg Landrum 03.05.98

*********/
#ifndef _POLYHED_
#define _POLYHED_

#define DOT_TOL 0.000001
#define BOX_TOL 0.000001

/****

  this is used to maintain a linked list of active
  triangles.

****/
typedef struct triangle_list_def {
  triangle_type tri;
  point_type bmax, bmin;
  struct triangle_list_def *next;
} triangle_list_type;

/* check to see if a point falls in a bounding region */
#define POINT_IS_INSIDE(_p_, _max_, _min_)                                     \
  (_p_.x - _max_.x < BOX_TOL && _p_.x - _min_.x > BOX_TOL &&                   \
   _p_.y - _max_.y < BOX_TOL && _p_.y - _min_.y > BOX_TOL &&                   \
   _p_.z - _max_.z < BOX_TOL && _p_.z - _min_.z > BOX_TOL)

#endif
