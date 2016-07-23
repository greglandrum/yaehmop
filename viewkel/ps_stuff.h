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

/***
  Recent Edit History:

   08.09.98 gL:
     support for toggling atom and bond types in ps_options


***/
/******************************

  ps_stuff.h

  This file contains the stuff for dealing with postscript printing

  Created by gL June 1995

*******************************/

#ifndef _PS_STUFF_
#define _PS_STUFF_

#ifndef _MY_DEFINES_
#include "defines.h"
#endif

/*********

  defines for modes and the default fonts

*********/
#define PS_PRINT_BOTTOM 0
#define PS_PRINT_MIDDLE 1
#define PS_PRINT_TOP 2

#define DEF_PS_FONT "Times-Roman"
#define DEF_PS_FONT_SIZE 12
#define DEF_PS_SCALE 1.0

/* justification modes */
#define LEFT_JUST 0
#define CENTER_JUST 1
#define RIGHT_JUST 2
#define BOTTOM_JUST 3
#define TOP_JUST 4

/*************

  Postscript printing options

**************/
typedef struct {
  int where_to_print;
  char fontname[120];
  float fontsize;
  float printscale;
  char tight_b_box;
  point_type min_pt, max_pt;
  int atom_sphere_type;
  int bond_type;
  float atom_shade, bond_shade;
  float atom_color[3], bond_color[3];
} PS_options_type;

#endif
