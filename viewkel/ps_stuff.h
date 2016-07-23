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
#define PS_PRINT_TOP    2

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
typedef struct{
  int where_to_print;
  char fontname[120];
  float fontsize;
  float printscale;
  char tight_b_box;
  point_type min_pt,max_pt;
  int atom_sphere_type;
  int bond_type;
  float atom_shade,bond_shade;
  float atom_color[3],bond_color[3];
} PS_options_type;



#endif
