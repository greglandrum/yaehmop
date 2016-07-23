/*******************************************************
*      Copyright (C) 1996, 1999 Greg Landrum
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


/******************************

  labels.h

  This file contains the stuff for labels
   needed by viewkel.

  Created by gL October 1996

*******************************/

/***
  Recent Edit History:

  24.09.98 gL:
    various additions to label_type
    
***/
#ifndef _LABELS_

#define _LABELS_

#ifndef _MY_DEFINES_
#include "defines.h"
#endif



/********

  used for labels

*********/
typedef struct{
  int num_lines;
  char text[MAX_LINES][NORMAL_STR_LEN];

  /*******

    this stuff is so we can draw lines to atoms
     from the label

  *******/
  int show_lines,show_arcs;
  int justification;
  int num_atoms_labelled;
  atom_type **atoms_to_label;

} label_type;

#endif
