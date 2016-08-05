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

#include "viewkel.h"


/*******************

  this has got the stuff for dealing with labels on graphs

  created by greg Landrum  Jun 1995

*******************/


/***
  Recent Edit History:

  24.09.98 gL:
    bounding boxes evaluated in draw_label to allow
    selection of labels.
    labels are now drawn using a call to g_label
  26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)
  22.05.99 gL:
     changes to make adding a label work.  yay.

***/

#define MAX_LINES_IN_LABEL 20

/****************************************************************************
 *
 *                   Procedure new_label
 *
 * Arguments: none
 *
 * Returns: none
 *
 * Action: sets up a new label for the graph
 *
 ****************************************************************************/
void new_label(char *text)
{
  int num_lines,i;
  object_type *new_label;
  char *string_to_change[MAX_LINES];

  num_lines = 0;

  /*********

    set up a new object to hold the label...
    to keep things convenient, we are going to
    store the labels in a different object tree

  *********/
  makenewobject();

  new_label = head->obj;

  /* now build the label primitive */
  new_label->prim = (prim_type *)D_CALLOC(1,sizeof(prim_type));
  if( !new_label->prim )fatal("Can't get space for label primitive.");
  new_label->prim->which = LABEL;

  new_label->prim->label = (label_type *)D_CALLOC(1,sizeof(label_type));
  if( !new_label->prim->label )
    fatal("Can't get space for label.");

  if( !text ){
    string_to_change[0] = new_label->prim->label->text[0];
    readstringparm( "label",string_to_change);
    new_label->prim->label->num_lines = 1;
  }else{
    strcpy(new_label->prim->label->text[0],text);
    new_label->prim->label->num_lines=1;
  }

  new_label->scale.x=new_label->scale.y=new_label->scale.z=2.0;
  new_label->cent.x=180;new_label->cent.y=180;
  new_label->cent.z=0;
  new_label->trans.x=0;new_label->trans.y=0;
  new_label->trans.z=0;

}




/****************************************************************************
 *
 *                   Procedure draw_label
 *
 * Arguments: prim: pointer to primitive_type
 *             obj: pointer to object_type
 *
 * Returns: none
 *
 * Action: Draws in a label
 *
 ****************************************************************************/
void draw_label(prim,obj)
  prim_type *prim;
  object_type *obj;
{
  label_type *label;

  label = prim->label;
  if( !label ) FATAL_BUG("Bogus label passed to draw_label");

  /* this just works for a single line of text right now */
  localmax.x = 0;  localmax.y = 0;
  localmin.x=10000;  localmin.y=10000;
  g_label(obj->cent.x,obj->cent.y,label);
  obj->bmax.x = localmax.x;  obj->bmax.y = localmax.y;
  obj->bmin.x = localmin.x;  obj->bmin.y = localmin.y;

}

/****************************************************************************
 *
 *                   Procedure adjust_label
 *
 * Arguments: label: pointer to label_type
 *
 * Returns: none
 *
 * Action: Adjusts parameters for a label
 *
 ****************************************************************************/
void adjust_label(label_type *label)
{
  if(label->num_lines){
    readintparm("show lines",&(label->show_lines));
  }
}
