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
   06.04.98 gL:
     Added some error checking to drawit and instantiate.
      previously these were not so good about checking for
      null pointers.... now they are better.

   26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)
***/

#include "viewkel.h"

extern double atan2 PROTO((double,double));


/****************************************************************************
 *
 *                   Procedure addprojectmat
 *
 * Arguments: none
 * Returns: none
 *
 * Action: Gets the stack set up for the projection.  This means that the
 *     projection matrix and scaling matrix are premultiplied onto the
 *     stack in this procedure.
 *
 ****************************************************************************/
void addprojectmat(void)
{
  matrix_type *stacksave;
  static matrix_type *temp_mat = 0;
  point_type temp_pt;

  if( redo_projection ){
    /* in order to get this all properly built it is necessary to make a new
       matrix stack.  The old stack is saved by pointing stacksave at it
       */
    stacksave=stack;
    stack=0;

    /* load up the identity matrix */
    loadmatrix(ident);
    pushmatrix(0);

    if( !temp_mat ){
      temp_mat = (matrix_type *)calloc(1,sizeof(matrix_type));
      if( !temp_mat )fatal("Memory allocation." );
    }


    /* save the view vector as a point */
    temp_pt.x = camera->la.x - camera->lf.x;
    temp_pt.y = camera->la.y - camera->lf.y;
    temp_pt.z = camera->la.z - camera->lf.z;

    /* put on the main projection matrix */
    temp_mat->matrix[0][0]=1;
    temp_mat->matrix[1][1]=1;
    temp_mat->matrix[2][2]=1;
    temp_mat->matrix[3][2]=1.0/camera->foclength;
    multmatrix(temp_mat);

    scale((float)g_xmax,(float)g_xmax,1.0);
    translate( 0.0, 0.0, -(camera->foclength));

    /* translate lookfrom to the origin */
    translate( -(camera->lf.x), -(camera->lf.y), -(camera->lf.z) );

    /* get a copy of the matrix */
    getmatrix(mainortho);
    popmatrix();

    /* restore things */
    stack = stacksave;
    redo_projection = 0;
  }
  /* premultiply by the projection matrix */
  premultmatrix(mainortho);
}



/****************************************************************************
 *
 *                   Procedure drawit
 *
 * Arguments: prim: a pointer to a primitive type
 *             obj: a pointer to object_type
 * Returns: none
 *
 * Action: transforms the points of the primitive and then displays it
 *   in the three D window and in the ortho view windows (if ortho views are
 *   currently active.
 *
 ****************************************************************************/
void drawit( prim_type *prim,object_type *obj )
{
  if (!prim || !obj) error("null pointer in drawit");

  if( projviewon ){
    /* now draw the damn thing */
    switch(prim->which){
    case MOLECULE:
      /* load the perspective view matrix */
      pushmatrix(0);
      addprojectmat();
      draw_3D_objects(prim,obj);
      /* clean up the stack */
      popmatrix();
      break;
    case GRAPH:
      draw_graph(prim,obj);
      break;
    case CONT_PLOT:
      draw_cont_plot(prim,obj);
      break;
    case PROP_GRAPH:
      draw_graph(prim,obj);
      break;
    case BAND_GRAPH:
      draw_band_graph(prim,obj);
      break;
    case WALSH_GRAPH:
      draw_walsh_graph(prim,obj);
      break;
    case FMO_DIAGRAM:
      draw_FMO_diagram(prim,obj);
      break;
    case MO_SURF:
      /* load the perspective view matrix */
      pushmatrix(0);
      addprojectmat();
      draw_3D_objects(prim,obj);
      /* clean up the stack */
      popmatrix();

      break;
    case LABEL:
      draw_label(prim,obj);
      break;

    default:
      FATAL_BUG("Bogus primitive in drawit.");
      break;
    }
  }
}

/****************************************************************************
 *
 *                   Procedure instantiate
 *
 * Arguments: object: pointer to object_type
 * Returns: none
 *
 * Action: recusively descend the current object, keeping the matrix stack
 *  updated and drawing each primitive as it is hit. This also updates the
 *  bounding box for each primitive.
 *
 ****************************************************************************/
void instantiate( object_type *object, char do_labels )
{
  int i;

  if( !object ) return;

  /*******

    stack operations only need to be performed for 3D objects

  ********/
  if( object->prim && (object->prim->which == MOLECULE
                       || object->prim->which == MO_SURF)){
    /* update the current matrix */
    pushmatrix(0);

    /* undo the translation */
    translate(object->trans.x/20.0,-object->trans.y/20.0,
              object->trans.z/20.0);

    /* the scaling */
    scale(object->scale.x,object->scale.y,object->scale.z);

    /* the rotation */
    xrot(object->rot.x);
    yrot(object->rot.y);
    zrot(object->rot.z);
  }

  /* if there is a primitive at this node then draw it now */
  if( !do_labels ){
    if( object->prim && !object->prim->label ){
      drawit( object->prim, object );
    }
  } else if(object->prim && object->prim->label){
    drawit( object->prim, object );
  }


  /* update the global bounding box */
  if( object->bmin.x<globalmin.x)globalmin.x=object->bmin.x;
  if( object->bmax.x>globalmax.x)globalmax.x=object->bmax.x;
  if( object->bmin.y<globalmin.y)globalmin.y=object->bmin.y;
  if( object->bmax.y>globalmax.y)globalmax.y=object->bmax.y;

  /* now go through the children of the node */
  i=0;
  while( object->children[i] ){
    /* instantiate this child */
    instantiate(object->children[i],do_labels);

    /* move to the next child */
    i++;
  }
  if( object->prim && (object->prim->which == MOLECULE ||
                       object->prim->which == MO_SURF)){
    /* clean up the stack */
    popmatrix();
  }
}

/****************************************************************************
 *
 *                   Procedure redrawgraph
 *
 * Arguments: none
 * Returns: none
 *
 * Action: goes through the objects and redraws them
 *
 ****************************************************************************/
void redrawgraph(void)
{
  head_type *temphead;

  g_clear_screen();

  /* set the global bounding box to zero */
  globalmin.x = globalmin.y = 1000.0;
  globalmax.x = globalmax.y = -1000.0;

  /********
    the stack is expected at this point to have only the main
     identity matrix atop it
  ********/
  temphead = head;
  /* first go through the head list and draw objects */
  while( temphead ){
    instantiate(temphead->obj,0);
    temphead = temphead->next;
  }
  temphead = head;
  /* go through again and draw labels */
  while( temphead ){
    instantiate(temphead->obj,1);
    temphead = temphead->next;
  }

  g_switch_buffers();

}

/****************************************************************************
 *
 *                   Procedure redraw
 *
 * Arguments: none
 * Returns: none
 *
 * Action: goes through the objects and redraws everything
 *
 ****************************************************************************/
void redraw(void)
{
#ifdef X_GRAPHICS
  if( doing_X ){
    /* redraw the buttons */
    g_draw_all_buttons(button_wins);
  }
#endif

  redrawgraph();
}


