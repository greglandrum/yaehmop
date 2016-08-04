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

/****************************************************************************
 *
 *                   Procedure traverseobj
 *
 * Arguments:       obj: pointer to type object
 *            cor1,cor2: integers
 *                  win: a window
 * Returns:  object_type
 *
 * Action: recursively descends the object tree looking for the bounding
 *         box which encloses the point cor1,cor2 in the window win
 *
 ****************************************************************************/
object_type *traverseobj(object_type *obj,int cor1,int cor2,Window win)
{
    int i, found;
    point_type *min, *max;
#ifdef X_GRAPHICS
    XEvent event;
    object_type *tempobj;

    tempobj=0;
    min = &(obj->bmin);
    max = &(obj->bmax);

    found=0;

    /* the base case for the recursion */
    if( !obj )return(0);
    /* check to see if the bounding box at this node contains the point */
    if(win==xywin){
      if(cor1>=min->x && cor1<=max->x && cor2>=min->y && cor2<=max->y)
          found = 1;
    }
    else if(win==yzwin){
        if(cor1>=min->z && cor1<=max->z && cor2>=min->y && cor2<=max->y)
          found = 1;
    }
    else if(win==xzwin){
       if(cor1>=min->x && cor1<=max->x && cor2>=min->z && cor2<=max->z)
          found = 1;
    }
    else printf( "Bogus window passed to bb routine.\n." );

    /* if the point was within a bounding box then draw in the box in the ortho
       windows. */
    if( found ){
        drawboundbox(obj);
        /* prompt the user */
        display("This one?" );
        /* now wait for a button press */
        while( found ){
            XNextEvent(disp,&event);
            if(event.type==ButtonPress){
                switch(event.xbutton.button){
                    case 1: display("Got it.");
                      return(obj);
                      break;
                    case 2:
                    case 3: found=0; break;
                }
            }
        }
    }
    /* if we've made it to this point, then the object for this node is
       either not selected or is not the one desired. we must continue
       down the tree by recursively descending the children */
    i=0;
    while( !tempobj && obj->children[i] ){
        tempobj=traverseobj(obj->children[i],cor1,cor2,win);
        /* if the object is a child of this node keep a pointer to here */
        if( obj->children[i]==tempobj ){
            parentobj=obj;
            parenthead=0;
        }
        i++;
    }
    return( tempobj );
#endif
}

/****************************************************************************
 *
 *                   Procedure bbsearch
 *
 * Arguments: cor1,cor2: integers
 *                  win: a window
 * Returns: none
 *
 * Action: finds the bounding box which encloses the point selected by the
 *         user in ortho view window win
 *
 ****************************************************************************/
void bbsearch( int cor1, int cor2, Window win )
{
    head_type *temphead;
    object_type *obj;

    obj=0;
    temphead=head;
    while( temphead && !obj ){
        /* go through the object off of this node */
        obj = traverseobj(temphead->obj,cor1,cor2,win);
        /* see if the new object is the child of this node */
        if( temphead->obj == obj ){
             parenthead=temphead;
             parentobj=0;
        }
        temphead = temphead->next;
    }
    /* at this point there should be some object selected.  If not, the user
       either clicked in a bad spot or went through all options at the point
       where they clicked  */
    if( obj ) whichobj=obj;
    else display( "No object selected." );
}

