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

