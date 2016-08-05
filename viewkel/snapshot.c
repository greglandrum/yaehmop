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
 *                   Procedure readobj
 *
 * Arguments:  snapptr: a pointer to snap_type
 *             obj: a pointer to object_type
 * Returns: a pointer to snap_type
 *
 * Action: reads the information in the snap node pointed to by snapptr
 *       into object obj,  then recurses over the children of obj and fills
 *       them in.  The value returned is the next unused snapshot
 *
 ****************************************************************************/
snap_type *readobj(snap_type *snapptr,object_type *obj)
{
  int i;
  snap_type *tempsnap;

  /* first read in the information from the snapshot */
  obj->cent.x=snapptr->xc;obj->cent.y=snapptr->yc;obj->cent.z=snapptr->zc;
  obj->trans.x=snapptr->xt;obj->trans.y=snapptr->yt;obj->trans.z=snapptr->zt;
  obj->scale.x=snapptr->xs;obj->scale.y=snapptr->ys;obj->scale.z=snapptr->zs;
  obj->rot.x=snapptr->xr;obj->rot.y=snapptr->yr;obj->rot.z=snapptr->zr;
  tempsnap = snapptr->next;

  /* now recurse over the children */
  i=0;
  while( obj->children[i] ){
    tempsnap = readobj( tempsnap, obj->children[i] );
    i++;
  }

  /* finally return the last unused snap shot */
  return( tempsnap );
}

/****************************************************************************
 *
 *                   Procedure readsnap
 *
 * Arguments:  which: an integer
 *             shead: a pointer to shead_type
 * Returns: none
 *
 * Action: reads snapshot which into the object hierarchy
 *
 ****************************************************************************/
void readsnap( int which, shead_type *shead )
{
  char string1[80], string2[80];
  shead_type *tempshead;
  snap_type *snapptr;
  head_type *temphead;

  tempshead=shead;

  /* find the proper keyframe */
  while( tempshead && tempshead->which < which )tempshead=tempshead->next;
  /* make sure that the keyframe actually exists */
  if(!tempshead || tempshead->which != which){
    /* display an error message and return */
    sprintf(string1,"%d",which);
    strcpy(string2, " <-no such keyframe");
    strcat(string1,string2);
    display( string1 );
    return;
  }
  /* set up the camera parms */
  camera->vup.x=tempshead->vup.x;camera->vup.y=tempshead->vup.y;
  camera->vup.z=tempshead->vup.z;

  camera->la.x=tempshead->la.x;camera->la.y=tempshead->la.y;
  camera->la.z=tempshead->la.z;

  camera->lf.x=tempshead->lf.x;camera->lf.y=tempshead->lf.y;
  camera->lf.z=tempshead->lf.z;

  camera->foclength=tempshead->foclength;

  /* now fill in the hierarchy tree */
  snapptr = tempshead->snap;
  temphead=head;
  while(snapptr){
    snapptr=readobj(snapptr,temphead->obj);
    /* no error checking is necessary here because the hierarchy tree
       not change.  If it does: oh well, you lose */
    temphead=temphead->next;
  }

  redrawgraph();
}

/****************************************************************************
 *
 *                   Procedure snapobj
 *
 * Arguments: obj: a pointer to object_type
 * Returns: a pointer to snap_type
 *
 * Action: puts the parameters for the current object into the next snapshot
 *  object and then recurses over the children of the object.  The pointer
 *  returned is to the first snap object that this fills (it points to all
 *  of the others encountered in the process of the recursion).
 *
 ****************************************************************************/
snap_type *snapobj( object_type *obj )
{
  int i;
  snap_type *newsnap, *tempsnap;

  newsnap=(snap_type *)D_CALLOC(1,sizeof(snap_type));
  if(!newsnap)fatal("Memory Allocation.");

  /* fill in the parameters for this object */
  newsnap->xc=obj->cent.x;newsnap->yc=obj->cent.y;newsnap->zc=obj->cent.z;
  newsnap->xt=obj->trans.x;newsnap->yt=obj->trans.y;newsnap->zt=obj->trans.z;
  newsnap->xs=obj->scale.x;newsnap->ys=obj->scale.y;newsnap->zs=obj->scale.z;
  newsnap->xr=obj->rot.x;newsnap->yr=obj->rot.y;newsnap->zr=obj->rot.z;

  /* now recurse over the children */
  i=0;
  tempsnap=newsnap;
  while(obj->children[i]){
    tempsnap->next=snapobj(obj->children[i]);
    /* now move to the end of the snapshot list returned */
    while(tempsnap->next)tempsnap=tempsnap->next;
    /* check the next child */
    i++;
  }
  /* return a pointer to the new snapshot */
  return( newsnap );
}

/****************************************************************************
 *
 *                   Procedure takesnap
 *
 * Arguments: which: an integer
 * Returns: none
 *
 * Action: Saves the current state into keyframe which
 *
 ****************************************************************************/
void takesnap( int which )
{
  shead_type *tempshead,*lastshead;
  head_type *temphead;

  /* find the correct location for this keyframe in the list */
  tempshead = shead;
  if( !tempshead ){
    /* this is the trivial case */
    shead = (shead_type *)D_CALLOC(1,sizeof(shead_type));
    if(!shead) fatal("Memory allocation");
    tempshead = shead;
  }
  else{
    while(tempshead && tempshead->which < which){
      lastshead = tempshead;
      tempshead = tempshead->next;
    }
    /* if this is supposed to be the first member of the list */
    if( tempshead == shead ){
      tempshead = (shead_type *)D_CALLOC(1,sizeof(shead_type));
      if(!tempshead)fatal("Memory Allocation" );
      /* make sure that the first element is not supposed to just be
         replaced */
      if( shead->which == which ){
        D_FREE(shead);
        shead = 0;
      }
      tempshead->next=shead;
      shead=tempshead;
    }
    /* otherwise we are in the middle somewhere */
    else{
      tempshead = (shead_type *)D_CALLOC(1,sizeof(shead_type));
      if(!tempshead)fatal("Memory Allocation" );
      /* check to see if replacement of the snap is needed */
      if( lastshead->next && lastshead->next->which == which ){
        tempshead->next = lastshead->next->next;
        D_FREE(lastshead->next);
        lastshead->next=tempshead;
      }
      else{
        tempshead->next = lastshead->next;
        lastshead->next = tempshead;
      }
    }
  }


  /* now that the new snap node is inserted, fill the camera parms */
  tempshead->vup.x=camera->vup.x;tempshead->vup.y=camera->vup.y;
  tempshead->vup.z=camera->vup.z;

  tempshead->la.x=camera->la.x;tempshead->la.y=camera->la.y;
  tempshead->la.z=camera->la.z;

  tempshead->lf.x=camera->lf.x;tempshead->lf.y=camera->lf.y;
  tempshead->lf.z=camera->lf.z;

  tempshead->foclength=camera->foclength;
  tempshead->which = which;

  /* that takes care of setting up the node.  time to take the snapshot */
  temphead = head;
  /* go through the whole forest of objects */
  while(temphead){
    tempshead->snap=snapobj(temphead->obj);
    temphead=temphead->next;
  }
}


/****************************************************************************
 *
 *                   Procedure animateit
 *
 * Arguments: from, to: integers
 * Returns: none
 *
 * Action: does the way cool linearly interpolated animation from keyframe from
 *       through keyframe to, using every keyframe in between.
 *
 ****************************************************************************/
void animateit(int from,int to)
{
  int start, end, i;
  float gap;
  char string1[80],string2[80];
  shead_type *frame1, *frame2, interp;
  snap_type *snap1, *snap2, *tempsnap;

  /* make sure that the two keyframes are in proper order */
  if( to < from ){
    start = to;
    to = from;
    from = start;
  }

  /* set up the interpolated snapshot head */
  interp.which=0;
  interp.next=0;

  /* now build the interpolation keyframe list.  This is just as long as
     any of the keyframe lists, so build it by traversing one */
  tempsnap = (snap_type *)D_CALLOC(1,sizeof(snap_type));
  if(!tempsnap)fatal("memory allocation");
  interp.snap = tempsnap;
  snap1 = tempsnap;
  snap2 = shead->snap;
  while(snap2){
    /* if there is another node after this one then add another node to
       the interpolation list */
    if(snap2->next){
      tempsnap = (snap_type *)D_CALLOC(1,sizeof(snap_type));
      if(!tempsnap)fatal("memory allocation");
      snap1->next=tempsnap;
      snap1=snap1->next;
    }
    snap2 = snap2->next;
  }
  start = from;
  /* find the first keyframe */
  frame1 = shead;
  while(frame1 && frame1->which < start ) frame1 = frame1->next;
  /* make sure that the keyframe actually was there */
  if( !frame1 || frame1->which != start ){
    /* display an error message and return */
    sprintf(string1,"%d",from);
    strcpy(string2, " <-no such keyframe");
    strcat(string1,string2);
    display( string1 );
    return;
  }

  /* loop until the first keyframe in the interval is not equal to to */
  while( start < to ){
    /* get the next keyframe */
    frame2 = frame1->next;
    /* make sure that the next frame actually exists */
    if( !frame2 ){
      /* display an error message and return */
      strcpy(string2, "Ran out of keyframes" );
      display( string2 );
      return;
    }
    end = frame2->which;

    /* now loop between the two keyframes */
    for(i=start;i<=end;i++){
      gap = (float)(i-start)/(float)(end-start);
      /* set the snapshot pointers back to the beginning of the lists */
      snap1=frame1->snap;
      snap2=frame2->snap;
      tempsnap = interp.snap;

      /* interpolate the camera parameters */
      interp.lf.x=(1-gap)*frame1->lf.x+gap*frame2->lf.x;
      interp.lf.y=(1-gap)*frame1->lf.y+gap*frame2->lf.y;
      interp.lf.z=(1-gap)*frame1->lf.z+gap*frame2->lf.z;

      interp.la.x=(1-gap)*frame1->la.x+gap*frame2->la.x;
      interp.la.y=(1-gap)*frame1->la.y+gap*frame2->la.y;
      interp.la.z=(1-gap)*frame1->la.z+gap*frame2->la.z;

      interp.vup.x=(1-gap)*frame1->vup.x+gap*frame2->vup.x;
      interp.vup.y=(1-gap)*frame1->vup.y+gap*frame2->vup.y;
      interp.vup.z=(1-gap)*frame1->vup.z+gap*frame2->vup.z;

      interp.foclength=(1-gap)*frame1->foclength+gap*frame2->foclength;

      /* fill in the interpolated snap list */
      while( snap1 ){
        /* do the object parameters */
        tempsnap->xt=(1-gap)*snap1->xt+gap*snap2->xt;
        tempsnap->yt=(1-gap)*snap1->yt+gap*snap2->yt;
        tempsnap->zt=(1-gap)*snap1->zt+gap*snap2->zt;

        tempsnap->xs=(1-gap)*snap1->xs+gap*snap2->xs;
        tempsnap->ys=(1-gap)*snap1->ys+gap*snap2->ys;
        tempsnap->zs=(1-gap)*snap1->zs+gap*snap2->zs;

        tempsnap->xc=(1-gap)*snap1->xc+gap*snap2->xc;
        tempsnap->yc=(1-gap)*snap1->yc+gap*snap2->yc;
        tempsnap->zc=(1-gap)*snap1->zc+gap*snap2->zc;

        tempsnap->xr=(1-gap)*snap1->xr+gap*snap2->xr;
        tempsnap->yr=(1-gap)*snap1->yr+gap*snap2->yr;
        tempsnap->zr=(1-gap)*snap1->zr+gap*snap2->zr;

        snap1 = snap1->next;
        snap2 = snap2->next;
        tempsnap = tempsnap->next;
      }
      /* the interpolated list is done, send it off */
      readsnap(0,&interp);
    }
    /* make this keyframe the new begin frame */
    frame1 = frame2;
    start = end;
  }
}

/****************************************************************************
 *
 *                   Procedure doanimation
 *
 * Arguments: none
 * Returns: none
 *
 * Action: reads in two integers from the user and then animates between
 *   those two keyframes
 *
 ****************************************************************************/
void doanimation(void)
{
  int start, end;

  display("Look in the xterm...");
  printf( "Start animation at which keyframe: ");
  scanf( "%d", &start );
  printf( "End at which keyframe: " );
  scanf( "%d", &end );
  animateit( start, end );
}

/****************************************************************************
 *
 *                   Procedure getkeyframe
 *
 * Arguments: none
 * Returns: none
 *
 * Action: gets an integer from the user and then reads in that keyframe
 *
 ****************************************************************************/
void getkeyframe(void)
{
  int which;

  display("Look in the xterm...");
  printf( "Read in which keyframe: ");
  scanf( "%d", &which );
  readsnap( which, shead );
}

/****************************************************************************
 *
 *                   Procedure takekeyframe
 *
 * Arguments: none
 * Returns: none
 *
 * Action: gets an integer from the user and then writes that keyframe
 *
 ****************************************************************************/
void takekeyframe(void)
{
  int which;

  display("Look in the xterm...");
  printf( "Save which keyframe: ");
  scanf( "%d", &which );
  takesnap( which );
}

