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
 *                   Procedure draw_param_surf
 *
 * Arguments: surf: pointer to param_surf_type
 *          center: pointer to point_type
 *            mode: char
 *
 * Returns: none
 *
 * Action: This draws a parametric surface.
 *      if mode == ALL then the entire surface will be drawn.
 *      if mode == BACK then all the surface points will be transformed,
 *            but only those with a z coordinate > 'center->z will be drawn,
 *            the rest are stored.
 *      if mode == FRONT then the stored surface points in front of the center
 *            will be drawn.
 *
 ****************************************************************************/
void draw_param_surf(surf,center,mode)
  param_surf_type *surf;
  point_type *center;
  char mode;
{
  static int num_p_allocated=0;
  static point_type *temp_points=0;
  static XPoint *back_xpoints1=0,*back_xpoints2=0;
  static XPoint *front_xpoints1=0,*front_xpoints2=0;
  static XSegment *xsegs=0;
  static front_num1,front_num2;

  int i,j;
  int num_p,back_num1,back_num2,num_segs;
  point_type *the_point;

  num_p = surf->samples1*surf->samples2;
  /* Check to see if we need to get space for the temporary points array */
  if( num_p_allocated < num_p ){
    /*******

      if we already got space, free it up
      I'm doing this instead of using D_REALLOC because I don't trust
      D_REALLOC.

      *******/
    if( temp_points ){
      D_FREE(temp_points);
      D_FREE(front_xpoints1);
      D_FREE(front_xpoints2);
      D_FREE(back_xpoints1);
      D_FREE(back_xpoints2);
      D_FREE(xsegs);
    }

    temp_points = (point_type *)D_CALLOC(num_p,sizeof(point_type));
    front_xpoints1 = (XPoint *)D_CALLOC(num_p,sizeof(XPoint));
    front_xpoints2 = (XPoint *)D_CALLOC(num_p,sizeof(XPoint));
    back_xpoints1 = (XPoint *)D_CALLOC(num_p,sizeof(XPoint));
    back_xpoints2 = (XPoint *)D_CALLOC(num_p,sizeof(XPoint));
    xsegs = (XSegment *)D_CALLOC(2000,sizeof(XSegment));
    if( !temp_points || !back_xpoints2 || !xsegs ){
      fatal("Can't allocate temporary storage in draw_param_surf.");
    }
    num_p_allocated = num_p;
  }

  /* okay, copy the points for this surface so that we can transform them */
  bcopy(surf->points,temp_points,num_p*sizeof(point_type));

  /***************

    Now loop through all the points and transform and draw them

    Since we know the total number of points, and they are in cartesian
    coordinates at this point, this can just be done by treating the point
    array as 1 dimensional.

    ***************/
  the_point = temp_points;
  num_segs = back_num1 = back_num2 = 0;

  if( mode != FRONT ){
    front_num1 = front_num2 = 0;
    for(i=0;i<num_p;i++){
      transform(the_point);

      /* now do the homogeneous coordinate transform */
      the_point->x /= the_point->z;
      the_point->y /= the_point->z;

      /******

        check the color.
        To add more colors, just change this case
        statement and get more xpoints arrays.

      *******/
      switch(mode){
      case ALL:
        switch(surf->colors[i]){
        case 1:
          front_xpoints1[front_num1].x = (int)the_point->x;
          front_xpoints1[front_num1].y = (int)the_point->y;
          front_num1++;
          break;
        case -1:
          front_xpoints2[front_num2].x = (int)the_point->x;
          front_xpoints2[front_num2].y = (int)the_point->y;
          front_num2++;
          break;
        default:
          error("Bogus color in draw_param_surf.");
        }
        break;
      case FRONT:
      case BACK:
        switch(surf->colors[i]){
        case 1:
          if( the_point->z > center->z ){
            back_xpoints1[back_num1].x = (int)the_point->x;
            back_xpoints1[back_num1].y = (int)the_point->y;
            back_num1++;
          }else{
            front_xpoints1[front_num1].x = (int)the_point->x;
            front_xpoints1[front_num1].y = (int)the_point->y;
            front_num1++;
          }
          break;
        case -1:
          if( the_point->z > center->z ){
            back_xpoints2[back_num2].x = (int)the_point->x;
            back_xpoints2[back_num2].y = (int)the_point->y;
            back_num2++;
          }else{
            front_xpoints2[front_num2].x = (int)the_point->x;
            front_xpoints2[front_num2].y = (int)the_point->y;
            front_num2++;
          }
          break;
        default:
          error("Bogus color in draw_param_surf.");
        }
        break;
      default:
        error("Bogus mode passed to draw_param_surf. This is a bug.");
      }
      the_point++;
    }
  }

  /* draw the surface as mere points for a while */
  if( mode == FRONT || mode == ALL ){
    if( front_num1 ){
      if( !draw_connectors ){
        XDrawPoints(disp,gpix,graphgc,front_xpoints1,front_num1,CoordModeOrigin);
      }
      else{
        XDrawLines(disp,gpix,graphgc,front_xpoints1,front_num1,CoordModeOrigin);
      }
    }
    if( front_num2 ){
      if( !draw_connectors ){
        XDrawPoints(disp,gpix,colorgc,front_xpoints2,front_num2,CoordModeOrigin);
      }
      else{
        XDrawLines(disp,gpix,colorgc,front_xpoints2,front_num2,CoordModeOrigin);
      }
    }
  }
  else{
    if( back_num1 ){
      if( !draw_connectors ){
        XDrawPoints(disp,gpix,graphgc,back_xpoints1,back_num1,CoordModeOrigin);
      }
      else{
        XDrawLines(disp,gpix,graphgc,back_xpoints1,back_num1,CoordModeOrigin);
      }
    }
    if( back_num2 ){
      if( !draw_connectors ){
        XDrawPoints(disp,gpix,colorgc,back_xpoints2,back_num2,CoordModeOrigin);
      }
      else{
        XDrawLines(disp,gpix,colorgc,back_xpoints2,back_num2,CoordModeOrigin);
      }
    }
  }
}

/****
  A hack for drawing surfaces with atoms with the painter's algorithm:
    Call draw_param_surf, WITH the atomic z coordinate
     have it draw the parts of the surface that are behind the atom, and save
      the rest in an array.
    Draw the atom
    Call draw_param_surf, and have it draw the parts in front of the atom.

  This should work!
****/

