/*******************************************************
*      Copyright (C) 1998, 1999 Greg Landrum
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
typedef struct triangle_list_def{
  triangle_type tri;
  point_type bmax,bmin;
  struct triangle_list_def *next;
} triangle_list_type;


/* check to see if a point falls in a bounding region */
#define POINT_IS_INSIDE(_p_,_max_,_min_)\
(_p_.x - _max_.x < BOX_TOL && \
_p_.x - _min_.x > BOX_TOL && \
_p_.y - _max_.y < BOX_TOL && \
_p_.y - _min_.y > BOX_TOL  && \
_p_.z - _max_.z < BOX_TOL && \
_p_.z - _min_.z > BOX_TOL  )

#endif
