/*******************************************************
*      Copyright (C) 1995, 1998,1999  Greg Landrum
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

  basic_objects.h

  This file contains the basic objects (points, lines, etc.)
   needed by viewkel.

  Created by gL June 1995

*******************************/

/*******
  Recent modification history

********/

#ifndef _BASIC_OBJECTS_

#define _BASIC_OBJECTS_

#ifndef _MY_DEFINES_
#include "defines.h"
#endif

/******** 

  this is the node for the matrix stack

********/
typedef struct matrix_type_def{
  float matrix[DIM][DIM];
  struct matrix_type_def *next;
} matrix_type;




/**********

  a 3D point

**********/
typedef struct {
  float x,y,z;
} point_type;



/*********** 

  a 2D point 

*************/
typedef struct {
  float x,y;
} point_type2D;
  

/************

  used to draw lines between atoms

************/
typedef struct{
  long int end1,end2;
  char custom,drawn,tube,breaking,arrow,dashed;
  float length;
  int thickness;
  int type;
} line_type;


/************

  a triangle with vertices as tabs (used to draw polygonal surfaces)
  

*************/
typedef struct{
  int vertices[3];
  point_type center,normal;
  char perp;
  int color;
} triangle_type;


/***********

  triangle vertex, with position and normal

************/
typedef struct vertex {		   
    point_type position, normal;	
} vertex_type;

/************

  a triangle with real vertices specified (used to draw polygonal surfaces)  

*************/
typedef struct{
  point_type vertices[3];
  point_type center,normal;
  int color;
} explicit_triangle_type;



/* a particle */
typedef struct particle_type_def{
  char on_surf;
  char phase;
  int num_neighbors;
  int num;
  float rad;
  float val;
  point_type loc;
  struct particle_type_def *next;
  struct particle_type_def *neighbors[6];
} particle_type;


#endif
