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
/*******************************************************
*
*  This file is part of yaehmop.
*
*  There are, of course, no warranties at all on this program.
*
********************************************************************/
#include "viewkel.h"
#include "matrix_ops.h"

/****************************

  These are operations on matrices and vectors adapted from:
     2d and 3d Vector C Library
     by Andrew Glassner
     from "Graphics Gems", Academic Press, 1990

*****************************/

/***
  Recent Edit History

   26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)

***/


/* returns squared length of input vector */
double V3SquaredLength(Vector3 *a)
{
  return((a->x * a->x)+(a->y * a->y)+(a->z * a->z));
}

/* returns length of input vector */
double V3Length(Vector3 *a)
{
  return(sqrt(V3SquaredLength(a)));
}

/* negates the input vector and returns it */
Vector3 *V3Negate(Vector3 *v)
{
  v->x = -v->x;  v->y = -v->y;  v->z = -v->z;
  return(v);
}

/* normalizes the input vector and returns it */
Vector3 *V3Normalize(  Vector3 *v)
{
  double len = V3Length(v);
  if (len != 0.0) { v->x /= len;v->y /= len; v->z /= len; }
  return(v);
}

/* multiplies the input vector by the constant */
Vector3 *V3ConstantScale(Vector3 *v, double scale)
{
  v->x *= scale;  v->y *= scale;  v->z *= scale;
  return(v);
}

/* scales the input vector to the new length and returns it */
Vector3 *V3Scale(Vector3 *v, double newlen)
{
  double len = V3Length(v);
  if (len != 0.0) {
    v->x *= newlen/len;   v->y *= newlen/len;  v->z *= newlen/len;
  }
  return(v);
}


/* return vector sum c = a+b */
Vector3 *V3Add(Vector3 *a, Vector3 *b, Vector3 *c)
{
  c->x = a->x+b->x;  c->y = a->y+b->y;  c->z = a->z+b->z;
  return(c);
}

/* return vector difference c = a-b */
Vector3 *V3Sub(Vector3 *a, Vector3 *b,Vector3 *c)
{
  c->x = a->x-b->x;  c->y = a->y-b->y;  c->z = a->z-b->z;
  return(c);
}

/* return the dot product of vectors a and b */
double V3Dot(Vector3 *a,Vector3 *b)
{
  return((a->x*b->x)+(a->y*b->y)+(a->z*b->z));
}

/* linearly interpolate between vectors by an amount alpha */
/* and return the resulting vector. */
/* When alpha=0, result=lo.  When alpha=1, result=hi. */
Vector3 *V3Lerp(Vector3 *lo, Vector3 *hi, double alpha, Vector3 *result)
{
  result->x = LERP(alpha, lo->x, hi->x);
  result->y = LERP(alpha, lo->y, hi->y);
  result->z = LERP(alpha, lo->z, hi->z);
  return(result);
}

/* make a linear combination of two vectors and return the result. */
/* result = (a * ascl) + (b * bscl) */
Vector3 *V3Combine (Vector3 *a, Vector3 *b, Vector3 *result,
                    double ascl,double bscl)
{
  result->x = (ascl * a->x) + (bscl * b->x);
  result->y = (ascl * a->y) + (bscl * b->y);
  result->z = (ascl * a->z) + (bscl * b->z);
  return(result);
}


/* multiply two vectors together component-wise and return the result */
Vector3 *V3Mul (Vector3 *a, Vector3 *b, Vector3 *result)
{
  result->x = a->x * b->x;
  result->y = a->y * b->y;
  result->z = a->z * b->z;
  return(result);
}

/* return the distance between two points */
double V3DistanceBetween2Points(Point3 *a, Point3 *b)
{
  double dx = a->x - b->x;
  double dy = a->y - b->y;
  double dz = a->z - b->z;
  return(sqrt((dx*dx)+(dy*dy)+(dz*dz)));
}

/* return the angle between 3 points (added by gL) */
double V3AngleBetween3Points(Point3 *a, Point3 *b, Point3 *c)
{
  double dotprod,lenprod;
  Vector3 V1,V2;
  V1.x = a->x-b->x;V1.y = a->y-b->y;V1.z = a->z-b->z;
  V2.x = c->x-b->x;V2.y = c->y-b->y;V2.z = c->z-b->z;

  /*********

    we're using:

      q = acos( (a . b) / (|a| |b|) )

  **********/
  lenprod = V3Length(&V1)*V3Length(&V2);
  dotprod = V3Dot(&V1,&V2);
  return(acos((dotprod/lenprod)));
}



/* return the dihedral angle between 4 points (added by gL) */
double V3DihedralAngle(Point3 *a, Point3 *b, Point3 *c, Point3 *d)
{
  double dotprod;
  Vector3 V1,V2,V3,V1xV2,V3xV2;
  V1.x = a->x-b->x;V1.y = a->y-b->y;V1.z = a->z-b->z;
  V2.x = c->x-b->x;V2.y = c->y-b->y;V2.z = c->z-b->z;
  V3.x = d->x-c->x;V3.y = d->y-c->y;V3.z = d->z-c->z;

  /* take the cross products */
  V3Cross(&V1,&V2,&V1xV2);
  V3Cross(&V3,&V2,&V3xV2);

  /* normalize them to make life easier */
  V3Normalize(&V1xV2);
  V3Normalize(&V3xV2);

  /*********
    now figure out the angle between them
    we're using:

      q = acos( (a . b) )
    since a and b are normalized

  **********/
  dotprod = V3Dot(&V1xV2,&V3xV2);
  return(acos(dotprod));
}

/* return the cross product c = a cross b */
Vector3 *V3Cross(Vector3 *a, Vector3 *b, Vector3 *c)
{
  c->x = (a->y*b->z) - (a->z*b->y);
  c->y = (a->z*b->x) - (a->x*b->z);
  c->z = (a->x*b->y) - (a->y*b->x);
  return(c);
}


void *V3MulPointByProjMatrix(Point3 *pin,Point3 *pout)
{
  double w;
  pout->x = (pin->x * stack->matrix[0][0]) +
    (pin->y * stack->matrix[1][0]) +
      (pin->z * stack->matrix[2][0]) + stack->matrix[3][0];
  pout->y = (pin->x * stack->matrix[0][1]) +
    (pin->y * stack->matrix[1][1]) +
      (pin->z * stack->matrix[2][1]) + stack->matrix[3][1];
  pout->z = (pin->x * stack->matrix[0][2]) +
    (pin->y * stack->matrix[1][2]) +
      (pin->z * stack->matrix[2][2]) + stack->matrix[3][2];
  w =    (pin->x * stack->matrix[0][3]) +
    (pin->y * stack->matrix[1][3]) +
      (pin->z * stack->matrix[2][3]) + stack->matrix[3][3];
  if (w != 0.0) { pout->x /= w;  pout->y /= w; }
  /*
     fprintf(stderr,"  (% 6.4lf, % 6.4lf, % 6.4lf)  % 6.4lf ->\
     \n\t (% 6.4lf, % 6.4lf, % 6.4lf)\n",
     pin->x, pin->y, pin->z, w,
     pout->x, pout->y,
     pout->z);
     */

}
