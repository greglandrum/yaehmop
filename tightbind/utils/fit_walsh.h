/************************************************************************
  This is the include file for the property fitting programs

   Created by greg Landrum June 1994
************************************************************************/
#include <stdio.h>
#include <math.h>
#include <string.h>


/******
  These are used by the fileio routines.
*******/
#define FATAL 0
#define ERROR 1
#define IGNORE 2


/******
  this is the data type that will be used to store the walsh data
  within the program.
******/
typedef struct walsh_point_type_def{
  real energy;
  real *symmetries;
  struct walsh_point_type_def *next;
} walsh_point_type;
  
