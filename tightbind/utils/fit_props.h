/************************************************************************
  This is the include file for the property fitting programs

   Created by greg Landrum March 1994
************************************************************************/
#include <stdio.h>
#include <math.h>

#ifdef USING_THE_MAC
#include <SIOUX.h>
extern FILE *choose_mac_file(char *,char);
#include "Mac_Fopen.h"
#endif

/******
  These are used by the fileio routines.
*******/
#define FATAL 0
#define ERROR 1
#define IGNORE 2

/********
  The program will compile to use doubles instead of floats to represent everything.
  If for some reason you want to use floats, then you must put -DUSE_FLOATS in the
  CFLAGS field of the makefile
*********/
#ifndef USE_FLOATS
typedef double real;
#else
typedef float real;
#endif

#ifndef USE_BZERO
#define bzero(a,b) ( memset((void *)(a),0,(b)) )
#define bcopy(a,b,c) ( memcpy((void *)(b),(const void *)(a),(c)) )
#endif

/* this is the step size used to generate the smoothed data (in eV) */
#define ENERGY_STEP .01

/* this is the broadening factor */
#define BROADENING .01

#define PI 3.141592653589793
#ifndef M_PI
#define M_PI 3.141592653589793
#endif

#define MAX_STR_LEN 2048

/******
  this is the data type that will be used to store the properties data
  within the program.
******/
typedef struct{
  real height;
  real energy;
} point_type;
  
