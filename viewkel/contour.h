/******************************

  contour.h

  This file contains the definitions required to draw contour plots

  Nabbed from gnuplot 3.5 by gL, June 1996

*******************************/

/***
  Recent modification history
 
  26.09.98 gL:
    added a prev pointer to iso_curve_type to make
      that a doubly linked list (to enable use of
      symmetry without sorting in the evaluation of
      MO planes).
***/
#ifndef _CONTOUR_
#define _CONTOUR_

#ifndef _MY_DEFINES_
#include "defines.h"
#endif

#ifndef _BASIC_OBJECTS_
#include "basic_objects.h"
#endif


struct curve_points {
	struct curve_points *next_cp;	/* pointer to next plot in linked list */
 	int p_max;					/* how many points are allocated */
	int p_count;					/* count of points in points */
	point_type *points;
};

struct gnuplot_contours {
	struct gnuplot_contours *next;
	point_type *coords;
 	char isNewLevel;
 	char label[12];
	int num_pts;
};

typedef struct gnuplot_contours gnuplot_contour_type;

typedef struct iso_curve {
	struct iso_curve *next,*prev;
 	int p_max;					/* how many points are allocated */
	int p_count;					/* count of points in points */
	point_type *points;
} iso_curve_type;

#endif
