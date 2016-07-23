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



/****************************************************************************
 *
 *          help.c
 *
 * This file contains the stuff needed for the program to be able
 *  to provide some help for it's use.
 *
 *  written by greg Landrum September 1994
****************************************************************************/
#include "viewkel.h"

/******
  this is the help text... add to it as you will, just make sure that
  the last entry has an & as the first character.
*******/
char help_lines[][80] = {"\tHelp for Viewkel","  ",
"Commands are not case sensitive.",
" ",
"File Commands:",
"  Graph: read in some graph data",
"  Props: read in some average properties data (DOS or COOP)",
"  Bands: read in band structure data",
"  Walsh: read in data for a Walsh diagram",
"    FMO: read in data for an FMO diagram",
"Format Commands:",
"  Xtics: toggles tic marks along the X axis.",
"  Ytics: toggles tic marks along the Y axis.",
"  Xlegend: enter a legend for the X axis",
"  Ylegend: enter a legend for the Y axis",
"  Title:   enter a title for the graph",
"  Xmin:  enter a new value for the minimum value displayed on the X axis",
"  Xmax:  enter a new value for the maximum value displayed on the X axis",
"  Ymin:  enter a new value for the minimum value displayed on the Y axis",
"  Ymax:  enter a new value for the maximum value displayed on the Y axis",
"  Xscale:  enter a scaling factor for the X direction",
"  Yscale:  enter a scaling factor for the Y direction",
"  Xtrans:  Move the active graph along the X direction",
"  Ytrans:  Move the active graph along the Y direction",
"  Curve n: toggles the display of curve n in the data.",
"           For example:  entering",
"             viewkel>  curve 3",
"	    would toggle the display status of curve 3.",
"  Integ n: similar to the curve command, toggles integration curves.",
"  Integscale:  Use the integration to provide the X scale of the graph.",
"  Fermi:  toggles display of the Fermi Energy on the graph",
"  Change Fermi: enter a new value for the Fermi Energy",
"  Fill: toggles filling of projected DOS's (not shown on screen)",
"  Linestyle n d: Sets the linestyle of curve n to style d",
"  Integstyle n d: Sets the linestyle of integration n to style d",
"FMO commands:",
"  Show all: show all electrons in the FMO diagram.",
"  Show homo: show only electrons occupying the HOMO in the FMO diagram.",
"  Show none: show no electrons in the FMO diagram.",
"  Left frag: change the number of the fragment shown at the left of the",
"             FMO diagram.  Set this to -1 to not show any left fragment.",
" Right frag: change the number of the fragment shown at the right of the",
"             FMO diagram.  Set this to -1 to not show any right fragment.",
"Misc. Commands:",
"  Kill:  kills the current graph",
"  Purge:  kills all open graphs",
"  Quit:  quit the program",
" ",
" ",
"&"};



/****************************************************************************
 *
 *                   Procedure show_help
 *
 * Arguments: none
 *            
 * Returns: none
 *
 * Action:  Prints out the help message contained in the variable 'help_strings
 *  
 ****************************************************************************/
void show_help(void)
{
  int i;

  for(i=0;help_lines[i][0] != '&'; i++){
    printf("%s\n",help_lines[i]);
    /* only print out 20 lines at a time, then wait for a carriage return */
    if( i && !(i%20) ){
      printf(">>>> hit enter to continue <<<<<");
      fflush(stdin);fgetc(stdin); fflush(stdin);
    }
  }
}
