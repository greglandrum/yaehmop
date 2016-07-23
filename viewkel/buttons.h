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

#ifndef _BUTTONS_

#define _BUTTONS_

#ifndef _MY_DEFINES_
#include "defines.h"
#endif

#ifdef MAC_GRAPHICS
#ifndef _MAC_DEFINES_
#include "Mac_defines.h"
#endif
#endif

/*
        this file contains all the data structures and defines
        needed to do my kind of buttons
*/

#define TOGGLE 0          /* associated with booleans (on or off) */
#define FLOATPARM 1       /* buttons that read in floating point parameters */
#define FUNCTION 2        /* calls a function */
#define MODE 3            /* for modal type things */
#define STRINGPARM 4      /* buttons that read in strings */
#define INTPARM 5         /* buttons that read in integer parameters */
#define CURVETOGGLE 6     /* for toggling curves and drawing styles */
#define MULTISTRINGPARM 7 /* multiple strings */

typedef struct button_def button_type;
typedef union guts_def button_guts;

/* this is all the spoo needed for mode buttons */
typedef struct {
  int *controlling_variable;
  char names[MAX_MODES][80];
  int num_modes;
} mode_type;

/* this is for having multi-line strings */
typedef struct {
  char *lines[MAX_LINES];
  char num_lines;
} string_parm_type;

/* this is to allow the function buttons to have arguments */
typedef struct {
  /* the function itself */
  void (*func)();

  /************
    The arguments themselves.
    NOTE: to make this general, all arguments should be cast as char *.
       they can be converted back to the appropriate type in the function
       which is called.
  ************/
  int num_args;
  char *args[MAX_ARGS];
} func_struct_type;

typedef struct {
  char *boolean;
  char *style;
  char *fill;
} curve_toggle_type;

union guts_def {
  char *boolean;
  curve_toggle_type curve;
  mode_type mode;
  float *floatparm;
  int *intparm;
  string_parm_type stringparm;
  func_struct_type function;
  button_type *children;
};

struct button_def {
  char type;
  char string[80];
  int xdim, ydim;
  button_guts guts;
  button_type *next;
#ifdef X_GRAPHICS
  /* multiple pixmaps in case the button has multiple settings */
  Pixmap pix1, pix2, pix3;
#endif

#ifdef MAC_GRAPHICS
  MenuHandle whichMenu;
  short menuID;
  short itemID;
#endif
};

#endif
