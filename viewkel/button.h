/*******************************************************

Copyright (C) 1999 Greg Landrum
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
#define TOGGLE 0     /* associated with booleans (on or off) */
#define FLOATPARM 1  /* buttons that read in floating point parameters */
#define FUNCTION 2   /* calls a function */
#define MODE 3       /* for modal type things */
#define STRINGPARM 4 /* buttons that read in strings */
#define RADIO 10     /* has child windows */

#define MAX_MODES 10 /* any more than this is too many */
#define MAX_LINES 10 /* the maximum number of lines in a multi-line title */
#define MAX_ARGS 10  /* the maximum number of arguments a function can take */

#define TEXTOFF 10
#define BUTTONOFF 10

/* these are the maps for all of the buttons */
#include "floatparm.xbm"
#include "function.xbm"
#include "toggleon.xbm"
#include "toggleoff.xbm"
#include "modal.xbm"

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

union guts_def {
  char *boolean;
  mode_type mode;
  float *floatparm;
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
  Pixmap pix1, pix2; /* two pixmaps in case the button has two settings */
};
