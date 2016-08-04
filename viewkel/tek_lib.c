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

/*************************************************

  This file contains stuff for dealing with Tektronix terminal emulators

    It was liberally hacked from the gnuplot vttek terminal driver.

      Created by greg Landrum  September, 1994

**************************************************/
#include "viewkel.h"

extern FILE *Tek_file;

#define max(a,b) ( a>b ? a : b )

static char *vt_linetype = "`a`abcdhijkl" ;
static int last_vt_linetype = 0;


#define TEK_XLAST (TEK_XMAX - 1)
#define TEK_YLAST (TEK_YMAX - 1)

#define HX 0x20                /* bit pattern to OR over 5-bit data */
#define HY 0x20
#define LX 0x40
#define LY 0x60

#define LOWER5 31
#define UPPER5 (31<<5)

void TEK_init(void)
{
  fprintf(Tek_file,"\033[?38h");
  fflush(Tek_file);
  /* sleep 1 second to allow screen time to clear on some terminals */
  sleep(1);
}

void TEK_reset(void)
{
  fprintf(Tek_file,"\033[?38l");
  fflush(Tek_file);
  /* sleep 1 second to allow screen time to clear on some terminals */
  sleep(1);
}


void TEK_graphics(void)
{
  fprintf(Tek_file,"\033\014");
  /*                   1
                       1. clear screen
                       */
  (void) fflush(Tek_file);
  sleep(1);
  /* sleep 1 second to allow screen time to clear on real
     tektronix terminals */
}

void TEK_text(void)
{
  TEK_move(0,12);
  fprintf(Tek_file,"\037");
  /*                   1
                       1. into alphanumerics
                       */
}


void TEK_move( unsigned int x,unsigned int y)
{
  (void) putc('\035', Tek_file);        /* into graphics */
  TEK_vector(x,y);
}


void TEK_vector(unsigned int x,unsigned int y)
{
  (void) putc((HY | (y & UPPER5)>>5), Tek_file);
  (void) putc((LY | (y & LOWER5)), Tek_file);
  (void) putc((HX | (x & UPPER5)>>5), Tek_file);
  (void) putc((LX | (x & LOWER5)), Tek_file);
}



/******************
  These are the linetypes for VT-type terminals in tektronix emulator mode:

  `=solid
  a=fine dots
  b=short dashes
  c=dash dot
  d=long dash dot
  h=bold solid
  i=bold fine dots
  j=bold short dashes,
  k=bold dash dot
  l=bold long dash dot
  ********************/
void TEK_linetype(int linetype)
{
  if (linetype >= 10)
    linetype %= 10;
  fprintf(Tek_file,"\033%c",vt_linetype[linetype+2]);
  last_vt_linetype = linetype;
}

void TEK_put_text(unsigned int x,unsigned int y,char *str)
{
  int linetype;
  linetype = last_vt_linetype;
  TEK_linetype(0);
  TEK_move(x,y-11);
  fprintf(Tek_file,"\037%s\n",str);
  TEK_linetype(linetype);
}

void TEK_center_text(unsigned int x,unsigned int y,char *str)
{
  TEK_put_text(x-strlen(str)*TEK_HCHAR/2,y,str);
}

void TEK_right_text(unsigned int x,unsigned int y,char *str)
{
  TEK_put_text(x-strlen(str)*TEK_HCHAR,y,str);
}
