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

/***
  Recent Edit History:

  04.09.98 gL:
    added ButtonReleaseMask to the list of events taken by gwin.
    this is needed for click and drag selection to work properly.

   26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)

***/

#include "viewkel.h"
#ifdef X_GRAPHICS
#define NAME "Viewkel!!!"
#define SMALLFONTNAME "6x10"
#define BIGFONTNAME "9x15"

#include "gray0.xbm"
#include "gray1.xbm"
#include "gray2.xbm"
#include "gray3.xbm"
#include "gray4.xbm"
#include "gray5.xbm"
#include "gray6.xbm"
#include "gray7.xbm"

/****************************************************************************
 *
 *  Procedure X_initgraphics
 *
 *   Arguments: xsize,ysize: ints
 *   Return Value: none
 *
 *  Actions: This does whatever is needed to prepare the system for display
 *      of graphics (opening windows, declaring fonts, etc)
 *
 *
 *****************************************************************************/
void X_initgraphics(int xsize,int ysize)
{
  int width, height;
  XGCValues xgcv;
  XSetWindowAttributes xswa;
  char *name;
  XColor thecolor,foocolor;
  Colormap colormap;


  /* if the windows have not already been opened, open them now */
  /* set the program's name */
  name = (char *)D_CALLOC(80,sizeof(char));
  if(!name)fatal("Memory Allocation.");
  (void)strcpy( name, NAME );

  /* first open the display (and make sure that it is open )*/
  if( !(disp = XOpenDisplay(getenv("DISPLAY"))))
    fatal( "Cannot open display" );

  /* find some information about the display */
  screen = DefaultScreen(disp);
  root = RootWindow( disp, screen );

#ifdef DEBUG
  XSynchronize(disp,True);
#endif
  /* get the foreground and background colors */
  bcolor = WhitePixel( disp, screen );
  fcolor = BlackPixel( disp, screen );

  /* check to see if there are parameters being passed in */

  if(!xsize || !ysize){
    width = GWINWIDTH+BUTWIDTH;
    height = GWINHEIGHT;
    g_xmax = GWINWIDTH;
    g_ymax = GWINHEIGHT;

  }
  else{
    width = xsize+BUTWIDTH;
    height = ysize;
    g_xmax = width;
    g_ymax = height;
  }
  g_ymax = height;


  /* Create the windows, and make sure they were created */
  graph_win = XCreateSimpleWindow(disp,root,0,0,
                                  width,height,1,fcolor,bcolor );
  /* These are the events that the window will have to deal with */
  xswa.event_mask = ExposureMask | StructureNotifyMask;
  /* Change the event mask for the window */
  XChangeWindowAttributes( disp, graph_win, CWEventMask, &xswa );
  XSetStandardProperties( disp,graph_win,name,name,None,
                         0,0,0);
  /* change the cursor */
  XDefineCursor(disp,graph_win,XCreateFontCursor(disp,34));

  /* set up the graphics contexts */
  xgcv.background = bcolor;
  xgcv.foreground = fcolor;
  xgcv.function = GXcopy;
  graphgc = XCreateGC( disp, graph_win,
                      GCForeground|GCBackground|GCFunction, &xgcv );
  XSetGraphicsExposures( disp, graphgc, 1 );
  XSetLineAttributes( disp, graphgc, 0, LineSolid, CapRound, JoinRound);


  /* try and get another color (ooooooo) */
  screen_depth = DefaultDepth(disp,screen);
#if 0
  colormap = DefaultColormap(disp,screen);
  if( screen_depth > 1 ){
    thecolor.red = 50000;
    thecolor.green = 50000;
    thecolor.blue = 50000;
    if(!XAllocColor(disp,colormap,&thecolor))fatal("Can't get a color");
    xgcv.foreground = thecolor.pixel;
    xgcv.background = bcolor;
    colorgc = XCreateGC( disp, graph_win,
                      GCForeground|GCBackground|GCFunction, &xgcv );
    thecolor.red = 25000;
    thecolor.green = 25000;
    thecolor.blue = 25000;
    if(!XAllocColor(disp,colormap,&thecolor))fatal("Can't get a color");
    xgcv.foreground = thecolor.pixel;
    xgcv.background = bcolor;
    graygc[0] = XCreateGC( disp, graph_win,
                      GCForeground|GCBackground|GCFunction, &xgcv );
    thecolor.red = 50000;
    thecolor.green = 50000;
    thecolor.blue = 50000;
    if(!XAllocColor(disp,colormap,&thecolor))fatal("Can't get a color");
    xgcv.foreground = thecolor.pixel;
    xgcv.background = bcolor;
    graygc[1] = XCreateGC( disp, graph_win,
                      GCForeground|GCBackground|GCFunction, &xgcv );
    thecolor.red = 75000;
    thecolor.green = 75000;
    thecolor.blue = 75000;
    if(!XAllocColor(disp,colormap,&thecolor))fatal("Can't get a color");
    xgcv.foreground = thecolor.pixel;
    xgcv.background = bcolor;
    graygc[2] = XCreateGC( disp, graph_win,
                      GCForeground|GCBackground|GCFunction, &xgcv );
    thecolor.red = 95000;
    thecolor.green = 95000;
    thecolor.blue = 95000;
    if(!XAllocColor(disp,colormap,&thecolor)){
      error("Can't get a color");
      graygc[3] = graygc[2];
    }
    else{
      xgcv.foreground = thecolor.pixel;
      xgcv.background = bcolor;
      graygc[3] = XCreateGC( disp, graph_win,
                            GCForeground|GCBackground|GCFunction, &xgcv );
    }
  }
  else{
    printf("No color here.....");
    colorgc = XCreateGC( disp, graph_win,
                      GCForeground|GCBackground|GCFunction, &xgcv );

  }
#endif


  /* set up the stipple gc's to use to fill the atoms */
  xgcv.background = bcolor;
  xgcv.foreground = fcolor;

  xgcv.fill_style = FillOpaqueStippled;
  xgcv.stipple = XCreateBitmapFromData(disp,graph_win,gray0_bits,
                                             gray0_width,gray0_height);
  graygc[0] = XCreateGC( disp, graph_win,
                        GCForeground|GCBackground|GCFunction|GCStipple|GCFillStyle,
                        &xgcv );

  xgcv.stipple = XCreatePixmapFromBitmapData(disp,graph_win,gray1_bits,
                                             gray1_width,gray1_height,fcolor,
                                             bcolor,1);
  graygc[1] = XCreateGC( disp, graph_win,
                        GCForeground|GCBackground|GCFunction|GCStipple|GCFillStyle,
                        &xgcv );

  xgcv.stipple = XCreatePixmapFromBitmapData(disp,graph_win,gray2_bits,
                                             gray2_width,gray2_height,fcolor,
                                             bcolor,1);
  graygc[2] = XCreateGC( disp, graph_win,
                        GCForeground|GCBackground|GCFunction|GCStipple|GCFillStyle,
                        &xgcv );

  xgcv.stipple = XCreatePixmapFromBitmapData(disp,graph_win,gray3_bits,
                                             gray3_width,gray3_height,fcolor,
                                             bcolor,1);
  graygc[3] = XCreateGC( disp, graph_win,
                        GCForeground|GCBackground|GCFunction|GCStipple|GCFillStyle,
                        &xgcv );

  xgcv.stipple = XCreatePixmapFromBitmapData(disp,graph_win,gray4_bits,
                                             gray4_width,gray4_height,fcolor,
                                             bcolor,1);
  graygc[4] = XCreateGC( disp, graph_win,
                        GCForeground|GCBackground|GCFunction|GCStipple|GCFillStyle,
                        &xgcv );

  xgcv.stipple = XCreatePixmapFromBitmapData(disp,graph_win,gray5_bits,
                                             gray5_width,gray5_height,fcolor,
                                             bcolor,1);
  graygc[5] = XCreateGC( disp, graph_win,
                        GCForeground|GCBackground|GCFunction|GCStipple|GCFillStyle,
                        &xgcv );

  xgcv.stipple = XCreatePixmapFromBitmapData(disp,graph_win,gray6_bits,
                                             gray6_width,gray6_height,fcolor,
                                             bcolor,1);
  graygc[6] = XCreateGC( disp, graph_win,
                        GCForeground|GCBackground|GCFunction|GCStipple|GCFillStyle,
                        &xgcv );

  xgcv.stipple = XCreatePixmapFromBitmapData(disp,graph_win,gray7_bits,
                                             gray7_width,gray7_height,fcolor,
                                             bcolor,1);
  graygc[7] = XCreateGC( disp, graph_win,
                        GCForeground|GCBackground|GCFunction|GCStipple|GCFillStyle,
                        &xgcv );

  xgcv.background=fcolor;
  xgcv.foreground=bcolor;
  blackgc = XCreateGC( disp, graph_win,
                      GCForeground|GCBackground|GCFunction, &xgcv );

  widegc = XCreateGC( disp, graph_win,
                      GCForeground|GCBackground|GCFunction, &xgcv );
  XSetLineAttributes( disp, widegc, 5, LineSolid, CapRound, JoinRound);

colormap = DefaultColormap(disp,screen);
if(XAllocNamedColor(disp,colormap,"red",&thecolor,&foocolor)){
xgcv.foreground = thecolor.pixel;
xgcv.background = bcolor;
} else{
xgcv.foreground = fcolor;
xgcv.background = bcolor;
}
colorgc = XCreateGC( disp, graph_win,
                      GCForeground|GCBackground|GCFunction, &xgcv );
if(XAllocNamedColor(disp,colormap,"blue",&thecolor,&foocolor)){
xgcv.foreground = thecolor.pixel;
xgcv.background = bcolor;
} else{
xgcv.foreground = fcolor;
xgcv.background = bcolor;
}

colorgc2 = XCreateGC( disp, graph_win,
                      GCForeground|GCBackground|GCFunction, &xgcv );

  /* now create the children */
#if 0
  gwin = XCreateSimpleWindow(disp,graph_win,BUTWIDTH,0,
                             GWINWIDTH,GWINHEIGHT,2,fcolor,bcolor);
  butwin = XCreateSimpleWindow(disp,graph_win,0,0,
                               BUTWIDTH,BUTHEIGHT,1,fcolor,bcolor);
#endif
  gwin = graph_win;
  /* set up the font gc's */
  xgcv.background = bcolor;
  xgcv.foreground=fcolor;
  bigtextgc = XCreateGC( disp, graph_win,
                        GCForeground|GCBackground|GCFunction, &xgcv );
  smalltextgc = XCreateGC( disp, graph_win,
                          GCForeground|GCBackground|GCFunction, &xgcv );
  XSelectInput(disp,gwin,KeyPressMask|ButtonPressMask|ButtonReleaseMask|
               ButtonMotionMask|ExposureMask|StructureNotifyMask);
#if 0
  XSelectInput(disp,butwin,KeyPressMask|ButtonPressMask|ButtonMotionMask|ExposureMask);


  /* change the cursor in the button window */
  XDefineCursor(disp,butwin,XCreateFontCursor(disp,60));
#endif

  /* load in the fonts */
  big_font = XLoadQueryFont( disp, BIGFONTNAME );
  small_font = XLoadQueryFont( disp, SMALLFONTNAME );
  if( !small_font || !big_font ) error( "Fonts not found." );

  /* make sure that the fonts can be created */
  XSetFont( disp, bigtextgc, big_font->fid );
  XSetFont( disp, smalltextgc, small_font->fid);

  /* create the pixmaps */
  gpix=XCreatePixmap(disp,gwin,GWINWIDTH,GWINHEIGHT,screen_depth);

  if( !gpix )fatal("Can't get main pixmap!");

  /*now expose the windows.*/
#if 0
  XMapSubwindows( disp, graph_win );
#endif
  XMapWindow( disp, graph_win );
/*  XMapWindow( disp, butwin ); */

  XFlush(disp);


  /* set the linestyle to the appropriate values */
  XSetLineAttributes(disp,graphgc,1,LineSolid,CapButt,JoinRound);
  doing_X = 1;
  doing_tek = 0;
  D_FREE(name);
}
#endif


