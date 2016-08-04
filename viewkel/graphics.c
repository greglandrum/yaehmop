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

   03.05.98 gL:
     created new breaking lines/tubes routines which stop the
     outer white line at a specified point to allow lines to
     stop at the edges of atoms.
   18.05.98 gL:
     changes to g_shaded_polygon
   21.05.98 gL:
     stop_tube calls now invoke stoptube in the PS output.
   08.09.98 gL:
     support for colored and/or shaded atoms
   24.09.98 gL:
     support for colored and/or shaded atoms updated
     bounding box update in g_center_text (for labels)
     changes to g_label to allow nicer looking labels.
   26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
      under gcc (yeah yeah... it's anal)
     linecap set for the white line in stop_lines in PS output.
      before this was ending up being round, which resulted in
      some penetration of atoms below... this was wrong.
   01.10.98 gL:
     crashing bug in draw_label fixed.
   18.10.98 gL:
     modifications to filled_polygon and open_polygon routines to allow
     them to deal with things with more than 3 sides.  This was added as
     part of the fatbands stuff, but it should generalize them significantly.
     NOTE: the Mac version of the filled polygon thing still only deals with
     three sides... this needs to be updated.  The tektronix stuff has not
     been touched at all... not that that matters.  Secretly, g_open_polygon
     has been broken in X mode for as long as it has existed... hahah.
     luckily it has never been used in a situation where this would be
     noticeable.
***/


/***************

  the following commands are to allow dependancies on the graphics
  system to be abstracted away as much as possible

***************/

#include "viewkel.h"

#ifdef MAC_GRAPHICS
#include "Mac_defines.h"
#include "Mac_protos.h"
#include "EasyApp.h"
#include <Quickdraw.h>
#endif

/* these are the linestyles used */
#define NUM_PS_LINESTYLES 4
char *PSlinestyles[NUM_PS_LINESTYLES] = {"[1 0]","[2 2]","[6 6]","[4 4 2 4]"};


#ifdef X_GRAPHICS

#define NUM_X_LINESTYLES 4
char Xstyle1[] = {1,1};
char Xstyle2[] = {4,4};
char Xstyle3[] = {8,8};
char Xstyle4[] = {8,8,4,8};
char *Xlinestyles[NUM_X_LINESTYLES] = { Xstyle1, Xstyle2, Xstyle3, Xstyle4};
char Xstylelengths[NUM_X_LINESTYLES] = {2,2,2,4};

GC *active_gc;

char need_colormap=1;
Colormap colormap;
XColor thecolor;
XGCValues xgcv;

#endif

#ifdef MAC_GRAPHICS
#define NUM_MAC_LINESTYLES 4
RGBColor Maclinestyles[NUM_MAC_LINESTYLES] =
       {{0x0000,0x0000,0x0000},
            {0xFFFF,0x0000,0x0000},
            {0x0000,0xFFFF,0x0000},
            {0x0000,0x0000,0xFFFF}};
#endif


/* local global variables used to cache values to save time */
float curr_line_width=-13;
int curr_color=-13;
int curr_line_style=-13;

/*****

  a line

*****/
void g_line(float x1,float y1,float x2,float y2)
{

  if( doing_ps ){
    fprintf(psfile,"%lf %lf M %lf %lf L SNP\n",
            x1,y1,x2,y2);
  }

#ifdef X_GRAPHICS
  if( doing_X ){
    XDrawLine(disp,gpix,graphgc,(short)x1,(short)y1,(short)x2,(short)y2);
  }
#endif

#ifdef TEK_GRAPHICS
  if( doing_tek ){
    TEK_move((unsigned int)x1,g_ymax - (unsigned int)y1);
    TEK_vector((unsigned int)x2,g_ymax - (unsigned int)y2);
  }
#endif

#ifdef MAC_GRAPHICS
        if(doing_Mac){
                MoveTo((short)(GRAPHICS_SCALE*x1),
                       (short)(GRAPHICS_SCALE*y1));
                LineTo((short)(GRAPHICS_SCALE*x2),
                       (short)(GRAPHICS_SCALE*y2));
        }
#endif
}

/*****

  Multiple lines

  'fill_it is used to indicate whether or not the region bounded by the line
  should be filled.

  *****/
void g_lines(XPoint *Ipoints,point_type2D *Fpoints,
             int num_points,char fill_it)
{
  int i;
#ifdef MAC_GRAPHICS
  RgnHandle theRegion;
#endif

  if( doing_ps ){
    fprintf(psfile,"2 setlinejoin\n");
    fprintf(psfile,"%6.2lf %6.2lf M\n",Fpoints[0].x,Fpoints[0].y);
    for(i=1;i<num_points;i++){
#if 0
      /* this is a space saving check to avoid overly long paths */
      if( i>num_points-1 ||
         !((Fpoints[i].x == Fpoints[i-1].x &&
            Fpoints[i].x == Fpoints[i+1].x) ||
           (Fpoints[i].y == Fpoints[i-1].y &&
            Fpoints[i].y == Fpoints[i+1].y) ) ){
        fprintf(psfile,"%lf %lf L\n",Fpoints[i].x,Fpoints[i].y);
      }
#endif
      fprintf(psfile,"%6.2lf %6.2lf L\n",Fpoints[i].x,Fpoints[i].y);
    }
    if( fill_it ){
      fprintf(psfile,"gsave\n");
      fprintf(psfile,"pat3 8 1 30 72 300 32 div div setpattern\n");
      fprintf(psfile,"fill grestore\n");
    }

    fprintf(psfile,"stroke\n");
  }

#ifdef X_GRAPHICS
  if( doing_X ){
    if( fill_it ){
      XFillPolygon(disp,gpix,graygc[1],Ipoints,num_points,Complex,
                   CoordModeOrigin);
    }
    XDrawLines(disp,gpix,graphgc,Ipoints,num_points,CoordModeOrigin);

  }
#endif

#ifdef TEK_GRAPHICS
  if( doing_tek ){
    TEK_move((unsigned int)Ipoints[0].x,g_ymax - (unsigned int)Ipoints[0].y);
    for(i=0;i<num_points;i++){
      TEK_vector((unsigned int)Ipoints[i].x,g_ymax - (unsigned int)Ipoints[i].y);
    }
  }
#endif

#ifdef MAC_GRAPHICS
  if(doing_Mac){
    if(fill_it){
      theRegion = NewRgn();
    }
    MoveTo((short)(GRAPHICS_SCALE*Ipoints[0].x),
           (short)(GRAPHICS_SCALE*Ipoints[0].y));
    for(i=1;i<num_points;i++){
      LineTo((short)(GRAPHICS_SCALE*Ipoints[i].x),
             (short)(GRAPHICS_SCALE*Ipoints[i].y));
    }

    if(fill_it){
      LineTo((short)(GRAPHICS_SCALE*Ipoints[0].x),
             (short)(GRAPHICS_SCALE*Ipoints[0].y));

      CloseRgn(theRegion);
      /*
         PenPat(&qd.gray);
         PaintRgn(theRegion);
         PenPat(&qd.black);
      */
      FillRgn(theRegion,&qd.black);
      FrameRgn(theRegion);
      DisposeRgn(theRegion);
    }
  }
#endif

}

#ifdef DEBUG_HIDDEN_LINE
void g_clines(XPoint *Ipoints,point_type2D *Fpoints,
              int num_points,char fill_it)
{
  int i;
#ifdef MAC_GRAPHICS
  RgnHandle theRegion;
#endif

#ifdef X_GRAPHICS
  if( doing_X ){
    XDrawLines(disp,gpix,colorgc2,Ipoints,num_points,CoordModeOrigin);
  }
#endif
}
#endif

/*******

  centered text

  ********/
void g_center_text(float x,float y,char *string)
{
  int textloc;
  int twidth;
  int theight;

  if( doing_ps ){
    fprintf(psfile,"%.1lf %.1lf textsize add M\n",x,y);
    ENHPS_put_text((int)x,(int)y,string,CENTER_JUST);
  }
#ifdef  X_GRAPHICS
  if( doing_X ){
    twidth = XTextWidth(big_font,string,strlen(string));
    theight = big_font->ascent;
    textloc = (int)x - twidth/2;
    XDrawString(disp,gpix,bigtextgc,textloc,
                (int)y+big_font->ascent,string,strlen(string));
    if( x + twidth/2 > localmax.x )
      localmax.x = x+twidth/2;
    if( x - twidth/2 < localmin.x )
      localmin.x = x-twidth/2;
    if( y+theight > localmax.y ) localmax.y = y+theight;
    if( y < localmin.y ) localmin.y = y;
  }
#endif

#ifdef TEK_GRAPHICS
  if( doing_tek ){
    TEK_center_text((unsigned int)x,g_ymax - (unsigned int)y,string);
    if( x+(TEK_HCHAR*strlen(string))/2 > localmax.x )
      localmax.x = x+(TEK_HCHAR*strlen(string))/2;
    if( x-(TEK_HCHAR*strlen(string))/2 < localmin.x )
      localmin.x = x-(TEK_HCHAR*strlen(string))/2;
    if( y > localmax.y ) localmax.y = y;
    if( y < localmin.y ) localmin.y = y;
  }

#endif

#ifdef MAC_GRAPHICS
  /* !!! needs to be updated */
  if(doing_Mac){
    textloc = (int)(GRAPHICS_SCALE*x) -
      (int)(TextWidth(string,0,strlen(string))/2);
    MoveTo((short)textloc,(short)(GRAPHICS_SCALE*y +
                                  Mac_globals.fontSize));
    DrawText(string,0,strlen(string));
  }
#endif

}


/**********

  Right Justified text

  ***********/
void g_right_text(float x,float y,char *string)
{
  int textloc;

  if( doing_ps ){
    fprintf(psfile,"%.1lf %.1lf M\n",x,y);
    ENHPS_put_text((int)x,(int)y,string,RIGHT_JUST);
  }

#ifdef  X_GRAPHICS
  if( doing_X ){
    textloc = (int)x - XTextWidth(big_font,string,strlen(string));
    XDrawString(disp,gpix,bigtextgc,textloc,
                (int)y+big_font->ascent/2,string,strlen(string));
    if( x > localmax.x ) localmax.x = x;
    if( textloc < localmin.x ) localmin.x = x;
    if( y > localmax.y ) localmax.y = y;
    if( y < localmin.y ) localmin.y = y;

  }
#endif

#ifdef TEK_GRAPHICS
  if( doing_tek ){
    TEK_right_text((unsigned int)x,g_ymax - (unsigned int)(y),string);
    if( x > localmax.x )
      localmax.x = x;
    if( x-(TEK_HCHAR*strlen(string)) < localmin.x )
      localmin.x = x-(TEK_HCHAR*strlen(string));
    if( y > localmax.y ) localmax.y = y;
    if( y < localmin.y ) localmin.y = y;
  }

#endif

#ifdef MAC_GRAPHICS
  if(doing_Mac){
    textloc = (int)(GRAPHICS_SCALE*x) -
      (int)TextWidth(string,0,strlen(string));
    MoveTo((short)textloc,(short)(GRAPHICS_SCALE*y +
                                  Mac_globals.fontSize/2));
    DrawText(string,0,strlen(string));
  }
#endif

}



/**********

  Left Justified text

  ***********/
void g_left_text(float x,float y,char *string)
{
  int textloc;

  if( doing_ps ){
    fprintf(psfile,"%.1lf %.1lf M\n",x,y);
    ENHPS_put_text((int)x,(int)y,string,LEFT_JUST);
  }

#ifdef  X_GRAPHICS
  if( doing_X ){
    textloc = (int)x;
    XDrawString(disp,gpix,bigtextgc,textloc,
                (int)y+big_font->ascent/2,string,strlen(string));
  }
#endif

#ifdef TEK_GRAPHICS
  if( doing_tek ){
    TEK_put_text((unsigned int)x,g_ymax - (unsigned int)(y),string);
  }

#endif
#ifdef MAC_GRAPHICS
  if(doing_Mac){
    textloc = (int)(GRAPHICS_SCALE*x);
    MoveTo((short)textloc,(short)(GRAPHICS_SCALE*y +
                                  Mac_globals.fontSize/2));
    DrawText(string,0,strlen(string));
  }
#endif

}




/************

  Changing linestyle

  ************/
void g_change_linestyle(int style)
{

/*  if( curr_line_style == style ) return;*/
  curr_line_style = style;

  if( doing_ps ){
    fprintf(psfile,"stroke newpath\n");
    fprintf(psfile,"%s centerdash\n",
            PSlinestyles[style%NUM_PS_LINESTYLES]);
  }

#ifdef X_GRAPHICS
  if( doing_X ){
    if( style%NUM_X_LINESTYLES){
      XSetDashes(disp,graphgc,0,Xlinestyles[style%NUM_X_LINESTYLES],
                 Xstylelengths[style%NUM_X_LINESTYLES]);
      XSetLineAttributes( disp, graphgc, 1, LineOnOffDash,
                         CapRound, JoinRound);
    }
    else{
      XSetLineAttributes( disp, graphgc, 0, LineSolid,
                         CapRound, JoinRound);
    }
  }
#endif

#ifdef TEK_GRAPHICS
  if( doing_tek ){
    TEK_linetype((unsigned int)style);
  }
#endif

#ifdef MAC_GRAPHICS
  if( doing_Mac ){
    RGBForeColor(&(Maclinestyles[style%NUM_MAC_LINESTYLES]));
  }
#endif
}

/************

  Changing the linewidth

  ************/
void g_change_linewidth(float thickness)
{

/*  if(thickness == curr_line_width) return;*/

  curr_line_width = thickness;

  if( doing_ps ){
    fprintf(psfile,"stroke newpath\n");
    fprintf(psfile,"%lf setlinewidth\n",thickness);
  }

#ifdef X_GRAPHICS
  if( doing_X ){
    if( thickness != 1 ){
      XSetLineAttributes( disp, graphgc, (short)thickness, LineSolid,
                         CapRound, JoinRound);
    }else{
      XSetLineAttributes( disp, graphgc, (short)0, LineSolid,
                         CapRound, JoinRound);
    }
  }
#endif

#ifdef TEK_GRAPHICS
  if( doing_tek ){
    /****

      I don't know how to do this in Tek (except an ugly way, which I don't
      feel like dealing with right now.)
      *****/
  }
#endif

#ifdef MAC_GRAPHICS
  if(doing_Mac){
    PenSize((short)thickness,(short)thickness);
  }
#endif
}


/*******

  a line which "breaks" other lines (for 3-D stuff)

  *********/
void g_draw_breaking_line(float x1,float y1,float x2,float y2,
                          int width)
{
#ifdef MAC_GRAPHICS
  RGBColor the_color;
  float slope;
#endif

/*  curr_line_width = (float)width;*/

  if( doing_ps ){
    fprintf(psfile,"%d setlinewidth 1 setgray\n",width+4);
    fprintf(psfile,"%lf %lf M %lf %lf L stroke\n",
            x1,y1,x2,y2);
    fprintf(psfile,"%d setlinewidth 0 setgray\n",width);
    fprintf(psfile,"%lf %lf M %lf %lf L stroke\n",
            x1,y1,x2,y2);
  }

#ifdef X_GRAPHICS
  if( doing_X ){

    XSetLineAttributes( disp, widegc, width+4, LineSolid,
                       CapNotLast, JoinRound);
    XDrawLine(disp,gpix,widegc,(short)x1,(short)y1,
              (short)x2,(short)y2);
    XSetLineAttributes( disp, graphgc, width, LineSolid,
                       CapRound, JoinRound);
    XDrawLine(disp,gpix,graphgc,(short)x1,(short)y1,(short)x2,(short)y2);

  }
#endif

#ifdef TEK_GRAPHICS
  if( doing_tek ){
    TEK_move((unsigned int)x1,g_ymax - (unsigned int)y1);
    TEK_vector((unsigned int)x2,g_ymax - (unsigned int)y2);
  }
#endif
#ifdef MAC_GRAPHICS
  if(doing_Mac){

    /***************

      this is ugly  on the mac because of the fact that
      quickdraw doesn't define line widths about the center of a line,
      but rather from the "top" of the line.  Since this is dumb and
      it would be a ton of work to get it right, the mac version won't
      always show the self-breaking lines correctly on screen.
      The generated postscript will still show the breaking lines
      properly.

      ****************/

    /* switch to a white line */
    the_color.red = 0xFFFF;
    the_color.green = 0xFFFF;
    the_color.blue = 0xFFFF;
    RGBForeColor(&the_color);

    /* use the slope of the line to figure out how to draw the 2 side lines */
    slope = (y1-y2)/(x1-x2);
    slope = fabs(slope);
    if( slope >= 1){
      PenSize((short)width+1,(short)width+1);
      MoveTo((short)(GRAPHICS_SCALE*x1+width),(short)(GRAPHICS_SCALE*y1));
      LineTo((short)(GRAPHICS_SCALE*x2+width),(short)(GRAPHICS_SCALE*y2));
      MoveTo((short)(GRAPHICS_SCALE*x1-width),(short)(GRAPHICS_SCALE*y1));
      LineTo((short)(GRAPHICS_SCALE*x2-width),(short)(GRAPHICS_SCALE*y2));
    }else{
      PenSize((short)width+1,(short)width+1);
      MoveTo((short)(GRAPHICS_SCALE*x1),(short)(GRAPHICS_SCALE*y1-width));
      LineTo((short)(GRAPHICS_SCALE*x2),(short)(GRAPHICS_SCALE*y2-width));
      MoveTo((short)(GRAPHICS_SCALE*x1),(short)(GRAPHICS_SCALE*y1+width));
      LineTo((short)(GRAPHICS_SCALE*x2),(short)(GRAPHICS_SCALE*y2+width));
    }


    the_color.red = 0x0000;
    the_color.green = 0x0000;
    the_color.blue = 0x0000;
    RGBForeColor(&the_color);

    PenSize((short)width+1,(short)width+1);
    MoveTo((short)(GRAPHICS_SCALE*x1),(short)(GRAPHICS_SCALE*y1));
    LineTo((short)(GRAPHICS_SCALE*x2),(short)(GRAPHICS_SCALE*y2));


  }
#endif
}
/*******

  a line which "breaks" other lines (for 3-D stuff)

  *********/
void g_draw_stop_line(float x1,float y1,float x2,float y2,
                          int width,float xs, float ys)
{
#ifdef MAC_GRAPHICS
  RGBColor the_color;
  float slope;
#endif

/*  curr_line_width = (float)width;*/

  if( doing_ps ){
    fprintf(psfile,"%d setlinewidth 1 setgray 0 setlinecap\n",width+4);
    fprintf(psfile,"%lf %lf M %lf %lf L stroke\n",
            xs,ys,x2,y2);
    fprintf(psfile,"%d setlinewidth 0 setgray 1 setlinecap\n",width);
    fprintf(psfile,"%lf %lf M %lf %lf L stroke\n",
            x1,y1,x2,y2);
  }

#ifdef X_GRAPHICS
  if( doing_X ){

    XSetLineAttributes( disp, widegc, width+4, LineSolid,
                       CapNotLast, JoinRound);
    XDrawLine(disp,gpix,widegc,(short)xs,(short)ys,
              (short)x2,(short)y2);
    XSetLineAttributes( disp, graphgc, width, LineSolid,
                       CapRound, JoinRound);
    XDrawLine(disp,gpix,graphgc,(short)x1,(short)y1,(short)x2,(short)y2);

  }
#endif

#ifdef TEK_GRAPHICS
  if( doing_tek ){
    TEK_move((unsigned int)x1,g_ymax - (unsigned int)y1);
    TEK_vector((unsigned int)x2,g_ymax - (unsigned int)y2);
  }
#endif
#ifdef MAC_GRAPHICS
  if(doing_Mac){

    /***************

      this is ugly  on the mac because of the fact that
      quickdraw doesn't define line widths about the center of a line,
      but rather from the "top" of the line.  Since this is dumb and
      it would be a ton of work to get it right, the mac version won't
      always show the self-breaking lines correctly on screen.
      The generated postscript will still show the breaking lines
      properly.

      ****************/

    /* switch to a white line */
    the_color.red = 0xFFFF;
    the_color.green = 0xFFFF;
    the_color.blue = 0xFFFF;
    RGBForeColor(&the_color);

    /* use the slope of the line to figure out how to draw the 2 side lines */
    slope = (y1-y2)/(x1-x2);
    slope = fabs(slope);
    if( slope >= 1){
      PenSize((short)width+1,(short)width+1);
      MoveTo((short)(GRAPHICS_SCALE*xs+width),(short)(GRAPHICS_SCALE*ys));
      LineTo((short)(GRAPHICS_SCALE*x2+width),(short)(GRAPHICS_SCALE*y2));
      MoveTo((short)(GRAPHICS_SCALE*xs-width),(short)(GRAPHICS_SCALE*ys));
      LineTo((short)(GRAPHICS_SCALE*x2-width),(short)(GRAPHICS_SCALE*y2));
    }else{
      PenSize((short)width+1,(short)width+1);
      MoveTo((short)(GRAPHICS_SCALE*xs),(short)(GRAPHICS_SCALE*ys-width));
      LineTo((short)(GRAPHICS_SCALE*x2),(short)(GRAPHICS_SCALE*y2-width));
      MoveTo((short)(GRAPHICS_SCALE*xs),(short)(GRAPHICS_SCALE*ys+width));
      LineTo((short)(GRAPHICS_SCALE*x2),(short)(GRAPHICS_SCALE*y2+width));
    }


    the_color.red = 0x0000;
    the_color.green = 0x0000;
    the_color.blue = 0x0000;
    RGBForeColor(&the_color);

    PenSize((short)width+1,(short)width+1);
    MoveTo((short)(GRAPHICS_SCALE*x1),(short)(GRAPHICS_SCALE*y1));
    LineTo((short)(GRAPHICS_SCALE*x2),(short)(GRAPHICS_SCALE*y2));


  }
#endif
}

/*******

  a "tube" which "breaks" other lines (for 3-D stuff)

*********/
void g_draw_tube_line(float x1,float y1,float x2,float y2,
                      int width)
{
#ifdef MAC_GRAPHICS
  RGBColor the_color;
  float slope;
#endif


/*  curr_line_width = width;*/

  if( doing_ps ){
    fprintf(psfile,"%d setlinewidth 1 setgray 0 setlinecap\n",width+8);
    fprintf(psfile,"%lf %lf M %lf %lf L stroke\n",
            x1,y1,x2,y2);
    fprintf(psfile,"%d setlinewidth 0 setgray 1 setlinecap\n",width+4);
    fprintf(psfile,"%lf %lf M %lf %lf L stroke\n",
            x1,y1,x2,y2);
    fprintf(psfile,"%d setlinewidth 1 setgray 1 setlinecap\n",width);
    fprintf(psfile,"%lf %lf M %lf %lf L stroke 0 setlinecap\n",
            x1,y1,x2,y2);
  }

#ifdef X_GRAPHICS
  if( doing_X ){
    XSetLineAttributes( disp, widegc, width+8, LineSolid,
                       CapNotLast, JoinRound);
    XDrawLine(disp,gpix,widegc,(short)x1,(short)y1,(short)x2,(short)y2);
    XSetLineAttributes( disp, graphgc, width+4, LineSolid,
                       CapRound,JoinRound);
    XDrawLine(disp,gpix,graphgc,(short)x1,(short)y1,(short)x2,(short)y2);
    XSetLineAttributes( disp, widegc, width, LineSolid, CapRound, JoinRound);
    XDrawLine(disp,gpix,widegc,(short)x1,(short)y1,(short)x2,(short)y2);


  }
#endif

#ifdef TEK_GRAPHICS
  if( doing_tek ){
    TEK_move((unsigned int)x1,g_ymax - (unsigned int)y1);
    TEK_vector((unsigned int)x2,g_ymax - (unsigned int)y2);
  }
#endif

#ifdef MAC_GRAPHICS
  if(doing_Mac){

    /***************

      this is ugly  on the mac because of the fact that
      quickdraw doesn't define line widths about the center of a line,
      but rather from the "top" of the line.  Since this is dumb and
      it would be a ton of work to get it right, the mac version won't
      always show the self-breaking lines correctly on screen.
      The generated postscript will still show the breaking lines
      properly.

      ****************/

    /* switch to a white line */
    the_color.red = 0xFFFF;
    the_color.green = 0xFFFF;
    the_color.blue = 0xFFFF;
    RGBForeColor(&the_color);

    /* use the slope of the line to figure out how to draw the 2 side lines */
    slope = (y1-y2)/(x1-x2);
    slope = fabs(slope);
#define OUTER_OFFSET 4
#define INNER_OFFSET 2
    if( slope >= 1){
      PenSize((short)width+OUTER_OFFSET,(short)width+OUTER_OFFSET);
      MoveTo((short)(GRAPHICS_SCALE*x1+width),
             (short)(GRAPHICS_SCALE*y1));
      LineTo((short)(GRAPHICS_SCALE*x2+width),
             (short)(GRAPHICS_SCALE*y2));
      MoveTo((short)(GRAPHICS_SCALE*x1-width-(OUTER_OFFSET-1)),
             (short)(GRAPHICS_SCALE*y1));
      LineTo((short)(GRAPHICS_SCALE*x2-width-(OUTER_OFFSET-1)),
             (short)(GRAPHICS_SCALE*y2));
    }else{
      PenSize((short)width+3,(short)width+3);
      MoveTo((short)(GRAPHICS_SCALE*x1),
             (short)(GRAPHICS_SCALE*y1-width-(OUTER_OFFSET-1)));
      LineTo((short)(GRAPHICS_SCALE*x2),
             (short)(GRAPHICS_SCALE*y2-width-(OUTER_OFFSET-1)));
      MoveTo((short)(GRAPHICS_SCALE*x1),
             (short)(GRAPHICS_SCALE*y1+width));
      LineTo((short)(GRAPHICS_SCALE*x2),
             (short)(GRAPHICS_SCALE*y2+width));
    }


    the_color.red = 0x0000;
    the_color.green = 0x0000;
    the_color.blue = 0x0000;
    RGBForeColor(&the_color);

    if( slope >= 1){
      PenSize((short)width+INNER_OFFSET,(short)width+INNER_OFFSET);
      MoveTo((short)(GRAPHICS_SCALE*x1+width),
             (short)(GRAPHICS_SCALE*y1));
      LineTo((short)(GRAPHICS_SCALE*x2+width),
             (short)(GRAPHICS_SCALE*y2));
      MoveTo((short)(GRAPHICS_SCALE*x1-width-(INNER_OFFSET-1)),
             (short)(GRAPHICS_SCALE*y1));
      LineTo((short)(GRAPHICS_SCALE*x2-width-(INNER_OFFSET-1)),
             (short)(GRAPHICS_SCALE*y2));
    }else{
      PenSize((short)width+2,(short)width+2);
      MoveTo((short)(GRAPHICS_SCALE*x1),
             (short)(GRAPHICS_SCALE*y1-width-(INNER_OFFSET-1)));
      LineTo((short)(GRAPHICS_SCALE*x2),
             (short)(GRAPHICS_SCALE*y2-width-(INNER_OFFSET-1)));
      MoveTo((short)(GRAPHICS_SCALE*x1),
             (short)(GRAPHICS_SCALE*y1+width));
      LineTo((short)(GRAPHICS_SCALE*x2),
             (short)(GRAPHICS_SCALE*y2+width));
    }


    /* switch to a white line */
    the_color.red = 0xFFFF;
    the_color.green = 0xFFFF;
    the_color.blue = 0xFFFF;
    RGBForeColor(&the_color);

    PenSize((short)width+1,(short)width+1);
    MoveTo((short)(GRAPHICS_SCALE*x1),(short)(GRAPHICS_SCALE*y1));
    LineTo((short)(GRAPHICS_SCALE*x2),(short)(GRAPHICS_SCALE*y2));

#undef OUTER_OFFSET
#undef INNER_OFFSET


  }
#endif
}



/*******

  a "tube" which "breaks" other lines (for 3-D stuff)

*********/
void g_draw_stop_tube(float x1,float y1,float x2,float y2,
                      int width,float xs,float ys)
{
#ifdef MAC_GRAPHICS
  RGBColor the_color;
  float slope;
#endif


/*  curr_line_width = width;*/

  if( doing_ps ){
#if 0
    fprintf(psfile,"%d setlinewidth 1 setgray 0 setlinecap\n",width+8);
    fprintf(psfile,"%lf %lf M %lf %lf L stroke\n",
            xs,ys,x2,y2);
    fprintf(psfile,"%d setlinewidth 0 setgray 1 setlinecap\n",width+4);
    fprintf(psfile,"%lf %lf M %lf %lf L stroke\n",
            x1,y1,x2,y2);
    fprintf(psfile,"%d setlinewidth 1 setgray 1 setlinecap\n",width);
    fprintf(psfile,"%lf %lf M %lf %lf L stroke 0 setlinecap\n",
            x1,y1,x2,y2);
#endif
    switch(PS_options.bond_type){
    case BOND_PLAIN:
      fprintf(psfile,"[%lf %d %lf %lf %lf %lf %lf %lf] stoptube\n",
              (float)1.0,width,x1,y1,x2,y2,xs,ys);
      break;
    case BOND_SHADE:
      fprintf(psfile,"[%lf %d %lf %lf %lf %lf %lf %lf] stoptube\n",
              (float)PS_options.bond_shade,width,x1,y1,x2,y2,xs,ys);
      break;
    case BOND_CSHADE:
      fprintf(psfile,"[%.4lf %.4lf %.4lf %d %lf %lf %lf %lf %lf %lf] Cstoptube\n",
              PS_options.bond_color[0],PS_options.bond_color[1],PS_options.bond_color[2],
              width,x1,y1,x2,y2,xs,ys);
      break;
    }

  }

#ifdef X_GRAPHICS
  if( doing_X ){
    XSetLineAttributes( disp, widegc, width+8, LineSolid,
                       CapNotLast, JoinRound);
    XDrawLine(disp,gpix,widegc,(short)xs,(short)ys,(short)x2,(short)y2);
    XSetLineAttributes( disp, graphgc, width+4, LineSolid,
                       CapRound,JoinRound);
    XDrawLine(disp,gpix,graphgc,(short)x1,(short)y1,(short)x2,(short)y2);
    XSetLineAttributes( disp, widegc, width, LineSolid, CapRound, JoinRound);
    XDrawLine(disp,gpix,widegc,(short)x1,(short)y1,(short)x2,(short)y2);


  }
#endif

#ifdef TEK_GRAPHICS
  if( doing_tek ){
    TEK_move((unsigned int)x1,g_ymax - (unsigned int)y1);
    TEK_vector((unsigned int)x2,g_ymax - (unsigned int)y2);
  }
#endif

#ifdef MAC_GRAPHICS
  if(doing_Mac){

    /***************

      this is ugly  on the mac because of the fact that
      quickdraw doesn't define line widths about the center of a line,
      but rather from the "top" of the line.  Since this is dumb and
      it would be a ton of work to get it right, the mac version won't
      always show the self-breaking lines correctly on screen.
      The generated postscript will still show the breaking lines
      properly.

      ****************/

    /* switch to a white line */
    the_color.red = 0xFFFF;
    the_color.green = 0xFFFF;
    the_color.blue = 0xFFFF;
    RGBForeColor(&the_color);

    /* use the slope of the line to figure out how to draw the 2 side lines */
    slope = (y1-y2)/(x1-x2);
    slope = fabs(slope);
#define OUTER_OFFSET 4
#define INNER_OFFSET 2
    if( slope >= 1){
      PenSize((short)width+OUTER_OFFSET,(short)width+OUTER_OFFSET);
      MoveTo((short)(GRAPHICS_SCALE*xs+width),
             (short)(GRAPHICS_SCALE*ys));
      LineTo((short)(GRAPHICS_SCALE*x2+width),
             (short)(GRAPHICS_SCALE*y2));
      MoveTo((short)(GRAPHICS_SCALE*xs-width-(OUTER_OFFSET-1)),
             (short)(GRAPHICS_SCALE*ys));
      LineTo((short)(GRAPHICS_SCALE*x2-width-(OUTER_OFFSET-1)),
             (short)(GRAPHICS_SCALE*y2));
    }else{
      PenSize((short)width+3,(short)width+3);
      MoveTo((short)(GRAPHICS_SCALE*xs),
             (short)(GRAPHICS_SCALE*ys-width-(OUTER_OFFSET-1)));
      LineTo((short)(GRAPHICS_SCALE*x2),
             (short)(GRAPHICS_SCALE*y2-width-(OUTER_OFFSET-1)));
      MoveTo((short)(GRAPHICS_SCALE*xs),
             (short)(GRAPHICS_SCALE*ys+width));
      LineTo((short)(GRAPHICS_SCALE*x2),
             (short)(GRAPHICS_SCALE*y2+width));
    }


    the_color.red = 0x0000;
    the_color.green = 0x0000;
    the_color.blue = 0x0000;
    RGBForeColor(&the_color);

    if( slope >= 1){
      PenSize((short)width+INNER_OFFSET,(short)width+INNER_OFFSET);
      MoveTo((short)(GRAPHICS_SCALE*x1+width),
             (short)(GRAPHICS_SCALE*y1));
      LineTo((short)(GRAPHICS_SCALE*x2+width),
             (short)(GRAPHICS_SCALE*y2));
      MoveTo((short)(GRAPHICS_SCALE*x1-width-(INNER_OFFSET-1)),
             (short)(GRAPHICS_SCALE*y1));
      LineTo((short)(GRAPHICS_SCALE*x2-width-(INNER_OFFSET-1)),
             (short)(GRAPHICS_SCALE*y2));
    }else{
      PenSize((short)width+2,(short)width+2);
      MoveTo((short)(GRAPHICS_SCALE*x1),
             (short)(GRAPHICS_SCALE*y1-width-(INNER_OFFSET-1)));
      LineTo((short)(GRAPHICS_SCALE*x2),
             (short)(GRAPHICS_SCALE*y2-width-(INNER_OFFSET-1)));
      MoveTo((short)(GRAPHICS_SCALE*x1),
             (short)(GRAPHICS_SCALE*y1+width));
      LineTo((short)(GRAPHICS_SCALE*x2),
             (short)(GRAPHICS_SCALE*y2+width));
    }


    /* switch to a white line */
    the_color.red = 0xFFFF;
    the_color.green = 0xFFFF;
    the_color.blue = 0xFFFF;
    RGBForeColor(&the_color);

    PenSize((short)width+1,(short)width+1);
    MoveTo((short)(GRAPHICS_SCALE*x1),(short)(GRAPHICS_SCALE*y1));
    LineTo((short)(GRAPHICS_SCALE*x2),(short)(GRAPHICS_SCALE*y2));

#undef OUTER_OFFSET
#undef INNER_OFFSET


  }
#endif
}




/*******

  a line which "breaks" other lines (for 3-D stuff)

  *********/
void g_draw_dashed_breaking_line(float x1,float y1,float x2,float y2,
                          int width,int pattern)
{
#ifdef MAC_GRAPHICS
  RGBColor the_color;
  float slope;
#endif

/*  curr_line_width = (float)width;*/

  if( doing_ps ){
    fprintf(psfile,"%d setlinewidth 1 setgray\n",width+4);
    fprintf(psfile,"%lf %lf M %lf %lf L stroke\n",
            x1,y1,x2,y2);
    fprintf(psfile,"%d setlinewidth 0 setgray\n",width);
    fprintf(psfile,"%s centerdash\n",
            PSlinestyles[pattern%NUM_PS_LINESTYLES]);
    fprintf(psfile,"%lf %lf M %lf %lf L stroke\n",
            x1,y1,x2,y2);
    fprintf(psfile,"%s centerdash\n",
            PSlinestyles[0]);
  }

#ifdef X_GRAPHICS
  if( doing_X ){

    XSetLineAttributes( disp, widegc, width+4, LineSolid,
                       CapNotLast, JoinRound);
    XDrawLine(disp,gpix,widegc,(short)x1,(short)y1,
              (short)x2,(short)y2);
    XSetDashes(disp,graphgc,0,Xlinestyles[pattern%NUM_X_LINESTYLES],
               Xstylelengths[pattern%NUM_X_LINESTYLES]);
    XSetLineAttributes( disp, graphgc, width, LineOnOffDash,
                        CapRound, JoinRound);
    XDrawLine(disp,gpix,graphgc,(short)x1,(short)y1,(short)x2,(short)y2);

  }
#endif

#ifdef TEK_GRAPHICS
  if( doing_tek ){
    TEK_move((unsigned int)x1,g_ymax - (unsigned int)y1);
    TEK_vector((unsigned int)x2,g_ymax - (unsigned int)y2);
  }
#endif
#ifdef MAC_GRAPHICS
  if(doing_Mac){

    /***************

      this is ugly  on the mac because of the fact that
      quickdraw doesn't define line widths about the center of a line,
      but rather from the "top" of the line.  Since this is dumb and
      it would be a ton of work to get it right, the mac version won't
      always show the self-breaking lines correctly on screen.
      The generated postscript will still show the breaking lines
      properly.

      ****************/

    /* switch to a white line */
    the_color.red = 0xFFFF;
    the_color.green = 0xFFFF;
    the_color.blue = 0xFFFF;
    RGBForeColor(&the_color);

    /* use the slope of the line to figure out how to draw the 2 side lines */
    slope = (y1-y2)/(x1-x2);
    slope = fabs(slope);
    if( slope >= 1){
      PenSize((short)width+1,(short)width+1);
      MoveTo((short)(GRAPHICS_SCALE*x1+width),(short)(GRAPHICS_SCALE*y1));
      LineTo((short)(GRAPHICS_SCALE*x2+width),(short)(GRAPHICS_SCALE*y2));
      MoveTo((short)(GRAPHICS_SCALE*x1-width),(short)(GRAPHICS_SCALE*y1));
      LineTo((short)(GRAPHICS_SCALE*x2-width),(short)(GRAPHICS_SCALE*y2));
    }else{
      PenSize((short)width+1,(short)width+1);
      MoveTo((short)(GRAPHICS_SCALE*x1),(short)(GRAPHICS_SCALE*y1-width));
      LineTo((short)(GRAPHICS_SCALE*x2),(short)(GRAPHICS_SCALE*y2-width));
      MoveTo((short)(GRAPHICS_SCALE*x1),(short)(GRAPHICS_SCALE*y1+width));
      LineTo((short)(GRAPHICS_SCALE*x2),(short)(GRAPHICS_SCALE*y2+width));
    }


    the_color.red = 0x0000;
    the_color.green = 0x0000;
    the_color.blue = 0x0000;
    RGBForeColor(&the_color);

    PenSize((short)width+1,(short)width+1);
    MoveTo((short)(GRAPHICS_SCALE*x1),(short)(GRAPHICS_SCALE*y1));
    LineTo((short)(GRAPHICS_SCALE*x2),(short)(GRAPHICS_SCALE*y2));


  }
#endif
}

/*******

  a "tube" which "breaks" other lines (for 3-D stuff)

*********/
void g_draw_dashed_tube_line(float x1,float y1,float x2,float y2,
                      int width,int pattern)
{
#ifdef MAC_GRAPHICS
  RGBColor the_color;
  float slope;
#endif


/*  curr_line_width = width;*/

  if( doing_ps ){
    fprintf(psfile,"%d setlinewidth 1 setgray 0 setlinecap\n",width+8);
    fprintf(psfile,"%lf %lf M %lf %lf L stroke\n",
            x1,y1,x2,y2);
    fprintf(psfile,"%d setlinewidth 0 setgray 1 setlinecap\n",width+4);
    fprintf(psfile,"%lf %lf M %lf %lf L stroke\n",
            x1,y1,x2,y2);
    fprintf(psfile,"%d setlinewidth 1 setgray 1 setlinecap\n",width);
    fprintf(psfile,"%s centerdash\n",
            PSlinestyles[pattern%NUM_PS_LINESTYLES]);
    fprintf(psfile,"%lf %lf M %lf %lf L stroke 0 setlinecap\n",
            x1,y1,x2,y2);
    fprintf(psfile,"%s centerdash\n",
            PSlinestyles[0]);
  }

#ifdef X_GRAPHICS
  if( doing_X ){
    XSetLineAttributes( disp, widegc, width+8, LineSolid,
                       CapNotLast, JoinRound);
    XDrawLine(disp,gpix,widegc,(short)x1,(short)y1,(short)x2,(short)y2);
    XSetLineAttributes( disp, graphgc, width+4, LineSolid,
                       CapRound,JoinRound);
    XDrawLine(disp,gpix,graphgc,(short)x1,(short)y1,(short)x2,(short)y2);

    XSetDashes(disp,graphgc,0,Xlinestyles[pattern%NUM_X_LINESTYLES],
               Xstylelengths[pattern%NUM_X_LINESTYLES]);
    XSetLineAttributes( disp, widegc, width, LineOnOffDash, CapRound, JoinRound);
    XDrawLine(disp,gpix,widegc,(short)x1,(short)y1,(short)x2,(short)y2);


  }
#endif

#ifdef TEK_GRAPHICS
  if( doing_tek ){
    TEK_move((unsigned int)x1,g_ymax - (unsigned int)y1);
    TEK_vector((unsigned int)x2,g_ymax - (unsigned int)y2);
  }
#endif

#ifdef MAC_GRAPHICS
  if(doing_Mac){

    /***************

      this is ugly  on the mac because of the fact that
      quickdraw doesn't define line widths about the center of a line,
      but rather from the "top" of the line.  Since this is dumb and
      it would be a ton of work to get it right, the mac version won't
      always show the self-breaking lines correctly on screen.
      The generated postscript will still show the breaking lines
      properly.

      ****************/

    /* switch to a white line */
    the_color.red = 0xFFFF;
    the_color.green = 0xFFFF;
    the_color.blue = 0xFFFF;
    RGBForeColor(&the_color);

    /* use the slope of the line to figure out how to draw the 2 side lines */
    slope = (y1-y2)/(x1-x2);
    slope = fabs(slope);
#define OUTER_OFFSET 4
#define INNER_OFFSET 2
    if( slope >= 1){
      PenSize((short)width+OUTER_OFFSET,(short)width+OUTER_OFFSET);
      MoveTo((short)(GRAPHICS_SCALE*x1+width),
             (short)(GRAPHICS_SCALE*y1));
      LineTo((short)(GRAPHICS_SCALE*x2+width),
             (short)(GRAPHICS_SCALE*y2));
      MoveTo((short)(GRAPHICS_SCALE*x1-width-(OUTER_OFFSET-1)),
             (short)(GRAPHICS_SCALE*y1));
      LineTo((short)(GRAPHICS_SCALE*x2-width-(OUTER_OFFSET-1)),
             (short)(GRAPHICS_SCALE*y2));
    }else{
      PenSize((short)width+3,(short)width+3);
      MoveTo((short)(GRAPHICS_SCALE*x1),
             (short)(GRAPHICS_SCALE*y1-width-(OUTER_OFFSET-1)));
      LineTo((short)(GRAPHICS_SCALE*x2),
             (short)(GRAPHICS_SCALE*y2-width-(OUTER_OFFSET-1)));
      MoveTo((short)(GRAPHICS_SCALE*x1),
             (short)(GRAPHICS_SCALE*y1+width));
      LineTo((short)(GRAPHICS_SCALE*x2),
             (short)(GRAPHICS_SCALE*y2+width));
    }


    the_color.red = 0x0000;
    the_color.green = 0x0000;
    the_color.blue = 0x0000;
    RGBForeColor(&the_color);

    if( slope >= 1){
      PenSize((short)width+INNER_OFFSET,(short)width+INNER_OFFSET);
      MoveTo((short)(GRAPHICS_SCALE*x1+width),
             (short)(GRAPHICS_SCALE*y1));
      LineTo((short)(GRAPHICS_SCALE*x2+width),
             (short)(GRAPHICS_SCALE*y2));
      MoveTo((short)(GRAPHICS_SCALE*x1-width-(INNER_OFFSET-1)),
             (short)(GRAPHICS_SCALE*y1));
      LineTo((short)(GRAPHICS_SCALE*x2-width-(INNER_OFFSET-1)),
             (short)(GRAPHICS_SCALE*y2));
    }else{
      PenSize((short)width+2,(short)width+2);
      MoveTo((short)(GRAPHICS_SCALE*x1),
             (short)(GRAPHICS_SCALE*y1-width-(INNER_OFFSET-1)));
      LineTo((short)(GRAPHICS_SCALE*x2),
             (short)(GRAPHICS_SCALE*y2-width-(INNER_OFFSET-1)));
      MoveTo((short)(GRAPHICS_SCALE*x1),
             (short)(GRAPHICS_SCALE*y1+width));
      LineTo((short)(GRAPHICS_SCALE*x2),
             (short)(GRAPHICS_SCALE*y2+width));
    }


    /* switch to a white line */
    the_color.red = 0xFFFF;
    the_color.green = 0xFFFF;
    the_color.blue = 0xFFFF;
    RGBForeColor(&the_color);

    PenSize((short)width+1,(short)width+1);
    MoveTo((short)(GRAPHICS_SCALE*x1),(short)(GRAPHICS_SCALE*y1));
    LineTo((short)(GRAPHICS_SCALE*x2),(short)(GRAPHICS_SCALE*y2));

#undef OUTER_OFFSET
#undef INNER_OFFSET


  }
#endif
}

/*******

  a line which "breaks" other lines (for 3-D stuff)

  *********/
void g_draw_dashed_stop_line(float x1,float y1,float x2,float y2,
                             int width,int pattern,float xs,float ys)
{
#ifdef MAC_GRAPHICS
  RGBColor the_color;
  float slope;
#endif

/*  curr_line_width = (float)width;*/

  if( doing_ps ){
    fprintf(psfile,"%d setlinewidth 1 setgray\n",width+4);
    fprintf(psfile,"%lf %lf M %lf %lf L stroke\n",
            xs,ys,x2,y2);
    fprintf(psfile,"%d setlinewidth 0 setgray\n",width);
    fprintf(psfile,"%s centerdash\n",
            PSlinestyles[pattern%NUM_PS_LINESTYLES]);
    fprintf(psfile,"%lf %lf M %lf %lf L stroke\n",
            x1,y1,x2,y2);
    fprintf(psfile,"%s centerdash\n",
            PSlinestyles[0]);
  }

#ifdef X_GRAPHICS
  if( doing_X ){

    XSetLineAttributes( disp, widegc, width+4, LineSolid,
                       CapNotLast, JoinRound);
    XDrawLine(disp,gpix,widegc,(short)xs,(short)ys,
              (short)x2,(short)y2);
    XSetDashes(disp,graphgc,0,Xlinestyles[pattern%NUM_X_LINESTYLES],
               Xstylelengths[pattern%NUM_X_LINESTYLES]);
    XSetLineAttributes( disp, graphgc, width, LineOnOffDash,
                        CapRound, JoinRound);
    XDrawLine(disp,gpix,graphgc,(short)x1,(short)y1,(short)x2,(short)y2);

  }
#endif

#ifdef TEK_GRAPHICS
  if( doing_tek ){
    TEK_move((unsigned int)x1,g_ymax - (unsigned int)y1);
    TEK_vector((unsigned int)x2,g_ymax - (unsigned int)y2);
  }
#endif
#ifdef MAC_GRAPHICS
  if(doing_Mac){

    /***************

      this is ugly  on the mac because of the fact that
      quickdraw doesn't define line widths about the center of a line,
      but rather from the "top" of the line.  Since this is dumb and
      it would be a ton of work to get it right, the mac version won't
      always show the self-breaking lines correctly on screen.
      The generated postscript will still show the breaking lines
      properly.

      ****************/

    /* switch to a white line */
    the_color.red = 0xFFFF;
    the_color.green = 0xFFFF;
    the_color.blue = 0xFFFF;
    RGBForeColor(&the_color);

    /* use the slope of the line to figure out how to draw the 2 side lines */
    slope = (y1-y2)/(x1-x2);
    slope = fabs(slope);
    if( slope >= 1){
      PenSize((short)width+1,(short)width+1);
      MoveTo((short)(GRAPHICS_SCALE*xs+width),(short)(GRAPHICS_SCALE*ys));
      LineTo((short)(GRAPHICS_SCALE*x2+width),(short)(GRAPHICS_SCALE*y2));
      MoveTo((short)(GRAPHICS_SCALE*xs-width),(short)(GRAPHICS_SCALE*ys));
      LineTo((short)(GRAPHICS_SCALE*x2-width),(short)(GRAPHICS_SCALE*y2));
    }else{
      PenSize((short)width+1,(short)width+1);
      MoveTo((short)(GRAPHICS_SCALE*xs),(short)(GRAPHICS_SCALE*ys-width));
      LineTo((short)(GRAPHICS_SCALE*x2),(short)(GRAPHICS_SCALE*y2-width));
      MoveTo((short)(GRAPHICS_SCALE*xs),(short)(GRAPHICS_SCALE*ys+width));
      LineTo((short)(GRAPHICS_SCALE*x2),(short)(GRAPHICS_SCALE*y2+width));
    }


    the_color.red = 0x0000;
    the_color.green = 0x0000;
    the_color.blue = 0x0000;
    RGBForeColor(&the_color);

    PenSize((short)width+1,(short)width+1);
    MoveTo((short)(GRAPHICS_SCALE*x1),(short)(GRAPHICS_SCALE*y1));
    LineTo((short)(GRAPHICS_SCALE*x2),(short)(GRAPHICS_SCALE*y2));


  }
#endif
}

/*******

  a dashed "tube" which "breaks" other lines (for 3-D stuff)

*********/
void g_draw_dashed_stop_tube(float x1,float y1,float x2,float y2,
                      int width,int pattern,float xs,float ys)
{
#ifdef MAC_GRAPHICS
  RGBColor the_color;
  float slope;
#endif


/*  curr_line_width = width;*/

  if( doing_ps ){
    fprintf(psfile,"%d setlinewidth 1 setgray 0 setlinecap\n",width+8);
    fprintf(psfile,"%lf %lf M %lf %lf L stroke\n",
            xs,ys,x2,y2);
    fprintf(psfile,"%d setlinewidth 0 setgray 1 setlinecap\n",width+4);
    fprintf(psfile,"%lf %lf M %lf %lf L stroke\n",
            x1,y1,x2,y2);
    fprintf(psfile,"%d setlinewidth 1 setgray 1 setlinecap\n",width);
    fprintf(psfile,"%s centerdash\n",
            PSlinestyles[pattern%NUM_PS_LINESTYLES]);
    fprintf(psfile,"%lf %lf M %lf %lf L stroke 0 setlinecap\n",
            x1,y1,x2,y2);
    fprintf(psfile,"%s centerdash\n",
            PSlinestyles[0]);
  }

#ifdef X_GRAPHICS
  if( doing_X ){
    XSetLineAttributes( disp, widegc, width+8, LineSolid,
                       CapNotLast, JoinRound);
    XDrawLine(disp,gpix,widegc,(short)xs,(short)ys,(short)x2,(short)y2);
    XSetLineAttributes( disp, graphgc, width+4, LineSolid,
                       CapRound,JoinRound);
    XDrawLine(disp,gpix,graphgc,(short)x1,(short)y1,(short)x2,(short)y2);

    XSetDashes(disp,graphgc,0,Xlinestyles[pattern%NUM_X_LINESTYLES],
               Xstylelengths[pattern%NUM_X_LINESTYLES]);
    XSetLineAttributes( disp, widegc, width, LineOnOffDash, CapRound, JoinRound);
    XDrawLine(disp,gpix,widegc,(short)x1,(short)y1,(short)x2,(short)y2);


  }
#endif

#ifdef TEK_GRAPHICS
  if( doing_tek ){
    TEK_move((unsigned int)x1,g_ymax - (unsigned int)y1);
    TEK_vector((unsigned int)x2,g_ymax - (unsigned int)y2);
  }
#endif

#ifdef MAC_GRAPHICS
  if(doing_Mac){

    /***************

      this is ugly  on the mac because of the fact that
      quickdraw doesn't define line widths about the center of a line,
      but rather from the "top" of the line.  Since this is dumb and
      it would be a ton of work to get it right, the mac version won't
      always show the self-breaking lines correctly on screen.
      The generated postscript will still show the breaking lines
      properly.

      ****************/

    /* switch to a white line */
    the_color.red = 0xFFFF;
    the_color.green = 0xFFFF;
    the_color.blue = 0xFFFF;
    RGBForeColor(&the_color);

    /* use the slope of the line to figure out how to draw the 2 side lines */
    slope = (y1-y2)/(x1-x2);
    slope = fabs(slope);
#define OUTER_OFFSET 4
#define INNER_OFFSET 2
    if( slope >= 1){
      PenSize((short)width+OUTER_OFFSET,(short)width+OUTER_OFFSET);
      MoveTo((short)(GRAPHICS_SCALE*xs+width),
             (short)(GRAPHICS_SCALE*ys));
      LineTo((short)(GRAPHICS_SCALE*x2+width),
             (short)(GRAPHICS_SCALE*y2));
      MoveTo((short)(GRAPHICS_SCALE*xs-width-(OUTER_OFFSET-1)),
             (short)(GRAPHICS_SCALE*ys));
      LineTo((short)(GRAPHICS_SCALE*x2-width-(OUTER_OFFSET-1)),
             (short)(GRAPHICS_SCALE*y2));
    }else{
      PenSize((short)width+3,(short)width+3);
      MoveTo((short)(GRAPHICS_SCALE*xs),
             (short)(GRAPHICS_SCALE*ys-width-(OUTER_OFFSET-1)));
      LineTo((short)(GRAPHICS_SCALE*x2),
             (short)(GRAPHICS_SCALE*y2-width-(OUTER_OFFSET-1)));
      MoveTo((short)(GRAPHICS_SCALE*xs),
             (short)(GRAPHICS_SCALE*ys+width));
      LineTo((short)(GRAPHICS_SCALE*x2),
             (short)(GRAPHICS_SCALE*y2+width));
    }


    the_color.red = 0x0000;
    the_color.green = 0x0000;
    the_color.blue = 0x0000;
    RGBForeColor(&the_color);

    if( slope >= 1){
      PenSize((short)width+INNER_OFFSET,(short)width+INNER_OFFSET);
      MoveTo((short)(GRAPHICS_SCALE*x1+width),
             (short)(GRAPHICS_SCALE*y1));
      LineTo((short)(GRAPHICS_SCALE*x2+width),
             (short)(GRAPHICS_SCALE*y2));
      MoveTo((short)(GRAPHICS_SCALE*x1-width-(INNER_OFFSET-1)),
             (short)(GRAPHICS_SCALE*y1));
      LineTo((short)(GRAPHICS_SCALE*x2-width-(INNER_OFFSET-1)),
             (short)(GRAPHICS_SCALE*y2));
    }else{
      PenSize((short)width+2,(short)width+2);
      MoveTo((short)(GRAPHICS_SCALE*x1),
             (short)(GRAPHICS_SCALE*y1-width-(INNER_OFFSET-1)));
      LineTo((short)(GRAPHICS_SCALE*x2),
             (short)(GRAPHICS_SCALE*y2-width-(INNER_OFFSET-1)));
      MoveTo((short)(GRAPHICS_SCALE*x1),
             (short)(GRAPHICS_SCALE*y1+width));
      LineTo((short)(GRAPHICS_SCALE*x2),
             (short)(GRAPHICS_SCALE*y2+width));
    }


    /* switch to a white line */
    the_color.red = 0xFFFF;
    the_color.green = 0xFFFF;
    the_color.blue = 0xFFFF;
    RGBForeColor(&the_color);

    PenSize((short)width+1,(short)width+1);
    MoveTo((short)(GRAPHICS_SCALE*x1),(short)(GRAPHICS_SCALE*y1));
    LineTo((short)(GRAPHICS_SCALE*x2),(short)(GRAPHICS_SCALE*y2));

#undef OUTER_OFFSET
#undef INNER_OFFSET


  }
#endif
}


#ifdef FORCED_PERSPECTIVE
/*******

  a "tube" with "perpective scaling" which "breaks" other lines (for 3-D stuff)

*********/
void  g_draw_trapezoid_bond(float P1x,float P1y,float P1z,
                            float P2x,float P2y,float P2z,
                            float bond_rad,float radius1,float radius2,
                            int width)
{
  float diffx,diffy;
  float dx1,dx2,dy1,dy2,slope,inv_slope;
  float prefact;
  float rad1,rad2;
  point_type2D P1t,P1b,P2t,P2b;
#ifdef MAC_GRAPHICS
  Rect the_rect;
#endif

  /* figure out the "radius" of the bond at each end */
/*  rad1 = bond_rad*rad_scale*(1-(P1z+camera->lf.z));
  rad2 = bond_rad*rad_scale*(1-(P2z+camera->lf.z));
  rad2 = rad1 + rad1 * (P1z - P2z);*/
  rad1 = radius1;
  rad2 = radius2;
/*printf("rad1: %lf (%lf)\t rad2:  %lf (%lf)\n",rad1,P1z,rad2,P2z);*/

  /* now figure out where the corners of the trapezoid are */
  diffx = P1x - P2x;
  diffy = P1y - P2y;
  if( fabs(diffy) > 1e-2 ){
    if( fabs(diffx) > 1e-2 ){
      slope = (P1y - P2y)/(P1x - P2x);

      inv_slope = -1.0/slope;
      prefact = 1/sqrt(inv_slope*inv_slope + 1);
      dx1 = prefact * rad1;
      dy1 = inv_slope * dx1;
      dx2 = prefact * rad2;
      dy2 = inv_slope * dx2;
      P1t.x = P1x + dx1; P2t.x = P2x + dx2;
      P1t.y = P1y + dy1; P2t.y = P2y + dy2;
      P1b.x = P1x - dx1; P2b.x = P2x - dx2;
      P1b.y = P1y - dy1; P2b.y = P2y - dy2;
    }else{
      /* no change in X, slope must be infinite */
      P1t.x = P1x + rad1; P2t.x = P2x + rad2;
      P1t.y = P1y; P2t.y = P2y;
      P1b.x = P1x - rad1; P2b.x = P2x - rad2;
      P1b.y = P1y; P2b.y = P2y;
    }
  }else{
    /* no change in Y, slope is zero */
    P1t.x = P1x; P2t.x = P2x;
    P1t.y = P1y + rad1; P2t.y = P2y + rad2;
    P1b.x = P1x; P2b.x = P2x;
    P1b.y = P1y - rad1; P2b.y = P2y - rad2;
  }


  /*************


    Here's the drawing order:
     filled circles at the end to "cap" the tubes
     breaking line from P1t->P2t, and a breaking line from
      P1b->P2b
     a filled trapezoid

  **************/
  if(doing_ps){
    fprintf(psfile,"1 setgray\n");
    fprintf(psfile,"%lf %lf %lf 0 360 arc fill stroke\n",
            P1x,P1y,rad1);
    fprintf(psfile,"%lf %lf %lf 0 360 arc fill stroke\n",
            P2x,P2y,rad2);
  }
#ifdef X_GRAPHICS
  if( doing_X ){
    XFillArc(disp,gpix,*active_gc,(int)(P1x-rad1),(int)(P1y-rad1),
             2*(int)rad1,2*(int)rad1,0,23040);
    XFillArc(disp,gpix,*active_gc,(int)(P2x-rad2),(int)(P2y-rad2),
             2*(int)rad2,2*(int)rad2,0,23040);
  }
#endif
  g_draw_breaking_line(P1t.x,P1t.y,P2t.x,P2t.y,width);
  g_draw_breaking_line(P1b.x,P1b.y,P2b.x,P2b.y,width);

}
#endif

/********

  Rectangles

  ********/
void g_rectangle(float x,float y,float dimx,float dimy)
{
#ifdef MAC_GRAPHICS
  Rect the_rect;
#endif
  if( doing_ps ){
    fprintf(psfile,"%lf %lf M %lf %lf L %lf %lf L\n",
            x,y,x+dimx,y,x+dimx,y-dimy);
    fprintf(psfile,"%lf %lf L %lf %lf L\n",
            x,y-dimy,x,y);
    fprintf(psfile,"stroke\n");
  }

#ifdef X_GRAPHICS
  if(doing_X){
    XDrawRectangle(disp,gpix,graphgc,(int)x,(int)(y-dimy),
                   (int)dimx,(int)dimy);
  }
#endif

#ifdef TEK_GRAPHICS
  if(doing_tek){
    TEK_move((unsigned int)x,g_ymax - (unsigned int)y);
    TEK_vector((unsigned int)(x+dimx),g_ymax - (unsigned int)y);
    TEK_vector((unsigned int)(x+dimx),g_ymax - (unsigned int)(y-dimy));
    TEK_vector((unsigned int)(x),g_ymax - (unsigned int)(y-dimy));
    TEK_vector((unsigned int)x,g_ymax - (unsigned int)y);
  }
#endif

#ifdef MAC_GRAPHICS
  if(doing_Mac){
    the_rect.top = (short)(GRAPHICS_SCALE*(y-dimy));
    the_rect.left = (short)(GRAPHICS_SCALE*x);
    /*****

      We have to put these little 1 pixel shifts in
      because of the way the Mac draws rectangles.

      ******/
    the_rect.bottom = (short)(GRAPHICS_SCALE*y)+1;
    the_rect.right = (short)(GRAPHICS_SCALE*(x+dimx))+1;
    FrameRect(&the_rect);
  }
#endif
}


/***************

  changing colors

  This only affects operations that use color

  ****************/
void g_change_color(int color)
{
  float graylevel;
#ifdef MAC_GRAPHICS
  RGBColor the_color;
#endif

/*  if( curr_color == color ) return;*/
  curr_color = color;

  if( doing_ps ){
    if( color != 0 ){
      graylevel = .9 - (float)(NUM_COLORS-color)/(float)NUM_COLORS;
    } else{
      graylevel = 0.0;
    }
    fprintf(psfile,"%4.2lf SG\n",graylevel);
  }

#ifdef X_GRAPHICS
  if(doing_X){
    active_gc = &(graygc[color%NUM_COLORS]);
  }
#endif

#ifdef MAC_GRAPHICS
  if(doing_Mac){
    the_color.red = the_color.green = the_color.blue =
      (short)((float)0xFFFF*(float)color/(float)NUM_COLORS);
    RGBForeColor(&the_color);
  }
#endif
}


/***************

  an open circle

  ****************/
void g_open_circle(float x,float y,float radius)
{
#ifdef MAC_GRAPHICS
  Rect the_rect;
#endif

  if(doing_ps){
    fprintf(psfile,"arcwid setlinewidth\n");
    fprintf(psfile,"SNP %lf %lf %lf 0 360 arc stroke\n",
            x,y,radius);
  }


#ifdef X_GRAPHICS
  if( doing_X ){
    XDrawArc(disp,gpix,graphgc,(int)(x-radius),(int)(y-radius),
             2*(int)radius,2*(int)radius,0,23040);
  }
#endif
#ifdef MAC_GRAPHICS
  if( doing_Mac ){
    /* figure out the bounding rectangle */
    the_rect.top = (short)((GRAPHICS_SCALE*y)-radius);
    the_rect.bottom = (short)((GRAPHICS_SCALE*y)+radius);
    the_rect.left = (short)((GRAPHICS_SCALE*x)-radius);
    the_rect.right = (short)((GRAPHICS_SCALE*x)+radius);
    FrameOval(&the_rect);
  }
#endif

}

/***************

  g_crossed_circle: a circle with a faked cross
   (to look 3D)

  ****************/
void g_crossed_circle(float x,float y,float radius)
{
#ifdef MAC_GRAPHICS
  Rect the_rect;
#endif

  if(doing_ps){
    fprintf(psfile,"arcwid setlinewidth\n");
    fprintf(psfile,"SNP %lf %lf %lf crossedcirc\n",x,y,radius);
  }

#ifdef X_GRAPHICS
  if( doing_X ){
    XDrawArc(disp,gpix,graphgc,(int)(x-radius),(int)(y-radius),
             2*(int)radius,2*(int)radius,0,23040);

    XDrawArc(disp,gpix,graphgc,(int)(x-radius),(int)(y-0.5*radius),
             2*(int)radius,(int)radius,0,-180*64);
    XDrawArc(disp,gpix,graphgc,(int)(x-0.5*radius),(int)(y-radius),
             (int)radius,2*(int)radius,90*64,180*64);
  }
#endif
#ifdef MAC_GRAPHICS
  if( doing_Mac ){
    /* figure out the bounding rectangle */
    the_rect.top = (short)((GRAPHICS_SCALE*y)-radius);
    the_rect.bottom = (short)((GRAPHICS_SCALE*y)+radius);
    the_rect.left = (short)((GRAPHICS_SCALE*x)-radius);
    the_rect.right = (short)((GRAPHICS_SCALE*x)+radius);
    FrameOval(&the_rect);
     /* now do the cross */
    the_rect.top = (short)((GRAPHICS_SCALE*y)-radius);
    the_rect.bottom = (short)((GRAPHICS_SCALE*y)+radius);
    the_rect.left = (short)((GRAPHICS_SCALE*x)-0.5*radius);
    the_rect.right = (short)((GRAPHICS_SCALE*x)+0.5*radius);
    FrameArc(&the_rect,0,-180);
    the_rect.top = (short)((GRAPHICS_SCALE*y)-0.5*radius);
    the_rect.bottom = (short)((GRAPHICS_SCALE*y)+0.5*radius);
    the_rect.left = (short)((GRAPHICS_SCALE*x)-radius);
    the_rect.right = (short)((GRAPHICS_SCALE*x)+radius);
    FrameArc(&the_rect,90,180);
  }
#endif

}



/***************

  a filled circle

  ****************/
void g_filled_circle(float x,float y,float radius,float shade,float color[3],
                     long int *Gpixel_val, long int *Cpixel_val)
{
  int i;
  float cur_rad,xoffset,yoffset;
#ifdef MAC_GRAPHICS
  Rect the_rect;
#endif
  if(doing_ps){
    switch(PS_options.atom_sphere_type){
    case ATOM_PLAIN_FILL:
      fprintf(psfile,"SNP %lf %lf %lf 0 360 arc fill stroke\n",
              x,y,radius);
      break;
    case ATOM_SHADE_FILL:
      fprintf(psfile,"[%4.4lf %lf %lf %lf] shadecirc\n",
              shade,x,y,radius);
      break;
    case ATOM_COLOR_FILL:
      fprintf(psfile,"%4.4lf %4.4lf %4.4lf setrgbcolor\n",
              color[0],color[1],color[2]);
      fprintf(psfile,"SNP %lf %lf %lf 0 360 arc fill stroke\n",
              x,y,radius);
      fprintf(psfile,"0 0 0 setrgbcolor\n");
      break;
    case ATOM_CSHADE_FILL:
      fprintf(psfile,"[%4.4lf %4.4lf %4.4lf %lf %lf %lf] Cshadecirc\n",
              color[0],color[1],color[2],x,y,radius);
      break;
    }
  }
#ifdef X_GRAPHICS
  if( doing_X ){
#ifdef SUPPORT_COLOR_X
    switch(PS_options.atom_sphere_type){
    case ATOM_SHADE_FILL:
      if( need_colormap ){
        colormap = DefaultColormap(disp,screen);
        need_colormap = 0;
      }
      if( Gpixel_val[0] == -1 || refresh_all_colormaps ){
        for(i=0;i<NUM_X_SHADES;i++){
          thecolor.red = (short)(65000*(shade+(1-shade)*i/NUM_X_SHADES));
          thecolor.green = (short)(65000*(shade+(1-shade)*i/NUM_X_SHADES));
          thecolor.blue = (short)(65000*(shade+(1-shade)*i/NUM_X_SHADES));
          if(!XAllocColor(disp,colormap,&thecolor))
            error("can't allocate color.");
          Gpixel_val[i] = thecolor.pixel;
        }
      }
      for(i=0;i<NUM_X_SHADES;i++){
        cur_rad = radius * (1 - (float)i*.08);
        xoffset = x - cur_rad * (1+(float)i*0.05);
        yoffset = y - cur_rad * (1- (float)i*0.05);
        xgcv.foreground = Gpixel_val[i];
        XChangeGC(disp,colorgc,GCForeground,&xgcv);
        XFillArc(disp,gpix,colorgc,(int)(xoffset),(int)(yoffset),
                 2*(int)cur_rad,2*(int)cur_rad,0,23040);
      }
      break;

    case ATOM_COLOR_FILL:
      if( need_colormap ){
        colormap = DefaultColormap(disp,screen);
        need_colormap = 0;
      }
      if( Cpixel_val[0]== -1 || refresh_all_colormaps ){
        for(i=0;i<NUM_X_SHADES;i++){
          thecolor.red = (short)(65000*(color[0]+
                                        (1-color[0])*i/NUM_X_SHADES));
          thecolor.green = (short)(65000*(color[1]+
                                        (1-color[1])*i/NUM_X_SHADES));
          thecolor.blue = (short)(65000*(color[2]+
                                        (1-color[2])*i/NUM_X_SHADES));
          if(!XAllocColor(disp,colormap,&thecolor))
            error("can't allocate color.");
          Cpixel_val[i] = thecolor.pixel;
        }
      }
      xgcv.foreground = Cpixel_val[0];
      XChangeGC(disp,colorgc,GCForeground,&xgcv);
      XFillArc(disp,gpix,colorgc,(int)(x-radius),(int)(y-radius),
               2*(int)radius,2*(int)radius,0,23040);
      break;
    case ATOM_CSHADE_FILL:
      if( need_colormap ){
        colormap = DefaultColormap(disp,screen);
        need_colormap = 0;
      }
      if( *Cpixel_val== -1 || refresh_all_colormaps ){
        for(i=0;i<NUM_X_SHADES;i++){
          thecolor.red = (short)(65000*(color[0]+
                                        (1-color[0])*i/NUM_X_SHADES));
          thecolor.green = (short)(65000*(color[1]+
                                        (1-color[1])*i/NUM_X_SHADES));
          thecolor.blue = (short)(65000*(color[2]+
                                        (1-color[2])*i/NUM_X_SHADES));
          if(!XAllocColor(disp,colormap,&thecolor))
            error("can't allocate color.");
          Cpixel_val[i] = thecolor.pixel;
        }
      }
      for(i=0;i<NUM_X_SHADES;i++){
        cur_rad = radius * (1 - (float)i*.08);
        xoffset = x - cur_rad * (1+(float)i*0.05);
        yoffset = y - cur_rad * (1- (float)i*0.05);
        xgcv.foreground = Cpixel_val[i];
        XChangeGC(disp,colorgc,GCForeground,&xgcv);
        XFillArc(disp,gpix,colorgc,(int)(xoffset),(int)(yoffset),
                 2*(int)cur_rad,2*(int)cur_rad,0,23040);
      }
      break;

    default:
      XFillArc(disp,gpix,*active_gc,(int)(x-radius),(int)(y-radius),
               2*(int)radius,2*(int)radius,0,23040);
      break;
    }
#else
    XFillArc(disp,gpix,*active_gc,(int)(x-radius),(int)(y-radius),
             2*(int)radius,2*(int)radius,0,23040);
#endif
  }
#endif
#ifdef MAC_GRAPHICS
  if( doing_Mac ){
    /* figure out the bounding rectangle */
    the_rect.top = (short)((GRAPHICS_SCALE*y)-radius);
    the_rect.bottom = (short)((GRAPHICS_SCALE*y)+radius);
    the_rect.left = (short)((GRAPHICS_SCALE*x)-radius);
    the_rect.right = (short)((GRAPHICS_SCALE*x)+radius);
    PaintOval(&the_rect);
  }
#endif
}

/***************

  a filled white circle

  ****************/
void g_white_circle(float x,float y,float radius)
{
#ifdef MAC_GRAPHICS
  Rect the_rect;
#endif
  if(doing_ps){
    fprintf(psfile,"SNP 1 SG %lf %lf %lf 0 360 arc fill stroke\n",
            x,y,radius);
  }
#ifdef X_GRAPHICS
  if( doing_X ){
    XFillArc(disp,gpix,widegc,(int)(x-radius),(int)(y-radius),
             2*(int)radius,2*(int)radius,0,23040);
  }
#endif
#ifdef MAC_GRAPHICS
  if( doing_Mac ){
    /* figure out the bounding rectangle */
    the_rect.top = (short)((GRAPHICS_SCALE*y)-radius);
    the_rect.bottom = (short)((GRAPHICS_SCALE*y)+radius);
    the_rect.left = (short)((GRAPHICS_SCALE*x)-radius);
    the_rect.right = (short)((GRAPHICS_SCALE*x)+radius);
    PaintOval(&the_rect);
  }
#endif
}


/**************

  A Filled polygon (convex only please)

  ***************/
void g_filled_polygon(XPoint *points,int num_points)
{
  int i;
#ifdef MAC_GRAPHICS
  static PolyHandle        thePoly;
  static RgnHandle theRgn;
  static char firstCall = 1;
  short max_x,min_x,max_y,min_y;
#endif

  if( doing_ps ){
    if( num_points == 3 ){
      fprintf(psfile,"%d %d %d %d %d %d %d %d fT\n",
              points[0].x,points[0].y,points[1].x,points[1].y,
              points[2].x,points[2].y,points[0].x,points[0].y);
    } else{
      fprintf(psfile,"%d %d M ",points[0].x,points[0].y);
      for(i=1;i<num_points;i++){
        fprintf(psfile,"%d %d L ",points[i].x,points[i].y);
      }
      fprintf(psfile,"%d %d L fill stroke\n",points[0].x,points[0].y);
    }
  }

#ifdef X_GRAPHICS
  if( doing_X ){
    XFillPolygon(disp,gpix,*active_gc,points,num_points,
                 Convex,CoordModeOrigin);
  }
#endif

#ifdef MAC_GRAPHICS
  if( doing_Mac ){
#if 0
    if(num_points > 3) return;
    if(firstCall){
      thePoly = OpenPoly( );        /* start recording */
      MoveTo((short)(GRAPHICS_SCALE*points[0].x),
             (short)(GRAPHICS_SCALE*points[0].y));
      LineTo((short)(GRAPHICS_SCALE*points[1].x),
             (short)(GRAPHICS_SCALE*points[1].y));
      LineTo((short)(GRAPHICS_SCALE*points[2].x),
             (short)(GRAPHICS_SCALE*points[2].y));
      LineTo((short)(GRAPHICS_SCALE*points[0].x),
             (short)(GRAPHICS_SCALE*points[0].y));
      ClosePoly();                /* stop recording */
      firstCall = 0;
    }else{
      /* just directly fill the polyPoints array in the poly record */
      (*thePoly)->polyPoints[0].h = (short)(GRAPHICS_SCALE*points[0].x);
      (*thePoly)->polyPoints[0].v = (short)(GRAPHICS_SCALE*points[0].y);
      (*thePoly)->polyPoints[1].h = (short)(GRAPHICS_SCALE*points[1].x);
      (*thePoly)->polyPoints[1].v = (short)(GRAPHICS_SCALE*points[1].y);
      (*thePoly)->polyPoints[2].h = (short)(GRAPHICS_SCALE*points[2].x);
      (*thePoly)->polyPoints[2].v = (short)(GRAPHICS_SCALE*points[2].y);
      /* set the bounding box */
      max_x = (short)(GRAPHICS_SCALE*points[0].x);
      min_x = (short)(GRAPHICS_SCALE*points[0].x);
      max_y = (short)(GRAPHICS_SCALE*points[0].y);
      min_y = (short)(GRAPHICS_SCALE*points[0].y);

      if((short)(GRAPHICS_SCALE*points[1].x) > max_x)
        max_x = (short)(GRAPHICS_SCALE*points[1].x);
      else if((short)(GRAPHICS_SCALE*points[1].x) < min_x)
        min_x = (short)(GRAPHICS_SCALE*points[1].x);
      if((short)(GRAPHICS_SCALE*points[1].y) > max_y)
        max_y = (short)(GRAPHICS_SCALE*points[1].y);
      else if((short)(GRAPHICS_SCALE*points[1].y) < min_y)
        min_y = (short)(GRAPHICS_SCALE*points[1].y);

      if((short)(GRAPHICS_SCALE*points[2].x) > max_x)
        max_x = (short)(GRAPHICS_SCALE*points[2].x);
      else if((short)(GRAPHICS_SCALE*points[2].x) < min_x)
        min_x = (short)(GRAPHICS_SCALE*points[2].x);
      if((short)(GRAPHICS_SCALE*points[2].y) > max_y)
        max_y = (short)(GRAPHICS_SCALE*points[2].y);
      else if((short)(GRAPHICS_SCALE*points[2].y) < min_y)
        min_y = (short)(GRAPHICS_SCALE*points[2].y);

      (*thePoly)->polyBBox.top = min_y;
      (*thePoly)->polyBBox.bottom = max_y;
      (*thePoly)->polyBBox.left = min_x;
      (*thePoly)->polyBBox.right = max_x;

    }
    PaintPoly( thePoly);        /* paint the inside */
#endif
#if 0
    thePoly = OpenPoly( );        /* start recording */
    MoveTo((short)(GRAPHICS_SCALE*points[0].x),
           (short)(GRAPHICS_SCALE*points[0].y) );
    LineTo((short)(GRAPHICS_SCALE*points[1].x),
           (short)(GRAPHICS_SCALE*points[1].y) );
    LineTo((short)(GRAPHICS_SCALE*points[2].x),
           (short)(GRAPHICS_SCALE*points[2].y) );
    ClosePoly();                /* stop recording */
    PaintPoly( thePoly);
    KillPoly(thePoly);
#endif
    if( firstCall ){
      theRgn = NewRgn();
      firstCall = 0;
    }
    OpenRgn();
    MoveTo((short)(GRAPHICS_SCALE*points[0].x),
           (short)(GRAPHICS_SCALE*points[0].y) );
    LineTo((short)(GRAPHICS_SCALE*points[1].x),
           (short)(GRAPHICS_SCALE*points[1].y) );
    LineTo((short)(GRAPHICS_SCALE*points[2].x),
           (short)(GRAPHICS_SCALE*points[2].y) );
    LineTo((short)(GRAPHICS_SCALE*points[0].x),
           (short)(GRAPHICS_SCALE*points[0].y) );
    CloseRgn(theRgn);                /* stop recording */
    PaintRgn( theRgn);
    /*DisposeRgn(theRgn);*/

  }
  RGBForeColor(&(Maclinestyles[0]));
#endif


}



/**************

  A shaded polygon (convex only please)

***************/
void g_shaded_polygon(XPoint *points,int num_points,point_type *normal)
{
  int i;
#ifdef MAC_GRAPHICS
  static PolyHandle        thePoly;
  static RgnHandle theRgn;
  static char firstCall = 1;
  short max_x,min_x,max_y,min_y;
  RGBColor the_color;
#endif
  float colorscale;

#ifdef CULL_POLYGONS
  colorscale = fabs(normal->z);
#else
  colorscale = -(normal->z);
#endif


  if( doing_ps && colorscale >= 0){
    if( num_points == 3 ){
      fprintf(psfile,"%f SNG %d %d %d %d %d %d %d %d fT\n",
              colorscale,
              points[0].x,points[0].y,points[1].x,points[1].y,
              points[2].x,points[2].y,points[0].x,points[0].y);
    }else{
      fprintf(psfile,"%f SNG %d %d M ",colorscale,
              points[0].x,points[0].y);
      for(i=1;i<num_points;i++){
        fprintf(psfile,"%d %d L ",points[i].x,points[i].y);
      }
      fprintf(psfile,"%d %d L fill stroke\n",points[0].x,points[0].y);
    }
  }

#ifdef X_GRAPHICS
  if( doing_X ){
    if( need_colormap ){
      colormap = DefaultColormap(disp,screen);
      need_colormap = 0;
    }

#ifdef CULL_POLYGONS
    colorscale = colorscale / 1.5 + .16;
    thecolor.red = (short)(65000 * colorscale);
    thecolor.green = (short)(65000 * colorscale);
    thecolor.blue = (short)(65000 * colorscale);
#else
    if( colorscale >= 0.0 ){
      thecolor.red = (short)(65000 * colorscale);
      thecolor.blue = (short)(65000 * colorscale);
      thecolor.green = 0;
    }else{
      thecolor.blue =  0;
      thecolor.green =(short)(-65000 * colorscale);
      thecolor.red =(short)(-65000 * colorscale);
    }
#endif
    if(!XAllocColor(disp,colormap,&thecolor))
      error("can't allocate color.");
    xgcv.foreground = thecolor.pixel;
    XChangeGC(disp,colorgc,GCForeground,&xgcv);
    XFillPolygon(disp,gpix,colorgc,points,num_points,
                 Convex,CoordModeOrigin);
  }
#endif

#ifdef MAC_GRAPHICS
  if( doing_Mac ){
#if 0
    if(num_points > 3)return;
    if(firstCall){
      thePoly = OpenPoly( );        /* start recording */
      MoveTo((short)(GRAPHICS_SCALE*points[0].x),
             (short)(GRAPHICS_SCALE*points[0].y));
      LineTo((short)(GRAPHICS_SCALE*points[1].x),
             (short)(GRAPHICS_SCALE*points[1].y));
      LineTo((short)(GRAPHICS_SCALE*points[2].x),
             (short)(GRAPHICS_SCALE*points[2].y));
      LineTo((short)(GRAPHICS_SCALE*points[0].x),
             (short)(GRAPHICS_SCALE*points[0].y));
      ClosePoly();                /* stop recording */
      firstCall = 0;
    }else{
      /* just directly fill the polyPoints array in the poly record */
      (*thePoly)->polyPoints[0].h = (short)(GRAPHICS_SCALE*points[0].x);
      (*thePoly)->polyPoints[0].v = (short)(GRAPHICS_SCALE*points[0].y);
      (*thePoly)->polyPoints[1].h = (short)(GRAPHICS_SCALE*points[1].x);
      (*thePoly)->polyPoints[1].v = (short)(GRAPHICS_SCALE*points[1].y);
      (*thePoly)->polyPoints[2].h = (short)(GRAPHICS_SCALE*points[2].x);
      (*thePoly)->polyPoints[2].v = (short)(GRAPHICS_SCALE*points[2].y);
      /* set the bounding box */
      max_x = (short)(GRAPHICS_SCALE*points[0].x);
      min_x = (short)(GRAPHICS_SCALE*points[0].x);
      max_y = (short)(GRAPHICS_SCALE*points[0].y);
      min_y = (short)(GRAPHICS_SCALE*points[0].y);

      if((short)(GRAPHICS_SCALE*points[1].x) > max_x)
        max_x = (short)(GRAPHICS_SCALE*points[1].x);
      else if((short)(GRAPHICS_SCALE*points[1].x) < min_x)
        min_x = (short)(GRAPHICS_SCALE*points[1].x);
      if((short)(GRAPHICS_SCALE*points[1].y) > max_y)
        max_y = (short)(GRAPHICS_SCALE*points[1].y);
      else if((short)(GRAPHICS_SCALE*points[1].y) < min_y)
        min_y = (short)(GRAPHICS_SCALE*points[1].y);

      if((short)(GRAPHICS_SCALE*points[2].x) > max_x)
        max_x = (short)(GRAPHICS_SCALE*points[2].x);
      else if((short)(GRAPHICS_SCALE*points[2].x) < min_x)
        min_x = (short)(GRAPHICS_SCALE*points[2].x);
      if((short)(GRAPHICS_SCALE*points[2].y) > max_y)
        max_y = (short)(GRAPHICS_SCALE*points[2].y);
      else if((short)(GRAPHICS_SCALE*points[2].y) < min_y)
        min_y = (short)(GRAPHICS_SCALE*points[2].y);

      (*thePoly)->polyBBox.top = min_y;
      (*thePoly)->polyBBox.bottom = max_y;
      (*thePoly)->polyBBox.left = min_x;
      (*thePoly)->polyBBox.right = max_x;

    }
    PaintPoly( thePoly);        /* paint the inside */
#endif
#if 0
    thePoly = OpenPoly( );        /* start recording */
    MoveTo((short)(GRAPHICS_SCALE*points[0].x),
           (short)(GRAPHICS_SCALE*points[0].y) );
    LineTo((short)(GRAPHICS_SCALE*points[1].x),
           (short)(GRAPHICS_SCALE*points[1].y) );
    LineTo((short)(GRAPHICS_SCALE*points[2].x),
           (short)(GRAPHICS_SCALE*points[2].y) );
    ClosePoly();                /* stop recording */
    PaintPoly( thePoly);
    KillPoly(thePoly);
#endif
    if( firstCall ){
      theRgn = NewRgn();
      firstCall = 0;
    }
    OpenRgn();
    MoveTo((short)(GRAPHICS_SCALE*points[0].x),
           (short)(GRAPHICS_SCALE*points[0].y) );
    LineTo((short)(GRAPHICS_SCALE*points[1].x),
           (short)(GRAPHICS_SCALE*points[1].y) );
    LineTo((short)(GRAPHICS_SCALE*points[2].x),
           (short)(GRAPHICS_SCALE*points[2].y) );
    LineTo((short)(GRAPHICS_SCALE*points[0].x),
           (short)(GRAPHICS_SCALE*points[0].y) );
    CloseRgn(theRgn);                /* stop recording */
        if(colorscale > 0 ){
    the_color.red = the_color.green = the_color.blue =
      (short)((float)0xFFFF*(float)colorscale);
      } else{
          the_color.red = the_color.green = the_color.blue = 0;
          }


    RGBForeColor(&the_color);


    PaintRgn( theRgn);
    /*DisposeRgn(theRgn);*/

  }
  RGBForeColor(&(Maclinestyles[0]));
#endif


}


/**************

  An open polygon (convex only please)

  for simplicities sake, these are drawn with lines.
  to make this work, 'points should be 'num_points+1
  elements long, with the first and last elements
  identical.

  ***************/
void g_open_polygon(XPoint *points,int num_points)
{
  int i;

  if( doing_ps ){
    if( num_points == 3 ){
      fprintf(psfile,"%d %d %d %d %d %d %d %d dT\n",
              points[0].x,points[0].y,points[1].x,points[1].y,
              points[2].x,points[2].y,points[3].x,points[3].y);
    } else{
      fprintf(psfile,"%d %d M ",points[0].x,points[0].y);
      for(i=1;i<num_points;i++){
        fprintf(psfile,"%d %d L ",points[i].x,points[i].y);
      }
      fprintf(psfile,"%d %d L stroke\n",points[0].x,points[0].y);
    }
  }

#ifdef X_GRAPHICS
  if( doing_X ){
    XDrawLines(disp,gpix,graphgc,points,num_points+1,CoordModeOrigin);
  }
#endif

#ifdef MAC_GRAPHICS
  if( doing_Mac ){
    MoveTo((short)(GRAPHICS_SCALE*points[0].x),
           (short)(GRAPHICS_SCALE*points[0].y));
    for(i=1;i<num_points;i++){
      LineTo((short)(GRAPHICS_SCALE*points[i].x),
             (short)(GRAPHICS_SCALE*points[i].y));
    }
    LineTo((short)(GRAPHICS_SCALE*points[0].x),
           (short)(GRAPHICS_SCALE*points[0].y));
  }
#endif
}



/**************

  The X legend

  ***************/
void g_xlegend(float x,float y,char *text)
{
  int textloc;
  int yloc;

  if( doing_ps ){
    fprintf(psfile,"%.1lf %.1lf textsize 2.5 mul add M\n",x,y);
    ENHPS_put_text((int)x,(int)y,text,CENTER_JUST);
  }



#ifdef X_GRAPHICS
  if( doing_X ){
    textloc = (int)x - XTextWidth(big_font,text,strlen(text))/2;
    yloc = (int)y + (int)3.0*big_font->ascent;
    XDrawString(disp,gpix,bigtextgc,textloc,
                yloc,text,strlen(text));
    if( yloc > localmax.y ) localmax.y = y;
  }
#endif

#ifdef TEK_GRAPHICS
  if(doing_tek){
    yloc = (int)(y+2.0*TEK_VCHAR);
    TEK_center_text((unsigned int)x,(unsigned int)(g_ymax - yloc),text);
    if( yloc > localmax.y ) localmax.y = y;
  }
#endif

#ifdef MAC_GRAPHICS
  if(doing_Mac){
    textloc = (int)(GRAPHICS_SCALE*x) -
      (int)(TextWidth(text,0,strlen(text))/2);
    MoveTo((short)textloc,(short)(GRAPHICS_SCALE*y +
                                  3*Mac_globals.fontSize));
    DrawText(text,0,strlen(text));
  }
#endif

}

/**************

  The Y legend

  ***************/
void g_ylegend(float x,float y,char *text)
{
  int textloc;

  if( doing_ps ){
    fprintf(psfile,"%.1lf %.1lf M\n",x,y);
    fprintf(psfile,"-90 rotit\n");
    fprintf(psfile,"0 -2 textsize mul M\n");
    ENHPS_put_text((int)x,(int)y,text,CENTER_JUST);
    fprintf(psfile,"unrot\n");
  }

#ifdef X_GRAPHICS
  if( doing_X ){
    textloc = (int)x - XTextWidth(big_font,text,strlen(text));
    XDrawString(disp,gpix,bigtextgc,textloc,(int)y,text,strlen(text));
    if( textloc < localmin.x ) localmin.x = textloc;
  }
#endif

#ifdef TEK_GRAPHICS
  if(doing_tek){
    TEK_right_text((unsigned int)x,g_ymax - (unsigned int)y,text);
  }
#endif

#ifdef MAC_GRAPHICS
  if(doing_Mac){
    textloc = (int)(GRAPHICS_SCALE*x) -
      (int)(TextWidth(text,0,strlen(text)));
    MoveTo((short)textloc,(short)(GRAPHICS_SCALE*y +
                                  Mac_globals.fontSize));
    DrawText(text,0,strlen(text));
  }
#endif


}


/********

  The title

  ********/
void g_title(float x,float y,char title[NUM_TITLE_LINES][NORMAL_STR_LEN])
{
  int i,lineshown;
  int xloc,yloc;

  if( doing_ps ){
    lineshown = 0;
    for(i=NUM_TITLE_LINES-1;i>=0;i--){
      if( title[i][0] != 0 ){
        lineshown++;
        fprintf(psfile,"%.1lf %.1lf 1.8 textsize mul %.1lf mul sub M\n",
                x,y,(float)lineshown);
        ENHPS_put_text((int)x,(int)y,title[i],CENTER_JUST);
      }
    }
  }

#ifdef X_GRAPHICS
  if( doing_X ){
    yloc = (int)y;
    for(i=NUM_TITLE_LINES-1;i>=0;i--){
      if( title[i][0] != 0 ){
        yloc -= big_font->ascent;
        xloc = (int)x - XTextWidth(big_font,title[i],strlen(title[i])/2);
        XDrawString(disp,gpix,bigtextgc,xloc,
                    yloc,title[i],strlen(title[i]));
        if(yloc < localmin.y) localmin.y = yloc;
      }
    }
  }
#endif

#ifdef TEK_GRAPHICS
  if( doing_tek ){
    yloc = (int)y;
    for(i=NUM_TITLE_LINES-1;i>=0;i--){
      if( title[i][0] != 0 ){
        yloc -= TEK_VCHAR;
        TEK_center_text((unsigned int)x,g_ymax -
                        (unsigned int)yloc,title[i]);
      }
    }
  }
#endif


#ifdef MAC_GRAPHICS
  if(doing_Mac){
    yloc = (int)(GRAPHICS_SCALE*y);
    for(i=NUM_TITLE_LINES-1;i>=0;i--){
      if(title[i][0] != 0){
        xloc = (int)(GRAPHICS_SCALE*x) -
          (int)(TextWidth(title[i],0,strlen(title[i]))/2);
        yloc -= Mac_globals.fontSize;
        MoveTo((short)xloc,(short)yloc);
        DrawText(title[i],0,strlen(title[i]));
        if(yloc < localmin.y) localmin.y = yloc;
      }
    }
  }
#endif

}


/********

  Labels

  ********/
void g_label(float x,float y,label_type *label)
{
  int i;
  char *string;
  float textloc;
  float twidth,theight;
  float xp1,xp2,yp1,yp2;

  string = label->text[0];

  if( doing_ps ){
    fprintf(psfile,"%.1lf %.1lf textsize -.5 thesize div mul add M\n",x,y);
    ENHPS_put_text((int)x,(int)y,string,CENTER_JUST);
    if( label->show_lines ){
      g_change_linestyle(1);
      for(i=0;i<label->num_atoms_labelled;i++){
        xp1 = x;
        yp1 = y;
        xp2 = label->atoms_to_label[i]->screen_loc.x;
        yp2 = label->atoms_to_label[i]->screen_loc.y;
        if( yp2 > yp1 ){
          fprintf(psfile,"%.1lf %.1lf M %.1f %.1lf L SNP\n",
                  xp1,yp1,xp2,yp2);
        } else{
          fprintf(psfile,
                  "%.1lf %.1lf textsize -1 thesize div mul add M %.1f %.1lf L SNP\n",
                  xp1,yp1,xp2,yp2);
        }
      }
      g_change_linestyle(0);
    }

  }

#ifdef  X_GRAPHICS
  if( doing_X ){
    twidth = XTextWidth(big_font,string,strlen(string));
    theight = big_font->ascent;
    textloc = x - twidth/2;
    XDrawString(disp,gpix,bigtextgc,(int)textloc,
                (int)y+big_font->ascent,string,strlen(string));
    if( x + twidth/2 > localmax.x )
      localmax.x = x+twidth/2;
    if( x - twidth/2 < localmin.x )
      localmin.x = x-twidth/2;
    if( y+theight > localmax.y ) localmax.y = y+theight;
    if( y < localmin.y ) localmin.y = y;

    if( label->show_lines ){
      g_change_linestyle(1);
      for(i=0;i<label->num_atoms_labelled;i++){
        xp1 = x;
        yp1 = y;
        xp2 = label->atoms_to_label[i]->screen_loc.x;
        yp2 = label->atoms_to_label[i]->screen_loc.y;
        if( yp2 > yp1 ) yp1 = localmax.y;
        XDrawLine(disp,gpix,graphgc,(short)xp1,(short)yp1,
                  (short)xp2,(short)yp2);
      }
      g_change_linestyle(0);
    }
  }
#endif

#ifdef TEK_GRAPHICS
  if( doing_tek ){
    TEK_center_text((unsigned int)x,g_ymax - (unsigned int)y,string);
    if( x+(TEK_HCHAR*strlen(string))/2 > localmax.x )
      localmax.x = x+(TEK_HCHAR*strlen(string))/2;
    if( x-(TEK_HCHAR*strlen(string))/2 < localmin.x )
      localmin.x = x-(TEK_HCHAR*strlen(string))/2;
    if( y > localmax.y ) localmax.y = y;
    if( y < localmin.y ) localmin.y = y;
  }

#endif

#ifdef MAC_GRAPHICS
  /* !!! needs to be updated */
  if(doing_Mac){
    textloc = (int)(GRAPHICS_SCALE*x) -
      (int)(TextWidth(string,0,strlen(string))/2);
    MoveTo((short)textloc,(short)(GRAPHICS_SCALE*y +
                                  Mac_globals.fontSize));
    DrawText(string,0,strlen(string));
  }
#endif

}




void g_initgraphics(int xsize,int ysize)
{

  /* set up the PS options */
  strcpy(PS_options.fontname,DEF_PS_FONT);
  PS_options.where_to_print = PS_PRINT_BOTTOM;
  PS_options.fontsize = DEF_PS_FONT_SIZE;
  PS_options.printscale = DEF_PS_SCALE;
#ifdef X_GRAPHICS
  if( doing_X ){
    X_initgraphics(xsize,ysize);
  }
#endif
#ifdef TEK_GRAPHICS
  if( doing_tek ){
    Tek_file = stdout;
    g_xmax = TEK_XMAX;
    g_ymax = TEK_YMAX;
  }
#endif
#ifdef MAC_GRAPHICS
  if( doing_Mac ){
    /**********

      use the EasyApp stuff to open the windows and fire up the
      menus

      **********/
    Mac_initgraphics();
  }
#endif


}
/*********

  this clears the graphics screen

  *********/
void g_clear_screen(void)
{

#ifdef X_GRAPHICS
  if( doing_X ){
    XFillRectangle(disp,gpix,blackgc,0,0,g_xmax,g_ymax);
  }
#endif

#ifdef TEK_GRAPHICS
  if( doing_tek ){
    TEK_init();
    TEK_graphics();
  }
#endif

#ifdef MAC_GRAPHICS
  /*************

    on the Mac, we need to do some magic to make things work.
    this is because of the fact that we're doing the drawing
    to an offscreen graphics world.

    It's more convenient to have this definition in one
    of the Mac specific files.

    ***************/
  if(doing_Mac){
    Mac_ClearScreen();
  }
#endif
}

/*********

  This is for double buffered drawing

  *********/
void g_switch_buffers(void)
{
#ifdef X_GRAPHICS
  if( doing_X ){
    XCopyArea(disp,gpix,gwin,graphgc,0,0,g_xmax,g_ymax,0,0);
    refresh_all_colormaps = 0;
  }
#endif

#ifdef MAC_GRAPHICS
  if( doing_Mac ){
    /* again, some macintosh magic (TM) is required to make this work */
    Mac_CopyGWorld();
  }
#endif

}


/*********

  This is for putting comments in output files to make them
  easier to hand edit

  *********/
void g_insert_comment(char *string)
{
  if(doing_ps){
    fprintf(psfile,"%% %s\n",string);
  }

}

