
/***************

  This is adapted from the file enhpost.trm for gnuplot v3.5

  It contains stuff for doing fancy things with postscript.

    This file is modified and distributed with the permission of
    David Denholm and Matt Heffron.

    If there are problems, please report them to me.

    Adapted by greg Landrum, July 1995

    Changes include:
      -removal of rotated text support
        This was taken out for user interface issues....
         as soon as those are worked out, putting the rotated
         text back in would be a damn good idea.
      -ANSI'ification of the code, including prototypes in
         "prototypes.h"

****************/


/*
 * Copyright (C) 1990 - 1993
 *
 * Permission to use, copy, and distribute this software and its
 * documentation for any purpose with or without fee is hereby granted,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.
 *
 * Permission to modify the software is granted, but not the right to
 * distribute the modified code.  Modifications are to be distributed
 * as patches to released version.
 *
 * This software  is provided "as is" without express or implied warranty.
 *
 * This file is included by ../term.c.
 *
 * This terminal driver supports:
 *     enhpost                Enhanced PostScript
 * It is an extension to the postscript terminal driver and depends on
 * that being included (in ../term.c) before this is included.
 *
 * AUTHORS
 *  Russell Lang
 *  David Denholm
 *  Matt Heffron
 *
 * DIV WAS HERE - we use different text routines that understand ^, _
 * and {<fontname> ... }
 * ^, _ and @ take only one character, unless you enclose stuff in {...}
 * Its up to the user to put in two \'s, so that we can use \ to specify
 * control charcters (eg "{/Symbol \245}"  gives infinity)
 * DIV was here AGAIN ! - teach it some more point types.
 * Div's latest addition - @{text} gives a 'phantom box' containing text, but
 * taking up no space, so @{text}more text will _overwrite_ text with
 * more text - hopefully useful for hats and such things..?
 * It's all done with recursion, so you can change font inside the phantom
 * box, etc.
 *
 * send your comments or suggestions to (info-gnuplot@dartmouth.edu).
 *
 * ALMOST the same as 'postscript' except for string handling.
 */

/***
  Recent Edit History:
    03.05.98 gL:
      fixes to the tight bounding box code.
    18.05.98 gL:
      changed dT (draw triangle) definition
      in ps header to use roundcapped lines... this
      generally looks better.

      Added SNG definition to PS header to allow easy
      adjustment of max/min values of shading

    21.05.98 gL:
      added Cshadecirc, shadecirc, Cstoptube, and stoptube
      to PS header to clean up the PS a little bit, and
      add the possibility for shaded circles.

    29.05.98 gL:
      fixed typo in shadecirc def

    05.09.98 gL:
      less circles in shadecirc and Cshadecirc defs
   08.09.98 gL:
     support for colored and/or shaded atoms
   16.09.98 gL:
     updated arcwid printing to be smaller with shading.
   26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)
     There was also a stupid mistake with bond colors... that's fixed now.
***/

#include "viewkel.h"

#define PS_XOFF        50        /* page offset in pts */
#define PS_YOFF        50

#define PS_XMAX 7200
#define PS_YMAX 5040

#define PS_SC 10

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

char *ENHPS_header[] = {
"/M {moveto} bind def\n",
"/L {lineto} bind def\n",
"/R {rmoveto} bind def\n",
"/V {rlineto} bind def\n",
"/vpt2 vpt 2 mul def\n",
"/hpt2 hpt 2 mul def\n",

/* For MFshow and MFwidth the tos is an array with the string and font info:  */
/*        [<fontname (a string)> <fontsize> <vertical offset> <width significant?> <text string>]  */

"/MFshow {{dup dup 0 get findfont exch 1 get scalefont setfont\n",
"     [ currentpoint ] exch dup 2 get 0 exch rmoveto dup 4 get show dup\n",
"     3 get {2 get neg 0 exch rmoveto pop} {pop aload pop moveto}ifelse} forall} bind def\n",
"/MFwidth {0 exch {dup 3 get{dup dup 0 get findfont exch 1 get scalefont setfont\n",
"      4 get stringwidth pop add}\n",
"    {pop} ifelse} forall} bind def\n",

/* flush left show */
"/Lshow { currentpoint stroke M\n",
"  0 exch R MFshow } bind def\n",

/* flush right show */
"/Rshow { currentpoint stroke M\n",
"  exch dup MFwidth neg 3 -1 roll R MFshow } def\n",

/* centred show */
"/Cshow { currentpoint stroke M\n",
"  exch dup MFwidth -2 div 3 -1 roll R MFshow } def\n",

NULL
};

/* add any extra PS definitions here... make sure NULL is the last
  thing in the array */
char *myPSdefns[]={
"% This draws a shaded circle with a highlight.  It is ugly, yet\n",
"%  functional.  Usage:  [r g b x y rad] Cshadecirc\n",
"/Cshadecirc { \n",
"/thearray 1 index def\n",
"/initoffset 0.3 def\n",
"/maxoffset 0.25 def\n",
"/minradscale 0.2 def\n",
"/lightang 150 def\n",
"/numcirc 10 def\n",
"/rv thearray 0 get def\n",
"/gv thearray 1 get def\n",
"/bv thearray 2 get def\n",
"/xp thearray 3 get def\n",
"/yp thearray 4 get def\n",
"/radv thearray 5 get def\n",
"/rvd 1 rv sub numcirc 1 add div def\n",
"/gvd 1 gv sub numcirc 1 add div def\n",
"/bvd 1 bv sub numcirc 1 add div def\n",
"/radvd radv radv minradscale mul sub numcirc div def\n",
"/xpd radv maxoffset mul lightang cos mul numcirc div def\n",
"/ypd radv maxoffset mul lightang sin mul numcirc div def\n",
"/sxp xp initoffset radv mul lightang cos mul add def\n",
"/syp yp initoffset radv mul lightang sin mul add def\n",
"% uncomment this to get black outlines around shaded circles\n",
"%1 setlinewidth 0 setgray xp yp radv 0 360 arc stroke\n",
"rv gv bv setrgbcolor xp yp radv 0 360 arc fill\n",
"gsave\n",
"newpath xp yp radv 0 360 arc clip\n",
"newpath\n",
"/inc 0 def\n",
"1 1 numcirc 1 sub {\n",
"pop \n",
"/inc 1 inc add def\n",
"newpath\n",
"rv rvd inc mul add gv gvd inc mul add bv bvd inc mul add setrgbcolor \n",
"sxp xpd inc mul add syp ypd inc mul add \n",
"radv radvd inc mul sub\n",
"0 360 arc fill\n",
"} for\n",
"grestore\n",
"pop\n",
"} def\n",
"\n",
"/shadecirc { \n",
"/thearray 1 index def\n",
"/initoffset 0.3 def\n",
"/maxoffset 0.25 def\n",
"/minradscale 0.2 def\n",
"/lightang 150 def\n",
"/numcirc 10 def\n",
"/grayv thearray 0 get def\n",
"/xp thearray 1 get def\n",
"/yp thearray 2 get def\n",
"/radv thearray 3 get def\n",
"/grayvd 1 grayv sub numcirc 1 add div def\n",
"/radvd radv radv minradscale mul sub numcirc div def\n",
"/xpd radv maxoffset mul lightang cos mul numcirc div def\n",
"/ypd radv maxoffset mul lightang sin mul numcirc div def\n",
"/sxp xp initoffset radv mul lightang cos mul add def\n",
"/syp yp initoffset radv mul lightang sin mul add def\n",
"% uncomment this to get black outlines around shaded circles\n",
"%1 setlinewidth 0 setgray xp yp radv 0 360 arc stroke\n",
"grayv setgray xp yp radv 0 360 arc fill\n",
"gsave\n",
"newpath xp yp radv 0 360 arc clip\n",
"newpath\n",
"/inc 0 def\n",
"1 1 numcirc 1 sub {\n",
"pop \n",
"/inc 1 inc add def\n",
"newpath\n",
"grayv grayvd inc mul add setgray\n",
"sxp xpd inc mul add syp ypd inc mul add \n",
"radv radvd inc mul sub\n",
"0 360 arc fill\n",
"} for\n",
"grestore\n",
"pop\n",
"} def\n",
"\n",
"% [r g b w x1 y1 x2 y2 xs ys] Cstoptube\n",
"/Cstoptube {\n",
"/thearray 1 index def\n",
"/rv thearray 0 get  def\n",
"/gv thearray 1 get def\n",
"/bv thearray 2 get def\n",
"/thewidth thearray 3 get def\n",
"/startx thearray 4 get def\n",
"/starty thearray 5 get def\n",
"/endx thearray 6 get def\n",
"/endy thearray 7 get def\n",
"/stopx thearray 8 get def\n",
"/stopy thearray 9 get def\n",
"thewidth 8 add setlinewidth\n",
"0 0 0 setrgbcolor\n",
"1 setgray\n",
"0 setlinecap\n",
"stopx stopy moveto endx endy lineto stroke\n",
"thewidth 4 add setlinewidth\n",
"0 setgray\n",
"1 setlinecap\n",
"startx starty moveto endx endy lineto stroke\n",
"thewidth setlinewidth\n",
"rv gv bv setrgbcolor\n",
"1 setlinecap\n",
"startx starty moveto endx endy lineto stroke\n",
"0 0 0 setrgbcolor\n",
"0 setgray\n",
"pop\n",
"} def\n",
"\n",
"% [gray w x1 y1 x2 y2 xs ys] stoptube\n",
"/stoptube {\n",
"/thearray 1 index def\n",
"/grayv thearray 0 get  def\n",
"/thewidth thearray 1 get def\n",
"/startx thearray 2 get def\n",
"/starty thearray 3 get def\n",
"/endx thearray 4 get def\n",
"/endy thearray 5 get def\n",
"/stopx thearray 6 get def\n",
"/stopy thearray 7 get def\n",
"thewidth 8 add setlinewidth\n",
"1 setgray\n",
"0 setlinecap\n",
"stopx stopy moveto endx endy lineto stroke\n",
"thewidth 4 add setlinewidth\n",
"0 setgray\n",
"1 setlinecap\n",
"startx starty moveto endx endy lineto stroke\n",
"thewidth setlinewidth\n",
"grayv setgray\n",
"1 setlinecap\n",
"startx starty moveto endx endy lineto stroke\n",
"0 0 0 setrgbcolor\n",
"0 setgray\n",
"pop\n",
"} def\n",

NULL};



/* added by Matt Heffron <heffron@falstaff.css.beckman.com> */
struct ENHPS_FontName {
        char *name;
        struct ENHPS_FontName *next;
} *ENHPS_DocFonts = NULL;

char ps_font[120];


char *ENHPS_RememberFont(char *fname)
{
  struct ENHPS_FontName *fnp;

  for (fnp=ENHPS_DocFonts; fnp && strcmp(fnp->name, fname); fnp = fnp->next);
  if (fnp)
    return fnp->name;        /* we must have found it in the list */

  fnp = (struct ENHPS_FontName *)D_MALLOC(sizeof(struct ENHPS_FontName));
  if( !fnp ) fatal("Can't get fnp memory");
  fnp->name = D_CALLOC(1+strlen(fname),sizeof(char));
  if( !fnp->name ) fatal("Can't get fnp->name");
  strcpy(fnp->name, fname);
  fnp->next = ENHPS_DocFonts;
  ENHPS_DocFonts = fnp;
  return fnp->name;
}


void ENHPS_init(void)
{
  strcpy(ps_font,PS_options.fontname);

  (void)ENHPS_RememberFont(ps_font);

}


void ENHPS_reset(void)
{
  fprintf(psfile,"%%%%Trailer\n");
  fprintf(psfile,"%%%%DocumentFonts: ");
  while (ENHPS_DocFonts) {
    struct ENHPS_FontName *fnp;
    fnp = ENHPS_DocFonts->next;
    fprintf(psfile, "%s%s", ENHPS_DocFonts->name, fnp ? ", " : "\n");
    D_FREE(ENHPS_DocFonts->name);
    D_FREE(ENHPS_DocFonts);
    ENHPS_DocFonts=fnp;
  }

  fprintf(psfile,"%%%%Pages: %d\n",1);
}


/* take out comments to output debugging info */

#define ENHPS_DEBUG(x)  /*printf x;*/

static char ENHps_opened_string;  /* try to cut out empty ()'s */

/* used in determining height of processed text */

static float ENHps_max_height, ENHps_min_height;


/*{{{  static char *ENHPS_recurse(char *p, int brace, ... )*/
/* process a bit of string, and return the last character used.
 * p is start of string
 * brace is TRUE to keep processing to }, FALSE for do one character
 * fontname & fontsize are obvious
 * base is the current baseline
 * widthflag is TRUE if the width of this should count,
 *              FALSE for zero width boxes
 */



char *ENHPS_recurse(char *p,char brace, char *fontname,
                           float fontsize, float base, char widthflag)
{

  /* close a postscript string if it has been opened */
#define ENHPS_FLUSH      \
{        if (ENHps_opened_string)  \
    {        fputs(")]\n", psfile);   \
      ENHps_opened_string = FALSE; \
}                         \
}

#define ENHPS_OPEN        \
{        if (!ENHps_opened_string) \
    { fprintf(psfile, "[(%s) %.1f %.1f %s (",  \
              fontname, fontsize, base, \
              widthflag ? "true" : "false");  \
                ENHps_opened_string = TRUE; \
}        \
}

  ENHPS_DEBUG(("RECURSE WITH [%p] %s, %d %s %.1f %.1f %d\n", p, p, brace, fontname, fontsize, base, widthflag))

    /* Start each recursion with a clean string */
    ENHPS_FLUSH

      if (base + fontsize > ENHps_max_height)
          {        ENHps_max_height = base + fontsize;
                ENHPS_DEBUG(("Setting max height to %.1f\n", ENHps_max_height));
              }

  if (base < ENHps_min_height)
      {        ENHps_min_height = base;
        ENHPS_DEBUG(("Setting min height to %.1f\n", ENHps_min_height));
      }

  for ( ; *p; ++p)
      {        float shift;
        float f=0;                /* used for getting new font size */
        char *localfontname, ch;

        /*{{{  look for 'special' characters - process then 'continue'*/
        switch (*p)
            {
            case '}'  :
              /*{{{  deal with it*/
              if (brace)
                return (p);

              fprintf(stderr, "enhpost printer driver - spurious }\n");
              break;
              /*}}}*/

            case '_'  :
            case '^'  :
              /*{{{  deal with super/sub script*/

              shift = (*p == '^') ? 0.5 : -0.3;

              ENHPS_FLUSH

                p = ENHPS_recurse(p+1, FALSE, fontname, fontsize*0.8, base+shift*fontsize, widthflag);

              break;
              /*}}}*/

            case '{'  :
              /*{{{  recurse (possibly with a new font) */
              ENHPS_DEBUG(("Dealing with {\n"))

                if (*++p == '/')
                    {                /* then parse a fontname, optional fontsize */
                      while (*++p == ' ');
                      localfontname = p;
                      while ((ch = *p) > ' ' && ch != '=')
                        ++p;
                      if (ch == '=')
                          {
                            *p++ = '\0';
                            /*{{{  get optional font size*/
                            ENHPS_DEBUG(("Calling strtod(%s) ...", p))
                              f = (float)strtod(p, &p);
                            ENHPS_DEBUG(("Retured %.1f and %s\n", f, p))

                              if (f)
                                f *= PS_SC; /* remember the scaling */
                              else
                                f = fontsize;

                            ENHPS_DEBUG(("Font size %.1f\n", f))
                              /*}}}*/
                          }
                      else
                          {
                            *p++ = '\0';
                            f = fontsize;
                          }

                      while (*p == ' ')
                        ++p;
                      if (*localfontname)
                        localfontname = ENHPS_RememberFont(localfontname);
                      else
                        localfontname = fontname;
                    }
                else
                    {
                      localfontname = fontname;
                      f = fontsize;
                    }
              /*}}}*/

              ENHPS_DEBUG(("Before recursing, we are at [%p] %s\n", p, p))

                p = ENHPS_recurse(p, TRUE, localfontname, f, base, widthflag);

              ENHPS_DEBUG(("BACK WITH %s\n", p));

              ENHPS_FLUSH

                break;
              /*}}}*/

            case '@' :
              /*{{{  phantom box - prints next 'char', then restores currentpoint */

              ENHPS_FLUSH

                p = ENHPS_recurse(++p, FALSE, fontname, fontsize, base, FALSE);

              break;
              /*}}}*/

            case '('  :
            case ')'  :
              /*{{{  an escape and print it */
              /* special cases */
              ENHPS_OPEN
                fputc('\\', psfile);
              fputc(*p, psfile);
              break;
              /*}}}*/

            case '\\'  :
              /*{{{  is it an escape */
              /* special cases */

              if (p[1]=='\\' || p[1]=='(' || p[1]==')')
                  {
                    ENHPS_OPEN
                      fputc('\\', psfile);
                  }
              else if ((ch = p[1]) >= '0' && ch <= '7')
                  {
                    /* up to 3 octal digits */
                    ENHPS_OPEN
                      fputc('\\', psfile);
                    fputc(ch, psfile);
                    ++p;
                    if ((ch = p[1]) >= '0' && ch <= '7')
                        {
                          fputc(ch, psfile);
                          ++p;
                          if ((ch = p[1]) >= '0' && ch <= '7')
                              {
                                fputc(ch, psfile);
                                ++p;
                              }
                        }
                    break;
                  }
              else if ((ch = p[1]) == 'A' && p[2] == 'A'){
                /* it was \AA, put in an Angstrom symbol */
                ENHPS_OPEN
                fputs("\\305",psfile);
                p+=3;
              }

              ++p;
              /* just go and print it (fall into the 'default' case) */

              /*}}}*/
            default:
              /*{{{  print it */
              ENHPS_OPEN

                fputc(*p, psfile);

              /*}}}*/

            }
        /*}}}*/

        /* like TeX, we only do one character in a recursion, unless it's
         * in braces
         */

        if (!brace)
            {
              ENHPS_FLUSH
                return(p);        /* the ++p in the outer copy will increment us */
            }

      }
  ENHPS_FLUSH
    return p;
}
/*}}}*/

/*{{{  ENHPS_put_text(x,y,str) */
int ENHPS_put_text(unsigned int x, unsigned int y,char *str,
                    char horiz_justify)
{
  /* flush any pending graphics (all the XShow routines do this...) */

  if (!strlen(str))
    return 0;

  fputs("[ ",psfile);

  /* set up the globals */

  ENHps_opened_string = FALSE;
  ENHps_max_height = -1000;
  ENHps_min_height = 1000;

  while (*(str = ENHPS_recurse(str, TRUE, ps_font,
                               (float)PS_options.fontsize,
                               (float)0.0, TRUE)));

  ENHps_max_height += ENHps_min_height;

  fprintf(psfile, "] %.1f ", -ENHps_max_height/3);

  switch(horiz_justify)
      {
      case LEFT_JUST : fprintf(psfile, "Lshow\n");
        break;
      case CENTER_JUST : fprintf(psfile, "Cshow\n");
        break;
      case RIGHT_JUST : fprintf(psfile, "Rshow\n");
        break;
      }

}
/*}}}*/


/****************************************************************************
 *
 *                Procedure open_ps_button_window
 *
 * Arguments:  none
 *
 * Returns: none
 *
 * Action:  opens up the ps options button window
 *
 ****************************************************************************/
void open_ps_button_window(void)
{
  build_PS_options_button_window(&button_wins);
}


/****************************************************************************
 *
 *                   Procedure do_ps_output
 *
 * Arguments: none
 * Returns: none
 *
 * Action:  dumps whatever is being currently displayed to a ps file.
 *
 *
 ****************************************************************************/
void do_ps_output(void)
{
  int i;
  char *stringarr[3],colorstring[240];
  char filename[80];
  char *theinline;
  float bbx1,bby1,bbx2,bby2;
  float translation;


  /* find the bounding box for the image */
  bbx1 = (float)globalmin.x;
  bbx2 = (float)globalmax.x;
  bby2 = (float)globalmin.y;
  bby1 = (float)globalmax.y;

  display("Look in xterm!");
#ifndef USE_READLINE
  printf("\nEnter postscript file name: ");
  scanf("%s\n",filename);
#else
    theinline= readline("Enter the postscript file name: ");
    add_history(theinline);
    if( theinline ){
      sscanf(theinline,"%s",filename);
      free(theinline);
    } else {
      error("Bad file name");
      filename[0] = 0;
    }
#endif

  psfile = fopen(filename,"w+");
  if(!psfile){
    error("Can't open output file.");
    return;
  }

  if( !PS_options.tight_b_box ){
    switch(PS_options.where_to_print){
    case PS_PRINT_TOP:
      translation = -1.8*(float)g_ymax;
      bby1 = 360.0;
      bby2 = 745.0;
      break;
    case PS_PRINT_MIDDLE:
      translation = -1.5*(float)g_ymax;
      bby1 = 240.0;
      bby2 = 625.0;
      break;
    case PS_PRINT_BOTTOM:
      translation = -(float)g_ymax;
      bby1 = 40.0;
      bby2 = 425.0;
    }
    bbx1 = 0;
    bbx2 = 640;
  } else{
#if 0
    bby1 -= (float)g_ymax;
    bby2 -= (float)g_ymax;
#endif
    /* scale the bounding box */
    bbx1 *= 8.5*72.0*PS_options.printscale/(float)g_xmax;
    bbx2 *= 8.5*72.0*PS_options.printscale/(float)g_xmax;
    bby1 *= -8.5*72.0*PS_options.printscale/(float)g_xmax;
    bby2 *= -8.5*72.0*PS_options.printscale/(float)g_xmax;
    /* now increase its size by 10% in each direction */
    printf("box: %6.2lf %6.2lf %6.2lf %6.2lf\n",
           bbx1,bby1,bbx2,bby2);
    bbx1 -= .1*fabs(bbx2-bbx1);
    bbx2 += .1*fabs(bbx2-bbx1);
    bby1 -= .1*fabs(bby2-bby1);
    bby2 += .1*fabs(bby2-bby1);
    translation = 0.0;
  }
  /* write out some header information */
  fprintf(psfile,"%%!PS-Adobe-2.0 EPSF-1.2\n");
  fprintf(psfile,"%%%%DocumentFonts: (atend) \n");
  fprintf(psfile,"%%%%Creator: Viewkel\n");
  fprintf(psfile,"%%%%Pages: 1\n");
  fprintf(psfile,"%%%%BoundingBox: %6.2lf %6.2lf %6.2lf %6.2lf\n",
          bbx1,bby1,bbx2,bby2);
  fprintf(psfile,"%%%%EndComments\n\n");


  ENHPS_init();

  fprintf(psfile,"%% This is hacked straight out of the Blue book (from Adobe) \n");
  fprintf(psfile,"/centerdash\n");
  fprintf(psfile,"  { /pattern exch def\n");
  fprintf(psfile,"    /pathlen pathlength def\n");
  fprintf(psfile,"    /patternlength 0 def\n");
  fprintf(psfile,"    pattern\n");
  fprintf(psfile,"      { patternlength add /patternlength exch def\n");
  fprintf(psfile,"      } forall\n");
  fprintf(psfile,"    pattern length 2 mod 0 ne\n");
  fprintf(psfile,"      { /patternlength patternlength 2 mul def } if\n");
  fprintf(psfile,"    /first pattern 0 get def\n");
  fprintf(psfile,"    /last patternlength first sub def\n");
  fprintf(psfile,"    /n pathlen last sub cvi patternlength idiv def\n");
  fprintf(psfile,"    /endpart pathlen patternlength n mul sub\n");
  fprintf(psfile,"       last sub 2 div def\n");
  fprintf(psfile,"    /offset first endpart sub def\n");
  fprintf(psfile,"    pattern offset setdash\n");
  fprintf(psfile,"  } def\n");
  fprintf(psfile,"\n");
  fprintf(psfile,"/pathlength\n");
  fprintf(psfile,"    { flattenpath\n");
  fprintf(psfile,"      /dist 0 def\n");
  fprintf(psfile,"      \n");
  fprintf(psfile,"      { /yfirst exch def /xfirst exch def\n");
  fprintf(psfile,"        /ymoveto yfirst def /xmoveto xfirst def }\n");
  fprintf(psfile,"      { /ynext exch def /xnext exch def\n");
  fprintf(psfile,"        /dist dist ynext yfirst sub dup mul\n");
  fprintf(psfile,"          xnext xfirst sub dup mul add sqrt add def\n");
  fprintf(psfile,"        /yfirst ynext def /xfirst xnext def }\n");
  fprintf(psfile,"      {}\n");
  fprintf(psfile,"      \n");
  fprintf(psfile,"      { /ynext ymoveto def /xnext xmoveto def\n");
  fprintf(psfile,"        /dist dist ynext yfirst sub dup mul\n");
  fprintf(psfile,"          xnext xfirst sub dup mul add sqrt add def\n");
  fprintf(psfile,"        /yfirst ynext def /xfirst xnext def }\n");
  fprintf(psfile,"      pathforall\n");
  fprintf(psfile,"      dist\n");
  fprintf(psfile,"    } def\n");

  fprintf(psfile,"%% This is the stuff for filling paths with lined patterns.\n");
  fprintf(psfile,"%%   It too borrows liberally from the code provided in Adobe's\n");
  fprintf(psfile,"%%   Blue Book.\n");
  fprintf(psfile,"/setuserscreendict 22 dict def\n");
  fprintf(psfile,"setuserscreendict begin\n");
  fprintf(psfile,"  /tempctm matrix def\n");
  fprintf(psfile,"  /temprot matrix def\n");
  fprintf(psfile,"  /tempscale matrix def\n");
  fprintf(psfile,"    \n");
  fprintf(psfile,"  /concatprocs\n");
  fprintf(psfile,"    { /proc2 exch cvlit def\n");
  fprintf(psfile,"      /proc1 exch cvlit def\n");
  fprintf(psfile,"      /newproc proc1 length proc2 length add\n");
  fprintf(psfile,"        array def\n");
  fprintf(psfile,"      newproc 0 proc1 putinterval\n");
  fprintf(psfile,"      newproc proc1 length proc2 putinterval\n");
  fprintf(psfile,"      newproc cvx\n");
  fprintf(psfile,"    } def\n");
  fprintf(psfile,"    \n");
  fprintf(psfile,"  /resmatrix matrix def\n");
  fprintf(psfile,"  /findresolution\n");
  fprintf(psfile,"    { 72 0 resmatrix defaultmatrix dtransform\n");
  fprintf(psfile,"      /yres exch def /xres exch def\n");
  fprintf(psfile,"      xres dup mul yres dup mul add sqrt\n");
  fprintf(psfile,"    } def\n");
  fprintf(psfile,"  end\n");
  fprintf(psfile,"\n");
  fprintf(psfile,"/setuserscreen\n");
  fprintf(psfile,"  { setuserscreendict begin\n");
  fprintf(psfile,"      /spotfunction exch def\n");
  fprintf(psfile,"      /screenangle exch def\n");
  fprintf(psfile,"      /cellsize exch def\n");
  fprintf(psfile,"      \n");
  fprintf(psfile,"      /m tempctm currentmatrix def\n");
  fprintf(psfile,"      /rm screenangle temprot rotate def\n");
  fprintf(psfile,"      /sm cellsize dup tempscale scale def\n");
  fprintf(psfile,"      \n");
  fprintf(psfile,"      sm rm m m concatmatrix m concatmatrix pop\n");
  fprintf(psfile,"      \n");
  fprintf(psfile,"      1 0 m dtransform /y1 exch def /x1 exch def\n");
  fprintf(psfile,"      \n");
  fprintf(psfile,"      /veclength x1 dup mul y1 dup mul add sqrt def\n");
  fprintf(psfile,"      /frequency findresolution veclength div def\n");
  fprintf(psfile,"      \n");
  fprintf(psfile,"      /newscreenangle y1 x1 atan def\n");
  fprintf(psfile,"      \n");
  fprintf(psfile,"      m 2 get m 1 get mul m 0 get m 3 get mul sub\n");
  fprintf(psfile,"        0 gt\n");
  fprintf(psfile,"        { {neg} /spotfunction load concatprocs\n");
  fprintf(psfile,"            /spotfunction exch def\n");
  fprintf(psfile,"        } if\n");
  fprintf(psfile,"        \n");
  fprintf(psfile,"      frequency newscreenangle /spotfunction load\n");
  fprintf(psfile,"        setscreen\n");
  fprintf(psfile,"    end\n");
  fprintf(psfile,"  } def\n");
  fprintf(psfile,"  \n");
  fprintf(psfile,"/setpatterndict 18 dict def\n");
  fprintf(psfile,"setpatterndict begin\n");
  fprintf(psfile,"  /bitison\n");
  fprintf(psfile,"    { /ybit exch def /xbit exch def\n");
  fprintf(psfile,"      /bytevalue bstring ybit bwidth mul xbit 8 idiv\n");
  fprintf(psfile,"        add get def\n");
  fprintf(psfile,"      \n");
  fprintf(psfile,"      /mask 1 7 xbit 8 mod sub bitshift def\n");
  fprintf(psfile,"      bytevalue mask and 0 ne\n");
  fprintf(psfile,"    } def\n");
  fprintf(psfile,"  end\n");
  fprintf(psfile,"  \n");
  fprintf(psfile,"/bitpatternspotfunction\n");
  fprintf(psfile,"  { setpatterndict begin\n");
  fprintf(psfile,"    /y exch def /x exch def\n");
  fprintf(psfile,"    \n");
  fprintf(psfile,"    /xindex x 1 add 2 div bpside mul cvi def\n");
  fprintf(psfile,"    /yindex y 1 add 2 div bpside mul cvi def\n");
  fprintf(psfile,"    \n");
  fprintf(psfile,"    xindex yindex bitison\n");
  fprintf(psfile,"      { /onbits onbits 1 add def 1 }\n");
  fprintf(psfile,"      { /offbits offbits 1 add def 0 }\n");
  fprintf(psfile,"      ifelse\n");
  fprintf(psfile,"    end\n");
  fprintf(psfile,"  } def\n");
  fprintf(psfile,"\n");
  fprintf(psfile,"/setpattern\n");
  fprintf(psfile,"  { setpatterndict begin\n");
  fprintf(psfile,"    /cellsz exch def\n");
  fprintf(psfile,"    /angle exch def\n");
  fprintf(psfile,"    /bwidth exch def\n");
  fprintf(psfile,"    /bpside exch def\n");
  fprintf(psfile,"    /bstring exch def\n");
  fprintf(psfile,"    \n");
  fprintf(psfile,"    /onbits 0 def /offbits 0 def\n");
  fprintf(psfile,"    cellsz angle /bitpatternspotfunction load\n");
  fprintf(psfile,"      setuserscreen\n");
  fprintf(psfile,"    { } settransfer\n");
  fprintf(psfile,"    offbits offbits onbits add div setgray\n");
  fprintf(psfile,"  end\n");
  fprintf(psfile,"} def\n");
  fprintf(psfile,"\n");
  fprintf(psfile,"%% these are the patterns used for lining\n");
  fprintf(psfile,"/pat3 <8888888888888888> def\n");
  fprintf(psfile,"/pat4 <aaaaaaaaaaaaaaaa> def\n");

  switch(PS_options.atom_sphere_type){
  case ATOM_PLAIN_FILL:
    fprintf(psfile,"/arcwid 3.0 def\n");
    break;
  default:
    fprintf(psfile,"/arcwid 1.0 def\n");
  }

  fprintf(psfile,"%% Re-encode the fonts to use ISOLatin1 so we can get extended characters\n");
  fprintf(psfile,"/%s findfont\n",PS_options.fontname);
  fprintf(psfile,"dup length dict begin\n");
  fprintf(psfile,"{1 index /FID ne {def} {pop pop} ifelse} forall\n");
  fprintf(psfile,"        /Encoding ISOLatin1Encoding def\n");



  fprintf(psfile,"          currentdict\n");
  fprintf(psfile,"end\n");
  fprintf(psfile,"/%s exch definefont pop\n",PS_options.fontname);


  fprintf(psfile,"/textsize %lf def\n",PS_options.fontsize);
  fprintf(psfile,
          "/normalfont {/%s findfont textsize scalefont setfont} def\n",
          PS_options.fontname);
  fprintf(psfile,
          "/symbolfont {/Symbol findfont textsize scalefont setfont} def\n");

  fprintf(psfile,"normalfont\n");
  fprintf(psfile,"/thesize %lf def\n",
          8.5*72.0*PS_options.printscale/(float)g_xmax);
  fprintf(psfile,"/scaleit {thesize -1 thesize mul scale} def\n");
  fprintf(psfile,"/iscaleit {1 thesize div -1 thesize div scale} def\n");
  fprintf(psfile,"/rotit {currentpoint gsave translate rotate 0 0 moveto} def\n");
  fprintf(psfile,"/unrot {grestore} def\n");
  fprintf(psfile,"/crossedcirc {newpath 2 index 2 index translate\n");
  fprintf(psfile,"              0 0 2 index 0 360 arc stroke\n");
  fprintf(psfile,"              1 -0.5 scale 0 0 2 index 0 180 arcn stroke\n");
  fprintf(psfile,"              0.5 2.0 scale 0 0 2 index 90 270 arc stroke\n");
  fprintf(psfile,"              2.0 -1.0 scale pop 1 index -1 mul 1 index -1 mul translate pop pop} def\n");

  /* some definitions which will be used for displaying text */
  fprintf(psfile,"/SNP {stroke newpath} def\n");
  fprintf(psfile,"/vshift -1 textsize mul def\n");
  fprintf(psfile,"/M {moveto} def\n");
  fprintf(psfile,"/L {lineto} def\n");
  fprintf(psfile,"/R {rmoveto} def\n");
  fprintf(psfile,"/SG {setgray} def\n");

  fprintf(psfile,"/MFshow {{dup dup 0 get findfont exch 1 get scalefont setfont\n");
  fprintf(psfile,"     [ currentpoint ] exch dup 2 get 0 exch rmoveto dup 4 get show dup\n");
  fprintf(psfile,"     3 get {2 get neg 0 exch rmoveto pop} {pop aload pop moveto}ifelse} forall} bind def\n");
  fprintf(psfile,"/MFwidth {0 exch {dup 3 get{dup dup 0 get findfont exch 1 get scalefont setfont\n");
  fprintf(psfile,"      4 get stringwidth pop add}\n");
  fprintf(psfile,"    {pop} ifelse} forall} bind def\n");
  fprintf(psfile,"/Lshow { currentpoint stroke M\n");
  fprintf(psfile,"  0 exch iscaleit R MFshow stroke scaleit} bind def\n");
  fprintf(psfile,"/Rshow { currentpoint stroke M\n");
  fprintf(psfile,"  exch dup iscaleit MFwidth neg 3 -1 roll R MFshow stroke scaleit } def\n");
  fprintf(psfile,"/Cshow { currentpoint stroke M\n");
  fprintf(psfile,"  exch dup iscaleit MFwidth -2 div 3 -1 roll R MFshow stroke scaleit} def\n");

  fprintf(psfile,"%%\n%% These next three are for triangles\n%%\n");
  fprintf(psfile,"%%    here you can change the shading of triangles\n");
  fprintf(psfile,"%%    increase 1 to decrease the range of colors and\n");
  fprintf(psfile,"%%    increase 0 to set the minimum darknes.\n");
  fprintf(psfile,"%%    remember: 0 is black and 1 is white.\n");
  fprintf(psfile,"%%    for example:\n");
  fprintf(psfile,"%%    /SNG { 2 div 0.25 add setgray } def\n");
  fprintf(psfile,"%%    gives a smaller range that tend to be dark, but\n");
  fprintf(psfile,"%%    doesn't make it all the way to black.\n");
  fprintf(psfile,"%%\n");
  fprintf(psfile,"%%    If you really want, you can even get color.\n");
  fprintf(psfile,"%%    uncomment this bit, then comment out SNG def below\n");
  fprintf(psfile,"%%    /SNG { dup\n");
  fprintf(psfile,"%%    1.5 div 0.5 add %% the red component\n");
  fprintf(psfile,"%%    dup 2 index \n");
  fprintf(psfile,"%%    2 div 0.5 add %% the green component\n");
  fprintf(psfile,"%%    2 index \n");
  fprintf(psfile,"%%    1.5 div 0.5 add %% the blue component\n");
  fprintf(psfile,"%%    setrgbcolor } def\n");
  fprintf(psfile,"%%    Gives a lovely faintly reddish aqua.  If you\n");
  fprintf(psfile,"%%    play with the divisors and offsets (0.5 in each\n");
  fprintf(psfile,"%%    of the above examples), you should get the idea\n");
  fprintf(psfile,"%%    pretty quickly.\n");
  fprintf(psfile,"/SNG { 1 div 0 add setgray } def\n");
  fprintf(psfile,"/fT { M L L L fill stroke} def\n");
  fprintf(psfile,"%% if you don't like the outlines, substitute the\n");
  fprintf(psfile,"%%    second definition by commenting out this line\n");
  fprintf(psfile,"%%    and uncommenting the next.\n");
  fprintf(psfile,"/dT { 1 setlinejoin M L L L stroke} def\n");
  fprintf(psfile,"%%/dT { pop pop pop pop pop pop pop pop} def\n");

  i=0;
  while(myPSdefns[i]!=0){
    fprintf(psfile,"%s",myPSdefns[i]);
    i++;
  }

  fprintf(psfile,"\n%%%%Page: 1 1\n\n");
  fprintf(psfile,"scaleit\n");

  fprintf(psfile,"0.0 %lf translate\n",translation);


  /*****
    prompt the user for some data, if required
  ******/
  switch(PS_options.bond_type){
  case BOND_SHADE:
    if( PS_options.bond_shade == 0 ){
      PS_options.bond_shade = 0.5;
    }
    readfloatparm("bond shading (between 0 and 1)",&PS_options.bond_shade);
    if( PS_options.bond_shade < 0 )
      PS_options.bond_shade = 0;
    if( PS_options.bond_shade > 1 )
      PS_options.bond_shade = 1;
    break;
  case BOND_CSHADE:
    if(colorstring[0] == 0 ){
      strcpy(colorstring,"0.0 0.0 0.5");
    }
    stringarr[0] = colorstring;
    readstringparm("bond color (RGB triple)",stringarr);
    sscanf(colorstring,"%lf %lf %lf",&PS_options.bond_color[0],
           &PS_options.bond_color[1],&PS_options.bond_color[2]);
    if( PS_options.bond_color[0] < 0 ) PS_options.bond_color[0] = 0;
    if( PS_options.bond_color[0] > 1 ) PS_options.bond_color[0] = 1;
    if( PS_options.bond_color[1] < 0 ) PS_options.bond_color[1] = 0;
    if( PS_options.bond_color[1] > 1 ) PS_options.bond_color[1] = 1;
    if( PS_options.bond_color[2] < 0 ) PS_options.bond_color[2] = 0;
    if( PS_options.bond_color[2] > 1 ) PS_options.bond_color[2] = 1;
    break;
  }


  /* set a toggle to let the redraw routines know to dump to ps */
  doing_ps = 1;

  redrawgraph();

  doing_ps = 0;

  fprintf(psfile,"showpage\n");
  ENHPS_reset();


  fclose(psfile);
}

