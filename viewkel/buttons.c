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

   22.05.98 gL:
     Modifications to curve_buttons to deal with individual
     shading toggles.
   16.06.98 gL:
     fixed the problem with toggle_buttons caused by the above addition.
   26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)
  18.10.98 gL:
     support for fatbands added.
  26.01.99 gL:
     support for global reads.
  29.01.99 gL:
     support for box toggling on molecules.

***/

#include "viewkel.h"

#ifndef _BUTTONS_
#include "buttons.h"
#endif

#ifdef X_GRAPHICS
/* these are the maps for all of the buttons */
#include "floatparm.xbm"
#include "function.xbm"
#include "toggleon.xbm"
#include "toggleoff.xbm"
#include "toggleonfill.xbm"
#include "modal.xbm"
#endif


#ifdef MAC_GRAPHICS
MenuHandle activeMenu;
short nextItemID;
enum{
function_width,function_height,
modal_width,modal_height,
toggleon_width,toggleon_height,
floatparm_width,floatparm_height} dumb_enum;
#endif


/**********************

  this is all of the code for displaying buttons and building the button
   list

   Buttons are handled differently depending upon which graphics
    system we are using.

   In X the buttons are actually buttons in separate windows
     (one window per object open)
   On the Mac, the buttons are entries in menus that are
     stuck in the menu bar.

***********************/


button_type *buttonlist;

/************

   the routines to draw buttons are not needed on the Mac, since
    the menubar updating stuff takes care of it for us.

*************/

#ifdef X_GRAPHICS
/****************************************************************************
*
*  Procedure g_floatbutton
*
*   Arguments:
*                    win: Window
*                 button: pointer to button_type
*              xpos,ypos: integers
*   Return Value: none
*
*  Actions: this draws the button at xpos, ypos in the window 'win.
*
*****************************************************************************/
void g_floatbutton(Window win,button_type *button,int xpos,int ypos)
{
  char string[80], floatstr[30];
  int middle;

  /* convert the float value into a string */
  sprintf(floatstr,"%.5lf",*(button->guts.floatparm));

  /* now form the string that is going to be put into the button */
  strcpy(string,button->string);
  strcat(string,": ");
  strcat(string,floatstr);

  /* set the gc correctly */
  XSetFunction(disp,smalltextgc,GXcopy);

  middle = xpos + button->xdim/2-
    XTextWidth(small_font,string,strlen(string))/2;
  XFillRectangle(disp,win,blackgc,xpos,ypos,button->xdim,
                 button->ydim);

  XCopyArea(disp,button->pix1,win,graphgc,0,0,
            floatparm_width,floatparm_height,xpos,ypos);
  XDrawString(disp,win,smalltextgc,middle,ypos+TEXTOFF,
              string,strlen(string));
}

/****************************************************************************
*
*  Procedure g_intbutton
*
*   Arguments:
*                    win: Window
*                 button: pointer to button_type
*              xpos,ypos: integers
*   Return Value: none
*
*  Actions: this draws the button at xpos, ypos in the window 'win.
*
*****************************************************************************/
void g_intbutton(Window win,button_type *button,int xpos,int ypos)
{
  char string[80], floatstr[30];
  int middle;

  /* convert the int value into a string */
  sprintf(floatstr,"%d",*(button->guts.intparm));

  /* now form the string that is going to be put into the button */
  strcpy(string,button->string);
  strcat(string,": ");
  strcat(string,floatstr);

  /* set the gc correctly */
  XSetFunction(disp,smalltextgc,GXcopy);

  middle = xpos + button->xdim/2-
    XTextWidth(small_font,string,strlen(string))/2;
  XFillRectangle(disp,win,blackgc,xpos,ypos,button->xdim,
                 button->ydim);
  XCopyArea(disp,button->pix1,win,graphgc,0,0,
            floatparm_width,floatparm_height,xpos,ypos);
  XDrawString(disp,win,smalltextgc,middle,ypos+TEXTOFF,
              string,strlen(string));
}


/****************************************************************************
*
*  Procedure g_stringbutton
*
*   Arguments:
*                    win: Window
*                 button: pointer to button_type
*              xpos,ypos: integers
*   Return Value: none
*
*  Actions: this draws the button at xpos, ypos in the window 'win.
*
*****************************************************************************/
void g_stringbutton(Window win,button_type *button,int xpos,int ypos)
{
  char string[80];
  int middle;

  /* now form the string that is going to be put into the button */
  strcpy(string,button->string);
  strcat(string,": ");
  strcat(string,button->guts.stringparm.lines[0]);

  /* set the gc correctly */
  XSetFunction(disp,smalltextgc,GXcopy);

  middle = xpos + button->xdim/2-
    XTextWidth(small_font,string,strlen(string))/2;
  XFillRectangle(disp,win,blackgc,xpos,ypos,button->xdim,
                 button->ydim);

  XCopyArea(disp,button->pix1,win,graphgc,0,0,
            floatparm_width,floatparm_height,xpos,ypos);
  XDrawString(disp,win,smalltextgc,middle,ypos+TEXTOFF,
              string,strlen(string));
}


/****************************************************************************
*
*  Procedure g_multistringbutton
*
*   Arguments:
*                    win: Window
*                 button: pointer to button_type
*              xpos,ypos: integers
*   Return Value: none
*
*  Actions: this draws the button at xpos, ypos in the window 'win.
*
*****************************************************************************/
void g_multistringbutton(Window win,button_type *button,int xpos,int ypos)
{
  char string[80];
  int middle;

  /* now form the string that is going to be put into the button */
  strcpy(string,button->string);

  /* set the gc correctly */
  XSetFunction(disp,smalltextgc,GXcopy);

  middle = xpos + button->xdim/2-
    XTextWidth(small_font,string,strlen(string))/2;
  XFillRectangle(disp,win,blackgc,xpos,ypos,button->xdim,
                 button->ydim);

  XCopyArea(disp,button->pix1,win,graphgc,0,0,
            floatparm_width,floatparm_height,xpos,ypos);
  XDrawString(disp,win,smalltextgc,middle,ypos+TEXTOFF,
              string,strlen(string));
}



/****************************************************************************
*
*  Procedure g_functionbutton
*
*   Arguments:
*                    win: Window
*                 button: pointer to button_type
*              xpos,ypos: integers
*   Return Value: none
*
*  Actions: this draws the button at xpos, ypos
*
*****************************************************************************/
void g_functionbutton(Window win,button_type *button,int xpos,int ypos)
{
  int middle;

  /* set the gc correctly */
  XSetFunction(disp,smalltextgc,GXcopy);

  middle = xpos + button->xdim/2-
    XTextWidth(small_font,button->string,strlen(button->string))/2;
  XFillRectangle(disp,win,blackgc,xpos,ypos,button->xdim,
                 button->ydim);
  XCopyArea(disp,button->pix1,win,graphgc,0,0,
            function_width,function_height,xpos,ypos);
  XDrawString(disp,win,smalltextgc,middle,ypos+TEXTOFF,
              button->string,strlen(button->string));
}

/****************************************************************************
*
*  Procedure g_togglebutton
*
*   Arguments:
*                    win: Window
*                 button: pointer to button_type
*              xpos,ypos: integers
*   Return Value: none
*
*  Actions: draws a toggle button.  The button is white if the toggled variable
*       is true, black if it is false.
*
*****************************************************************************/
void g_togglebutton(Window win,button_type *button,int xpos,int ypos)
{
  Pixmap *pix;
  int middle;

  if( *(button->guts.boolean) ){
    pix = &button->pix1;
  }
  else{
    pix = &button->pix2;
  }
  middle = xpos +  (button->xdim-BUTTONOFF)/2-
    XTextWidth(small_font,button->string,strlen(button->string))/2;

  XCopyArea(disp,*pix,win,graphgc,0,0,toggleon_width,toggleon_height,
            xpos,ypos);
  XDrawString(disp,win,smalltextgc,middle,ypos+TEXTOFF,
              button->string,strlen(button->string));

}


/****************************************************************************
*
*  Procedure g_curvebutton
*
*   Arguments:
*                    win: Window
*                 button: pointer to button_type
*              xpos,ypos: integers
*   Return Value: none
*
*  Actions: draws a curve button.  This is a toggle button
*     that draws a sample of the line style adjacent to the button.
*
*****************************************************************************/
void g_curvebutton(Window win,button_type *button,int xpos,int ypos)
{
  Pixmap *pix;
  int middle;

  if( *button->guts.curve.boolean) {
    if( button->guts.curve.fill && *(button->guts.curve.fill) ){
      pix = &button->pix3;
    } else{
      pix = &button->pix1;
    }
  }
  else{
      pix = &button->pix2;
  }
  middle = xpos +  (button->xdim-BUTTONOFF)/2-
    XTextWidth(small_font,button->string,strlen(button->string))/2;

  XCopyArea(disp,*pix,win,graphgc,0,0,toggleon_width,toggleon_height,
            xpos,ypos);
  XDrawString(disp,win,smalltextgc,middle,ypos+TEXTOFF,
              button->string,strlen(button->string));

  /* now draw in the line sample */
  XFillRectangle(disp,win,blackgc,xpos+button->xdim+5,
                 ypos + button->ydim/2-2,xpos+button->xdim + CURVE_SEG_LEN,
                 ypos + button->ydim/2+2);
  g_change_linestyle(*(button->guts.curve.style));
  XDrawLine(disp,win,graphgc,xpos+button->xdim + 5,
            ypos + button->ydim/2,xpos+button->xdim + CURVE_SEG_LEN,
            ypos + button->ydim/2);
  g_change_linestyle(0);
}

/****************************************************************************
*
*  Procedure g_modebutton
*
*   Arguments:
*                    win: Window
*                 button: pointer to button_type
*              xpos,ypos: integers
*   Return Value: none
*
*  Actions: draws a mode choice button.  The current mode is displayed.
*
*****************************************************************************/
void g_modebutton(Window win,button_type *button,int xpos,int ypos)
{
  Pixmap *pix;
  int middle;
  char *which_string;


  pix = &button->pix1;

  /* figure out which string to choose */
  which_string =
    button->guts.mode.names[*(button->guts.mode.controlling_variable)];

  middle = xpos +  (button->xdim-BUTTONOFF)/2-
    XTextWidth(small_font,which_string,strlen(which_string))/2;

  XCopyArea(disp,*pix,win,graphgc,0,0,toggleon_width,toggleon_height,
            xpos,ypos);
  XDrawString(disp,win,smalltextgc,middle,ypos+TEXTOFF,
              which_string,strlen(which_string));

}

/****************************************************************************
*
*  Procedure g_drawbuttons
*
*   Arguments: button_win: pointer to button_win_type
*
*
*   Return Value: none
*
*  Actions: this draws in the buttons which are in the button window 'button_win
*
*
*****************************************************************************/
void g_drawbuttons(button_win_type *button_win)
{
  int xpos, ypos, offset;
  Window which_win;
  button_type *button;

  offset = BUTTONOFF;
  xpos = BUTTONOFF;
  ypos = offset;

  if( !button_win ) FATAL_BUG("Null list passed to g_drawbuttons.");

  which_win = button_win->which_win;
  button = button_win->button_list;
  while(button){
    switch(button->type){
    case TOGGLE:
      g_togglebutton(which_win,button,xpos,ypos);
      break;
    case FLOATPARM:
      g_floatbutton(which_win,button,xpos,ypos);
      break;
    case INTPARM:
      g_intbutton(which_win,button,xpos,ypos);
      break;
    case STRINGPARM:
      g_stringbutton(which_win,button,xpos,ypos);
      break;
    case MULTISTRINGPARM:
      g_multistringbutton(which_win,button,xpos,ypos);
      break;
    case FUNCTION:
      g_functionbutton(which_win,button,xpos,ypos);
      break;
    case MODE:
      g_modebutton(which_win,button,xpos,ypos);
      break;
    case CURVETOGGLE:
      g_curvebutton(which_win,button,xpos,ypos);
      break;
    }
    ypos += button->ydim + offset;
    button = button->next;
  }

  XFlush(disp);
}


/****************************************************************************
*
*  Procedure g_draw_all_buttons
*
*   Arguments: button_win: pointer to button_win_type
*
*
*   Return Value: none
*
*  Actions: this draws in all the buttons in the linked list 'button_win
*
*
*****************************************************************************/
void g_draw_all_buttons(button_win_type *button_win)
{
  button_win_type *which_button_win;

  if( !button_win ) return;

  which_button_win = button_win;

  while(which_button_win){
    g_drawbuttons(which_button_win);
    if( which_button_win->child) g_drawbuttons(which_button_win->child);
    which_button_win = which_button_win->next;
  }
}

#endif

/****************************************************************************
 *
 *  Function new_button
 *
 *   Arguments: none
 *   Return Value: button_type *
 *
 *  Actions: This gets space for a new button and makes sure that the
 *    space was allocated correctly
 *
 *****************************************************************************/
button_type *newbutton(void)
{
  button_type *ptr;

  ptr = (button_type *)D_CALLOC(1,sizeof(button_type));
  if(!ptr )fatal("Memory Allocation getting a button");

#ifdef MAC_GRAPHICS
        /****************
                insert the item into the menu bar with a bogus string,
                we'll fix the string later
  *****************/
  AppendMenu(activeMenu,"\pFoo item");
  ptr->itemID = nextItemID;
  nextItemID++;
#endif

  return(ptr);
}

/****************************************************************************
 *
 *  Procedure free_but_win
 *
 *   Arguments: but_win: pointer to button_win_type
 *             but_list: pointer to pointer to button_win_type
 *
 *   Return Value: none
 *
 *  Actions: Frees up the button window pointed to by 'but_win and removes
 *   it from the list of windows pointed to by 'but_list
 *
 *****************************************************************************/
void free_but_win(button_win_type *but_win,button_win_type **but_list)
{
  button_win_type *temp;
  button_type *button,*last;
  char found;

  /* first find the window in the list */
  temp = *but_list;

  /* check to see if it's at the top of the list */
  found = 0;
  if( temp && but_list && temp == but_win ){
    *but_list = (*but_list)->next;
    found = 1;
  } else {
    while( temp && !found){
      if( temp->next == but_win ){
        temp->next = temp->next->next;
        found = 1;
      }
      else temp = temp->next;
    }
  }

  /* make sure that we found it in the list */
  if( !found ){
    if(but_list) FATAL_BUG("Can't find a but_win in the but_win list.");
    else return;
  }

  /* free the list of buttons in the button window */
  button = but_win->button_list;
  while(button){
#ifdef X_GRAPHICS
    /* free up the pixmaps */
    if( button->pix1 ) XFreePixmap(disp,button->pix1);
    if( button->pix2 ) XFreePixmap(disp,button->pix2);
#endif
    last = button;
    button = button->next;
    D_FREE(last);
  }

#ifdef X_GRAPHICS
  /* now close the window */
  XDestroyWindow(disp,but_win->which_win);
#endif

#ifdef MAC_GRAPHICS
  /* remove it from the menu bar and free the handle */
        DeleteMenu(but_win->menuID);
        DisposeMenu(but_win->theMenu);
        /* redraw the menu bar so that menu goes away */
        DrawMenuBar();
#endif

  /* free the structure */
  D_FREE(but_win);
}

/****************************************************************************
 *
 *  Procedure free_child_but_win
 *
 *   Arguments: but_win: pointer to button_win_type
 *
 *
 *   Return Value: none
 *
 *  Actions: Free's up the child button window pointed to by 'but_win
 *
 *****************************************************************************/
void free_child_but_win(button_win_type *but_win)
{
  button_type *button,*last;

  /* free the list of buttons in the button window */
  button = but_win->button_list;
  while(button){
#ifdef X_GRAPHICS
    /* free up the pixmaps */
    if( button->pix1 ) XFreePixmap(disp,button->pix1);
    if( button->pix2 ) XFreePixmap(disp,button->pix2);
#endif
    last = button;
    button = button->next;
    D_FREE(last);
  }

#ifdef X_GRAPHICS
  /* now close the window */
  XDestroyWindow(disp,but_win->which_win);
#endif

#ifdef MAC_GRAPHICS
  /* remove it from the menu bar and free the handle */
        DeleteMenu(but_win->menuID);
        DisposeMenu(but_win->theMenu);
        /* redraw the menu bar so that menu goes away */
        DrawMenuBar();
#endif

  /* free the structure */
  D_FREE(but_win);
}



/****************************************************************************
 *
 *                   Procedure handlebutton
 *
 * Arguments:  button_win: pointer to button_win_type
 *             xpos,ypos: integers
 *          mouse_button: integer
 *             key_state: integer
 * Returns: none
 *
 * Action: Figures out which button the user clicked on (in the button box) and
 *   performs the appropriate action;
 *   'mouse_button indicates which mouse button was clicked.  This allows the
 *    routines to deal with multiple types of mouse clicks.
 *
 ****************************************************************************/
void handlebutton(button_win_type *button_win,int xpos,int ypos,int mouse_button,
                  int key_state)
{
  int xcor, ycor, xdim, ydim, offset;
  button_type *button;
  int *mode_ptr;

  /*******

    the Mac equivalent of this is the function ViewkelObjectMenuDispatch
    in the file Mac_Menus.c

    ********/

#ifdef X_GRAPHICS
  offset = BUTTONOFF;
  xcor = offset;
  ycor = offset;

  button = button_win->button_list;
  while(button){
    xdim = button->xdim;
    ydim = button->ydim;
    if(xpos>xcor && xpos<xcor+xdim && ypos>ycor && ypos<ycor+ydim){
      switch(mouse_button){
      case 1:
        switch(button->type){
        case FLOATPARM:
          readfloatparm(button->string,button->guts.floatparm);
          break;
        case INTPARM:
          readintparm(button->string,button->guts.intparm);
          break;
        case STRINGPARM:
          readstringparm(button->string,button->guts.stringparm.lines);
          break;
        case MULTISTRINGPARM:
          readmultistringparm(button->string,button->guts.stringparm.num_lines,
                              button->guts.stringparm.lines);
          break;
        case TOGGLE:
          *(button->guts.boolean) = !*(button->guts.boolean);
          break;
        case CURVETOGGLE:
          *(button->guts.curve.boolean) = !*(button->guts.boolean);
          break;
        case FUNCTION:
          (button->guts.function.func)(button->guts.function.num_args,
                                       button->guts.function.args);
          break;
        case MODE:
          mode_ptr = button->guts.mode.controlling_variable;
          *mode_ptr = (*mode_ptr+1) % button->guts.mode.num_modes;
          break;
        }
      case 2:
        break;
      case 3:
        switch(button->type){
        case CURVETOGGLE:
          if( key_state & ShiftMask && button->guts.curve.fill ){
            *(button->guts.curve.fill) = !(*(button->guts.curve.fill));
          } else{
            *(button->guts.curve.style) = *button->guts.curve.style + 1;
          }
          break;
        }
        break;
      }
    }
    ycor += ydim + offset;
    button = button->next;
  }

  g_drawbuttons(button_win);
#endif
}



/****************************************************************************
 *
 *                   Procedure find_button_win
 *
 * Arguments:  button_win: pointer to button_win_type
 *                    win: Window
 *              xpos,ypos: integers
 *           mouse_button:  integer
 *              key_state: integer
 *
 * Returns: none
 *
 * Action: Determines which button_window was clicked in, then calls the
 *   handlebutton function to process the click.
 *
 ****************************************************************************/
void find_button_win(button_win_type *button_win,Window win,
                     int xpos,int ypos,int mouse_button,int key_state)
{
  button_win_type *which;
  int found = 0;

  /*******

    the Mac equivalent of this is the function ViewkelObjectMenuDispatch
    in the file Mac_Menus.c

    ********/
#ifdef X_GRAPHICS
  which = button_win;

  while( which && !found ){
    if( which->which_win == win ){
      found = 1;
      /* process the click */
      handlebutton(which,xpos,ypos,mouse_button,key_state);
    }
    else if( which->child && which->child->which_win == win ){
      found = 1;
      /* process the click */
      handlebutton(which->child,xpos,ypos,mouse_button,key_state);
    }
    else{
      which = which->next;
    }
  }
#endif
}


/****************************************************************************
 *
 *  Function new_button_win
 *
 *   Arguments: button_win:  pointer to button_win_type
 *   Return Value: a pointer to button_win_type
 *
 *  Actions: Opens a new button window and gets all the space needed for it.
 *
 *****************************************************************************/
button_win_type *new_button_win(button_win_type *button_win,char *string)
{
#ifdef MAC_GRAPHICS
        static short nextMenuID = 1;
        char pstring[255];
#endif
  button_win_type *new_win;

  /* get space */
  new_win = (button_win_type *)D_CALLOC(1,sizeof(button_type));
  if(!new_win) fatal("Can't get space for a new button window.");

#ifdef X_GRAPHICS
  /* open the window */
  new_win->which_win = XCreateSimpleWindow(disp,root,0,0,BUTWIDTH,BUTHEIGHT,1,
                                           fcolor,bcolor);
  new_win->xdim = BUTWIDTH;
  new_win->ydim = BUTHEIGHT;

  /* ask for the events that we need */
  XSelectInput(disp,new_win->which_win,
               KeyPressMask|ButtonPressMask|ButtonMotionMask|ExposureMask);

  /* change the cursor in the button window */
  XDefineCursor(disp,new_win->which_win,XCreateFontCursor(disp,60));

        /* set the title */
  XSetStandardProperties( disp,new_win->which_win,string,string,None,
                         0,0,0);


  /* expose the window */
  XMapWindow(disp,new_win->which_win);
  XFlush(disp);
#endif

#ifdef MAC_GRAPHICS
        /* set the title (goofy pascal strings...) */
        strcpy(pstring,string);
        CtoPstr(pstring);
        /* allocate the handle and set the ID */
        new_win->theMenu = NewMenu(nextMenuID,(unsigned char *)pstring);
        new_win->menuID = nextMenuID;
        nextMenuID++;
        /* insert the menu into the menubar */
        InsertMenu(new_win->theMenu,0);

        nextItemID = 1;
        activeMenu = new_win->theMenu;
#endif
  /* put the window at the head of the list */
  new_win->next = button_win;

  return(new_win);

  /* that's it.... */
}

/****************************************************************************
 *
 *  Procedure set_button_pixmaps
 *
 *   Arguments: button_win: pointer to button_win_type
 *   Return Value: none
 *
 *  Actions: This procedure assigns the appropriate pixmaps to each
 *    button window
 *
 *****************************************************************************/
void set_button_pixmaps(button_win_type *button_win)
{
  button_type *b;

#ifdef X_GRAPHICS
  /* now loop through the button list and add in the pixmaps */
  b = button_win->button_list;
  while( b ){
    switch(b->type){
    case TOGGLE:
      b->pix1=XCreatePixmapFromBitmapData(disp,button_win->which_win,
                                          toggleon_bits,
                                          toggleon_width,toggleon_height,fcolor,
                                          bcolor,screen_depth);
      b->pix2=XCreatePixmapFromBitmapData(disp,button_win->which_win,
                                          toggleoff_bits,
                                          toggleoff_width,toggleoff_height,
                                          fcolor,bcolor,screen_depth);
      break;
    case CURVETOGGLE:
      b->pix1=XCreatePixmapFromBitmapData(disp,button_win->which_win,
                                          toggleon_bits,
                                          toggleon_width,toggleon_height,fcolor,
                                          bcolor,screen_depth);
      b->pix2=XCreatePixmapFromBitmapData(disp,button_win->which_win,
                                          toggleoff_bits,
                                          toggleoff_width,toggleoff_height,
                                          fcolor,bcolor,screen_depth);
      b->pix3=XCreatePixmapFromBitmapData(disp,button_win->which_win,
                                          toggleonfill_bits,
                                          toggleonfill_width,toggleonfill_height,
                                          fcolor,bcolor,screen_depth);
      break;
    case FUNCTION:
      b->pix1=XCreatePixmapFromBitmapData(disp,button_win->which_win,
                                          function_bits,
                                          function_width,function_height,
                                          fcolor,bcolor,screen_depth);
      break;
    case FLOATPARM:
      b->pix1=XCreatePixmapFromBitmapData(disp,button_win->which_win,
                                          floatparm_bits,
                                          floatparm_width,floatparm_height,
                                          fcolor,bcolor,screen_depth);
      break;
    case MULTISTRINGPARM:
      b->pix1=XCreatePixmapFromBitmapData(disp,button_win->which_win,
                                          floatparm_bits,
                                          floatparm_width,floatparm_height,
                                          fcolor,bcolor,screen_depth);
      break;
    case STRINGPARM:
      b->pix1=XCreatePixmapFromBitmapData(disp,button_win->which_win,
                                          floatparm_bits,
                                          floatparm_width,floatparm_height,
                                          fcolor,bcolor,screen_depth);
      break;
    case INTPARM:
      b->pix1=XCreatePixmapFromBitmapData(disp,button_win->which_win,
                                          floatparm_bits,
                                          floatparm_width,floatparm_height,
                                          fcolor,bcolor,screen_depth);
      break;
    case MODE:
      b->pix1=XCreatePixmapFromBitmapData(disp,button_win->which_win,
                                          modal_bits,
                                          modal_width,modal_height,
                                          fcolor,bcolor,screen_depth);
      break;
    default:
      FATAL_BUG("Can't do a pixmap for a button.");
      break;
    }
    b = b->next;
  }
#endif
}



/****************************************************************************
 *
 *  Procedure build_main_button_window
 *
 *   Arguments: button_winP: pointer to pointer to button_win_list;
 *   Return Value: none
 *
 *  Actions: This procedure builds the main button list.  The buttons will be
 *    displayed in the order that they are put into the list.
 *
 *****************************************************************************/
void build_main_button_window(button_win_type **button_winP)
{
  button_type *b,*lastbutton;
  button_win_type *button_win;
  char string[80];
  int xdim,ydim;

  xdim = 0;
  ydim = 0;

  strcpy(string,"Main Options");

  /* open the window and get some space */
  *button_winP = new_button_win(*button_winP,string);
  button_win = *button_winP;

#ifdef X_GRAPHICS
  button_win->button_list = b = newbutton();
  b->xdim = modal_width;
  b->ydim = modal_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Main Mode");
  b->type = MODE;
  b->guts.mode.controlling_variable = &mainmode;
  strcpy(b->guts.mode.names[NONE],"Waiting");
  strcpy(b->guts.mode.names[ROT],"Rotate");
  strcpy(b->guts.mode.names[SCALE],"Scale");
  strcpy(b->guts.mode.names[TRANS],"Translate");
  strcpy(b->guts.mode.names[CENT],"Center");
  strcpy(b->guts.mode.names[CHOOSE],"Choose");
  b->guts.mode.num_modes = 6;
  lastbutton = b;

  if( global_read_on ){
    b = newbutton();
    lastbutton->next = b;
    b->xdim = function_width;
    b->ydim = function_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Read File");
    b->type = FUNCTION;
    b->guts.function.func = general_read;
    lastbutton = b;
  } else {
    b = newbutton();
    lastbutton->next = b;
    b->xdim = function_width;
    b->ydim = function_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Read Molecule");
    b->type = FUNCTION;
    b->guts.function.func = new_molecule;
    lastbutton = b;

    b = newbutton();
    lastbutton->next = b;
    b->xdim = function_width;
    b->ydim = function_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Read MO");
    b->type = FUNCTION;
    b->guts.function.func = new_MO_surface;
    lastbutton = b;

#ifdef INCLUDE_ADF_PLOTS
    b = newbutton();
    lastbutton->next = b;
    b->xdim = function_width;
    b->ydim = function_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Read ADF MO");
    b->type = FUNCTION;
    b->guts.function.func = new_ADF_MO_surface;
    lastbutton = b;
#endif

#ifdef INCLUDE_ADF_PLOTS
    b = newbutton();
    lastbutton->next = b;
    b->xdim = function_width;
    b->ydim = function_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Read ADF Vibn");
    b->type = FUNCTION;
    b->guts.function.func = new_vibration;
    lastbutton = b;
#endif

    b = newbutton();
    lastbutton->next = b;
    b->xdim = function_width;
    b->ydim = function_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Read Graph");
    b->type = FUNCTION;
    b->guts.function.func = new_graph;
    lastbutton = b;

    b = newbutton();
    lastbutton->next = b;
    b->xdim = function_width;
    b->ydim = function_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Read Contours");
    b->type = FUNCTION;
    b->guts.function.func = new_cont_plot;
    lastbutton = b;

    b = newbutton();
    lastbutton->next = b;
    b->xdim = function_width;
    b->ydim = function_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Read Props");
    b->type = FUNCTION;
    b->guts.function.func = new_prop_graph;
    lastbutton = b;

    b = newbutton();
    lastbutton->next = b;
    b->xdim = function_width;
    b->ydim = function_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Read Bands");
    b->type = FUNCTION;
    b->guts.function.func = new_band_graph;
    lastbutton = b;

    b = newbutton();
    lastbutton->next = b;
    b->xdim = function_width;
    b->ydim = function_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Read Walsh");
    b->type = FUNCTION;
    b->guts.function.func = new_walsh_graph;
    lastbutton = b;

    b = newbutton();
    lastbutton->next = b;
    b->xdim = function_width;
    b->ydim = function_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Read FMO");
    b->type = FUNCTION;
    b->guts.function.func = new_FMO_diagram;
    lastbutton = b;
  }
#if 0
  b = newbutton();
  lastbutton->next = b;
  b->xdim = function_width;
  b->ydim = function_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Add Label");
  b->type = FUNCTION;
  b->guts.function.func = new_label;
  lastbutton = b;
#endif
#endif
  b = newbutton();
  lastbutton->next = b;
  b->xdim = function_width;
  b->ydim = function_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Printing Options");
  b->type = FUNCTION;
  b->guts.function.func = open_ps_button_window;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = function_width;
  b->ydim = function_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Print");
  b->type = FUNCTION;
  b->guts.function.func = do_ps_output;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Fill Proj.s");
  b->type = TOGGLE;
  b->guts.boolean = &fill_projections;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = function_width;
  b->ydim = function_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Save All");
  b->type = FUNCTION;
  b->guts.function.func = save_all;

  lastbutton = b;
  b = newbutton();
  lastbutton->next = b;
  b->xdim = function_width;
  b->ydim = function_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Read All");
  b->type = FUNCTION;
  b->guts.function.func = read_all;
  lastbutton = b;


  b = newbutton();
  lastbutton->next = b;
  b->xdim = function_width;
  b->ydim = function_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Purge!");
  b->type = FUNCTION;
  b->guts.function.func = clearall;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = function_width;
  b->ydim = function_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Reset Labels");
  b->type = FUNCTION;
  b->guts.function.func = clear_labels;
  lastbutton = b;

  lastbutton->next = 0;
#ifdef X_GRAPHICS
  set_button_pixmaps(button_win);


  /* resize the window to fit the number of buttons */
  xdim += BUTTONOFF*2;
  ydim += BUTTONOFF*2;
  XResizeWindow(disp,button_win->which_win,xdim,ydim);
#endif
}


/****************************************************************************
 *
 *  Procedure build_PS_options_button_window
 *
 *   Arguments: button_winP: pointer to pointer to button_win_list;
 *   Return Value: none
 *
 *  Actions: This procedure builds the postscript options button list.
 *    The buttons will be
 *    displayed in the order that they are put into the list.
 *
 *****************************************************************************/
void build_PS_options_button_window(button_win_type **button_winP)
{
  button_type *b,*lastbutton;
  button_win_type *button_win;
  char string[80];
  int xdim,ydim;

  xdim = 0;
  ydim = 0;

  /* open the window and get some space */
  strcpy(string,"PS Options");
  *button_winP = new_button_win(*button_winP,string);
  button_win = *button_winP;

  button_win->button_list = b = newbutton();
  b->xdim = modal_width;
  b->ydim = modal_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Where: ");
  b->type = MODE;
  b->guts.mode.controlling_variable = &PS_options.where_to_print;
  strcpy(b->guts.mode.names[PS_PRINT_TOP],"Top");
  strcpy(b->guts.mode.names[PS_PRINT_MIDDLE],"Middle");
  strcpy(b->guts.mode.names[PS_PRINT_BOTTOM],"Bottom");
  b->guts.mode.num_modes = 3;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Font");
  b->type = STRINGPARM;
  b->guts.stringparm.num_lines = 1;
  b->guts.stringparm.lines[0] = PS_options.fontname;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Fontsize");
  b->type = FLOATPARM;
  b->guts.floatparm = &(PS_options.fontsize);
  lastbutton = b;


  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Tight Bounds?");
  b->type = TOGGLE;
  b->guts.boolean = &PS_options.tight_b_box;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Scale");
  b->type = FLOATPARM;
  b->guts.floatparm = &(PS_options.printscale);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = modal_width;
  b->ydim = modal_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Atoms: ");
  b->type = MODE;
  b->guts.mode.controlling_variable = &PS_options.atom_sphere_type;
  strcpy(b->guts.mode.names[ATOM_PLAIN_FILL],"Plain Atoms");
  strcpy(b->guts.mode.names[ATOM_SHADE_FILL],"Shaded Atoms");
  strcpy(b->guts.mode.names[ATOM_COLOR_FILL],"Color Atoms");
  strcpy(b->guts.mode.names[ATOM_CSHADE_FILL],"C-Shaded Atoms");
  b->guts.mode.num_modes = 4;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = modal_width;
  b->ydim = modal_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Bonds: ");
  b->type = MODE;
  b->guts.mode.controlling_variable = &PS_options.bond_type;
  strcpy(b->guts.mode.names[BOND_PLAIN],"Plain Bonds");
  strcpy(b->guts.mode.names[BOND_SHADE],"Shaded Bonds");
  strcpy(b->guts.mode.names[BOND_CSHADE],"Color Bonds");
  b->guts.mode.num_modes = 3;
  lastbutton = b;


  lastbutton->next = 0;

#ifdef X_GRAPHICS
  set_button_pixmaps(button_win);

  /* resize the window to fit the number of buttons */
  xdim += BUTTONOFF*2;
  ydim += BUTTONOFF*2;
  XResizeWindow(disp,button_win->which_win,xdim,ydim);
#endif
}

/****************************************************************************
 *
 *  Procedure build_label_button_window
 *
 *   Arguments: button_winP: pointer to pointer to button_win_list;
 *               label: pointer to label_type
 *   Return Value: none
 *
 *  Actions: This procedure builds the label options button list.
 *    The buttons will be
 *    displayed in the order that they are put into the list.
 *
 *****************************************************************************/
void build_label_button_window(button_win_type **button_winP,label_type *label)
{
  button_type *b,*lastbutton;
  button_win_type *button_win;
  char string[80];
  int xdim,ydim,i;

  xdim = 0;
  ydim = 0;

  /* open the window and get some space */
  strcpy(string,"Label Options");
  *button_winP = new_button_win(*button_winP,string);
  button_win = *button_winP;

  button_win->button_list = b = newbutton();
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Text");
  b->type = MULTISTRINGPARM;
  b->guts.stringparm.num_lines = MAX_LINES;
  for(i=0;i<b->guts.stringparm.num_lines;i++){
    b->guts.stringparm.lines[i] = label->text[i];
  }
  lastbutton = b;

  lastbutton->next = 0;

#ifdef X_GRAPHICS
  set_button_pixmaps(button_win);

  /* resize the window to fit the number of buttons */
  xdim += BUTTONOFF*2;
  ydim += BUTTONOFF*2;
  XResizeWindow(disp,button_win->which_win,xdim,ydim);
#endif
}


/****************************************************************************
 *
 *  Procedure build_molec_button_window
 *
 *   Arguments: button_winP: pointer to pointer to button_win_list;
 *                    molec: pointer to molec_type
 *
 *   Return Value: none
 *
 *  Actions: This procedure builds the button list for a molecule.  The buttons will
 *    be displayed in the order that they are put into the list.
 *
 *****************************************************************************/
void build_molec_button_window(button_win_type **button_winP,molec_type *molec)
{
  button_type *b,*lastbutton;
  button_win_type *button_win;
  char string[80];
  int xdim,ydim;


  xdim = 0;
  ydim = 0;

  /* open the window and get some space */
#ifdef X_GRAPHICS
  strcpy(string,"Molecule Options");
#endif
#ifdef MAC_GRAPHICS
        strcpy(string,"Molec");
#endif
  *button_winP = new_button_win(*button_winP,string);
  button_win = *button_winP;

  button_win->button_list = b = newbutton();
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Hydrogens?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->hydrogens_on;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Dummies?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->dummies_on;
  lastbutton = b;


  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Polyhedron?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->draw_polyhed;
  lastbutton = b;


  b = newbutton();
  lastbutton->next = b;
  b->xdim = function_width;
  b->ydim = function_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Center");
  b->type = FUNCTION;
  b->guts.function.func = center_molecule;
  b->guts.function.num_args = 1;
  b->guts.function.args[0] = (char *)molec;
  lastbutton = b;


  b = newbutton();
  lastbutton->next = b;
  b->xdim = function_width;
  b->ydim = function_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Hide Atoms");
  b->type = FUNCTION;
  b->guts.function.func = hide_atoms;
  b->guts.function.num_args = 1;
  b->guts.function.args[0] = (char *)molec;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = function_width;
  b->ydim = function_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Show Atoms");
  b->type = FUNCTION;
  b->guts.function.func = show_atoms;
  b->guts.function.num_args = 1;
  b->guts.function.args[0] = (char *)molec;
  lastbutton = b;


  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Axes?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->axes_on;
  lastbutton = b;


  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Outlines?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->outlines_on;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Shading?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->shading_on;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Crosses?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->crosses_on;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Connectors?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->draw_connectors;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Fancy Lines?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->fancy_lines;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Tubes?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->tubes_on;
  lastbutton = b;
#ifdef FORCED_PERSPECTIVE
  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"PTubes?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->trapezoidal_bonds;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Bond_Rad");
  b->type = FLOATPARM;
  b->guts.intparm = &(molec->bond_rad);
  lastbutton = b;
#endif

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Break Lines?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->breaking_lines;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Numbers?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->numbers_on;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Symbols?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->symbols_on;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Line Width");
  b->type = INTPARM;
  b->guts.intparm = &(molec->line_width);
  lastbutton = b;



  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Bond Tol");
  b->type = FLOATPARM;
  b->guts.floatparm = &(molec->bond_tol);
  lastbutton = b;

#ifdef INCLUDE_BOND_VALENCE
  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"ValenceTol?");
  b->type = TOGGLE;
  b->guts.boolean = &(molec->valence_for_bonds);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Bond Tol2");
  b->type = FLOATPARM;
  b->guts.floatparm = &(molec->bond_tol2);
  lastbutton = b;
#endif

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Rad Scale");
  b->type = FLOATPARM;
  b->guts.floatparm = &(molec->rad_mult);
  lastbutton = b;

  /* if there is a movie, here, allow the frame to be selected */
  if( molec->num_frames > 1 ){
    b = newbutton();
    lastbutton->next = b;
    b->xdim = floatparm_width;
    b->ydim = floatparm_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Frame");
    b->type = INTPARM;
    b->guts.intparm = &(molec->current_frame);
    lastbutton = b;
  }


  /*******
    if this is a solid, then offer the option of growing the xtal.
  *******/
  if( molec->num_dim >= 1 ){
    b = newbutton();
    lastbutton->next = b;
    b->xdim = toggleon_width;
    b->ydim = toggleon_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Lattice?");
    b->type = TOGGLE;
    b->guts.boolean = &(molec->draw_lattice);
    lastbutton = b;

    b = newbutton();
    lastbutton->next = b;
    b->xdim = toggleon_width;
    b->ydim = toggleon_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Box?");
    b->type = TOGGLE;
    b->guts.boolean = &(molec->draw_box);
    lastbutton = b;

    b = newbutton();
    lastbutton->next = b;
    b->xdim = function_width;
    b->ydim = function_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Grow Xtal!");
    b->type = FUNCTION;
    b->guts.function.func = grow_solid;
    b->guts.function.num_args = 1;
    b->guts.function.args[0] = (char *)molec;
    lastbutton = b;
  }

#ifdef INCLUDE_ADF_PLOTS
  /*******

    if this is a vibration, allow the scaling factor to be altered

  *******/
  if( molec->num_vibrations >= 1 ){
    b = newbutton();
    lastbutton->next = b;
    b->xdim = floatparm_width;
    b->ydim = floatparm_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Which Vibn");
    b->type = INTPARM;
    b->guts.intparm = &(molec->active_vibn);
    lastbutton = b;

    b = newbutton();
    lastbutton->next = b;
    b->xdim = floatparm_width;
    b->ydim = floatparm_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Vibn Scale");
    b->type = FLOATPARM;
    b->guts.floatparm = &(molec->vibration_scale);
    lastbutton = b;
  }
#endif
  lastbutton->next = 0;

#ifdef X_GRAPHICS
  set_button_pixmaps(button_win);

  /* resize the window to fit the number of buttons */
  xdim += BUTTONOFF*2;
  ydim += BUTTONOFF*2;
  XResizeWindow(disp,button_win->which_win,xdim,ydim);
#endif

}


/****************************************************************************
 *
 *  Procedure build_graph_button_window
 *
 *   Arguments: button_winP: pointer to pointer to button_win_list;
 *                    graph: pointer to graph_type
 *
 *   Return Value: none
 *
 *  Actions: This procedure builds the button list for a graph.  The buttons will
 *    be displayed in the order that they are put into the list.
 *
 *****************************************************************************/
void build_graph_button_window(button_win_type **button_winP,graph_type *graph)
{
  button_type *b,*lastbutton;
  button_win_type *button_win;
  char string[80];
  int xdim,ydim;
  int i;

  xdim = 0;
  ydim = 0;

  /* open the window and get some space */
#ifdef X_GRAPHICS
  strcpy(string,"Graph Options");
#endif
#ifdef MAC_GRAPHICS
        strcpy(string,"Graph");
#endif
  *button_winP = new_button_win(*button_winP,string);

  button_win = *button_winP;

  button_win->button_list = b = newbutton();
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Curve 1");
  b->type = CURVETOGGLE;
  b->guts.curve.boolean = &(graph->curves_to_display[0]);
  b->guts.curve.style = &(graph->styles[0]);
  lastbutton = b;

  for(i=1;i<graph->num_curves;i++){
    b = newbutton();
    lastbutton->next = b;
    b->xdim = toggleon_width;
    b->ydim = toggleon_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    sprintf(b->string,"Curve %d",i+1);
    b->type = CURVETOGGLE;
    b->guts.curve.boolean = &(graph->curves_to_display[i]);
    b->guts.curve.style = &(graph->styles[i]);
    lastbutton = b;
  }

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"X Min");
  b->type = FLOATPARM;
  b->guts.floatparm = &(graph->min_x);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"X Max");
  b->type = FLOATPARM;
  b->guts.floatparm = &(graph->max_x);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Y Min");
  b->type = FLOATPARM;
  b->guts.floatparm = &(graph->min_y);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Y Max");
  b->type = FLOATPARM;
  b->guts.floatparm = &(graph->max_y);
  lastbutton = b;


  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"X tics");
  b->type = TOGGLE;
  b->guts.boolean = &(graph->do_x_tics);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Y tics");
  b->type = TOGGLE;
  b->guts.boolean = &(graph->do_y_tics);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Show Title");
  b->type = TOGGLE;
  b->guts.boolean = &(graph->do_title);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Title");
  b->type = MULTISTRINGPARM;
  b->guts.stringparm.num_lines = NUM_TITLE_LINES;
  for(i=0;i<b->guts.stringparm.num_lines;i++){
    b->guts.stringparm.lines[i] = graph->title[i];
  }
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"X Legend");
  b->type = STRINGPARM;
  b->guts.stringparm.num_lines = 1;
  b->guts.stringparm.lines[0] = graph->xlegend;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Y Legend");
  b->type = STRINGPARM;
  b->guts.stringparm.num_lines = 1;
  b->guts.stringparm.lines[0] = graph->ylegend;
  lastbutton = b;

  lastbutton->next = 0;

#ifdef X_GRAPHICS
  set_button_pixmaps(button_win);

  /* resize the window to fit the number of buttons */
  xdim += BUTTONOFF*2;
  ydim += BUTTONOFF*2;
  XResizeWindow(disp,button_win->which_win,xdim,ydim);
#endif

}



/****************************************************************************
 *
 *  Procedure build_walsh_button_window
 *
 *   Arguments: button_winP: pointer to pointer to button_win_list;
 *                    graph: pointer to graph_type
 *
 *   Return Value: none
 *
 *  Actions: This procedure builds the button list for a walsh_graph.  The buttons will
 *    be displayed in the order that they are put into the list.
 *
 *****************************************************************************/
void build_walsh_button_window(button_win_type **button_winP,
                               walsh_graph_type *graph)
{
  button_type *b,*lastbutton;
  button_win_type *button_win;
  char string[80];
  int xdim,ydim;
  int i;

  xdim = 0;
  ydim = 0;

  /* open the window and get some space */
#ifdef X_GRAPHICS
  strcpy(string,"Walsh Options");
#endif
#ifdef MAC_GRAPHICS
        strcpy(string,"Walsh");
#endif
  *button_winP = new_button_win(*button_winP,string);

  button_win = *button_winP;


  button_win->button_list = b = newbutton();
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"MO's");
  b->type = CURVETOGGLE;
  b->guts.curve.boolean = &(graph->the_data->curves_to_display[0]);
  b->guts.curve.style = &(graph->the_data->styles[0]);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Total E");
  b->type = CURVETOGGLE;
  b->guts.curve.boolean = &(graph->the_data->curves_to_display[1]);
  b->guts.curve.style = &(graph->the_data->styles[1]);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"X Min");
  b->type = FLOATPARM;
  b->guts.floatparm = &(graph->the_data->min_x);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"X Max");
  b->type = FLOATPARM;
  b->guts.floatparm = &(graph->the_data->max_x);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Y Min");
  b->type = FLOATPARM;
  b->guts.floatparm = &(graph->the_data->min_y);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Y Max");
  b->type = FLOATPARM;
  b->guts.floatparm = &(graph->the_data->max_y);
  lastbutton = b;

    b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Tot E Min");
  b->type = FLOATPARM;
  b->guts.floatparm = &(graph->total_E->min_y);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Tot E Max");
  b->type = FLOATPARM;
  b->guts.floatparm = &(graph->total_E->max_y);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"X tics");
  b->type = TOGGLE;
  b->guts.boolean = &(graph->the_data->do_x_tics);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"MO tics");
  b->type = TOGGLE;
  b->guts.boolean = &(graph->the_data->do_y_tics);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Tot E tics");
  b->type = TOGGLE;
  b->guts.boolean = &(graph->total_E->do_y_tics);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Show Title");
  b->type = TOGGLE;
  b->guts.boolean = &(graph->the_data->do_title);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Title");
  b->type = MULTISTRINGPARM;
  b->guts.stringparm.num_lines = NUM_TITLE_LINES;
  for(i=0;i<b->guts.stringparm.num_lines;i++){
    b->guts.stringparm.lines[i] = graph->the_data->title[i];
  }
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"X Legend");
  b->type = STRINGPARM;
  b->guts.stringparm.num_lines = 1;
  b->guts.stringparm.lines[0] = graph->the_data->xlegend;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Y Legend");
  b->type = STRINGPARM;
  b->guts.stringparm.num_lines = 1;
  b->guts.stringparm.lines[0] = graph->the_data->ylegend;
  lastbutton = b;

  lastbutton->next = 0;

#ifdef X_GRAPHICS
  set_button_pixmaps(button_win);

  /* resize the window to fit the number of buttons */
  xdim += BUTTONOFF*2;
  ydim += BUTTONOFF*2;
  XResizeWindow(disp,button_win->which_win,xdim,ydim);
#endif

}


/****************************************************************************
 *
 *  Procedure build_band_button_window
 *
 *   Arguments: button_winP: pointer to pointer to button_win_list;
 *               band_graph: pointer to band_graph_type
 *
 *   Return Value: none
 *
 *  Actions: This procedure builds the button list for a band graph.  The buttons will
 *    be displayed in the order that they are put into the list.
 *
 *****************************************************************************/
void build_band_button_window(button_win_type **button_winP,
                              band_graph_type *band_graph)
{
  button_type *b,*lastbutton;
  button_win_type *button_win;
  char string[80];
  int xdim,ydim;
  int i;

  xdim = 0;
  ydim = 0;

  /* open the window and get some space */
#ifdef X_GRAPHICS
  strcpy(string,"Band Options");
#endif
#ifdef MAC_GRAPHICS
        strcpy(string,"Bands");
#endif
  *button_winP = new_button_win(*button_winP,string);

  button_win = *button_winP;

  button_win->button_list = b = newbutton();
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Bands");
  b->type = CURVETOGGLE;
  b->guts.curve.boolean = &(band_graph->the_data->curves_to_display[0]);
  b->guts.curve.style = &(band_graph->the_data->styles[0]);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Y Min");
  b->type = FLOATPARM;
  b->guts.floatparm = &(band_graph->the_data->min_y);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Y Max");
  b->type = FLOATPARM;
  b->guts.floatparm = &(band_graph->the_data->max_y);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"X Labels");
  b->type = TOGGLE;
  b->guts.boolean = &(band_graph->the_data->do_x_tics);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Y tics");
  b->type = TOGGLE;
  b->guts.boolean = &(band_graph->the_data->do_y_tics);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Show Title");
  b->type = TOGGLE;
  b->guts.boolean = &(band_graph->the_data->do_title);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Title");
  b->type = MULTISTRINGPARM;
  b->guts.stringparm.num_lines = NUM_TITLE_LINES;
  for(i=0;i<b->guts.stringparm.num_lines;i++){
    b->guts.stringparm.lines[i] = band_graph->the_data->title[i];
  }
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"X Legend");
  b->type = STRINGPARM;
  b->guts.stringparm.num_lines = 1;
  b->guts.stringparm.lines[0] = band_graph->the_data->xlegend;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Y Legend");
  b->type = STRINGPARM;
  b->guts.stringparm.num_lines = 1;
  b->guts.stringparm.lines[0] = band_graph->the_data->ylegend;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Show Fermi");
  b->type = TOGGLE;
  b->guts.boolean = &(band_graph->show_fermi);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Fermi E");
  b->type = FLOATPARM;
  b->guts.floatparm = &(band_graph->Fermi_E);
  lastbutton = b;

#ifdef SUPPORT_FATBANDS
  if( band_graph->num_fatbands ){
    b = newbutton();
    lastbutton->next = b;
    b->xdim = floatparm_width;
    b->ydim = floatparm_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Fat Bands");
    b->type = TOGGLE;
    b->guts.boolean = &(band_graph->fatbands_on);
    lastbutton = b;

    b = newbutton();
    lastbutton->next = b;
    b->xdim = modal_width;
    b->ydim = modal_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Filling");
    b->type = MODE;
    b->guts.mode.controlling_variable = &(band_graph->fatband_fill);
    strcpy(b->guts.mode.names[FATBANDS_SHADE],"Shade");
    strcpy(b->guts.mode.names[FATBANDS_LINE],"Lines");
    b->guts.mode.num_modes = 2;
    lastbutton = b;

    b = newbutton();
    lastbutton->next = b;
    b->xdim = floatparm_width;
    b->ydim = floatparm_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Fat Scale");
    b->type = FLOATPARM;
    b->guts.floatparm = &(band_graph->fatband_scale);
    lastbutton = b;
  }
#endif

  lastbutton->next = 0;

#ifdef X_GRAPHICS
  set_button_pixmaps(button_win);

  /* resize the window to fit the number of buttons */
  xdim += BUTTONOFF*2;
  ydim += BUTTONOFF*2;
  XResizeWindow(disp,button_win->which_win,xdim,ydim);
#endif

}


/****************************************************************************
 *
 *  Procedure build_prop_button_window
 *
 *   Arguments: button_winP: pointer to pointer to button_win_list;
 *               prop_graph: pointer to prop_graph_type
 *
 *   Return Value: none
 *
 *  Actions: This procedure builds the button list for a properties graph.
 *
 *****************************************************************************/
void build_prop_graph_button_window(button_win_type **button_winP,
                                    prop_graph_type *prop_graph)
{
  button_type *b,*lastbutton;
  button_win_type *button_win;
  char string[80];
  int xdim,ydim;
  int i;

  xdim = 0;
  ydim = 0;

  /* open the window and get some space */
#ifdef X_GRAPHICS
  strcpy(string,"Property Options");
#endif
#ifdef MAC_GRAPHICS
        strcpy(string,"Props");
#endif
  *button_winP = new_button_win(*button_winP,string);

  button_win = *button_winP;


  button_win->button_list = b = newbutton();
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Curve 1");
  b->type = CURVETOGGLE;
  b->guts.curve.boolean = &(prop_graph->the_data->curves_to_display[0]);
  b->guts.curve.style = &(prop_graph->the_data->styles[0]);
  b->guts.curve.fill = &(prop_graph->the_data->fills[0]);
  lastbutton = b;

  if( prop_graph->the_integration && prop_graph->the_integration->num_curves != 0 ){
    b = newbutton();
    lastbutton->next = b;
    b->xdim = toggleon_width;
    b->ydim = toggleon_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Integration 1");
    b->type = CURVETOGGLE;
    b->guts.curve.boolean =
      &(prop_graph->the_integration->curves_to_display[0]);
    b->guts.curve.style =
      &(prop_graph->the_integration->styles[0]);
    lastbutton = b;
  }

  for(i=1;i<prop_graph->the_data->num_curves;i++){
    b = newbutton();
    lastbutton->next = b;
    b->xdim = toggleon_width;
    b->ydim = toggleon_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    sprintf(b->string,"Curve %d",i+1);
    b->type = CURVETOGGLE;
    b->guts.curve.boolean =
      &(prop_graph->the_data->curves_to_display[i]);
    b->guts.curve.style =
      &(prop_graph->the_data->styles[i]);
    b->guts.curve.fill = &(prop_graph->the_data->fills[i]);
    lastbutton = b;

    if( prop_graph->the_integration && prop_graph->the_integration->num_curves != 0 ){
      b = newbutton();
      lastbutton->next = b;
      b->xdim = toggleon_width;
      b->ydim = toggleon_height;
      if( b->xdim > xdim ) xdim = b->xdim;
      ydim += b->ydim + BUTTONOFF;
      sprintf(b->string,"Integration %d",i+1);
      b->type = CURVETOGGLE;
      b->guts.curve.boolean =
        &(prop_graph->the_integration->curves_to_display[i]);
      b->guts.curve.style =
        &(prop_graph->the_integration->styles[i]);
      lastbutton = b;
    }
  }

  if( prop_graph->the_integration && prop_graph->the_integration->num_curves != 0 ){
    b = newbutton();
    lastbutton->next = b;
    b->xdim = toggleon_width;
    b->ydim = toggleon_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Integ Scale");
    b->type = TOGGLE;
    b->guts.boolean = &(prop_graph->integs_for_tics);
    lastbutton = b;
  }


  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"X Min");
  b->type = FLOATPARM;
  b->guts.floatparm = &(prop_graph->min_x);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"X Max");
  b->type = FLOATPARM;
  b->guts.floatparm = &(prop_graph->max_x);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Y Min");
  b->type = FLOATPARM;
  b->guts.floatparm = &(prop_graph->min_y);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Y Max");
  b->type = FLOATPARM;
  b->guts.floatparm = &(prop_graph->max_y);
  lastbutton = b;


  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"X tics");
  b->type = TOGGLE;
  b->guts.boolean = &(prop_graph->the_data->do_x_tics);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Y tics");
  b->type = TOGGLE;
  b->guts.boolean = &(prop_graph->the_data->do_y_tics);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Show Fermi");
  b->type = TOGGLE;
  b->guts.boolean = &(prop_graph->show_fermi);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Fermi E");
  b->type = FLOATPARM;
  b->guts.floatparm = &(prop_graph->Fermi_E);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Show Title");
  b->type = TOGGLE;
  b->guts.boolean = &(prop_graph->the_data->do_title);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Title");
  b->type = MULTISTRINGPARM;
  b->guts.stringparm.num_lines = NUM_TITLE_LINES;
  for(i=0;i<b->guts.stringparm.num_lines;i++){
    b->guts.stringparm.lines[i] = prop_graph->the_data->title[i];
  }
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"X Legend");
  b->type = STRINGPARM;
  b->guts.stringparm.num_lines = 1;
  b->guts.stringparm.lines[0] = prop_graph->the_data->xlegend;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Y Legend");
  b->type = STRINGPARM;
  b->guts.stringparm.num_lines = 1;
  b->guts.stringparm.lines[0] = prop_graph->the_data->ylegend;
  lastbutton = b;

  lastbutton->next = 0;

#ifdef X_GRAPHICS
  set_button_pixmaps(button_win);

  /* resize the window to fit the number of buttons */
  xdim += BUTTONOFF*2;
  ydim += BUTTONOFF*2;
  XResizeWindow(disp,button_win->which_win,xdim,ydim);
#endif

}

/****************************************************************************
 *
 *  Procedure build_FMO_button_window
 *
 *   Arguments: button_winP: pointer to pointer to button_win_list;
 *                  diagram: pointer to FMO_diagram_type
 *
 *   Return Value: none
 *
 *  Actions: This procedure builds the button list for an FMO interaction
 *    diagram.  The buttons will be displayed in the order that they are
 *    put into the list.
 *
 *****************************************************************************/
void build_FMO_button_window(button_win_type **button_winP,
                             FMO_diagram_type *diagram)
{
  button_type *b,*lastbutton;
  button_win_type *button_win;
  char string[80];
  int xdim,ydim;
  int i;

  xdim = 0;
  ydim = 0;

  /* open the window and get some space */
#ifdef X_GRAPHICS
  strcpy(string,"FMO Options");
#endif
#ifdef MAC_GRAPHICS
        strcpy(string,"FMO");
#endif
  *button_winP = new_button_win(*button_winP,string);

  button_win = *button_winP;

  button_win->button_list = b = newbutton();
  b->xdim = modal_width;
  b->ydim = modal_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Electrons");
  b->type = MODE;
  b->guts.mode.controlling_variable = &(diagram->electron_filling_mode);
  strcpy(b->guts.mode.names[FMO_FILL_NONE],"None");
  strcpy(b->guts.mode.names[FMO_FILL_HOMO],"HOMO");
  strcpy(b->guts.mode.names[FMO_FILL_ALL],"All");
  b->guts.mode.num_modes = 3;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Levels");
  b->type = TOGGLE;
  b->guts.boolean = &(diagram->show_data);
  lastbutton = b;

  if( diagram->num_frags > 0 ){
    b = newbutton();
    lastbutton->next = b;
    b->xdim = toggleon_width;
    b->ydim = toggleon_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Connectors");
    b->type = TOGGLE;
    b->guts.boolean = &(diagram->show_connects);
    lastbutton = b;
  }

  b = newbutton();
  lastbutton->next = b;
  b->xdim = function_width;
  b->ydim = function_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Change Tol");
  b->type = FUNCTION;
  b->guts.function.func = change_FMO_degen_tol;
  b->guts.function.num_args = 1;
  b->guts.function.args[0] = (char *)diagram;
  lastbutton = b;

  if( diagram->num_frags > 0 ){
    b = newbutton();
    lastbutton->next = b;
    b->xdim = toggleon_width;
    b->ydim = toggleon_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Left Fragment");
    b->type = INTPARM;
    b->guts.intparm = &(diagram->left_fragment);
    lastbutton = b;

    b = newbutton();
    lastbutton->next = b;
    b->xdim = toggleon_width;
    b->ydim = toggleon_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Right Fragment");
    b->type = INTPARM;
    b->guts.intparm = &(diagram->right_fragment);
    lastbutton = b;
  }
  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Y min");
  b->type = FLOATPARM;
  b->guts.floatparm = &(diagram->min_y);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Y max");
  b->type = FLOATPARM;
  b->guts.floatparm = &(diagram->max_y);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Level Width");
  b->type = INTPARM;
  b->guts.intparm = &(diagram->level_width);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Thickness");
  b->type = INTPARM;
  b->guts.intparm = &(diagram->level_thickness);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Electron len");
  b->type = INTPARM;
  b->guts.intparm = &(diagram->electron_length);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Y tics");
  b->type = TOGGLE;
  b->guts.boolean = &(diagram->do_y_tics);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Frame");
  b->type = TOGGLE;
  b->guts.boolean = &(diagram->show_box);
  lastbutton = b;


  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Show Title");
  b->type = TOGGLE;
  b->guts.boolean = &(diagram->do_title);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Title");
  b->type = MULTISTRINGPARM;
  b->guts.stringparm.num_lines = NUM_TITLE_LINES;
  for(i=0;i<b->guts.stringparm.num_lines;i++){
    b->guts.stringparm.lines[i] = diagram->title[i];
  }
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Main Label");
  b->type = STRINGPARM;
  b->guts.stringparm.num_lines = 1;
  b->guts.stringparm.lines[0] = diagram->label;
  lastbutton = b;

  for(i=0;i<diagram->num_frags;i++){
    b = newbutton();
    lastbutton->next = b;
    b->xdim = floatparm_width;
    b->ydim = floatparm_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    sprintf(b->string,"Frag %d Label",i+1);
    b->type = STRINGPARM;
    b->guts.stringparm.num_lines = 1;
    b->guts.stringparm.lines[0] = diagram->frags[i].label;
    lastbutton = b;
  }


  lastbutton->next = 0;

#ifdef X_GRAPHICS
  set_button_pixmaps(button_win);

  /* resize the window to fit the number of buttons */
  xdim += BUTTONOFF*2;
  ydim += BUTTONOFF*2;
  XResizeWindow(disp,button_win->which_win,xdim,ydim);
#endif

}


/****************************************************************************
 *
 *  Procedure build_MO_surf_button_window
 *
 *   Arguments: button_winP: pointer to pointer to button_win_list;
 *                MO_surf: pointer to MO_surface_type
 *
 *   Return Value: none
 *
 *  Actions: This procedure builds the button list for an MO isosurface
 *    The buttons will be displayed in the order that they are put into the
 *     list.
 *
 *****************************************************************************/
void build_MO_surf_button_window(button_win_type **button_winP,
                                 MO_surface_type *MO_surf)
{
  button_type *b,*lastbutton;
  button_win_type *button_win;
  char string[80];
  int xdim,ydim;
  molec_type *molec;

  xdim = 0;
  ydim = 0;

  molec = MO_surf->molec;

  /* open the window and get some space */
#ifdef X_GRAPHICS
  strcpy(string,"MO Options");
#endif
#ifdef MAC_GRAPHICS
        strcpy(string,"MO");
#endif
  *button_winP = new_button_win(*button_winP,string);

  button_win = *button_winP;

  button_win->button_list = b = newbutton();
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Surface?");
  b->type = TOGGLE;
  b->guts.boolean = &MO_surf->display_surf;
  lastbutton = b;



  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Molecule?");
  b->type = TOGGLE;
  b->guts.boolean = &MO_surf->display_molec;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Isosurface");
  b->type = FLOATPARM;
  b->guts.floatparm = &(MO_surf->surface_value);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Voxel Size");
  b->type = FLOATPARM;
  b->guts.floatparm = &(MO_surf->voxel_size);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Search Rad");
  b->type = FLOATPARM;
  b->guts.floatparm = &(MO_surf->search_radius);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Slop");
  b->type = FLOATPARM;
  b->guts.floatparm = &(MO_surf->slop);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = function_width;
  b->ydim = function_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Exclude Atoms");
  b->type = FUNCTION;
  b->guts.function.func = exclude_atoms;
  b->guts.function.num_args = 1;
  b->guts.function.args[0] = (char *)MO_surf;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = function_width;
  b->ydim = function_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Include Atoms");
  b->type = FUNCTION;
  b->guts.function.func = include_atoms;
  b->guts.function.num_args = 1;
  b->guts.function.args[0] = (char *)MO_surf;
  lastbutton = b;


  b = newbutton();
  lastbutton->next = b;
  b->xdim = function_width;
  b->ydim = function_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Surf evolve");
  b->type = FUNCTION;
  b->guts.function.func = construct_MO_isosurface;
  b->guts.function.num_args = 1;
  b->guts.function.args[0] = (char *)MO_surf;
  lastbutton = b;
#if 0
  b = newbutton();
  lastbutton->next = b;
  b->xdim = function_width;
  b->ydim = function_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Save Surf");
  b->type = FUNCTION;
  b->guts.function.func = save_triangle_locs;
  b->guts.function.num_args = 1;
  b->guts.function.args[0] = (char *)MO_surf;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = function_width;
  b->ydim = function_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Read Surf");
  b->type = FUNCTION;
  b->guts.function.func = read_triangle_locs;
  b->guts.function.num_args = 1;
  b->guts.function.args[0] = (char *)MO_surf;
  lastbutton = b;
#endif
  b = newbutton();
  lastbutton->next = b;
  b->xdim = function_width;
  b->ydim = function_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Change MO");
  b->type = FUNCTION;
  b->guts.function.func = change_active_MO;
  b->guts.function.num_args = 1;
  b->guts.function.args[0] = (char *)MO_surf;
  lastbutton = b;


  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Tri_Shading?");
  b->type = TOGGLE;
  b->guts.boolean = &MO_surf->do_shading;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Tri_Outlines?");
  b->type = TOGGLE;
  b->guts.boolean = &MO_surf->do_lines;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Contours?");
  b->type = TOGGLE;
  b->guts.boolean = &MO_surf->display_conts;
  lastbutton = b;



  b = newbutton();
  lastbutton->next = b;
  b->xdim = modal_width;
  b->ydim = modal_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Kind");
  b->type = MODE;
  b->guts.mode.controlling_variable = &(MO_surf->levels_kind);
  strcpy(b->guts.mode.names[LEVELS_AUTO],"Auto");
  strcpy(b->guts.mode.names[LEVELS_INCREMENTAL],"Incremental");
  strcpy(b->guts.mode.names[LEVELS_DISCRETE],"Discrete");
  b->guts.mode.num_modes = 3;
  lastbutton = b;


  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Planes");
  b->type = STRINGPARM;
  b->guts.stringparm.num_lines = 1;
  b->guts.stringparm.lines[0] = MO_surf->plane_dirs;
  lastbutton = b;



  b = newbutton();
  lastbutton->next = b;
  b->xdim = function_width;
  b->ydim = function_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Contour it!");
  b->type = FUNCTION;
  b->guts.function.func = construct_MO_contours;
  b->guts.function.num_args = 1;
  b->guts.function.args[0] = (char *)MO_surf;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = function_width;
  b->ydim = function_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Contour Surf");
  b->type = FUNCTION;
  b->guts.function.func = MO_contour_surf;
  b->guts.function.num_args = 1;
  b->guts.function.args[0] = (char *)MO_surf;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Sort Contours?");
  b->type = TOGGLE;
  b->guts.boolean = &MO_surf->sort_conts;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Invert Phase?");
  b->type = TOGGLE;
  b->guts.boolean = &MO_surf->invert_phase;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Hidden?");
  b->type = TOGGLE;
  b->guts.boolean = &MO_surf->hidden_surf;
  lastbutton = b;


#if 0

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Hydrogens?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->hydrogens_on;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Dummies?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->dummies_on;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = function_width;
  b->ydim = function_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Hide Atoms");
  b->type = FUNCTION;
  b->guts.function.func = hide_atoms;
  b->guts.function.num_args = 1;
  b->guts.function.args[0] = (char *)molec;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = function_width;
  b->ydim = function_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Show Atoms");
  b->type = FUNCTION;
  b->guts.function.func = show_atoms;
  b->guts.function.num_args = 1;
  b->guts.function.args[0] = (char *)molec;
  lastbutton = b;


  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Axes?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->axes_on;
  lastbutton = b;



  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Outlines?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->outlines_on;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Shading?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->shading_on;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Crosses?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->crosses_on;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Connectors?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->draw_connectors;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Fancy Lines?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->fancy_lines;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Break Lines?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->breaking_lines;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Numbers?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->numbers_on;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Symbols?");
  b->type = TOGGLE;
  b->guts.boolean = &molec->symbols_on;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Line Width");
  b->type = INTPARM;
  b->guts.intparm = &(molec->line_width);
  lastbutton = b;


  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Bond Tol");
  b->type = FLOATPARM;
  b->guts.floatparm = &(molec->bond_tol);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Rad Scale");
  b->type = FLOATPARM;
  b->guts.floatparm = &(molec->rad_mult);
  lastbutton = b;
#endif

  /*******
    if this is a solid, then offer the option of growing the xtal.
  *******/
  if( molec->num_dim >= 1 ){
    b = newbutton();
    lastbutton->next = b;
    b->xdim = function_width;
    b->ydim = function_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    strcpy(b->string,"Grow Xtal!");
    b->type = FUNCTION;
    b->guts.function.func = grow_solid_with_surface;
    b->guts.function.num_args = 2;
    b->guts.function.args[0] = (char *)molec;
    b->guts.function.args[1] = (char *)MO_surf;
    lastbutton = b;
  }

  lastbutton->next = 0;

#ifdef X_GRAPHICS
  set_button_pixmaps(button_win);

  /* resize the window to fit the number of buttons */
  xdim += BUTTONOFF*2;
  ydim += BUTTONOFF*2;
  XResizeWindow(disp,button_win->which_win,xdim,ydim);
#endif

}


/****************************************************************************
 *
 *  Procedure build_cont_plot_button_window
 *
 *   Arguments: button_winP: pointer to pointer to button_win_list;
 *                    graph: pointer to graph_type
 *
 *   Return Value: none
 *
 *  Actions: This procedure builds the button list for a contour plot.
 *
 *****************************************************************************/
void build_cont_plot_button_window(button_win_type **button_winP,
                                   contour_plot_type *cont_plot)
{
  button_type *b,*lastbutton;
  button_win_type *button_win;
  char string[80];
  int xdim,ydim;
  int i;

  xdim = 0;
  ydim = 0;

  /* open the window and get some space */
#ifdef X_GRAPHICS
  strcpy(string,"Cont_Plot Options");
#endif
#ifdef MAC_GRAPHICS
        strcpy(string,"Cont_Plot");
#endif
  *button_winP = new_button_win(*button_winP,string);

  button_win = *button_winP;

  button_win->button_list = b = newbutton();
  b->xdim = modal_width;
  b->ydim = modal_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Interp");
  b->type = MODE;
  b->guts.mode.controlling_variable = &(cont_plot->interp_kind);
  strcpy(b->guts.mode.names[INTERP_NOTHING],"None");
  strcpy(b->guts.mode.names[INTERP_CUBIC],"Cubic");
  strcpy(b->guts.mode.names[APPROX_BSPLINE],"B spline");
  b->guts.mode.num_modes = 3;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = modal_width;
  b->ydim = modal_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Kind");
  b->type = MODE;
  b->guts.mode.controlling_variable = &(cont_plot->levels_kind);
  strcpy(b->guts.mode.names[LEVELS_AUTO],"Auto");
  strcpy(b->guts.mode.names[LEVELS_INCREMENTAL],"Incremental");
  strcpy(b->guts.mode.names[LEVELS_DISCRETE],"Discrete");
  b->guts.mode.num_modes = 3;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Num Levels");
  b->type = INTPARM;
  b->guts.intparm = &(cont_plot->num_levels);
  lastbutton = b;


  b = newbutton();
  lastbutton->next = b;
  b->xdim = function_width;
  b->ydim = function_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Contour!");
  b->type = FUNCTION;
  b->guts.function.func = contour_data_set;
  b->guts.function.num_args = 1;
  b->guts.function.args[0] = (char *)cont_plot;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Curve 1");
  b->type = CURVETOGGLE;
  b->guts.curve.boolean = &(cont_plot->curves_to_display[0]);
  b->guts.curve.style = &(cont_plot->styles[0]);
  lastbutton = b;

  for(i=1;i<cont_plot->num_curves;i++){
    b = newbutton();
    lastbutton->next = b;
    b->xdim = toggleon_width;
    b->ydim = toggleon_height;
    if( b->xdim > xdim ) xdim = b->xdim;
    ydim += b->ydim + BUTTONOFF;
    sprintf(b->string,"Curve %d",i+1);
    b->type = CURVETOGGLE;
    b->guts.curve.boolean = &(cont_plot->curves_to_display[i]);
    b->guts.curve.style = &(cont_plot->styles[i]);
    lastbutton = b;
  }

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Z Min");
  b->type = FLOATPARM;
  b->guts.floatparm = &(cont_plot->min_z);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Z Max");
  b->type = FLOATPARM;
  b->guts.floatparm = &(cont_plot->max_z);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"X Min");
  b->type = FLOATPARM;
  b->guts.floatparm = &(cont_plot->min_x);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"X Max");
  b->type = FLOATPARM;
  b->guts.floatparm = &(cont_plot->max_x);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Y Min");
  b->type = FLOATPARM;
  b->guts.floatparm = &(cont_plot->min_y);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Y Max");
  b->type = FLOATPARM;
  b->guts.floatparm = &(cont_plot->max_y);
  lastbutton = b;


  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"X tics");
  b->type = TOGGLE;
  b->guts.boolean = &(cont_plot->do_x_tics);
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = toggleon_width;
  b->ydim = toggleon_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Y tics");
  b->type = TOGGLE;
  b->guts.boolean = &(cont_plot->do_y_tics);
  lastbutton = b;


  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"X Legend");
  b->type = STRINGPARM;
  b->guts.stringparm.num_lines = 1;
  b->guts.stringparm.lines[0] = cont_plot->xlegend;
  lastbutton = b;

  b = newbutton();
  lastbutton->next = b;
  b->xdim = floatparm_width;
  b->ydim = floatparm_height;
  if( b->xdim > xdim ) xdim = b->xdim;
  ydim += b->ydim + BUTTONOFF;
  strcpy(b->string,"Y Legend");
  b->type = STRINGPARM;
  b->guts.stringparm.num_lines = 1;
  b->guts.stringparm.lines[0] = cont_plot->ylegend;
  lastbutton = b;

  lastbutton->next = 0;

#ifdef X_GRAPHICS
  set_button_pixmaps(button_win);

  /* resize the window to fit the number of buttons */
  xdim += BUTTONOFF*2;
  ydim += BUTTONOFF*2;
  XResizeWindow(disp,button_win->which_win,xdim,ydim);
#endif

}
