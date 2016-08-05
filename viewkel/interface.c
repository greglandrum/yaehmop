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
   26.04.98 gL:
     added "camera motion"  okay, so it's only translations
       along z and focal length adjustment, but that's better
       than nothing.
     added near plane clipping toggling

   22.05.98 gL:
     added event.xbutton.state to calls to button handling
     routines to allow shift/control clicks.

   30.05.98 gL:
     added toggling for polyhedron outlining.

   09.07.98 gL:
      switched '1' action to dump VRML instead of rayshade.

   18.08.98 gL:
     modified action upon hitting '1' in the main window.
     now this prompts for dumping a VRML file or a rayshade file.

   04.09.98 gL:
     added "support" for click and drag selection.
     This is implemented in an unpleasing fashion, but it
     works and I need to do my research proposals.

   20.09.98 gL:
     updated documentation

   26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)

   18.01.99 gL:
     added toggling of dump_grids_on

   17.02.99 gL:
     added proportionate scaling of contour plots.

   18.02.99 gL:
     fixed said functionality so that it's actually possible to
     get into scale mode with a contour plot (duh).

   01.03.99 gL:
     various modifications to make all of the selection stuff
     also work with MO surfaces.

  06.03.99 gL:
     support for printing curve names added.

  06.04.99 gL:
     fixed stupid core dumping bug when the window was clicked in
     without an active object.  Ahem!

  12.04.99 gL:
     modified action upon hitting '1' in the main window yet again.
     now it makes it clear that it is possible to cancel.

  26.06.1999 gL:
     Bogosity in the if statement for the '+' action in do_keypress fixed.
     It is now possible to unhide hidden atoms in MO surfaces.

     Tolerance for starting to go into draw_select mode increased to 10 pixels.
     Line width for the dragged out box is now explicitly set.

  28.06.1999 gL:
     the handling of MotionNotify events when whichobj was null caused
     seg faults.  this has been fixed.

  29.06.1999 gL:
     added support for dumping MO volumes (binary files which contain MO data)
***/


/****************************************************************************
 *
 *          interface.c
 * This file contains the event loop and user interface routines for the
 * program.
 *
 *
 *  To port this to another event driven, interactive system is going to
 *    require some work... I'm not doing it now  gL Sept. 1994
 *
 *
 *
 *
 ****************************************************************************/
#include "viewkel.h"
#include <time.h>

FILE *rot_pipe=0;
int first_rot_pipe_read = 1;


/****************************************************************************
 *
 *                   Procedure show_keys
 *
 * Arguments: none
 * Returns: none
 *
 * Action: This routine shows the valid keys and a short description
 *    of their function
 *
 ****************************************************************************/
void show_keys()
{
  printf("Here are the valid keys, case is not important.\n");
  printf("\t\tthe marker [3D] indicates that this key only works for\n\
\t\tthree-dimensional objects (molecules, surfaces, etc.)\n");
  printf("\ta: [3D] When an atom is selected: alter properties of the atom.\n \
\t\tWhen two atoms are selected: alter properties of the bond.\n \
\tc: Switch to center mode.\n \
\t\t[3D] In choose mode with one atom selected, this allows alteration\n \
\t\tof the color of the selected atom. \n \
\td: [3D] Print current molecular coordinates in window\n \
\tr: Switch to rotate mode.\n \
\ts: Switch to scale mode.\n \
\t[3D] In choose mode with one atom selected, this allows alteration\n \
\t\tof the gray shade of the selected atom.\n \
\t\tIf a contour plot is active, then the y scaling will be adjusted \n \
\t\tto match the x scaling.\n \
\tx: [3D] Set view down x axis.\n \
\ty: [3D] Set view down y axis.\n \
\tz: [3D] Set view down z axis.\n \
\th: [3D] Hide selected atoms.\n \
\t+: [3D] Show hidden atoms.\n \
\t`: [3D] Inverts selected atoms (i.e. selected <-> unselected) \n \
\t): [3D] If the selected molecule has animation associated with it\n \
\t\t (e.g. Vibrations, Trajectories, etc.), then move to the next frame.\n \
\t(: [3D] If the selected molecule has animation associated with it\n \
\t\t (e.g. Vibrations, Trajectories, etc.), then move to the previous frame.\n \
\t1: [3D] Create a rayshade or VRML input file for the currently selected object.\n \
\t2: [3D] Generates a coordination polyhedron for the selected atom(s).\n \
\t4: [3D] If there are animated objects, the animation will be shown/stopped.\n \
\t5: [3D] Generate the coordination environment of the selected atom.\n \
\to: [3D] Toggles outlining of polyhedra on screen.\n \
\tg: [3D] Toggles grid dumping when an MO is contoured.\n \
\t!: [3D] Toggles near plane clipping.\n \
\t<Delete>: Deletes active object.\n \
\t<Tab>: Insert a text label onscreen.\n \
\t?: Show this message.\n \
\tq: Quit the program now.\n");

}

/****************************************************************************
 *
 *                   Procedure do_keypress
 *
 * Arguments: event: an event structure
 * Returns: none
 *
 * Action: This routine processes the users key presses
 *
 *
 *  NOTE: on the Mac, we should *not* be in this procedure if the command
 *          key is down.  Deal with command key presses in the hooks.c part
 *          to keep this cross-system code as clean as possible.
 *
 ****************************************************************************/
#ifdef X_GRAPHICS
void do_keypress( XEvent event )
#endif
#ifdef MAC_GRAPHICS
void do_keypress(EventRecord *event)
#endif
{
  int i;
  char string[4],instring[80]="V";
  head_type *head1,*head2;
  string[0]=0;
  /* find out what the key was */
#ifdef X_GRAPHICS
  XLookupString(&event.xkey,string,1,0,0);
#endif
#ifdef MAC_GRAPHICS
  string[0] = event->message & charCodeMask;
#endif

  switch( string[0] ){
  case '\t':
    new_label(0);

  case 'A':
  case 'a':
      if( mainmode == CHOOSE ){
        if( num_selected > 0 && (whichobj->prim->molec ||
         (whichobj->prim->MO_surf && whichobj->prim->MO_surf->molec))){
          adjust_style(num_selected,whichobj);
          unselect_all_atoms(num_selected,whichobj);
          num_selected = 0;
        } else if (whichobj->prim->label){
          adjust_label(whichobj->prim->label);
        }
      }
      break;
  case 'C':
  case 'c':
    if( mainmode == CHOOSE && num_selected == 1 && (whichobj->prim->molec ||
        (whichobj->prim->MO_surf && whichobj->prim->MO_surf->molec))){
      adjust_color(num_selected,whichobj,ATOM_COLOR_FILL);
      unselect_all_atoms(num_selected,whichobj);
      num_selected = 0;
    }
    else {
      mainmode = CENT;
    }
#ifdef X_GRAPHICS
    g_draw_all_buttons(button_wins);
#endif
    break;
  case ' ':
    mainmode = CHOOSE;
#ifdef X_GRAPHICS
    g_draw_all_buttons(button_wins);
#endif
    break;
  case 'D':
  case 'd':
    /**** insert some error checking here ******/
    if( whichobj->prim->molec )
      dump_molecular_coords(whichobj->prim->molec);
    break;
  case 'i':
    if(mainmode==CENT) if(whichobj) whichobj->cent.y -= 1;
    if(mainmode==TRANS) if(whichobj) whichobj->trans.y += 1;
    if(mainmode==SCALE) if(whichobj) whichobj->scale.y+=.1;
    if(mainmode==ROT) if(whichobj) whichobj->rot.y += .05;
    break;
  case 'I':
    if(mainmode==CENT) if(whichobj) whichobj->cent.y -= 10;
    if(mainmode==TRANS) if(whichobj) whichobj->trans.y += 10;
    if(mainmode==SCALE) if(whichobj) whichobj->scale.y+=.5;
    if(mainmode==ROT) if(whichobj) whichobj->rot.y += .2;
    break;
  case 'j':
    if(mainmode==CENT) if(whichobj) whichobj->cent.x -= 1;
    if(mainmode==TRANS) if(whichobj) whichobj->trans.x -= 1;
    if(mainmode==SCALE) if(whichobj) whichobj->scale.x -=.1;
    if(mainmode==ROT) if(whichobj) whichobj->rot.x -= .05;
    break;
  case 'J':
    if(mainmode==CENT) if(whichobj) whichobj->cent.x -= 10;
    if(mainmode==TRANS) if(whichobj) whichobj->trans.x -= 10;
    if(mainmode==SCALE) if(whichobj) whichobj->scale.x -=.5;
    if(mainmode==ROT) if(whichobj) whichobj->rot.x -= .2;
    break;
  case 'k':
    if(mainmode==CENT) if(whichobj) whichobj->cent.y += 1;
    if(mainmode==TRANS) if(whichobj) whichobj->trans.y -= 1;
    if(mainmode==SCALE) if(whichobj) whichobj->scale.y -= .1;
    if(mainmode==ROT) if(whichobj) whichobj->rot.y -= .05;
    break;
  case 'K':
    if(mainmode==CENT) if(whichobj) whichobj->cent.y += 10;
    if(mainmode==TRANS) if(whichobj) whichobj->trans.y -= 10;
    if(mainmode==SCALE) if(whichobj) whichobj->scale.y -=.5;
    if(mainmode==ROT) if(whichobj) whichobj->rot.y -= .2;
    break;
  case 'l':
    if(mainmode==CENT) if(whichobj) whichobj->cent.x += 1;
    if(mainmode==TRANS) if(whichobj) whichobj->trans.x += 1;
    if(mainmode==SCALE) if(whichobj) whichobj->scale.x+=.1;
    if(mainmode==ROT) if(whichobj) whichobj->rot.x += .05;
    break;
  case 'L':
    if(mainmode==CENT) if(whichobj) whichobj->cent.x += 10;
    if(mainmode==TRANS) if(whichobj) whichobj->trans.x += 10;
    if(mainmode==SCALE) if(whichobj) whichobj->scale.x+=.5;
    if(mainmode==ROT) if(whichobj) whichobj->rot.x += .2;
    break;

  case 'p':
    if(mainmode==TRANS) if(whichobj) whichobj->trans.z += 1;
    if(mainmode==SCALE) if(whichobj) whichobj->scale.z+=.1;
    if(mainmode==ROT) if(whichobj) whichobj->rot.z += .05;
    break;
  case 'P':
    if(mainmode==TRANS) if(whichobj) whichobj->trans.z += 10;
    if(mainmode==SCALE) if(whichobj) whichobj->scale.z+=.5;
    if(mainmode==ROT) if(whichobj) whichobj->rot.z += .2;
    break;
  case ';':
    if(mainmode==TRANS) if(whichobj) whichobj->trans.z -= 1;
    if(mainmode==SCALE) if(whichobj) whichobj->scale.z -= .1;
    if(mainmode==ROT) if(whichobj) whichobj->rot.z -= .05;
    break;
  case ':':
    if(mainmode==TRANS) if(whichobj) whichobj->trans.z -= 10;
    if(mainmode==SCALE) if(whichobj) whichobj->scale.z -= .5;
    if(mainmode==ROT) if(whichobj) whichobj->rot.z -= .2;
    break;

  case 'r':
  case 'R':
    mainmode=ROT;
#ifdef X_GRAPHICS
    g_draw_all_buttons(button_wins);
#endif
    break;
  case 'S':
  case 's':
      if( mainmode == CHOOSE && num_selected == 1 && (whichobj->prim->molec ||
          (whichobj->prim->MO_surf && whichobj->prim->MO_surf->molec))){
        adjust_color(num_selected,whichobj,ATOM_SHADE_FILL);
        unselect_all_atoms(num_selected,whichobj);
        num_selected = 0;
      }
      else if(mainmode == SCALE && whichobj->prim->cont_plot){
        whichobj->scale.y = whichobj->scale.x *
          (whichobj->prim->cont_plot->max_y-whichobj->prim->cont_plot->min_y) /
          (whichobj->prim->cont_plot->max_x-whichobj->prim->cont_plot->min_x);
      } else{
        mainmode=SCALE;
      }
#ifdef X_GRAPHICS
    g_draw_all_buttons(button_wins);
#endif
    break;
#if 0
  case 'S':
    /* make the print out small */
    PS_options.fontsize = 8;
    PS_options.printscale = 0.6;
    if(whichobj->prim->molec) whichobj->prim->molec->rad_mult = 0.5;
    if(whichobj->prim->MO_surf && whichobj->prim->MO_surf->molec)
      whichobj->prim->MO_surf->molec->rad_mult = 0.5;
#endif
  case 't':
    mainmode=TRANS;
#ifdef X_GRAPHICS
    g_draw_all_buttons(button_wins);
#endif
    break;
  case 'z':
  case 'Z':
    whichobj->rot.x = 0.0;
    whichobj->rot.y = 0.0;
    whichobj->rot.z = 0.0;
    break;
  case 'x':
  case 'X':
    whichobj->rot.x = 0.0;
    whichobj->rot.y = 0.5*M_PI;
    whichobj->rot.z = 0.0;
    break;
  case 'y':
  case 'Y':
    whichobj->rot.x = 0.5*M_PI;
    whichobj->rot.y = 0.0;
    whichobj->rot.z = 0.0;
    break;

  case 'n':
  case 'N':
    if( whichobj->prim->prop_graph){
      printf("Curve names:\n");
      for(i=0;i<whichobj->prim->prop_graph->the_data->num_curves;i++){
        if( whichobj->prim->prop_graph->the_data->curve_names[i*NORMAL_STR_LEN] ){
          printf("\t%d: %s\n",i+1,
                 &(whichobj->prim->prop_graph->the_data->curve_names[i*NORMAL_STR_LEN]));
        }
      }
    }
    break;
  case 'h':
  case 'H':
    if( mainmode == CHOOSE ){
      if( num_selected > 0 && (whichobj->prim->molec ||
         (whichobj->prim->MO_surf && whichobj->prim->MO_surf->molec))){
        hide_selected_atoms(num_selected,whichobj);
        unselect_all_atoms(num_selected,whichobj);
        num_selected = 0;
      }
    }
    break;

  case '`':
    if(mainmode == CHOOSE && whichobj &&
       whichobj->prim &&
       (whichobj->prim->molec || (whichobj->prim->MO_surf &&
                                  whichobj->prim->MO_surf->molec))){
      invert_selected_atoms(&num_selected,whichobj);
    }
    break;
  case '+':
    if(whichobj && whichobj->prim &&
       (whichobj->prim->molec || (whichobj->prim->MO_surf &&
                                  whichobj->prim->MO_surf->molec))){
      show_all_atoms(whichobj);
    }
    break;

  case '':
    if(whichobj){
      free_obj(whichobj);
      if( head ){
        if( head->obj == whichobj ){
          D_FREE(head->obj);
          if( head->next ){
            head1 = head->next;
            D_FREE(head);
            head = head1;
          } else{
            D_FREE(head);
            head = 0;
          }
        } else {
          head1 = head->next;
          head2 = head;
          while(head1 && head1->obj != whichobj){
            head2 = head1;
            head1 = head1->next;
          }
          if( !head1 ) FATAL_BUG("headless object removed... ACK!");
          head2->next = head1->next;
          D_FREE(head1->obj);
          D_FREE(head1);
        }
      } else {
        FATAL_BUG("I just removed an object, but there's no head.");
      }
    }
    whichobj = 0;
    break;
  case ')':
    if(whichobj && whichobj->prim && whichobj->prim->molec &&
       whichobj->prim->molec->num_frames > 1){
      whichobj->prim->molec->current_frame++;
      if(whichobj->prim->molec->current_frame >=
         whichobj->prim->molec->num_frames){
        whichobj->prim->molec->current_frame =
          whichobj->prim->molec->current_frame %
          whichobj->prim->molec->num_frames;
      }
#ifdef X_GRAPHICS
      if( doing_X ){
        /* redraw the buttons */
        g_draw_all_buttons(button_wins);
      }
#endif
    }
    break;
  case '(':
    if(whichobj && whichobj->prim && whichobj->prim->molec &&
       whichobj->prim->molec->num_frames > 1){
      whichobj->prim->molec->current_frame--;
      if(whichobj->prim->molec->current_frame >=
         whichobj->prim->molec->num_frames){
        whichobj->prim->molec->current_frame =
          whichobj->prim->molec->current_frame %
          whichobj->prim->molec->num_frames;
      }
      if(whichobj->prim->molec->current_frame < 0)
        whichobj->prim->molec->current_frame =
          whichobj->prim->molec->num_frames-1;
#ifdef X_GRAPHICS
      if( doing_X ){
        /* redraw the buttons */
        g_draw_all_buttons(button_wins);
      }
#endif
    }
    break;
  case '2':
    if( whichobj && whichobj->prim && whichobj->prim->molec ){
      if( num_selected == 0 ){
        printf("You must select at least one atom first.\n");
      } else {
        gen_coord_polyhed(whichobj->prim->molec);
      }
    }
    break;
  case '5':
    if(whichobj->prim->molec && num_selected == 1){
      show_coord_env(num_selected,whichobj);
    }
    break;

  case '4':
    non_blocking_event_loop = !non_blocking_event_loop;
    if( non_blocking_event_loop && whichobj->prim->molec &&
        whichobj->prim->molec->num_frames >1 ){
      animating_molecule = 1;
    } else{
      animating_molecule = 0;
    }
    break;

  case '3':
    camera->foclength += 1;
    redo_projection = 1;
    break;

  case '.':
    camera->foclength -= 1;
    redo_projection = 1;
    break;

  case '9':
    camera->lf.z += 10;
    redo_projection = 1;
    break;
  case '6':
    camera->lf.z -= 10;
    redo_projection = 1;
    break;

  case '1':
#ifdef USE_LASSP_ROTATE
    if( whichobj && whichobj->prim && whichobj->prim->molec ){
      if( !rot_pipe ){
        sprintf(command_string,
                "/home/landrum/bin/aix3.x/rotate -initial %lf %lf %lf",
                whichobj->rot.z,whichobj->rot.y,whichobj->rot.x);
        rot_pipe = popen(command_string,"r");
        non_blocking_event_loop = 1;
        whichobj->prim->molec->rotate_tool_on = 1;
        whichobj->prim->molec->rot_matrix.matrix[0][0] = 1.0;
        whichobj->prim->molec->rot_matrix.matrix[1][1] = 1.0;
        whichobj->prim->molec->rot_matrix.matrix[2][2] = 1.0;
        whichobj->prim->molec->rot_matrix.matrix[3][3] = 1.0;
      } else{
        pclose(rot_pipe);
        rot_pipe = 0;
        whichobj->prim->molec->rotate_tool_on = 0;
        non_blocking_event_loop = 0;
      }
    }
#endif

    readcharparm("(V)rml or (R)ayshade (anything else to cancel)",instring);
    switch(instring[0]){
    case 'r':
    case 'R':
      dump_3D_objects(whichobj->prim,whichobj);
      break;
    case 'v':
    case 'V':
      dump_VRML(whichobj->prim,whichobj);
      break;
    default:
      fprintf(stderr,"invalid option: %c\n",instring[0]);
      break;
    }
    break;
  case '?':
    show_keys();
    break;
#ifdef SUPPORT_VOLUMES
  case 'v':
  case 'V':
    if( whichobj && whichobj->prim && whichobj->prim->MO_surf ){
      eval_MO_volume(whichobj->prim->MO_surf);
    }
    break;
#endif
  case '!':
    if(near_plane_clipping_on){
      printf("Near plane clipping off.\n");
      near_plane_clipping_on = 0;
    } else{
      printf("Near plane clipping on.\n");
      near_plane_clipping_on = 1;
    }
    break;

  case 'o':
  case 'O':
    if(outline_polyhed_on){
      printf("Polyhedron outlining off.\n");
      outline_polyhed_on = 0;
    } else{
      printf("Polyhedron outlining on.\n");
      outline_polyhed_on = 1;
    }
    break;

  case 'g':
  case 'G':
    if(dump_grids_on){
      printf("Grid dumping off.\n");
      dump_grids_on = 0;
    } else{
      printf("Grid dumping on.\n");
      dump_grids_on = 1;
    }
    break;

  case 'q':
  case 'Q': quit=1;break;
    }
  redrawgraph();
  }

#ifdef X_GRAPHICS

int click_x, click_y, end_x, end_y;
char drag_select=0;


/****************************************************************************
 *
 *                   Procedure do_button
 *
 * Arguments: event: an event structure
 * Returns: none
 *
 * Action: This routine processes the user's button clicks ;-)
 *
 ****************************************************************************/
void do_button( XEvent event )
{
  int xpos, ypos;

  xpos=event.xbutton.x;
  ypos=event.xbutton.y;
  click_x=event.xbutton.x;
  click_y=event.xbutton.y;

  /*
    printf( "Window: %d Button:  %d  (%d, %d)\n", (int)event.xany.window,
    event.xbutton.button, xpos, ypos);
    */
  switch( event.xbutton.button ){
  case 1:
    if( event.xany.window == gwin ){

      switch(mainmode){
#if 0
        whichobj->rot.x = 6.28*(g_xmax-xpos)/g_xmax;
        whichobj->rot.y = 6.28*(g_ymax-ypos)/g_ymax;
        redrawgraph();
#endif
      case ROT:
        click_x = xpos;
        click_y = ypos;
        break;
      case TRANS:
        click_x = xpos;
        click_y = ypos;
        break;
      case CENT:
        if( whichobj ){
          whichobj->cent.x = (float)(xpos);
          whichobj->cent.y = (float)(ypos);
          redrawgraph();
        }
        break;
      case CHOOSE:
        if( whichobj ){
          if( whichobj->prim->molec ){
            select_atom(whichobj->prim->molec,xpos,ypos);
            redraw();
          }
          else if( whichobj->prim->MO_surf && whichobj->prim->MO_surf->molec ){
            select_atom(whichobj->prim->MO_surf->molec,xpos,ypos);
            redraw();
          } else{
            select_object(xpos,ypos);
          }
        }
        break;
      }
    }
    else{
      find_button_win(button_wins,event.xany.window,xpos,ypos,1,
                      event.xbutton.state);
    }
    break;
  case 2:
    if( event.xany.window != gwin ){
      find_button_win(button_wins,event.xany.window,xpos,ypos,2,
                      event.xbutton.state);
    }
    redraw();
    break;
  case 3:
    if( event.xany.window == gwin ){
      if( !num_selected ){
        select_object(xpos,ypos);
      }
      else{
        if( whichobj && (whichobj->prim->molec ||
           (whichobj->prim->MO_surf && whichobj->prim->MO_surf->molec))){
          show_selected_data(num_selected,whichobj,xpos,ypos);
          unselect_all_atoms(num_selected,whichobj);
          num_selected = 0;
          redrawgraph();
        }
      }
    }
    else{
      find_button_win(button_wins,event.xany.window,xpos,ypos,3,
                      event.xbutton.state);
    }
    break;
  default:
    printf( "Foo Button %d %d \n",xpos,ypos);break;
  }
}
#endif

#ifdef X_GRAPHICS
/****************************************************************************
 *
 *                   Procedure do_events
 *
 * Arguments: none
 * Returns: none
 *
 * Action: This procedure implements the event loop for interactive use in an
 *  event driven environment.
 *
 ****************************************************************************/
void do_events(void)
{
  XEvent event;
  int xpos,ypos;
  int box_x,box_y,box_width,box_height;
  int num;

  g_draw_all_buttons(button_wins);
  XFlush(disp);
  quit = 0;
  drag_select = 0;

  /* This is the event loop */
  while( !quit ){
    if( non_blocking_event_loop ){
      /* check to see if there are any events waiting... */
      num = XPending(disp);

    }else num = 1;
    if(num){
      XNextEvent(disp, &event);
      switch(event.type){
      case KeyPress:
        do_keypress( event );
        break;

      case ButtonPress:
        do_button( event );
        break;

      case ButtonRelease:
        if( drag_select && event.xany.window==gwin &&
            event.xbutton.button == 1 &&
            mainmode == CHOOSE &&
            (whichobj->prim->molec || (whichobj->prim->MO_surf &&
                                       whichobj->prim->MO_surf->molec))){
          drag_select = 0;
          if( click_x < xpos ) box_x = click_x;
          else box_x = xpos;
          if( click_y < ypos ) box_y = click_y;
          else box_y = ypos;
          box_width = ABS(click_x-xpos);
          box_height = ABS(click_y-ypos);
          if( whichobj->prim->molec ){
            select_atoms_in_region(whichobj->prim->molec,box_x,box_y,
                                   box_x+box_width,box_y+box_height);
          } else {
            select_atoms_in_region(whichobj->prim->MO_surf->molec,box_x,box_y,
                                   box_x+box_width,box_y+box_height);
          }

          redrawgraph();
        }
        break;

      case MotionNotify:
        if( whichobj && event.xany.window==gwin &&
           event.xmotion.state & Button1Mask){
          num=0;
          /* make sure that we don't get trapped taking mouse motion events... */
          while(num<5){
            xpos = event.xmotion.x;
            ypos = event.xmotion.y;

            XPeekEvent(disp,&event);
            if(event.type == MotionNotify)XNextEvent(disp,&event);
            else break;
            num++;
          }

          switch(mainmode){
          case ROT:
            whichobj->rot.x -= (float)(ypos - click_y)/(100.0*M_PI);
            whichobj->rot.y += (float)(xpos - click_x)/(100.0*M_PI);
            click_x = xpos;
            click_y = ypos;
            redrawgraph();
            break;
          case CENT:
            /*            whichobj->cent.x = g_xmax/2 + (float)(xpos);
            whichobj->cent.y = (float)(ypos) - whichobj->trans.y;
            */
            click_x = xpos;
            click_y = ypos;
            redrawgraph();
            break;
          case TRANS:
            whichobj->trans.x += (float)(xpos - click_x)/10.0;
            whichobj->trans.y -= (float)(ypos - click_y)/10.0;
            click_x = xpos;
            click_y = ypos;
            redrawgraph();
            break;
          case CHOOSE:
            if( ABS(xpos - click_x) + ABS(ypos-click_y) > 10 ){
              drag_select = 1;
              redrawgraph();
              /***
                !!!memo to self...
                this is evil, please to be doing it properly RSN.
              ***/
              if( click_x < xpos ) box_x = click_x;
              else box_x = xpos;
              if( click_y < ypos ) box_y = click_y;
              else box_y = ypos;
              box_width = ABS(click_x-xpos);
              box_height = ABS(click_y-ypos);
              XSetLineAttributes( disp, graphgc, (short)0, LineSolid,
                                  CapRound, JoinRound);
              XDrawRectangle(disp,gwin,graphgc,box_x,box_y,
                             box_width,box_height);
            }

            break;

          }
        }
        break;


      case ConfigureNotify:
        /* must resize the window now.... */
        g_xmax = event.xconfigure.width;
        g_ymax = event.xconfigure.height;
        /* we need a new sized pixmap */
        XFreePixmap(disp,gpix);

        gpix=XCreatePixmap(disp,gwin,g_xmax,g_ymax,screen_depth);
        if(!gpix){
          fatal("Can't D_REALLOCate pixmap\n");
        }

        redo_projection = 1;
        camera->yaspect = g_ymax/g_xmax;
        redraw();
        break;

      case Expose:
        num=0;
        /* make sure that we don't get trapped taking exposure events... */
        while(XCheckMaskEvent(disp,ExposureMask,&event));
        redraw();
        break;

      default:
        break;
      }
    } else {
#ifdef USE_LASSP_ROTATE
      /********

        if the rotation window is open, then check to see
        if we need to read from that.

      ********/
      if( rot_pipe ){
        int nfds,readfds,writefds,exceptfds;
        char instring[80];
        struct timeval timeout;

        readfds = 1 << fileno(rot_pipe);
        nfds = fileno(rot_pipe) + 1;
        writefds = 0;
        exceptfds = 0;
        timeout.tv_sec = 0;
        timeout.tv_usec = 0;
        select(nfds,&readfds,&writefds,&exceptfds,&timeout);
        if( readfds ){
          if( first_rot_pipe_read ){
            first_rot_pipe_read = 0;
          } else{
            fprintf(stderr,"Got one!\n");
            if( whichobj && whichobj->prim && whichobj->prim->molec ){
              fgets(instring,80,rot_pipe);
              sscanf(instring,"%lf %lf %lf",
                     &(whichobj->prim->molec->rot_matrix.matrix[0][0]),
                     &(whichobj->prim->molec->rot_matrix.matrix[0][1]),
                     &(whichobj->prim->molec->rot_matrix.matrix[0][2]));
              fgets(instring,80,rot_pipe);
              sscanf(instring,"%lf %lf %lf",
                     &(whichobj->prim->molec->rot_matrix.matrix[1][0]),
                     &(whichobj->prim->molec->rot_matrix.matrix[1][1]),
                     &(whichobj->prim->molec->rot_matrix.matrix[1][2]));
              fgets(instring,80,rot_pipe);
              sscanf(instring,"%lf %lf %lf",
                     &(whichobj->prim->molec->rot_matrix.matrix[2][0]),
                     &(whichobj->prim->molec->rot_matrix.matrix[2][1]),
                     &(whichobj->prim->molec->rot_matrix.matrix[2][2]));
              fgets(instring,80,rot_pipe);
              redraw();
            } else{
              fgets(instring,80,rot_pipe);
              fprintf(stderr,"%s",instring);
              fgets(instring,80,rot_pipe);
              fprintf(stderr,"%s",instring);
              fgets(instring,80,rot_pipe);
              fprintf(stderr,"%s",instring);
              fgets(instring,80,rot_pipe);
              fprintf(stderr,"%s",instring);
            }
          }
        }
      }
#endif
      /******
         there are no events waiting, and we may need to be
         updating a surface
      *******/
      redrawgraph();
    }
  }
}
#endif
/****************************************************************************
 *
 *                   Procedure parse_commands
 *
 * Arguments: none
 * Returns: none
 *
 * Action: This is the parser for commands in a non-event driven (command driven)
 *   environment.
 *
 ****************************************************************************/
void parse_commands(void)
{
  int i,which_curve,which_style;
  graph_type *the_graph;
  prop_graph_type *prop_graph;
  prim_type *prim;
  char instring[80],tempstring[80];
  char *string_to_change[MAX_LINES];
  float *float_to_change;
  int *int_to_change;


  printf("\n\nViewkel: a program for displaying the results of extended\n");
  printf("      Hueckel calculations.\n\n");
  printf("Written by Greg Landrum\n");
  printf("    send comments/suggestions/bugs to landrum@xtended.tn.cornell.edu\n");
  printf("\n\nWelcome to the command line driven version of viewkel.\n");
  printf("  You can get a list of commands by typing 'help' or '?'\n\n");

  while(!quit){

    prop_graph = 0;
    if( whichobj && whichobj->prim){
      switch(whichobj->prim->which){
      case GRAPH:
        the_graph = whichobj->prim->graph;
        break;
      case PROP_GRAPH:
        the_graph =  whichobj->prim->prop_graph->the_data;
        prop_graph = whichobj->prim->prop_graph;
        break;
      case BAND_GRAPH:
        the_graph = whichobj->prim->band_graph->the_data;
        break;
      case WALSH_GRAPH:
        the_graph =  whichobj->prim->walsh_graph->the_data;
        break;
      }
    } else{
      the_graph = 0;
    }

    if( whichobj ) prim = whichobj->prim;

    printf("viewkel> ");
    fgets(instring,80,stdin);
    upcase(instring);

    if( !strncmp(instring,"QUIT",4) ){
      quit = 1;
    }
    else if(!strncmp(instring,"GRAP",4)){
      new_graph((char *)0);
      redraw();
      fgets(instring,80,stdin);
    }
    else if(!strncmp(instring,"BAND",4)){
      new_band_graph((char *)0);
      redraw();
      fgets(instring,80,stdin);
    }
    else if(!strncmp(instring,"PROP",4)){
      new_prop_graph((char *)0);
      redraw();
      fgets(instring,80,stdin);
    }
    else if(!strncmp(instring,"FMO",3)){
      new_FMO_diagram(NULL);
      redraw();
      fgets(instring,80,stdin);
    }
    else if(!strncmp(instring,"WAL",3)){
      new_walsh_graph(0);
      redraw();
      fgets(instring,80,stdin);
    }
    else if(!strncmp(instring,"RED",3)){
      redraw();
      fgets(instring,80,stdin);
    }
    else if(!strncmp(instring,"XTIC",4)){
      if(the_graph){
        the_graph->do_x_tics = !the_graph->do_x_tics;
        redraw();
        fgets(instring,80,stdin);
      }
    }
    else if(!strncmp(instring,"YTIC",4)){
      if(the_graph){
        the_graph->do_y_tics = !the_graph->do_y_tics;
        redraw();
        fgets(instring,80,stdin);
      }
    }
    else if(!strncmp(instring,"XLEG",4)){
      if(the_graph){
        string_to_change[0] = the_graph->xlegend;
        readstringparm("X legend",string_to_change);
        redraw();
        fgets(instring,80,stdin);
      }
    }
    else if(!strncmp(instring,"YLEG",4)){
      if(the_graph){
        string_to_change[0] = the_graph->ylegend;
        readstringparm("Y legend",string_to_change);
        redraw();
        fgets(instring,80,stdin);
      }
    }
    else if(!strncmp(instring,"XMAX",4)){
      if(the_graph){
        if( prop_graph ){
          float_to_change = &(prop_graph->max_x);
        } else{
          float_to_change = &(the_graph->max_x);
        }
        readfloatparm("x max",float_to_change);
        redraw();
        fgets(instring,80,stdin);
      }
    }
    else if(!strncmp(instring,"XMIN",4)){
      if(the_graph){
        if( prop_graph ){
          float_to_change = &(prop_graph->min_x);
        } else{
          float_to_change = &(the_graph->min_x);
        }
        readfloatparm("x min",float_to_change);
        redraw();
        fgets(instring,80,stdin);
      }
    }
    else if(!strncmp(instring,"YMAX",4)){
      if(the_graph){
        if( prop_graph ){
          float_to_change = &(prop_graph->max_y);
        } else{
          float_to_change = &(the_graph->max_y);
        }
        readfloatparm("y max",float_to_change);
        redraw();
        fgets(instring,80,stdin);
      }
      else if( prim->FMO_diagram ){
        float_to_change = &(prim->FMO_diagram->max_y);
        readfloatparm("y max",float_to_change);
        redraw();
        fgets(instring,80,stdin);
      }
    }
    else if(!strncmp(instring,"YMIN",4)){
      if(the_graph){
        if( prop_graph ){
          float_to_change = &(prop_graph->min_y);
        } else{
          float_to_change = &(the_graph->min_y);
        }
        readfloatparm("y min",float_to_change);
        redraw();
        fgets(instring,80,stdin);
      }
      else if( prim->FMO_diagram ){
        float_to_change = &(prim->FMO_diagram->min_y);
        readfloatparm("y min",float_to_change);
        redraw();
        fgets(instring,80,stdin);
      }

    }
    else if(!strncmp(instring,"XSCAL",5)){
      if( whichobj ){
        readfloatparm("X scaling",&(whichobj->scale.x));
      }
    }
    else if(!strncmp(instring,"YSCAL",5)){
      if( whichobj ){
        readfloatparm("Y scaling",&(whichobj->scale.y));
      }
    }
    else if(!strncmp(instring,"XTRANS",6)){
      if( whichobj ){
        readfloatparm("X translation",&(whichobj->trans.x));
      }
    }
    else if(!strncmp(instring,"YTRANS",6)){
      if( whichobj ){
        readfloatparm("Y translation",&(whichobj->trans.y));
      }
    }
    else if(!strncmp(instring,"YMAX",4)){
      if(the_graph){
        if( prop_graph ){
          float_to_change = &(prop_graph->max_y);
        } else{
          float_to_change = &(the_graph->max_y);
        }
        readfloatparm("y max",float_to_change);
        redraw();
        fgets(instring,80,stdin);
      }
    }
    else if(!strncmp(instring,"TITL",4)){
      if(the_graph){
        for(i=0;i<NUM_TITLE_LINES;i++){
          string_to_change[i] = the_graph->title[i];
        }
        readmultistringparm("Title",NUM_TITLE_LINES,string_to_change);
        the_graph->do_title = 1;
        redraw();
        fgets(instring,80,stdin);
      }
    }
    else if(!strncmp(instring,"PRIN",4)){
      do_ps_output();
    }
    else if(!strncmp(instring,"CURVE",5)){
      if( the_graph ){
        /* figure out which curve to use */
        sscanf(instring,"%s %d",tempstring,&which_curve);

        if( which_curve <= the_graph->num_curves ){
          the_graph->curves_to_display[which_curve-1] =
            !(the_graph->curves_to_display[which_curve-1]);
        }else{
          fprintf(stderr,"There are only %d curves in the data set!\n",
                  the_graph->num_curves);
        }
      }
    }
    else if(!strncmp(instring,"LINESTY",7)){
      if( the_graph ){
        /* figure out which curve and which style to use */
        sscanf(instring,"%s %d %d",tempstring,&which_curve,&which_style);

        if( which_curve <= the_graph->num_curves ){
          the_graph->styles[which_curve-1] = which_style;
        }else{
          fprintf(stderr,"There are only %d curves in the data set!\n",
                  the_graph->num_curves);
        }
      }
    }

    else if(!strncmp(instring,"INTEGSTY",8)){
      if( whichobj->prim->which != PROP_GRAPH ){
        fprintf(stderr,"This type of graph doesn't have an integration.\n");
      }
      else{
        the_graph = prop_graph->the_integration;

        if( the_graph ){
          /* figure out which curve to use */
          sscanf(instring,"%s %d %d",tempstring,&which_curve,&which_style);

          if( which_curve <= the_graph->num_curves ){
            the_graph->styles[which_curve-1] = which_style;
          }else{
            fprintf(stderr,"There are only %d curves in the data set!\n",
                    the_graph->num_curves);
          }
        }
      }
    }

    else if(!strncmp(instring,"INTEGSC",7)){
      if( whichobj->prim->which != PROP_GRAPH ){
        fprintf(stderr,"This type of graph doesn't have an integration.\n");
      }
      else{
        if(prop_graph->the_integration){
          prop_graph->integs_for_tics = !(prop_graph->integs_for_tics);
        }
      }
    }
    else if(!strncmp(instring,"INTEG",5)){
      if( whichobj->prim->which != PROP_GRAPH ){
        fprintf(stderr,"This type of graph doesn't have an integration.\n");
      }
      else{
        the_graph = prop_graph->the_integration;

        if( the_graph ){
          /* figure out which curve to use */
          sscanf(instring,"%s %d",tempstring,&which_curve);

          if( which_curve <= the_graph->num_curves ){
            the_graph->curves_to_display[which_curve-1] =
              !(the_graph->curves_to_display[which_curve-1]);
          }else{
            fprintf(stderr,"There are only %d curves in the data set!\n",
                    the_graph->num_curves);
          }
        }
      }
    }
    else if(!strncmp(instring,"KILL",4)){
      if(whichobj){
        free_obj(whichobj);
        D_FREE(whichobj);
        whichobj = 0;
      }
      if( head->next ){
        head = head->next;
      }
      else{
        head = 0;
      }
    }
    else if(!strncmp(instring,"PURGE",5)){
      clearall();
    }
    else if(!strncmp(instring,"FILL",4)){
      if( fill_projections ){
        fill_projections = 0;
        printf("Filling of DOS Projections is off.\n");
      }
      else{
        fill_projections = 1;
        printf("Filling of DOS Projections is on.\n");
      }
    }
    else if(!strncmp(instring,"FERMI",5)){
      if(whichobj && whichobj->prim ){
        switch(whichobj->prim->which){
        case BAND_GRAPH:
          whichobj->prim->band_graph->show_fermi =
            !(whichobj->prim->band_graph->show_fermi);
          redraw();
          break;
        case PROP_GRAPH:
          whichobj->prim->prop_graph->show_fermi =
            !(whichobj->prim->prop_graph->show_fermi);
          redraw();
          break;
        default:
          fprintf(stderr,"This kind of graph doesn't have a Fermi level.\n");
          break;
        }
      }
    }
    else if(!strncmp(instring,"CHANGE FER",10)){
      if(whichobj && whichobj->prim ){
        switch(whichobj->prim->which){
        case BAND_GRAPH:
          readfloatparm("Fermi Level",&(whichobj->prim->band_graph->Fermi_E));
          whichobj->prim->band_graph->show_fermi = 1;
          redraw();
          break;
        case PROP_GRAPH:
          readfloatparm("Fermi Level",&(whichobj->prim->prop_graph->Fermi_E));
          whichobj->prim->prop_graph->show_fermi = 1;
          redraw();
          break;
        default:
          fprintf(stderr,"This kind of graph doesn't have a Fermi level.\n");
          break;
        }
      }
    }
    else if(!strncmp(instring,"SHOW",4)){
      if(whichobj && prim ){
        if( prim->which == FMO_DIAGRAM ){
          if( strstr(instring,"HOMO") ){
            prim->FMO_diagram->electron_filling_mode = FMO_FILL_HOMO;
          }else if( strstr(instring,"ALL") ){
            prim->FMO_diagram->electron_filling_mode = FMO_FILL_ALL;
          }else if( strstr(instring,"NONE") ){
            prim->FMO_diagram->electron_filling_mode = FMO_FILL_NONE;
          }
        }
      }
    }
    else if(strstr(instring,"FRAG")){
      if(whichobj && prim ){
        if( prim->which == FMO_DIAGRAM ){
          if( strstr(instring,"LEFT") ){
            int_to_change = &(prim->FMO_diagram->left_fragment);
            readintparm("left fragment",int_to_change);
          } else if(strstr(instring,"RIGHT")){
            int_to_change = &(prim->FMO_diagram->right_fragment);
            readintparm("right fragment",int_to_change);
          }
        }
      }
    }
    else if(!strncmp(instring,"HELP",4)){
      show_help();
    }
    else if(!strncmp(instring,"?",1)){
      show_help();
    }
    else if(instring[0] == '\n');
    else fprintf(stderr,"Invalid command: %s\n",instring);
  }
}
