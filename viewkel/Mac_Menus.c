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

#include "viewkel.h"
#include "EasyApp.h"
#include "Mac_defines.h"
#include "Mac_protos.h"

void UpdateViewkelObjectMenus(WindowPtr macWindow);
Boolean ViewkelMenuDispatch(WindowPtr macWindow, long menuChoice);
void action_menu_hook(short menuItem);
void read_data_submenu_hook(short menuItem);
void mode_menu_hook(short menuItem);
Boolean ViewkelObjectMenuDispatch(short menuID, short menuItem);



/****************************************************************************
 *
 *                   Function SetupViewkelMenus
 *
 * Arguments: none
 *
 * Returns: OSErr
 *
 * Action: Does the initialization for the menus that will
 *      appear when viewkel is started up.
 *
 ****************************************************************************/
OSErr SetupViewkelMenus(void)
{
  OSErr error = noErr;

  /* get space for the menus and make sure we got it */
  Mac_globals.modeMenu = GetMenu(MODE_MENU_ID);
  Mac_globals.actionMenu = GetMenu(ACTION_MENU_ID);
  Mac_globals.read_dataMenu = GetMenu(READ_DATA_SUBMENU_ID);
  if(!Mac_globals.modeMenu || !Mac_globals.actionMenu ||
     !Mac_globals.read_dataMenu ){
    error = kNoResourceError;
    return error;
  }

  /* insert the menus into the menu bar now */
  InsertMenu(Mac_globals.actionMenu, 0 );
  InsertMenu(Mac_globals.read_dataMenu, -1 );
  InsertMenu(Mac_globals.modeMenu, 0 );

  return error;
}

/****************************************************************************
 *
 *                   Procedure UpdateViewkelMenus
 *
 * Arguments: macWindow: WindowPtr
 *
 * Returns: void
 *
 * Action: Deals with all the stuff needed to updated the menus
 *     so that they will be displayed properly.  This includes
 *     setting marks next to active items, etc.
 *
 ****************************************************************************/
void UpdateViewkelMenus(WindowPtr macWindow)
{
  int i;
  MenuHandle theModeMenu;

        EasyDisableMenu(g.fileMenu);
        EnableItem(g.fileMenu,kQuitItem);

  theModeMenu=GetMenuHandle(MODE_MENU_ID);

  /* turn off all the checks */
  for(i=MODE_NONE_ITEM_ID;i<MODE_LAST_ITEM_ID;i++){
    CheckItem(theModeMenu,i,FALSE);
  }

  /* update the mode menu to reflect the currently active mode */
  switch(mainmode){
  case NONE:
    CheckItem(theModeMenu,MODE_NONE_ITEM_ID,TRUE);break;
  case ROT:
    CheckItem(theModeMenu,MODE_ROTATE_ITEM_ID,TRUE);break;
  case SCALE:
    CheckItem(theModeMenu,MODE_SCALE_ITEM_ID,TRUE);break;
  case TRANS:
    CheckItem(theModeMenu,MODE_TRANSLATE_ITEM_ID,TRUE);break;
  case CENT:
    CheckItem(theModeMenu,MODE_CENTER_ITEM_ID,TRUE);break;
  }

  /* check (or not) the fill proj's item */
  CheckItem(theModeMenu,MODE_FILL_ITEM_ID,(Boolean)fill_projections);

        UpdateViewkelObjectMenus(macWindow);
}







/****************************************************************************
 *
 *                   Procedure UpdateViewkelObjectMenus
 *
 * Arguments: macWindow: WindowPtr
 *
 * Returns: void
 *
 * Action: Deals with all the stuff needed to update the object menus
 *     so that they will be displayed properly.  This includes
 *     setting marks next to active items, etc.
 *
 ****************************************************************************/
void UpdateViewkelObjectMenus(WindowPtr macWindow)
{
        char pstring[255];
        button_win_type *button_win;
        button_type *button;

        /* loop over all of the button_wins and their children */
        button_win = button_wins;
        while(button_win){
                button = button_win->button_list;
                while(button){
                        /* start by setting the text for this button */
                        strcpy(pstring,button->string);
                        CtoPstr(pstring);
                        SetItem(button_win->theMenu,button->itemID,(unsigned char *)pstring);

                        /* now make sure that the buttons are properly displayed */
                        switch(button->type){
                        case TOGGLE:
                                CheckItem(button_win->theMenu,button->itemID,
                                                                        (Boolean)(*button->guts.boolean));
                                break;
                        case CURVETOGGLE:
                                CheckItem(button_win->theMenu,button->itemID,
                                                                        (Boolean)(*button->guts.boolean));
                                break;

                        }
                        button = button->next;
                }

                /* update the children if there are any */
                if( button_win->child ){
                                button = button_win->child->button_list;
                while(button){
                        /* start by setting the text for this button */
                        strcpy(pstring,button->string);
                        CtoPstr(pstring);
                        SetItem(button_win->child->theMenu,button->itemID,(unsigned char *)pstring);

                        /* now make sure that the buttons are properly displayed */
                        switch(button->type){
                        case TOGGLE:
                                CheckItem(button_win->child->theMenu,button->itemID,
                                                                        (Boolean)(*button->guts.boolean));
                                break;
                        case CURVETOGGLE:
                                CheckItem(button_win->child->theMenu,button->itemID,
                                                                        (Boolean)(*button->guts.boolean));
                                break;

                        }
                        button = button->next;
                }
                }

                button_win = button_win->next;
        }
}


Boolean ViewkelMenuDispatch(WindowPtr macWindow, long menuChoice)
{
  short menuID,menuItem;
  Boolean handled=false;


    /* EasyApp didn't touch it, now we deal */
    menuID = HiWord(menuChoice);
    menuItem = LoWord(menuChoice);

    switch(menuID){
    case READ_DATA_SUBMENU_ID:
      read_data_submenu_hook(menuItem);
      handled = true;
      break;
          case ACTION_MENU_ID:
            action_menu_hook(menuItem);
            break;
    case MODE_MENU_ID:
      mode_menu_hook(menuItem);
      handled = true;
      break;
    default:
      handled = ViewkelObjectMenuDispatch(menuID,menuItem);
      break;
    }
         return handled;
}

/***************
  action_menu_hook
  takes a menu choice as an argument and calls
  the appropriate function to deal with what should be read in.
 ****************/
void action_menu_hook(short menuItem)
{

        switch(menuItem){
        case PURGE_ITEM_ID:
         clearall();
         break;
        case PRINT_ITEM_ID:
                do_ps_output();
                break;
        case CLEAR_LABEL_ITEM_ID:
                clear_labels();
                break;
        default:
                break;
        }

}

/***************
  read_data_submenu_hook
  takes a menu choice as an argument and calls
  the appropriate function to deal with what should be read in.
 ****************/
void read_data_submenu_hook(short menuItem)
{

  switch(menuItem){
  case READ_PROPS_ITEM_ID:
    new_prop_graph(NULL);
    break;
  case READ_BANDS_ITEM_ID:
    new_band_graph(NULL);
    break;
  case READ_CONT_ITEM_ID:
    new_cont_plot(NULL);
    break;
  case READ_FMO_ITEM_ID:
          new_FMO_diagram(NULL);
          break;
  case READ_WALSH_ITEM_ID:
          new_walsh_graph(0);
          break;
  case READ_MOLEC_ITEM_ID:
          new_molecule(NULL);
          break;
        case READ_MO_ITEM_ID:
                new_MO_surface(NULL);
  default:
    break;
  }
}


/***************
  mode_menu_hook
  takes a menu choice as an argument and sets the
  appropriate mode in the mode menu
****************/
void mode_menu_hook(short menuItem)
{

  switch(menuItem){
  case MODE_NONE_ITEM_ID:
    mainmode = NONE;
    break;
  case MODE_TRANSLATE_ITEM_ID:
    mainmode = TRANS;
    break;
  case MODE_ROTATE_ITEM_ID:
    mainmode = ROT;
    break;
  case MODE_CENTER_ITEM_ID:
    mainmode = CENT;
    break;
  case MODE_SCALE_ITEM_ID:
    mainmode = SCALE;
    break;
        case MODE_FILL_ITEM_ID:
                fill_projections = !fill_projections;
                break;
  default:
    break;
  }
}


/*******************************************

        ViewkelObjectMenuDispatch

        Arguments: menuID, menuItem: shorts
        Returns: Boolean

        This tracks through the menus and deals with calling the
                appropriate function if the user has selected
                an item in an ObjectMenu

********************************************/
Boolean ViewkelObjectMenuDispatch(short menuID, short menuItem)
{
  Boolean handled=0;
        button_win_type *button_win;
        button_type *button;
        int *mode_ptr;

        /* first find the appropriate button_win */
        button_win = button_wins;
        while(button_win && !handled){
                if(button_win->menuID == menuID){
                        /* find the item */
                        button = button_win->button_list;
                        while(button && !handled){
                                if(button->itemID == menuItem) handled = 1;
                                else button = button->next;
                        }
                }
                if(button_win->child && button_win->child->menuID == menuID){
                        /* find the item */
                        button = button_win->child->button_list;
                        while(button && !handled){
                                if(button->itemID == menuItem) handled = 1;
                                else button = button->next;
                        }
                }

                if( !handled ){
                        button_win = button_win->next;
                }
        }

        if( handled ){
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
        }

        return handled;
}

