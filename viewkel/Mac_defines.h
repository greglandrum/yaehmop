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

#ifndef _MAC_DEFINES_
#define _MAC_DEFINES_

#include <QDOffscreen.h>
#include "EasyApp.h"

/**********

 start off with resources and menu ID's

**********/

/* controls */
#define VSCROLL_CNTRL_ID 1001
#define HSCROLL_CNTRL_ID 1002

/* things on the menu bar and their menu items to make things tidy */
#define ACTION_MENU_ID 1003
#define READ_DATA_ITEM_ID 1
#define PURGE_ITEM_ID 2
#define PRINT_ITEM_ID 3
#define CLEAR_LABEL_ITEM_ID 4

#define MODE_MENU_ID 1004
#define MODE_NONE_ITEM_ID 1
#define MODE_ROTATE_ITEM_ID 2
#define MODE_SCALE_ITEM_ID 3
#define MODE_TRANSLATE_ITEM_ID 4
#define MODE_CENTER_ITEM_ID 5
/************

 this should be the highest numbered item_id for the possible
 main modes.  It corresponds to the separator underneath
 the main mode listing

 ***********/
#define MODE_LAST_ITEM_ID 6
#define MODE_FILL_ITEM_ID 7

/* submenus and their items */
#define READ_DATA_SUBMENU_ID 128
#define READ_PROPS_ITEM_ID 1
#define READ_BANDS_ITEM_ID 2
#define READ_CONT_ITEM_ID 3
/* separator */
#define READ_FMO_ITEM_ID 5
#define READ_WALSH_ITEM_ID 6
/* separator */
#define READ_MOLEC_ITEM_ID 8
#define READ_MO_ITEM_ID 9
/* separator */
#define READ_ADF_MO_ITEM_ID 11
#define READ_ADF_FREQ_ITEM_ID 12

/* the main window resource */
#define MAIN_WINDOW_ID 1001

/* dimensions of the full graphics map (the off screen graphics world) */
#define MAC_DOCUMENT_SIZE 500

#define APPCREATOR = 'Vwkl'
#define APPTYPE = 'Vwkl'

/****************
        To make everything easier and be consistant with the
        EasyApp stuff, store all of my resources and everything
        in a global struct
*****************/
typedef struct {
  MenuHandle actionMenu;
  MenuHandle modeMenu;
  MenuHandle read_dataMenu;
  ControlHandle vScrollBar, hScrollBar;
  short fontSize;

} Mac_global_struct;

extern Mac_global_struct Mac_globals;
extern EasyGlobal g;
