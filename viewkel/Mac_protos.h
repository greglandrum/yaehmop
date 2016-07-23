/*******************************************************

Copyright (C) 1996 Greg Landrum
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
/***************

  prototypes for the Mac specific functions

  Since we're guaranteed to be on a Mac here,
    we don't have to be as general as before.

  Created by greg Landrum  Jan 16, 1996

****************/

#ifndef MAC_PROTOS
#define MAC_PROTOS
extern OSErr InitViewkelApp(void);
extern OSErr InitViewkelWindow(WindowPtr, void *);
extern OSErr SetupViewkelMenus(void);
extern void UpdateViewkelMenus(WindowPtr);
extern Boolean ViewkelMenuDispatch(WindowPtr, long);
extern void action_menu_hook(short);
extern void read_data_submenu_hook(short);
extern void mode_menu_hook(short);
extern OSErr ClickInViewkel(WindowPtr macWindow, void *data);
extern OSErr UpdateViewkel(WindowPtr macWindow, void *data);
extern OSErr BackgroundViewkelClick(WindowPtr macWindow, void *data);
extern OSErr ViewkelResize(WindowPtr macWindow, void *data);
extern OSErr SizeViewkelWindow(WindowPtr macWindow, void *data);
extern FILE *choose_mac_file(char *file_name, char action);
extern void Mac_initgraphics();
extern void Mac_ClearScreen();
extern void Mac_CopyGWorld();
extern void SetViewkelWindowInit(EasyWindowPtr easyWindow, OSType theFile);

#endif
