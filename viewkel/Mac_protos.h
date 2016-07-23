/***************
  Copyright 1999, Greg Landrum

  prototypes for the Mac specific functions
  
  Since we're guaranteed to be on a Mac here,
    we don't have to be as general as before.
    
  Created by greg Landrum  Jan 16, 1996
    
****************/

#ifndef MAC_PROTOS
#define MAC_PROTOS
extern OSErr InitViewkelApp(void);
extern OSErr InitViewkelWindow(WindowPtr,void *);
extern OSErr SetupViewkelMenus (void);
extern void UpdateViewkelMenus(WindowPtr);
extern Boolean ViewkelMenuDispatch(WindowPtr,long);
extern void action_menu_hook(short);
extern void read_data_submenu_hook(short);
extern void mode_menu_hook(short);
extern OSErr ClickInViewkel (WindowPtr macWindow, void *data);
extern OSErr UpdateViewkel (WindowPtr macWindow, void *data);
extern OSErr BackgroundViewkelClick (WindowPtr macWindow, void *data);
extern OSErr ViewkelResize (WindowPtr macWindow, void *data);
extern OSErr SizeViewkelWindow (WindowPtr macWindow, void *data);
extern FILE *choose_mac_file(char *file_name,char action);
extern void Mac_initgraphics();
extern void Mac_ClearScreen();
extern void Mac_CopyGWorld();
extern void SetViewkelWindowInit (EasyWindowPtr easyWindow, OSType theFile);

#endif
