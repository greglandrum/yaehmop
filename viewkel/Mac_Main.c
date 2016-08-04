/*****************************************************************/
/*
  EasyApp

  The Macintosh Application Shell

  James E. Trudeau, Nebula, Inc.

  With thanks to Eric Shapiro, Bill Worzel, Brad Mohr, Bill
  Hofmann, and many more folks, including the fine people at
  Apple's Developer Technical Support and Developer University.

  Version 1.0  4/3/95

  This was hacked up for inclusion in viewkel by greg Landrum
  1/11/96

  */
/*****************************************************************/

#include        "EasyApp.h"
#include "viewkel.h"
#include "Mac_defines.h"
#include "Mac_protos.h"

/* define a single instance of our global variables */
EasyGlobal                        g;

/****************************************************************
  Mac_initgraphics
  ****************************************************************/
/*
  To achieve maximum flexibility using the EasyApp shell, most of the
  time the Mac_initgraphics() function calls "hook" routines to set up the app.
  See hooks.c for an explanation of what hooks are and how they're used.

  The tasks performed by Mac_initgraphics() include
  ¥        set up application memory
  ¥        initialize the Mac Toolbox
  ¥        examines the operating environement for features
  ¥        sets up cursors for use in the application
  ¥        sets up the menu bar
  ¥        installs operating system callback routines (handlers)
  ¥        calls the event loop

  Requires: a Macintosh computer
  Receives: nothing
  Changes:  who knows? You might change the world.
  Returns:  nothing
  */
void Mac_initgraphics (void)
{
  OSErr                error;


  /*         Set up application heap (memory) and initialize the Toolbox.         */
  MaxApplZone();                                                /* take all the memory we can */
  EasyGetMoreMasters(kNumberMasters);        /* allocate master pointer blocks */
  InitToolboxHook();                                        /* initialize the Toolbox */


  /*
    Set up the error dialog first, so we can display startup errors.
    This also involves some application global initializing.
    */
  error = EasyInitErrorDialog();
  if (error)
      {
        fatal("Can't open error Dialog\n");
        return;        /* g.error dialog set to NIL, we can't display error */
      }

  /*
    See what we've got in the operating environment, sets environment globals.
    Any error from here to the event loop is fatal and will cause the
    application to quit after displaying an error notice to the user.
    */

  error = CheckEnvironmentHook();
  if (error)
      {
        ShowErrorHook(error);
        return;
      }

  /* set up the rest of the application globals */
  error = InitAppHook();
  if (error)
      {
        ShowErrorHook(error);
        if(error != kNoPrinterError) /* no printer not fatal */
            {
              return;
            }
      }

  /* Set up any cursors you want to use in the application. */
  error = LoadCursorsHook();
  if (error)
      {
        ShowErrorHook(error);
        return;
      }

  /* set up the menus for the application */
  /*
    Note: instead of calling a hook directly, we call an EasyApp routine
    first. EasySetupMenus() does the Apple, File, and Edit menus, then
    calls SetupAppMenusHook() so your app can set up its own menus.
    */
  error = EasySetupMenus();
  if (error)
      {
        ShowErrorHook(error);
        return;
      }

  /* install Apple event handlers */
  error = InstallAEHandlersHook();
  if (error)
      {
        ShowErrorHook(error);
        return;
      }


#if 0
  /********

    open the graphics window and do the necessary initialization.

  *********/
  main_window = GetNewCWindow(MAIN_WINDOW_ID,NULL,(WindowPtr)-1);
  if( !main_window ){
    fatal("Can't open CWindow");
  }

#endif


  /* let's do it */
  EasyEventLoop();        /* where all the real fun happens */

  /* let's shut it down */
  error = ShutdownHook();

  if (error)
      {        /* display error dialog with correct message before we die or quit */
        ShowErrorHook(error);
      }

  /* we're all done */
}

/***********************************
  InitViewkelApp
  Handles initialization behaviour that differs from the EasyApp
  defaults
  ************************************/
OSErr InitViewkelApp(void)
{



}
