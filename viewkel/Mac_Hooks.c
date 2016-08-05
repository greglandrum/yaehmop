/*****************************************************************/
/*
  EasyApp

  The Macintosh Application Shell

  James E. Trudeau, Nebula, Inc.

  With thanks to Eric Shapiro, Bill Worzel, Brad Mohr, Bill
  Hoffman, and many more folks, including the fine people at
  Apple's Developer Technical Support and Developer University.

  Version 1.0 b2 3/8/95
  */
/*****************************************************************/

/*
  This file contains all the hooks called by the EasyApp shell.
  You modify these routines to call your own application routines
  to implement application-specific behavior.

  For the most part, each hook routine calls a corresponding EasyApp
  routine to get work done. For example, InitAppHook() calls
  EasyInitApp().

  In general there are three or four ways you can change the
  behavior of the EasyApp shell.

  First, you can replace the call currently in the hook with
  a call to your own routine. You are then responsible for all of
  the functionality implemented in the Easy routine. In other
  words, you take over completely.

  Second, in some cases you can leave the call to the EasyApp
  routine, and then follow or precede it with a call to your routine.
  In this way you do not replace default behavior, you augment it.

  Third, many of EasyApp's default routines send messages to windows.
  You can modify the window behavior by installing different functions
  in the window. For example, when it is appropriate to activate a
  window, EasyApp's default routines send a DoActivate message to the
  window. You can install your own DoActivate behavior in your own
  windows so they behave in a unique way.

  You are also free to modify the EasyApp routines themselves if you
  wish, but we don't recommend doing so. You'd be better off making
  a copy of the routine, renaming and making desired changes, and then
  calling the new routine.

  Applications can be sensitive things. If you start making changes in
  the core behavior of EasyApp, you may break something.

  Besides, as you can see, EasyApp is plenty flexible. You should be
  able to get it to do almost anything you want without having to
  change it.
  */

#include        "EasyApp.h"
#include "viewkel.h"
#include "Mac_defines.h"
#include "Mac_protos.h"



Mac_global_struct Mac_globals;

/* STARTUP & SHUTDOWN HOOKS */


/*****************************************************************/
/* InitToolboxHook */
/*****************************************************************/
/*
  EasyApp calls this function to initialize the Macintosh Toolbox.

  The default behavior should be sufficient for most applications.
  However, as new technology appears on the horizon, or if you want
  to use some other parts of the Toolbox, you may need to modify
  the behavior. You can do that easily by calling your own routine
  after the existing call to EasyInitToolbox()

  Requires: nothing
  Receives: nothing
  Changes:  initializes the Toolbox
  Returns:  nothing
  */

void        InitToolboxHook (void)
{
  EasyInitToolbox();

  /* call your own routine for additional initializing */


}



/*****************************************************************/
/* CheckEnvironmentHook */
/*****************************************************************/
/*
  EasyApp calls this function at launch after setting up the error
  dialog. We do that first in case an error happens here.

  The default behavior examines the operating environment to determine
  what system-level services are available.

  There are two kinds of things to look for: what you need, and what
  makes a difference. If something you need is missing, your application
  can't run. For example, EasyApp requires System 7. If System 7 is not
  available, EasyApp quits gracefully.

  Things that make a difference should be handled differently. For
  example, you may want to know whether color is available. If it
  is, you can do some things in color, if not you can do them in
  black and white.

  Once again, the best way to modify EasyApp is to augment it - call
  your own routine after calling EasyCheckEnvironment().If you are
  missing something you need, return an error. Otherwise, return no
  error. Any error is assumed to be fatal.

  Requires: nothing
  Receives: nothing
  Changes:  various flags in application global variables
  Returns:  error value, 0 => no error, assumes any error is fatal
  */

OSErr        CheckEnvironmentHook (void)
{
  OSErr        error = noErr;

  error = EasyCheckEnvironment();        /* see what's in the environment */
  if (error){
    return error;        /* any error is a fatal error         */
  }


  /* do your own checks here */


  return error;
}


/*****************************************************************/
/* InitAppHook */
/*****************************************************************/
/*
  EasyApp calls this routine at launch time, after it has set up
  memory, initialized the Toolbox, and examined the environment.
  So your app will know any environment-dependent facts you gleaned
  from examining the operating environment, like if you're running
  on a color-capable machine, what system version you're using,
  and any other details you want to look for.

  The default behavior sets up global values and a print record.
  See EasyInitApp() for how it's done. The best way to modify
  this behavior would be to augment it - call your own routine
  after calling EasyInitApp().

  Requires: Toolbox initialized
  Receives: nothing
  Changes:  nothing
  Returns:  error value, any error at this point is fatal
  */

OSErr        InitAppHook (void)
{
  OSErr        error = noErr;

  error = EasyInitApp();
  if(error){
    return error;        /* any error is fatal */
  }

  /* call your own routine here */
  InitViewkelApp();

  return error;
}


/*****************************************************************/
/* LoadCursorsHook */
/*****************************************************************/
/*
  EasyApp calls this hook after initializing the app. This
  hook allows you to set up your own cursors for application use.

  The default behavior sets up the arrow and watch cursors in the
  global variables.

  Once again, the best way to modify EasyApp is to augment it - call
  your own routine after calling EasyLoadCursors(). Any error
  returned is assumed to be fatal and causes the app to quit.

  Requires: nothing
  Receives: nothing
  Changes:  may change app globals pointing to cursors
  Returns:  error value 0 => no error, any error fatal
  */

OSErr        LoadCursorsHook (void)
{
  OSErr        error = noErr;

  error = EasyLoadCursors();
  if(error)
      {
        return error;        /* any error is fatal */
      }

  /* set up your own cursors here however you want */


  return error;
}



/*****************************************************************/
/* SetupAppMenusHook */
/*****************************************************************/

/*
  EasyApp calls this function at application startup to allow your
  app to add menus to the menu bar. EasyApp takes care of the Apple,
  File, and Edit menus before coming here. You can safely add any
  application-specific menus at this time.

  If you wish to modify the Apple, File, or Edit menus, just change
  the MENU resources. EasyApp handles installing them in the menu bar.

  EasyApp draws the menu bar after you're done.

  Requires: nothing
  Receives: nothing
  Changes:  nothing (app routine adds menus to menu bar)
  Returns:  error value, 0 => OK, any error fatal
  */

OSErr        SetupAppMenusHook (void)
{
  OSErr        error = noErr;

  /* call our menu setup function */
  error = SetupViewkelMenus();
  return error;
}


/*****************************************************************/
/* InstallAEHandlersHook */
/*****************************************************************/
/*
  EasyApp calls this routine at startup while the application is
  launching. It should install any handlers for apple events your
  application recognizes.

  The default behavior - EasyInstallAEHandlers() - installs the
  default handlers for the four required Apple events: open app,
  open doc, print doc, and quit. It also installs a generic
  Apple event handler that takes care of any unrecognized Apple
  event by displaying a dialog indicating that an unknown event
  has been received, and giving the type code for that event.

  You can either replace or augment this call. This is the last
  startup hook called before the application enters the main
  event loop.

  Requires: System 7, EasyApp tests before this routine is called
  Receives: nothing
  Changes:  calls routine to install handlers
  Returns:  error value
*/

OSErr InstallAEHandlersHook (void)
{
  return EasyInstallAEHandlers();
}


/*****************************************************************/
/* ShutdownHook */
/*****************************************************************/
/*
  EasyApp calls this function after the user quits the application.
  At the time this function is called, all windows have been closed and
  saved, and the application has left the main event loop. At this
  point the application will shut down, whatever happens.

  Typically you would release memory or perform any other cleanup
  required to gracefully terminate the application. While disposing
  of memory is not technically required, it's the right thing to do.

  EasyApp will display any error before quitting, but you cannot
  abort the quit process at this time. Anything that might happen
  while quitting that should abort the process should be handled
  by the QuitHook() function.

  Requires: nothing
  Receives: nothing
  Changes:  nothing
  Returns:  OSErr value, application will quit regardless, but you
  can display an error message.
  */

OSErr ShutdownHook (void)
{
  return EasyShutdown();
}




/*****************************************************************/
/* ShowErrorHook */
/*****************************************************************/
/*
  EasyApp calls this function when a reportable error occurs. The
  application may not be in the foreground at the time, so you still
  must exercise care before displaying an error dialog.

  You may replace or augment this with your own error management
  and reporting mechanism.

  Requires: nothing
  Receives: a short value, the error number
  Changes:  nothing directly, may result in error dialog appearing
  Returns:  nothing
  */

void ShowErrorHook (short error)
{

  EasyShowError(error);
}



/* EVENT HOOKS */


/*****************************************************************/
/* CheckMemoryHook */
/*****************************************************************/
/*
  EasyApp calls this function every time through the event loop.
  This gives you the opportunity to perform some simple memory
  management strategies. EasyApp uses an emergency memory reserve.
  The reserve is allocated at startup and preserved in a handle
  stored in g.memoryReserve.

  The default size of the reserve is 40K, set in the constant
  kMemoryReserveSize declared in userDeclare.h. You can set the
  size to anything you want by changing that constant.

  EasyCheckMemoryReserve() does two things. If there is a memory
  reserve, it looks to see if there is a minimum amount of memory
  available. If there is not, it D_FREEs the memory reserve and
  warns the user. If there is no memory reserve (because it has
  already been released), the routine tries to recreate the reserve.

  You may wish to implement a more comprehensive or flexible
  memory management strategy. If you do, call your own routine
  here.

  Requires: nothing
  Receives: pointer to event record
  Changes:  nothing
  Returns:  nothing
  */

void CheckMemoryHook (void)
{
  EasyCheckMemoryReserve();
}




/*****************************************************************/
/* PreprocessEventHook */
/*****************************************************************/
/*
  EasyApp calls this routine every time it gets an event. This allows
  you to pass the event to managers that might need to handle the
  event before you do (like QuickTime or PowerTalk). You can
  do any other app-specific event processing you want at this time,
  although typically you wouldn't have to do any.

  Make sure you leave the contents of the event record unchanged!

  If you have completely handled the event, return the value true.
  If you do, EasyApp will not process the event further.

  Requires: nothing
  Receives: pointer to event record
  Changes:  nothing directly
  Returns:  Boolean indication whether event has been handled
  */

Boolean        PreprocessEventHook (EventRecord *theEvent)
{
  Boolean handled = false;        /* assume not handled */


  return handled;
}


/*****************************************************************/
/* ModelessDialogHook */
/*****************************************************************/
/*
  EasyApp calls this function from inside the main event handler
  if an event is an event inside a movable dialog. Movable dialogs
  may be modal or modeless. The default behavior is a trivial
  implementation of a movable dialog handler, just to give you an
  idea how to manage one.

  To modify EasyApp behavior, you would typically replace this call
  with a call to your own handler.

  Requires: nothing
  Receives: pointer to event record for event in a modeless dialog
  Changes:  nothing directly
  Returns:  nothing
  */

void MovableDialogHook(DialogPtr theDialog, EventRecord *theEvent)
{
  EasyMovableFilter(theDialog, theEvent);
}


/*****************************************************************/
/* MouseDownEventHook */
/*****************************************************************/
/*
  EasyApp calls this function when a mouseDown event is received.

  The default behavior parses the mousedown event and implements
  most of the EasyApp behavior - select, resize, zoom, and close
  windows, menu selection, and more. You will probably not need
  to override this hook.

  Requires: nothing
  Receives: pointer to event record
  Changes:  nothing directly
  Returns:  nothing
  */

void MouseDownEventHook (EventRecord *theEvent)
{
  short whereinWin;
  WindowPtr whichWin;

  EasyMouseDown(theEvent);
#if 0
  /**********

    figure out which window they clicked in so that we
    can decide how to deal with this event.

    ***********/
  whereinWin = FindWindow(theEvent->where,&whichWin);
  if(whichWin == main_window ){
    /* bring this window to the front */
    SelectWindow(whichWin);
    /* now deal with the click if we're in the content region */
    switch(whereinWin){
    case inContent:
      redraw();
      break;
    default:
      break;
    }
  }
#endif
}



/*****************************************************************/
/* KeyDownEventHook */
/*****************************************************************/
/*
  EasyApp calls this function when a keyDown event is received.

  The default behavior checks if this is a menu selection (cmd-key)
  and if not, passes the key to the front window's DoKey behavior.
  You will probably not need to override this hook.

  Requires: nothing
  Receives: pointer to event record
  Changes:  nothing directly
  Returns:  nothing
  */

void KeyDownEventHook (EventRecord *theEvent)
{
  /* handles clicks in the easyApp window */
  EasyKeyDown(theEvent);

  /* mabye we need to do something in our window? */
  if( FrontWindow() == g.windowList &&
     !(theEvent->modifiers & cmdKey)){

    /* deal with the keypress! */
//    printf("Key: %c\n",theEvent->message & charCodeMask);
    do_keypress(theEvent);
  }
}



/*****************************************************************/
/* NullEventHook */
/*****************************************************************/
/*
  EasyApp calls this function when the application receives a null
  event. You may receive null events when in the background if you
  set the application to receive them. You set the "SIZE" flags for
  the app to tell the system you want background null events.

  EasyApp is capable of receiving background null events. However, the
  default EasyApp behavior sets the sleep time when in the background
  to the value -1, so it never receives background null events.

  The default behavior sends a DoIdle message to each open window.
  You can provide your windows with a DoIdle behavior. Or you can
  call your own function here to do other idle time processing.

  Requires: nothing
  Receives: nothing
  Changes:  nothing directly
  Returns:  nothing
  */

void NullEventHook (void)
{
  EasyNullEvent ();
}


/*****************************************************************/
/* MouseMovedHook */
/*****************************************************************/
/*
  EasyApp calls this function when the application receives a mouseMoved
  event. This happens if you specify a regoin handle in the global
  g.mouseRegion, used in the call to WaitNextEvent() in the main
  event loop. If the mouse moves outside the specified region, the
  event is generated.

  Typically you respond to that by updating the cursor. EasyApp handles
  cursor updates continuously through the AdjustCursorHook(), so there
  is no need to do that here. As a result, this hook simple calls the
  null-event routine.

  If you want to respond to mouseMoved events, you can call your own
  routine here.

  Requires: nothing
  Receives: pointer to event record
  Changes:  nothing directly
  Returns:  nothing
  */

void MouseMovedHook (EventRecord *theEvent)
{
  EasyNullEvent ();
}


/*****************************************************************/
/* ClipToPrivateHook */
/*****************************************************************/
/*
  EasyApp calls this function when it receives a resume event. This
  means the application is becoming active. Some other app may have
  put data on the clipboard. Your app may wish to use this data. If
  you maintain a private scrap, you must convert the scrap at this
  time.

  For example, if you are using TextEdit, this is a good place
  to convert the TextEdit internal scrap.

  Typically you would call an app-specific function to determine if
  there is data you can use on the clipboard. If there is, you can
  transfer it to your application's private scrap.

  You probably don't need the event record. It is passed in just in case.

  Requires: nothing
  Receives: pointer to resume event record
  Changes:  nothing
  Returns:  nothing
  */

void ClipToPrivateHook (EventRecord *theEvent)
{
  TEFromScrap();        /* ignore error */
}




/*****************************************************************/
/* PrivateToClipHook */
/*****************************************************************/
/*
  EasyApp calls this function when it receives a suspend event. This
  means the application is going inactive. If you use a private
  scrap, the user may want to use that data in some other app. You
  should put it on the public clipboard before going into the
  background.

  Typically you would call an app-specific function to convert any
  data on the private scrap into some form of data appropriate for
  the public clipboard (like data in PICT, TEXT, or other common
  form of data interchange).

  You probably don't need the event record. It is passed in just in case.

  Requires: nothing
  Receives: pointer to resume event record
  Changes:  nothing
  Returns:  nothing
  */

void PrivateToClipHook (EventRecord *theEvent)
{
  TEToScrap();        /* ignore error */
}



/* WINDOW HOOKS */


/*****************************************************************/
/* MakeWindowHook */
/*****************************************************************/
/*
  EasyApp calls this routine to make a new window and initialize
  the various fields in the window record. This routine does not
  add CONTENT to the window. See the EasyMakeWindow() function
  for an example of what must be done.

  You can augment or replace the default behavior. You may wish
  to add new behaviors that are specific to your application, or new
  fields to track the state of your application.

  There are a couple of ways you can do this. You may create your own
  window structure based on the EasyApp "EasyWindow" structure. Or you
  may use the window's refCon field to point to your own app's window-
  related data. However you do this, you MUST include all the EasyWindow
  fields. You may add more as you wish, including new WindowProcs to
  add entirely new behaviors for your windows. But you cannot remove
  any fields from the EasyWindow.

  Your window-making function must return an  error if a problem occurs.
  The caller is responsible for handling the error. In the default
  EasyApp, new windows are invisible. The caller is responsible for
  showing the window if successful.


  Requires: nothing
  Receives: pointer to WindowPtr to store result
  OSType for the file type associated with the window
  Changes:  creates window in memory
  Returns:  error value
  */


OSErr        MakeWindowHook (WindowPtr *macWindow, OSType fileType)
{
  return EasyMakeWindow(macWindow, fileType);
}



/*****************************************************************/
/* GetWindowResIDHook */
/*****************************************************************/
/*
  Called from the EasyMakeWindow() routine. You may be able to use
  that routine to build your own windows. However, if you want to
  use different window types, you need to specify the ID number for
  the resource you want to use to built the window. This function
  receives an OSType corresponding to the file type being created.

  The default hook simply returns kBaseResourceID, which is 1000.
  That's the ID number of the default EasyAPp window resource.

  To modify this behavior you can return the resource ID you want
  to use, or switch based on file type and return any one of
  a selection of resource ID numbers.

  Requires: nothing
  Receives: OSTYpe for the file type being created.
  Changes:  nothing
  Returns:  short value = to WIND resource ID
  */

short GetWindowResIDHook (OSType fileType)
{
  return kBaseResourceID;
}




/*****************************************************************/
/* SetWindowInitHook */
/*****************************************************************/
/*
  Called from the EasyMakeWindow() routine. You may be able to use
  that routine to build your own windows. However, you will want
  your window to have different behavior than the blank EasyApp
  window. The DoInitialize behavior of the window is what sets all
  the functions for the window. So the shell calls here to give
  you a chance to set your own unique initializing behavior for
  a standard EasyApp window.

  Your initializing behavior should be patterned after BaseInitialize().
  Set every field, but set the window behaviors to what you want. If you
  have more than one kind of window, you could use the fileType
  parameter to distinguish among them and set the correct
  DoInitialize behavior for each one.

  Requires: that an EasyWindow exist
  Receives: pointer to the EasyWindow record, file type for this window
  Changes:  DoInitialize field of EasyWindow
  Returns:  nothing

  */

void SetWindowInitHook (EasyWindowPtr easyWindow, OSType fileType)
{
  /*        EasySetWindowInit (easyWindow, fileType); */
  SetViewkelWindowInit(easyWindow,fileType);
}



/* MENU HOOKS */



/*****************************************************************/
/* UpdateMenusHook */
/*****************************************************************/
/*
  Called before MenuSelect() or MenuKey() to make sure menus are
  accurate based on state of application.

  NOTE: If a modal dialog is in front, we never end up here. Calling
  ModalDialog() turns off menus. We will come here if any application
  window is open, or if a modeless dialog is open. We can distinguish
  between them because app windows are in the window list, and
  modeless dialogs are not.

  EasyUpdateMenus() calls various menu hook routines to update
  individual menus, including UpdateAppMenusHook() so you can
  update the app's own menus. It then calls the front window's
  DoUpdateMenu behavior so the window can make the final
  adjustments to the menus.

  You can adjust this process in four distinct ways.

  First, you can replace the call to EasyUpdateMenus() in this
  routine with a call to your own routine, replacing the EasyApp
  behavior completely. You are then responsible for all menu updating.

  Second, if you want unique behavior for one of the three EasyApp
  menus (Apple, File, Edit) you can use the UpdateAppleMenuHook(),
  UpdateFileMenuHook(), or UpdateEditMenuHook() functions to update
  the menu correctly for your app. We do not recommend that you modify
  EasyUpdateMenus(). There should be no need to do so.

  Third, in the UpdateAppMenusHook() you should make a call to
  a routine that updates your application menus based on the state
  of the application.

  Fourth, when you have a window open you can substitute your own
  version of the window's DoUpdateMenu behavior. EasyApp default
  windows use BaseUpdateMenus() in windowProc.c

  Requires: nothing
  Receives: pointer to front window, may be nil
  Changes:  may change state of menu items
  Returns:  nothing
  */

void        UpdateMenusHook (WindowPtr macWindow)
{
  EasyUpdateMenus(macWindow);
}


/*****************************************************************/
/* UpdateAppleMenuHook */
/*****************************************************************/
/*
  Called to update the Apple menu.

  Requires: nothing
  Receives: pointer to front window, may be nil
  Changes:  state of menu items
  Returns:  nothing
  */

void UpdateAppleMenuHook (WindowPtr macWindow)
{
  EasyUpdateAppleMenu(macWindow);
}



/*****************************************************************/
/* UpdateFileMenuHook */
/*****************************************************************/
/*
  Called to update the File menu.

  Requires: nothing
  Receives: pointer to front window, may be nil
  Changes:  state of menu items
  Returns:  nothing
  */

void UpdateFileMenuHook (WindowPtr macWindow)
{
  EasyUpdateFileMenu(macWindow);
}




/*****************************************************************/
/* UpdateEditMenuHook */
/*****************************************************************/
/*
  Called to update the Edit menu.

  Requires: nothing
  Receives: pointer to front window, may be nil
  Changes:  state of menu items
  Returns:  nothing
  */

void UpdateEditMenuHook (WindowPtr macWindow)
{
  EasyUpdateEditMenu(macWindow);
}






/*****************************************************************/
/* UpdateAppMenusHook */
/*****************************************************************/
/*
  EasyApp calls this hook from the EasyUpdateMenus() routine so
  that your app can update its own menus beyond Apple, File, and Edit.
  See the comments for the UpdateMenusHook() for more details.

  The front window will have an opportunity to modify the menus
  after this call.

  Requires: nothing
  Receives: pointer to front window, may be nil
  Changes:  state of menu items
  Returns:  nothing
  */

void        UpdateAppMenusHook (WindowPtr macWindow)
{
  /* call our own handler */
  UpdateViewkelMenus(macWindow);
}


/*****************************************************************/
/* MenuDispatchHook */
/*****************************************************************/
/*
  EasyApp calls this function whenever a menu item is selected.
  Like menu updating, you can modify this process in several ways.

  First, you can replace the call to EasyMenuDispatch() entirely.
  Then you are responsible for all menu dispatching.

  Second, EasyMenuDispatch() calls the front window's DoHandleMenu
  behavior, so the window gets first shot at responding to a menu
  selection. You can substitute your own DoHandleMenu behavior
  unique to your own windows if you wish.

  Third, if you want unique behavior for one of the three EasyApp
  menus (Apple, File, Edit) you can use the AppleMenuHook(),
  FileMenuHook(), or EditMenuHook() functions to dispatch control
  correctly for your app. We do not recommend that you modify
  EasyMenuDispatch(). There should be no need to do so.

  Fourth, after the call to EasyMenuDispatch() below, if you add
  more menus specific to your app, you should call your own
  dispatch routine to handle them.

  Requires: nothing
  Receives: pointer to front window, could be nil or not my window
  short value for menuID, short value for item hit in menu
  Changes:  nothing directly
  Returns:  Boolean indicating whether the selection was handled
  */

Boolean        MenuDispatchHook (WindowPtr macWindow, long menuChoice)
{
  Boolean        handled;
  short menuID,menuItem;

  /* handles Apple, File, Edit menus, replace if you need to */
  handled = EasyMenuDispatch(macWindow, menuChoice);

  if( !handled ){
    /* call our own menu dispatcher */
    handled = ViewkelMenuDispatch(macWindow,menuChoice);
  }

  return handled;
}



/*****************************************************************/
/* AppleMenuHook */
/*****************************************************************/
/*
  EasyApp calls this function when the user selects an item in the
  Apple menu.

  Requires: nothing
  Receives: item number of item hit in Apple menu
  Changes:  nothing directly
  Returns:  nothing
  */

void AppleMenuHook (short menuItem)
{
  EasyHandleAppleMenu(menuItem);

}



/*****************************************************************/
/* AboutBoxHook */
/*****************************************************************/
/*
  EasyApp calls this function when the user selects the About item
  in the Apple menu. Typically you would create your own About box
  and call the routine that handles that from here.

  Requires: nothing
  Receives: nothing
  Changes:  nothing
  Returns:  nothing
  */

void AboutBoxHook (void)
{
  EasyAboutBox();
}




/*****************************************************************/
/* FileMenuHook */
/*****************************************************************/
/*
  EasyApp calls this function when the user chooses an item in the
  File menu.  If you wish to modify the menu, or handle this
  differently, call your own routine here.

  In the default behavior, whatever item the user chooses, EasyApp
  calls the corresponding hook or EasyAPp routine (e.g. NewHook(),
  OpenHook(), EasySave(), and so forth).

  Requires: nothing
  Receives: item number of item hit in File menu
  Changes:  nothing directly
  Returns:  nothing
  */

void FileMenuHook (short menuItem)
{
  EasyHandleFileMenu(menuItem);
}


/*****************************************************************/
/* NewHook */
/*****************************************************************/
/*
  EasyApp calls this function when the user selects New from the
  File menu.

  Caller is responsible for showing new window or error.

  Requires: nothing
  Receives: pointer to a WindowPtr, where result is stored
  OSType for the file type associated with the window
  Changes:  value pointed to by newWindow parameter
  routine called may allocate memory and create window
  Returns:  error value
  */

OSErr        NewHook (WindowPtr *macWindow, OSType fileType)
{
  return EasyNew (macWindow, fileType);
}


/*****************************************************************/
/* OpenHook */
/*****************************************************************/
/*
  EasyApp calls this function when the user selects the Open item
  in the File menu.

  The caller is responsible for showing the window if successful,
  or displaying an error message if not.

  Requires: nothing
  Receives: pointer to a WindowPointer to store result
  pointer to an FSSpec, which is either empty or points to
  the file to open
  Changes:  value pointed to by newWindow parameter
  routine called may allocate memory and create window
  Returns:  error value
  */

OSErr        OpenHook        (WindowPtr *newWindow, FSSpec *theFile)
{
  return EasyOpen(newWindow, theFile);
}



/*****************************************************************/
/* CloseHook */
/*****************************************************************/
/*
  EasyApp calls this function when the user chooses the Close item
  in the File menu or clicks the goAway box in a window. The window
  may be an EasyWindow, or a modeless dialog..

  Requires: an open window
  Receives: Window pointer, short value indicating whether to ask user
  Changes:  closes the window
  Returns:  error value, 0 => closed successfully
  */

OSErr        CloseHook        (WindowPtr macWindow, short saveOptions)
{

  return EasyClose(macWindow, saveOptions);
}



/*****************************************************************/
/* PageSetupHook */
/*****************************************************************/
/*
  EasyApp calls this function when the user selects the Page Setup
  item in the File menu. This routine is provided as a hook so that
  you can support other printing architectures like QuickDraw GX.

  When there is no window open, the default routine uses an
  application default print record. Otherwise, it tells the window
  to handle the Page Setup dialog.

  Requires: the window have a handle to a TPrint record
  Receives: pointer to the window in question
  Changes:  routine called may change contents of TPrint record
  Returns:  nothing
  */

OSErr        PageSetupHook (WindowPtr macWindow)
{
  return EasyPageSetup (macWindow);
}





/*****************************************************************/
/* QuitHook */
/*****************************************************************/
/*
  EasyApp calls this function when the user selects the Print item
  in the File menu.

  Typically you would close all windows (asking the user to save
  changes), perform any cleanup necessary for a graceful exit,
  and set g.done to true.

  Requires: nothing
  Receives: nothing
  Changes:  may shut application down
  Returns:  error value, 0 => we can quit
  */

OSErr        QuitHook (short saveOptions)
{
  return EasyQuit (saveOptions);
}



/*****************************************************************/
/* EditMenuHook */
/*****************************************************************/
/*
  EasyApp calls this function when the user selects an item in the
  Edit menu.  If you wish to modify the menu, or handle this
  differently, call your own routine here.

  The default behavior - EasyHandleEditMenu() - calls the hooks
  or EasyApp routines for the various individual menu items and
  displays any error.

  Requires: nothing
  Receives: item number of item hit in Apple menu
  Changes:  nothing directly
  Returns:  nothing
  */

void EditMenuHook (short menuItem)
{
  EasyHandleEditMenu(menuItem);

}



/*****************************************************************/
/* UndoHook */
/*****************************************************************/
/*
  EasyApp calls this function when the user selects the Undo item
  in the Edit menu.

  The default behavior calls the window's DoUndo behavior. However,
  you may want to undo something other than window actions. This
  hook gives you the flexibility to implement an application-wide
  undo strategy that does not depend on windows.

  Requires: some way to track previous action
  Receives: pointer to window
  Changes:  nothing directly
  Returns:  error value
  */

OSErr        UndoHook (WindowPtr macWindow)
{
  return EasyUndo (macWindow);
}
