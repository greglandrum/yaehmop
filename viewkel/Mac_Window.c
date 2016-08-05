/*****************************************************************/
/*

                                                Adapted from
                               Easy Puzzle

                   Using the EasyApp application shell

                     James E. Trudeau, Nebula, Inc.
                                  and
                            Bradley D. Mohr

        With thanks to Eric Shapiro, Bill Worzel, Bill Hofmann,
        and many more folks, including the fine people at Apple's
        Developer Technical Support and Developer University.

        version 1.0  4/3/95

        Adaptation for viewkel done by greg Landrum 01/15/96
*/
/*****************************************************************/

/*
        Routines for making and managing Mac windows.
*/


#include "viewkel.h"
#include "EasyApp.h"
#include "Mac_defines.h"
#include "Mac_protos.h"

extern EasyGlobal g;
/*extern        PuzzleGlobal puzzleG;*/


/*****************************************************************/
/* WINDOW PROCS FOR THE Viewkel WINDOW */
/*****************************************************************/

/****************************************************************/
/* InitViewkelWindow                                                                                        */
/****************************************************************/
/*
        There are two ways you can do this. In this approach, we call
        the BaseInitialize() routine, let it do all the setup work, and
        then redo some of the work. We set different windowProcs, dispose
        of the scroll bars, and so forth. This is less efficient in terms
        of performance, but more logical in the sense that you start from
        the base and build on it.

        Alternatively, you could simply replace BaseInitialize() entirely
        by duplicating it here and then substituting code where necessary.

        The data parameter is a pointer to the EasyWindow record.

        Requires: window structure exist in memory
        Receives: pointer to window in question, pointer to EasyWindow record
        Changes:  contents of the window structure
        Returns:  error value
*/

OSErr        InitViewkelWindow (WindowPtr macWindow, void *data)
{
        EasyWindowPtr        easyWindow                 = (EasyWindowPtr)data;
        OSErr                        error                         = noErr;
        Handle                        myProcsHandle;
        Handle                        myViewkelHandle;
        GWorldPtr  viewkelGWorld;
        Rect  GWorldRect;

        /* set up a "normal" window */
        error = BaseInitialize(macWindow, data);
        if(error)
        {
                return error;        /* all done, problem */
        }

#if 0
        /* puzzle has no grow box */
        easyWindow->hasGrow        = false;

        /* no scroll bars */
        DisposeControl(easyWindow->horizScroll);
        DisposeControl(easyWindow->vertScroll);
        easyWindow->horizScroll = 0;
        easyWindow->vertScroll = 0;
#endif

        easyWindow->DoContentClick                = (WindowProc)ClickInViewkel;
/*        easyWindow->DoBackgroundClick        = (WindowProc)BackgroundViewkelClick;
*/
        easyWindow->DoResize = (WindowProc)SizeViewkelWindow;
/*
        easyWindow->DoUpdate                    = (WindowProc)UpdateViewkel;
*/

#if 0
        /* allocate memory for special puzzle procs */
        error = EasyAllocateHandle( sizeof(MoreProcs), &myProcsHandle);
        if (error)
            {
              return error;
            }

        /* move it high and lock it in memory permanently, ignore return */
        EasyLockHandle(myProcsHandle);

        /* stuff a pointer to this in the moreProcs field */
        myProcs = (MoreProcsPtr)*myProcsHandle;
        easyWindow->moreProcs = (long) myProcs;

        /* set our own procs */
        myProcs->DoImport = (WindowProc)ImportPuzzlePICT;

        /* allocate memory for a PuzzleData structure */
        error = EasyAllocateHandle( sizeof(PuzzleData), &myPuzzleHandle);
        if (error)
            {
              return error;
            }

        /* move it high and lock it in memory permanently, ignore return */
        EasyLockHandle(myPuzzleHandle);

        /* stuff a pointer to this in the windowData field */
        myPuzzle = (PuzzlePtr) *myPuzzleHandle;
        easyWindow->windowData = (long) myPuzzle;

        /* set our own data */
        myPuzzle->frameColor = puzzleG.frameColor;        /* default color */

        /* identify drag and drop routines for later installation */
        /* track routine is the same, replacing only the receive routine */
        easyWindow->dragRoutine    = (ProcPtr)EasyDragTrackHandler;
        easyWindow->receiveRoutine = (ProcPtr)PuzzleReceiveHandler;
#endif
        easyWindow->docWidth  = MAC_DOCUMENT_SIZE;
        easyWindow->docHeight = MAC_DOCUMENT_SIZE;
        g_xmax = MAC_DOCUMENT_SIZE;
        g_ymax = MAC_DOCUMENT_SIZE;

        /* get the offscreen graphics world to draw in */
        GWorldRect.left = 0;
        GWorldRect.right = easyWindow->docWidth;
        GWorldRect.top = 0;
        GWorldRect.bottom = easyWindow->docHeight;

        error = NewGWorld(&viewkelGWorld,0,&GWorldRect,NULL,NULL,NULL);
        if( !viewkelGWorld || error != noErr ){
                fatal("Can't get GWorld... We're doomed!");
        }

        /* store the location of the GWorld in the window's data structure */
        easyWindow->windowData = (long)viewkelGWorld;

        /* store the font size */
        Mac_globals.fontSize = 12;
        TextSize(12);

        return error;
}



/*****************************************************************/
/* DrawViewkelPage */
/*****************************************************************/
/*
        This function should draw one page of your document. It should
        assume that the port has been set correctly on entry, because
        you might use this routine to draw for any number of purposes,
        not just to render the contents on screen. For example, you
        may print a page using this routine.

        EasyApp calls this function from the BasePrint() function.
        and from the BaseDrawWindow() function..

        The data parameter is a pointer to a short value holding the
        page number to draw. The number is one-based, that is the first
        page in the document is page 1 (as opposed to zero-based, where
        the pages would be counted from zero).

        For now, only one page is being used, this makes portability issues
         a little simpler.

        Requires: QuickDraw port set correctly
        Receives: pointer to window, pointer to a short value holding page #
        Changes:  should draw a page, change nothing
        Returns:  error value, not used

*/
OSErr        DrawViewkelPage (WindowPtr macWindow, void *data)
{
        EasyWindowPtr        easyWindow = GetEasyWindowPtr(macWindow);
#if 0
        PicHandle                thePicture = GetViewkelPicture(easyWindow);
        RgnHandle                displayRgn = GetViewkelDisplayRgn(easyWindow);
#endif
        Rect                        theFrame;
        RgnHandle                saveClip;

#if 0
        /**********


          set the display region as the clip region
          this is the area of the pict that is visible inside
          the main window.

          ***********/
        if (solutionRgn)
            {
              saveClip = NewRgn();
              GetClip(saveClip);
              SetClip(solutionRgn);
            }

        /* draw the picture */
        GetPuzzleFrame(easyWindow, &theFrame);
        DrawPicture(thePicture, &theFrame);

        /* restore clip */
        if (solutionRgn)
            {
              SetClip(saveClip);
              DisposeRgn(saveClip);
            }

        /* if there's a scramble dialog open, draw a grid */
        if (scramble)
            {
              DrawPuzzleGrid (macWindow, scramble);
            }

        if (thisPiece)
            {
              ForEveryPiece(thisPiece, DrawPiece, (void *) macWindow->visRgn);
            }

#endif
        return noErr;
}




/***************************************************************
  SizeViewkelWindow
***************************************************************/
/*
        We will use this as the window's DoResize behavior.

        We want to resize the window and update the g_xmax and
         g_ymax parameters so that the window can be drawn properly

         For the moment, this isn't going to do anything,
          we'll rely on the scroll bars to move shit around
          for us.

        data parameter is not used in this routine

        Requires: PICT data attached to window
        Receives: WindowPtr
        Changes:  size of window
        Returns:  error value, not used
*/
OSErr        SizeViewkelWindow (WindowPtr macWindow, void* data)
{
        EasyWindowPtr        easyWindow        = GetEasyWindowPtr(macWindow);
        Rect                        displayFrame;
        Rect                        portRect;
        short                        newWidth   = LoWord(*(long*)data);
        short                        newHeight  = HiWord(*(long*)data);
        short                        currentHeight;        /* of picture */
        short                        currentWidth;
        double                        scaleFactor;


#if 0
        /* display height will be 200 pixels */
        newHeight = kMaxHeight;        /* all puzzles will be 200 pixels high */

        /* calculate new width */
        scaleFactor =  kMaxHeight / currentHeight;

        newWidth = (short)(currentWidth * scaleFactor);

        if (newWidth < kMinWidth)        /* too narrow, stretch it to fit */
            {
              newWidth = kMinWidth;
            }
        else if (newWidth > kMaxWidth) /* too wide, squeeze it to fit */
            {
              newWidth = kMaxWidth;
            }

        /* set the display frame in the window */
        displayFrame.top    = kPuzzleMargin;
        displayFrame.left   = kPuzzleMargin;
        displayFrame.bottom = kPuzzleMargin + newHeight;
        displayFrame.right  = kPuzzleMargin + newWidth;
        SetPuzzleFrame(easyWindow, &displayFrame);

        /* now set window size */
        newHeight += 2 * kPuzzleMargin;
        newWidth  += 2 * kPuzzleMargin;

        SizeWindow(macWindow, newWidth, newHeight, true);
        easyWindow->docWidth  = newWidth;
        easyWindow->docHeight = newHeight;

        /* because this is a new picture, dispose existing regions */
        if (solutionRgn)
            {
              DisposeRgn(solutionRgn);
            }

        if (borderRgn)
            {
              DisposeRgn(borderRgn);
            }

        /* make a new solution region */
        solutionRgn = NewRgn();
        RectRgn(solutionRgn, &displayFrame);
        SetPuzzleSolutionRgn (easyWindow, solutionRgn);

        /* make a new border region */

        /* start with the window's portRect */
        portRect = macWindow->portRect;
        borderRgn = NewRgn();
        RectRgn(borderRgn, &portRect);

        /* subtract the solution area */
        DiffRgn(borderRgn, solutionRgn, borderRgn);
        SetPuzzleBorderRgn (easyWindow, borderRgn);
#endif
        SetPort(macWindow);
        if (easyWindow->hasGrow)        /* must ensure this gets invalidated */
        {
                Rect growBoxRect;

            SetRect(&growBoxRect, macWindow->portRect.right - 15,
                                                      macWindow->portRect.bottom - 15,
                                                      macWindow->portRect.right,
                                                      macWindow->portRect.bottom);
            InvalRect(&growBoxRect);
        }

        SizeWindow(macWindow, newWidth, newHeight, true);
/*
        easyWindow->docWidth  = newWidth-15;
        easyWindow->docHeight = newHeight-15;
*/
        /* resize scroll bars */
        easyWindow->DoSizeScrollBars(macWindow, nil);

#if 0
        g_xmax = newWidth-15;
        g_ymax = newHeight-15;
#endif
        return noErr;
}




/*****************************************************************/
/* ForcePuzzleUpdate */
/*****************************************************************/
/*
        Inalidate just the area of the puzzle proper, not the whole
        window.

        Requires: nothing
        Receives: WindowPtr
        Changes:  invalidates window
        Returns:  nothing
*/

void ForceViewkelUpdate (WindowPtr macWindow)
{
        EasyWindowPtr        easyWindow        = GetEasyWindowPtr(macWindow);
        GrafPtr                        oldPort;
        Rect                        theFrame;

#if 0
        GetPort(&oldPort);
        SetPort(macWindow);
        /* put an invalidate rectangle call in here */
        SetPort(oldPort);
#endif
}



/*****************************************************************/
/* WINDOW RELATED UTILITIES */
/*****************************************************************/

/*****************************************************************/
/* NewViewkel */
/*****************************************************************/
/*
        Called from the NewHook() routine. Caller shows window or error.
        All we have to do is create the window.

        This is copied almost completely from EasyNew(). We ignore the
        fileType, there's only one type in the Puzzle App.

        In the Puzzle App, we want to do two new things. When the user
        opens a "new" puzzle we'll load a pict from inside the app itself.
        They can always replace it with cut, copy, or paste. And we want
        to resize the window after loading the puzzle.

        Requires: nothing
        Receives: pointer to WindowPtr to store new WindowPtr
        Changes:  creates and initializes window
        Returns:  error value
*/

OSErr        NewViewkel (WindowPtr *newWindow, OSType fileType)
{
        OSErr                        error                 = noErr;
        WindowPtr                 macWindow;
        EasyWindowPtr         easyWindow;
        Handle                        thePicture;

#if 0

/* for now, do nothing.  we'll add content here later (mabye) */
        error = EasyNew(newWindow, fileType);
        if(error)
        {
                return error;
        }
#endif

        return error;
}


/****************************************************************/
/* SetPuzzleWindowInit                                                                                        */
/****************************************************************/
/*
        Called from the SetWindowInitHook().

        Set the correct initializing routine for a puzzle window.

        Requires: nothing
        Receives: pointer to EasyWindow record, fileType (not used)
        Changes: value of DoInitialize field of the easyWindow
        Returns:  nothing
*/

void SetViewkelWindowInit (EasyWindowPtr easyWindow, OSType theFile)
{
#if 0
        if (theFile == kTEXTType)
        {
                easyWindow->DoInitialize = (WindowProc)InitTextWindow;
        }
        else
        {
                easyWindow->DoInitialize = (WindowProc)InitViewkelWindow;
        }
#else
                easyWindow->DoInitialize = (WindowProc)InitViewkelWindow;
#endif
}






/***************************************************************
   DisposeViewkelContents

   clears out the objects in the window
***************************************************************/
void DisposeViewkelContents (WindowPtr macWindow)
{
  EasyWindowPtr        easyWindow   = GetEasyWindowPtr(macWindow);
  return;
}

/****************************************************************/
/* ClickInViewkel
/****************************************************************/
/*
        This is the puzzle's DoContentClick behavior. EasyApp sends that
        message to the window when the window is in front, is an EasyApp
        window, and the click is in the window content.

        We will look for clicks in puzzle pieces, starting at the back
        and going forward. Why? Because pieces are drawn from the front
        going back, so the last piece is on top.

        If the click is in a puzzle piece, and the user holds the mouse
        button down, we will drag it around with the mouse movement.

        When released, if it is near it's solution position, we will
        flash the piece, drop it into the solution, and remove it from
        the piece list.

        The additional data parameter is a pointer to the event record.
        The mouse location is in the where field, in global coordinates.

        Requires: nothing
        Receives: WindowPtr, pointer to event record for this click
        Changes:  may move a piece around
        Returns:  OSErr, not used
*/

OSErr        ClickInViewkel (WindowPtr macWindow, void *data)
{
        EasyWindowPtr        easyWindow = GetEasyWindowPtr(macWindow);
        EventRecord                *theEvent  = (EventRecord *)data;
        OSErr                        error           = noErr;
        Point                        localPoint;
        Boolean                        done = false;
        GrafPtr                        oldPort;
        int xpos,ypos;
        short vScrollPos,hScrollPos;
        int found;

        GetPort(&oldPort);
        SetPort(macWindow);

        /* start with the default behavior */
        BaseContentClick(macWindow,data);

        /* convert to local coordinates for this window */
        EasyGlobalToLocal(macWindow, theEvent->where, &localPoint);

        /************

         get the position of the scroll bars so that clicks are dealt with properly
           when the scroll bars have been moved.

        *************/
        hScrollPos = GetCtlValue(easyWindow->horizScroll);
        vScrollPos = GetCtlValue(easyWindow->vertScroll);

        xpos = (int)((localPoint.h+hScrollPos)/GRAPHICS_SCALE);
        ypos = (int)((localPoint.v+vScrollPos)/GRAPHICS_SCALE);

        /* check to see if anything was clicked in */
        switch(mainmode){
        case CHOOSE:
          if( whichobj->prim->molec ){
        found = select_atom(whichobj->prim->molec,xpos,ypos);
        if(!found){
        show_selected_data(num_selected,whichobj,xpos,ypos);
        unselect_all_atoms(num_selected,whichobj);
        num_selected = 0;
        }

    }
    else if( whichobj->prim->MO_surf && whichobj->prim->MO_surf->molec ){
        found = select_atom(whichobj->prim->MO_surf->molec,xpos,ypos);
         if(!found){
        show_selected_data(num_selected,whichobj,xpos,ypos);
        unselect_all_atoms(num_selected,whichobj);
        num_selected = 0;
        }

    } else{
        select_object(xpos,ypos);
    }
    break;
  default:
                select_object(xpos,ypos);
        }


        /* do the redraw  */
        redraw();
        SetPort(oldPort);

        return error;
}



/*****************************************************************/
/* BaseBackgroundClick */
/*****************************************************************/
/*
        Bring the window forward, set the port, and redraw

        Requires: nothing
        Receives: window pointer, pointer to event record in data
        Changes:  may result in drag beginning
        Returns:  error value, 0 => OK
*/

OSErr BackgroundViewkelClick (WindowPtr macWindow, void *data)
{
        OSErr                        error = noErr;
        EasyWindowPtr        easyWindow = GetEasyWindowPtr(macWindow);
        EventRecord                *theEvent = (EventRecord*)data;
        Boolean                        activate = true;        /* Bring window forward? Assume yes. */
        Point                        localPoint;
        GrafPtr                        oldPort;
        RgnHandle                solved;

        SetPort(macWindow);

        SelectWindow(macWindow);

        redraw();
        return error;
}



OSErr        UpdateViewkel (WindowPtr macWindow, void *data)
{
        OSErr                        error           = noErr;
#if 0
        EasyWindowPtr        easyWindow = GetEasyWindowPtr(macWindow);
        EventRecord                *theEvent  = (EventRecord *)data;

        Point                        localPoint;
        Boolean                        done = false;


        /* convert to local coordinates for this window */
        EasyGlobalToLocal(macWindow, theEvent->where, &localPoint);

        /*******

          for now do nothing, eventually we'll want to check if any objects
          have been selected.

        *******/

        /* set the active port and redraw */
        SetPort(macWindow);
        redraw();
#endif
        return error;
}

CGrafPtr origPort;
GDHandle origDev;


/*********************

        Mac_ClearScreen

        This clears out everything that's in the offscreen graphics
                world stored in the windowData field of the easyWindow,
                locks the pixmap, sets the active port to be the gworld,
                and basically does everything needed to get ready for a
                redraw.

        NOTE:  This function *MUST* be called at the beginning of
                 every call to redraw for things to get updated properly.

*******************/
void Mac_ClearScreen(void)
{
        WindowPtr macWindow;
        EasyWindowPtr easyWindow;
        OSErr theErr;
        Boolean locked;
        Rect the_rect;
        PixMapHandle thePixmap;
        GWorldPtr theGWorld;

        /* figure out where we are copying to */
        macWindow = g.windowList;
        easyWindow = GetEasyWindowPtr(macWindow);

        /* save the original port so that we can restore it later */
        GetGWorld(&origPort,&origDev);

        /* start off by setting the GWorld as the active port */
        theGWorld = (GWorldPtr)easyWindow->windowData;
        SetGWorld(theGWorld,0);

        /* get and lock the pixmap */
        thePixmap = GetGWorldPixMap(theGWorld);
        locked = LockPixels(thePixmap);
        if( !locked ) fatal("Can't lock pixmap... sorry.");

        /* clear out the pixmap */
        EraseRect(&(theGWorld->portRect));

        /*****

         that's it.  The redraw should now happen, followed by a
         call to g_switch_buffers (which calls Mac_CopyGWorld)

        ******/
}

/*********************

        Mac_CopyGWorld

        This copies the displayable region of the offscreen GWorld into
                the active window.

        NOTES:
          -This function *MUST* be called at the end of
                 every call to redraw for things to get updated properly.
          -The pixmap associated with the GWorld should already be locked
             by the time we get here.

*******************/
void Mac_CopyGWorld(void)
{
        WindowPtr macWindow;
        EasyWindowPtr easyWindow;
        Rect sourceRect,destRect;
        short width, height;
        short vScrollPos,hScrollPos;
        GWorldPtr theGWorld;

        /* figure out where we are copying to */
        macWindow = g.windowList;
        destRect = macWindow->portRect;
        easyWindow = GetEasyWindowPtr(macWindow);

        /* figure out where we are copying to */
        macWindow = g.windowList;
        destRect = macWindow->portRect;
        easyWindow = GetEasyWindowPtr(macWindow);

        /* remove the region covered by the scroll bars */
        destRect.bottom -= 15;
        destRect.right -= 15;

        /* figure out the width and height of the region */
        width = destRect.right - destRect.left;
        height = destRect.bottom - destRect.top;

        /* get the position of the scroll bars */
        hScrollPos = GetCtlValue(easyWindow->horizScroll);
        vScrollPos = GetCtlValue(easyWindow->vertScroll);

        /*****
         use the scroll bar positions to give us the top and left
         positions of the window.

         make sure that we don't copy from invalid areas!
        ******/
        if( (vScrollPos + height) <= (short)g_ymax ){
                sourceRect.top = vScrollPos;
                sourceRect.bottom = vScrollPos + height;
        }else{
                sourceRect.top = (short)g_ymax - height;
                sourceRect.bottom = (short)g_ymax;
        }
        if( (hScrollPos + width) <= (short)g_xmax ){
                sourceRect.left = hScrollPos;
                sourceRect.right = hScrollPos + width;
        }else{
                sourceRect.left = (short)g_xmax - width;
                sourceRect.right = (short)g_xmax;
        }

        /* get a pointer to our GWorld */
        theGWorld = (GWorldPtr)easyWindow->windowData;

        /* unlock the pixels associated with the graphics world now */
        UnlockPixels(GetGWorldPixMap(theGWorld));

        /* restore the original port */
        SetGWorld(origPort,origDev);


        /* okey dokey, we're set, so do the copy now */
        CopyBits(&(((GrafPtr)theGWorld)->portBits),
                         &(((GrafPtr)macWindow)->portBits),
                         &sourceRect,&destRect,srcCopy,0);

}





