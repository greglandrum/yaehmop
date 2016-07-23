
/*******************************************************
*      Copyright (C) 1995, 1998, 1999 Greg Landrum
*
*  This file is part of yaehmop.
*
*   This is free software.
* 
*  Permission is granted to modify, or otherwise fold, spindle, and mutilate this
*    code provided all copyright notices are left intact.
*
*  This code may be distributed to your heart's content, in whatever form,
*    provided no fee is charged for the distribution, all copyright notices are
*    left intact, and the source is distributed (without fee) along with any
*    binaries to anyone who requests it.
*
*  There are, of course, no warranties at all on this program.
*
********************************************************************/

/***
  Recent Edit History:
   26.04.98 gL:
     added initialization of near plane clipping toggle
   30.05.98 gL:
     added initialization of polyhedron outline toggle
   26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)
   18.01.99 gL:
     added initialization of dump_grid toggle
   29.01.99 gL:
     added initialization of global_read_on toggle and modified
     command line parsing to allow multiple arguments
   28.06.1999 gL:
     added an XFlush to the code to free up X resources.
   14.02.2004 gL:
     added a "10th anniversary" message :-)

***/

const char greetings[]="Welcome to the 10th Anniversary edition of YAeHMOP!\n";

#include "viewkel.h"

#ifdef MAC_GRAPHICS
#include <SIOUX.h>
#endif

void main(argc, argv)
int argc;
char **argv;
{
  int i;
  int xsize=0,ysize=0;
#ifdef X_GRAPHICS
  Display *local_disp;
#endif
#ifdef MAC_GRAPHICS
/****
	 set up stuff so that we can use SIOUX to handle console i/o but
	 we can still open our own windows and menus.  fun fun.
*****/
SIOUXSettings.standalone = FALSE;
SIOUXSettings.setupmenus = FALSE;
SIOUXSettings.initializeTB = FALSE;
SIOUXSettings.asktosaveonclose = FALSE;
#endif

#ifdef PARTICLES_FOR_SURF
  srand48(time());
#endif
  fprintf(stderr,greetings);

  /* check to see if we can open the display */
#ifdef X_GRAPHICS
  local_disp = XOpenDisplay(getenv("DISPLAY"));
  if(!local_disp){
    doing_X = 0;
    doing_tek = 1;
  }
  else{
    doing_X = 1;
    doing_tek = 0;
  }
#endif /* X_GRAPHICS */
  global_read_on = 1;
#ifdef X_GRAPHICS
  char *fileName=0;
  /* parse the command line arguments */
  if(argc > 1){
    /****
       when adding things here, make sure to add multiple letter
       arguments first!
    ****/
    for( i = 1; i<argc; i++){
      if(!strncmp(argv[i],"-gr",3)){
	global_read_on = 0;
      } else if(!strncmp(argv[i],"-g",2)){
	fprintf(stderr,"New size!\n");
	if( argc < i+2 ){
	  fprintf(stderr, "You tease, give me a real size!\n");
	  xsize = 0;
	  ysize = 0;
	}
	else{
	  sscanf(argv[i+1],"%d", &xsize);
	  sscanf(argv[i+2],"%d", &ysize);
	  i+=2;
	}
      } else if( !strncmp(argv[i],"-t",2) ){
	doing_X = 0;
	doing_tek = 1;
      } else if( !strncmp(argv[i],"-b",2) ){
	useButtons=0;
      }
      else if (!strncmp(argv[i],"-f",2) ){
	fileName = argv[i+1];
	i++;
      }
    }
  }
#endif

#ifdef MAC_GRAPHICS
	doing_Mac = 1;
#endif
  ident = (matrix_type *)D_CALLOC(1,sizeof(matrix_type));
  mainortho = (matrix_type *)D_CALLOC(1,sizeof(matrix_type));  
  if( !ident || !mainortho )fatal( "Memory Allocation." );
  
  /* initialize the identity matrix */
  for(i=0;i<DIM;i++) ident->matrix[i][i]=1.0;

  /* put the identity matrix on the stack */
  loadmatrix( ident );

  /* point the head at null */
  head = 0;

  projviewon = 1;
  
  near_plane_clipping_on = 1;
  outline_polyhed_on = 1;
  dump_grids_on = 0;
  mainmode = NONE;
  /* get some space for the camera storage node */
  camera = (camera_type *)D_CALLOC(1,sizeof(camera_type));
  if( !camera ) fatal("Memory Allocation." );
  /* initialize the camera */
  camera->lf.z=-159.75;
  camera->vup.y=1.0;
  camera->hsize=10.0;
  camera->yaspect=(float)g_ymax/GWINWIDTH;
  camera->foclength = fabs(camera->lf.z / camera->hsize);
/*  camera->foclength = 100.0;*/
  /*************************
	
    important note:
    this is where the event loop starts on the Mac,
    so we can't rely on the program coming back after this.
		  
  *************************/
  g_initgraphics(xsize,ysize);


#ifdef DEBUGGING_FROM_HOME
  new_molecule();
  redraw();
  whichobj->trans.z += 5.0;
  redraw();
  whichobj->trans.z += 5.0;
  redraw();
  exit();
#endif  

  button_wins = 0;
#ifdef X_GRAPHICS
  if( doing_X ){
    if(useButtons) build_main_button_window(&button_wins);
/*    build_PS_options_button_window(&button_wins);*/
    if(fileName) general_read(fileName);
    do_events();
  }
#endif

#ifdef TEK_GRAPHICS
  if( doing_tek ){
    parse_commands();
  }
#endif


#ifdef X_GRAPHICS
  if( doing_X ){
    XFlush(disp);
    XFreePixmap(disp,gpix);
    XFreeGC(disp,graphgc);
    XFreeGC(disp,blackgc);
    XFreeGC(disp,bigtextgc);
    XFreeGC(disp,smalltextgc);
    for(i=0;i<NUM_COLORS;i++){
      XFreeGC(disp,graygc[i]);
    }
  }
#endif

  clearall();

  if(useButtons) free_but_win(button_wins,&button_wins);

  if(head) D_FREE(head);
  D_FREE(camera);
  D_FREE(ident);
  D_FREE(mainortho);


#ifdef MEM_DEBUG
  d_check_core();
#endif  

}
