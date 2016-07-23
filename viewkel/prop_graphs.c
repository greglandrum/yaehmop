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



/********

  this has got the stuff for dealing with property graphs....

  Created by greg Landrum June 1994
*********/  

/********
  Recent Edit History

  06.03.99 gL:
     support for printing curve names added.

*********/

#include "viewkel.h"



/****************************************************************************
 *
 *                   Procedure new_prop_graph
 *
 * Arguments: filename (optional): pointer to type char
 *            
 * Returns: none
 *
 * Action: does everything to get space for and read in a new prop_graph
 *
 ****************************************************************************/
void new_prop_graph(char *filename)
{
  char instring[MAX_STR_LEN],file_name[80],foostring[80];
  char *theinline;
  char failed=0;
  graph_type *the_graph;
  FILE *infile;
  int eof_hit;

  /* set up a new object to hold the thing */
  makenewobject();
  whichobj = head->obj;

  /* now build the prop_graph primitive */
  whichobj->prim = (prim_type *)D_CALLOC(1,sizeof(prim_type));
  if( !whichobj->prim )fatal("Can't get space for prop_graph primitive.");
  whichobj->prim->which = PROP_GRAPH;
  
  whichobj->prim->prop_graph =
    (prop_graph_type *)D_CALLOC(1,sizeof(prop_graph_type));
  if( !whichobj->prim->prop_graph )
    fatal("Can't get space for prop_graph.");

  whichobj->prim->prop_graph->the_data =
    (graph_type *)D_CALLOC(1,sizeof(graph_type));
  if( !whichobj->prim->prop_graph->the_data )
    fatal("Can't get space for prop_graph->the_data.");


#ifndef USING_THE_MAC
  if( !filename ){
    display("Look in the xterm...");
    
#ifndef USE_READLINE
    printf("Enter the file name containing the property data: ");
    scanf("%s",file_name);
#else
    theinline= readline("Enter the file name containing the property data: ");
    add_history(theinline);
    if( theinline ){
      sscanf(theinline,"%s",file_name);
      D_FREE(theinline);
    } else {
      error("Bad file name");
      file_name[0] = 0;
    }
#endif

  } else{
    strcpy(file_name,filename);
  }

  /* open the file */
  infile = fopen(file_name,"r");
  if(!infile){
    printf("Problems opening file: %s\n",file_name);
    display("oooooops!");
    failed = 1;
  }
#else
	if(!filename){
		infile = choose_mac_file(file_name,MAC_FOPEN_OPEN_CD);
	} else{
		strcpy(file_name,filename);
		infile = fopen(file_name,"r");
	}
	if( !infile ){
    printf("Problems opening file: %s\n",file_name);
 	 	display("oooooops!");
   	failed = 1;
  }
#endif

  /* get the initial data */
  if( !failed ){
    strcpy(whichobj->prim->prop_graph->filename,file_name);
    /* check to see if this is a COOP or a DOS */
    skipcomments(infile,instring);
    upcase(instring);
    if( strstr(instring,"COOP") ){
      whichobj->prim->prop_graph->type = COOP_PROP;
    } else{
      whichobj->prim->prop_graph->type = DOS_PROP;
    } 
    read_graph_data(infile,whichobj->prim->prop_graph->the_data,1,1);
    preprocess_graph_data(whichobj->prim->prop_graph->the_data);
  }

  /* check to see if any curves were actually read in.... */
  if(!whichobj->prim->prop_graph->the_data->num_curves || failed){
    /* no... free the memory that we asked for */
    D_FREE(whichobj->prim->prop_graph->the_data);
    D_FREE(whichobj->prim->prop_graph);
    D_FREE(whichobj->prim);
    D_FREE(whichobj);
    whichobj=0;
    head->obj = 0;
    head = head->next;
    return;
  }
  else{
    whichobj->scale.x=whichobj->scale.y=whichobj->scale.z=2.0*GRAPHICS_SCALE;
    whichobj->cent.x=-180;whichobj->cent.y=180;
    whichobj->cent.z=0*GRAPHICS_SCALE;
    whichobj->trans.x=0;whichobj->trans.y=0;
    whichobj->trans.z=0;

    /* now fill in the legends... */
    the_graph = whichobj->prim->prop_graph->the_data;
    strcpy(the_graph->ylegend,"Energy (eV)");

    whichobj->prim->prop_graph->max_x = the_graph->max_x;
    whichobj->prim->prop_graph->max_y = the_graph->max_y;
    whichobj->prim->prop_graph->min_x = the_graph->min_x;
    whichobj->prim->prop_graph->min_y = the_graph->min_y;

    /***********
      find out whether or not there is integration data in the
      data file.
      ***********/
    if(skipcomments(infile,instring) >= 0 ){
      upcase(instring);
      if( instring[0] == '#' && strstr(instring,"INTEGRATION") ){
	
	/* there is, get space for it, then read in the data */
	whichobj->prim->prop_graph->the_integration =
	  (graph_type *)D_CALLOC(1,sizeof(graph_type));
	if( !whichobj->prim->prop_graph->the_integration )
	  fatal("Can't get space for prop_graph->the_integration.");
	
	read_graph_data(infile,whichobj->prim->prop_graph->the_integration,1,0);
	
	/*******
	  scale the integration data
	  *******/
	if( whichobj->prim->prop_graph->type == COOP_PROP ){
	  the_graph = whichobj->prim->prop_graph->the_integration;
	  the_graph->max_x = whichobj->prim->prop_graph->max_x;
	  the_graph->max_y = whichobj->prim->prop_graph->max_y;
	  the_graph->min_x = whichobj->prim->prop_graph->min_x;
	  the_graph->min_y = whichobj->prim->prop_graph->min_y;
	} else{
	  /* for DOS curves the integration data should start on a scale of 0 to 1 */
	  the_graph = whichobj->prim->prop_graph->the_integration;
	  the_graph->min_x = 0.0;
	  the_graph->max_x = 1.0;
	}
	    
	preprocess_graph_data(whichobj->prim->prop_graph->the_integration);
      }
    }
    
    /* find the fermi energy */
    rewind(infile);
    eof_hit = skipcomments(infile,instring);
    upcase(instring);
    
    while(eof_hit >= 0 && !strstr(instring,"FERMI_ENERGY")){
      eof_hit = skipcomments(infile,instring);
      upcase(instring);
    }
    
    if( eof_hit >= 0 ){
      sscanf(instring,"%s %lf",foostring,&(whichobj->prim->prop_graph->Fermi_E));
      whichobj->prim->prop_graph->show_fermi = 1;
    }
    
    /* make the button window */
#ifdef INTERACTIVE_USE
    build_prop_graph_button_window(&button_wins,whichobj->prim->prop_graph);
    whichobj->prim->but_win = button_wins;
#endif
  }
}
