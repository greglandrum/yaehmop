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

  03.05.98 gL:
    stupid update bug involving total_E x range not changing fixed.
    bounding boxes adjusted.
   26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)

  28.01.99 gL:
     read_walsh_data now deals with the beginning of the file having
     lines with '#'.
***/

/********

  this has got the stuff for dealing with Walsh Diagrams

  Created by greg Landrum June 1994
*********/
#include "viewkel.h"


/****************************************************************************
 *
 *                   Procedure read_walsh_data
 *
 * Arguments: infile: pointer to type FILE
 *               graph: a pointer to graph_type
 *               tot_E: a pointer to graph_type
 *
 * Returns: none
 *
 * Action:  Reads the data out of 'infile.
 *       This reads until it hits a line beginning with a # mark.
 *
 ****************************************************************************/
void read_walsh_data(FILE *infile,graph_type *graph,graph_type *tot_E)
{
  char instring[MAX_STR_LEN];
  int num_curves;
  int i,j;
  int num_p;
  float xval;
  float min_x,min_y,max_x,max_y;


  if( !graph )FATAL_BUG("No graph passed to read_walsh_data.");
  if( !tot_E )FATAL_BUG("No tot_E passed to read_walsh_data.");

  if( skipcomments(infile,instring) < 0){
    error("That file is empty!");
    display("Too Bad....");
    return;
  }

  while(instring[0] == '#' ) skipcomments(infile,instring);

  /* read out the number of orbitals and number of steps */
  sscanf(instring,"%d",&num_curves);

  skipcomments(infile,instring);
  sscanf(instring,"%d",&num_p);

  if(num_curves < 1){
    error("That data file has less than one data set in it... check it please.");
    display("Yikes!");
    return;
  }


  /*******
    get space for the data
  *******/
  graph->raw_data = (point_type2D *)D_CALLOC(num_p*num_curves,
                                           sizeof(point_type2D));
  tot_E->raw_data = (point_type2D *)D_CALLOC(num_p,sizeof(point_type2D));
  if( !graph->raw_data || !tot_E->raw_data){
    error("Can't get space for raw_graph data....");
    display("Oh well...");
    return;
  }

  /**********

    get space to keep track of which curves should be displayed
      we just need two slots here, one for the MO's and one for
      the total energy.

  ******/
  graph->curves_to_display = (char *)D_CALLOC(2,sizeof(char));
  if( !graph->curves_to_display )fatal("Can't get space for curves_to_display.");

  graph->curves_to_display[0] = 1;
  graph->curves_to_display[1] = 1;

  graph->styles = (char *)D_CALLOC(2,sizeof(char));
  if( !graph->styles )fatal("Can't get space for styles.");
  graph->styles[0] = 0;
  graph->styles[1] = 1;


  graph->do_x_tics = 1;
  graph->do_y_tics = 1;

  if( skipcomments(infile,instring) < 0 ){
    error("EOF hit while reading file.  read_walsh_data aborted\n");
    graph->num_curves = 0;
    return;
  }
  for( i=0; i<num_p; i++){
    /*******

      assume that the x coordinate is at the end of the line

      *******/
    sscanf((const char *)strtok(instring," "),"%lf",&(graph->raw_data[i*num_curves].y));
    for(j=1;j<num_curves;j++){
      sscanf((const char *)strtok(0," "),"%lf",&(graph->raw_data[i*num_curves+j].y));
    }

    /* read the total energy */
    sscanf((const char *)strtok(0," "),"%lf",&(tot_E->raw_data[i].y));

    /* now get and copy the x coordinate */
    sscanf((const char *)strtok(0," "),"%lf",&xval);
    for(j=0;j<num_curves;j++){
      graph->raw_data[i*num_curves+j].x = xval;
    }
    tot_E->raw_data[i].x = xval;

    if( skipcomments(infile,instring) < 0 && i != num_p - 1){
      error("EOF hit while reading file.  read_walsh_data aborted\n");
      graph->num_curves = 0;
      return;
    }

  }
  graph->num_p = num_p;
  graph->num_curves = num_curves;

  tot_E->num_p = num_p;
  tot_E->num_curves = 1;

  /* find the minimum and maximum values (to use in scaling) */
  min_y = min_x = 1e10;
  max_y = max_x = -1e10;
  for( i=0; i<num_p; i++ ){
    if(graph->raw_data[i*num_curves].x > max_x)
      max_x = graph->raw_data[i*num_curves].x;
    if(graph->raw_data[i*num_curves].x < min_x)
      min_x = graph->raw_data[i*num_curves].x;

    for(j=0;j<num_curves;j++){
      if(graph->raw_data[i*num_curves+j].y > max_y)
        max_y = graph->raw_data[i*num_curves+j].y;
      if(graph->raw_data[i*num_curves+j].y < min_y)
        min_y = graph->raw_data[i*num_curves+j].y;
    }
  }
  graph->max_x = max_x;
  graph->max_y = max_y;
  graph->min_x = min_x;
  graph->min_y = min_y;

  /* now find the max and min of the total energy curve */
  min_y = 1e10;
  max_y = -1e10;
  for( i=0; i<num_p; i++ ){
    if(tot_E->raw_data[i].y > max_y)
      max_y = tot_E->raw_data[i].y;
    if(tot_E->raw_data[i].y < min_y)
      min_y = tot_E->raw_data[i].y;
  }
  tot_E->max_x = max_x;
  tot_E->max_y = max_y;
  tot_E->min_x = min_x;
  tot_E->min_y = min_y;


  graph->data = (point_type2D *)D_CALLOC((unsigned)num_p*num_curves,sizeof(point_type2D));
  tot_E->data = (point_type2D *)D_CALLOC((unsigned)num_p,sizeof(point_type2D));
  if( !graph->data || !tot_E->data ) fatal("Can't get space for graph data.");

  /* that's it */
}


/****************************************************************************
 *
 *                   Procedure new_walsh_graph
 *
 * Arguments: none
 *
 * Returns: none
 *
 * Action: does everything to get space for and read in a new walsh_graph
 *
 ****************************************************************************/
void new_walsh_graph(char *filename)
{
  char failed = 0;
  char *theinline;
  char file_name[80];
  FILE *infile;

  /* set up a new object to hold the thing */
  makenewobject();
  whichobj = head->obj;

  /* now build the prop_graph primitive */
  whichobj->prim = (prim_type *)D_CALLOC(1,sizeof(prim_type));
  if( !whichobj->prim )fatal("Can't get space for walsh_graph primitive.");
  whichobj->prim->which = WALSH_GRAPH;

  whichobj->prim->walsh_graph =
    (walsh_graph_type *)D_CALLOC(1,sizeof(walsh_graph_type));
  if( !whichobj->prim->walsh_graph )
    fatal("Can't get space for walsh_graph.");

  whichobj->prim->walsh_graph->the_data =
    (graph_type *)D_CALLOC(1,sizeof(graph_type));
  if( !whichobj->prim->walsh_graph->the_data )
    fatal("Can't get space for walsh_graph->the_data.");

  whichobj->prim->walsh_graph->total_E =
    (graph_type *)D_CALLOC(1,sizeof(graph_type));
  if( !whichobj->prim->walsh_graph->total_E )
    fatal("Can't get space for walsh_graph->total_E.");


  display("Look in the xterm...");

#ifndef USING_THE_MAC
  if( !filename ){
#ifndef USE_READLINE
    printf("Enter the file name containing the Walsh data: ");
    scanf("%s\n",file_name);
#else
    theinline= readline("Enter the file name containing the Walsh data: ");
    add_history(theinline);
    if( theinline ){
      sscanf(theinline,"%s",file_name);
      free(theinline);
    } else {
      error("Bad file name");
      file_name[0] = 0;
    }
#endif

  }else{
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
    read_walsh_data(infile,whichobj->prim->walsh_graph->the_data,
                    whichobj->prim->walsh_graph->total_E);
  }

  /* check to see if any curves were actually read in.... */
  if(!whichobj->prim->walsh_graph->the_data->num_curves || failed){
    /* no... free the memory that we asked for */
    D_FREE(whichobj->prim->walsh_graph->the_data);
    D_FREE(whichobj->prim->walsh_graph);
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

#ifdef INTERACTIVE_USE
    /* make the button window */
    build_walsh_button_window(&button_wins,whichobj->prim->walsh_graph);
    whichobj->prim->but_win = button_wins;
#endif
  }
}



/****************************************************************************
 *
 *                   Procedure draw_walsh_graph
 *
 * Arguments: prim: pointer to primitive_type
 *             obj: pointer to object_type
 *
 * Returns: none
 *
 * Action: Draws in a Walsh diagram
 *   This takes care of the user specified scaling and translating.
 *
 ****************************************************************************/
void draw_walsh_graph(prim_type *prim,object_type *obj)
{
  static XPoint *points=0;
  static point_type2D *fpoints = 0;
  static int num_points=0;
  int i,j;
  graph_type *the_graph,*total_E;
  walsh_graph_type *walsh_graph;
  point_type2D origin,dim;
  char numstring[20];
  float xloc,yloc;
  float inv_xscale,inv_yscale;
  float xref,yref;
  float xval,yval;

  g_change_color(0);
  g_change_linewidth(1);
  walsh_graph = prim->walsh_graph;
  the_graph = walsh_graph->the_data;
  total_E = walsh_graph->total_E;

  /* check to see if we need to re-determine the location of tic marks */
  if( the_graph->old_max_x != the_graph->max_x ||
     the_graph->old_max_y != the_graph->max_y ||
     the_graph->old_min_x != the_graph->min_x ||
     the_graph->old_min_y != the_graph->min_y ){

    preprocess_graph_data(the_graph);
    the_graph->old_max_x = the_graph->max_x;
    the_graph->old_max_y = the_graph->max_y;
    the_graph->old_min_x = the_graph->min_x;
    the_graph->old_min_y = the_graph->min_y;
  }

  total_E->max_x = the_graph->max_x;
  total_E->min_x = the_graph->min_x;
  if( total_E->old_max_x != total_E->max_x ||
     total_E->old_max_y != total_E->max_y ||
     total_E->old_min_x != total_E->min_x ||
     total_E->old_min_y != total_E->min_y ){

    preprocess_graph_data(total_E);
    total_E->old_max_x = total_E->max_x;
    total_E->old_max_y = total_E->max_y;
    total_E->old_min_x = total_E->min_x;
    total_E->old_min_y = total_E->min_y;
  }

  /* check to see if we need memory for the Xpoints */
  if( !points || num_points < the_graph->num_p ){
    if( points ) free(points);
    num_points = the_graph->num_p;
    points = (XPoint *)calloc(num_points,sizeof(XPoint));
    if( !points ) fatal("Can't allocate memory for point storage in draw_graph.");
    if( fpoints ) free(fpoints);
    fpoints = (point_type2D *)calloc(num_points,sizeof(point_type2D));
    if( !fpoints ) fatal("Can't allocate memory for fpoint storage in draw_graph.");
  }

  /* inverse scaling terms */
  inv_xscale = (the_graph->max_x - the_graph->min_x) / DEF_GRAPH_X;
  inv_yscale = (the_graph->max_y - the_graph->min_y) / DEF_GRAPH_Y;

  /* determine the location of the origin on screen */
  origin.x = obj->cent.x + obj->trans.x + g_xmax / 2;
  origin.y = obj->cent.y - obj->trans.y + g_ymax / 2;
  dim.x = DEF_GRAPH_X * obj->scale.x;
  dim.y = DEF_GRAPH_Y * obj->scale.y;

  /* draw in a box surrounding the plot */
  g_change_linestyle(0);
  g_rectangle(origin.x,origin.y,dim.x,dim.y);

  /* find the point in data space which will appear at the origin */
  xref = the_graph->min_x*DEF_GRAPH_X/(the_graph->max_x-the_graph->min_x);
  yref = the_graph->min_y*DEF_GRAPH_Y/(the_graph->max_y-the_graph->min_y);

  /* determine the bounding box for this object */
  localmin.x = obj->bmin.x = origin.x;
  localmin.y = obj->bmin.y = origin.y - dim.y;
  localmax.x = obj->bmax.x = origin.x + dim.x;
  localmax.y = obj->bmax.y = origin.y;

  /* put in the tic marks now... */
  if( the_graph->do_x_tics ){
    for(i=0;i<(int)rint(the_graph->num_tics_x);i++){
      xloc = origin.x + obj->scale.x * (the_graph->tic_start_x +
                                        i * the_graph->tic_sep_x - xref);
      g_line(xloc,origin.y+obj->scale.y*TIC_DIM,xloc,origin.y);

      /*****

        do the labels

        *****/
      xval = (the_graph->tic_start_x + i*the_graph->tic_sep_x)*inv_xscale;
      if( fabs(xval) < 1e-12 ) xval = 0.0;
      sprintf(numstring,"%lg",xval);
      yloc = origin.y+obj->scale.y*TIC_DIM*(1.1);

      g_center_text(xloc,yloc,numstring);
    }
  }
  /* Y tics */
  if( the_graph->do_y_tics ){
    for(i=0;i<(int)rint(the_graph->num_tics_y);i++){
      yloc = origin.y + obj->scale.y * (yref - the_graph->tic_start_y -
                                        i * the_graph->tic_sep_y);
      g_line(origin.x-obj->scale.x*TIC_DIM,yloc,origin.x,yloc);

      yval = (the_graph->tic_start_y + i*the_graph->tic_sep_y)*inv_yscale;
      if( fabs(yval) < 1e-12 ) yval = 0.0;
      sprintf(numstring,"%lg",yval);
      xloc = origin.x-obj->scale.x*TIC_DIM*1.1;

      g_right_text(xloc,yloc,numstring);
    }
  }

  /* put in legends if they are needed */
  if( the_graph->xlegend[0] != 0 && the_graph->do_x_tics ){
    xloc = origin.x + dim.x/2;
    yloc = origin.y+obj->scale.y*TIC_DIM*(1.1);

    g_xlegend(xloc,yloc,the_graph->xlegend);
    obj->bmax.y += 2.1*PS_options.fontsize;

  }

  if( the_graph->ylegend[0] != 0 && the_graph->do_y_tics ){
    yloc = origin.y - dim.y/2;
    xloc = origin.x-obj->scale.x*TIC_DIM*5.0;

    g_ylegend(xloc,yloc,the_graph->ylegend);
    obj->bmin.x -= 3.0*PS_options.fontsize;
  }

  /* the title */
  if( the_graph->title[0] != 0 && the_graph->do_title){
    xloc = origin.x + dim.x/2;
    yloc = origin.y - dim.y;

    g_title(xloc,yloc,the_graph->title);
  }

  /* now do the plot itself */
  if( the_graph->curves_to_display[0] ){
    for(i=0;i<the_graph->num_curves;i++){
      g_change_linestyle(the_graph->styles[0]);

      /*******

        Bug fix: top band not being drawn
         05/22/95 gL

      ********/
      for(j=0;j<the_graph->num_p;j++){
        points[j].x = fpoints[j].x = origin.x +
          (the_graph->data[j*the_graph->num_curves + i].x - xref)*obj->scale.x;
        points[j].y = fpoints[j].y = origin.y +
          (yref - the_graph->data[j*the_graph->num_curves + i].y)*obj->scale.y;
      }
      g_lines(points,fpoints,the_graph->num_p,0);

      /* set the line style back to the default value */
      g_change_linestyle(0);
    }
  }

  /* draw in the total energy curve */
  if( the_graph->curves_to_display[1] ){


    /* find the point in data space which will appear at the origin */
    xref = total_E->min_x*DEF_GRAPH_X/(total_E->max_x-total_E->min_x);
    yref = total_E->min_y*DEF_GRAPH_Y/(total_E->max_y-total_E->min_y);
    inv_yscale = (total_E->max_y - total_E->min_y) / DEF_GRAPH_Y;
    /* Y tics */
    if( total_E->do_y_tics ){
      for(i=0;i<(int)rint(total_E->num_tics_y);i++){
        yloc = origin.y + obj->scale.y * (yref - total_E->tic_start_y -
                                          i * total_E->tic_sep_y);
        g_line(origin.x+dim.x+obj->scale.x*TIC_DIM,yloc,origin.x+dim.x,yloc);

        yval = (total_E->tic_start_y + i*total_E->tic_sep_y)*inv_yscale;
        if( fabs(yval) < 1e-12 ) yval = 0.0;
        sprintf(numstring,"%lg",yval);
        xloc = origin.x+dim.x+obj->scale.x*TIC_DIM*1.1;

        g_left_text(xloc,yloc,numstring);
        obj->bmax.x += 3.0*PS_options.fontsize;
      }
    }

    g_change_linestyle(the_graph->styles[1]);
    for(j=0;j<total_E->num_p;j++){
      points[j].x = fpoints[j].x = origin.x +
        (total_E->data[j].x - xref)*obj->scale.x;
      points[j].y = fpoints[j].y = origin.y +
        (yref - total_E->data[j].y)*obj->scale.y;
    }
    g_lines(points,fpoints,total_E->num_p,0);

    /* set the line style back to the default value */
    g_change_linestyle(0);
  }
}
