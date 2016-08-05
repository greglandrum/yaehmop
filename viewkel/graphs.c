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
    bounding boxes adjusted.
  26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)
  06.03.99 gL:
     support for printing curve names added.
  11.05.99 gL:
     added sarcastic comment to find_tic_sep
***/


/********

  this has got the stuff for dealing with graphs....

*********/
#include "viewkel.h"

/****************************************************************************

                  TODO list

****************************************************************************/

/****************************************************************************
 *
 *                   Procedure find_tic_sep
 *
 * Arguments: graph: pointer to graph_type
 *
 * Returns: none
 *
 * Action:  Determines a "reasonable" tic mark separation for the data set
 *    stored in 'graph.
 *
 *  NOTE:  This is all voodoo magic.
 *
 ****************************************************************************/
void find_tic_sep(graph_type *graph)
{
  float range,norm,tics,tic_sep,log_range;
  float xmax,xmin,ymax,ymin;

  /* first do the x tics */
  range = fabs(graph->min_x-graph->max_x);

  log_range = log10(range);

  norm = pow(10.0,log_range-(float)((log_range >= 0.0 ) ? (int)log_range :
                                    ((int)log_range-1)));
  if (norm <= 2)
    tics = 0.2;
  else if (norm <= 5)
    tics = 0.5;
  else tics = 1.0;
  tic_sep = tics * pow(10.0,(float)((log_range >= 0.0 ) ? (int)log_range :
                                    ((int)log_range-1)));

  graph->tic_sep_x = tic_sep;

  /* figure out the max and min values spanned by the tic marks */
  xmin = tic_sep * floor(graph->min_x/tic_sep);
  xmax = tic_sep * ceil(graph->max_x/tic_sep);
  graph->num_tics_x = 1 + (xmax - xmin) / tic_sep;

  /* if the first tic mark is outside the range the user specified, skip it! */
  if( (int)(10000.0*xmin) < (int)(10000.0*graph->min_x) ){
    graph->tic_start_x = xmin + tic_sep;
    xmin += tic_sep;
    graph->num_tics_x--;
  }
  else{
    graph->tic_start_x = xmin;
  }
  /* same deal with the last tic mark */
  if( (int)(10000.0*xmax) > (int)(10000.0*graph->max_x) ){
    graph->num_tics_x--;
    xmax -= tic_sep;
  }
  if( (int)(10000.0*xmax) > (int)(10000.0*graph->max_x) ){
    graph->num_tics_x--;
    xmax -= tic_sep;
  }

  /******

    repeat the process for the y data

  *******/
  range = fabs(graph->min_y-graph->max_y);
  log_range = log10(range);
  norm = pow(10.0,log_range-(float)((log_range >= 0.0 ) ? (int)log_range :
                                    ((int)log_range-1)));
  if (norm <= 2)
    tics = 0.2;
  else if (norm <= 5)
    tics = 0.5;
  else tics = 1.0;
  tic_sep = tics * pow(10.0,(float)((log_range >= 0.0 ) ? (int)log_range :
                                    ((int)log_range-1)));
  graph->tic_sep_y = tic_sep;
  ymin = tic_sep * floor(graph->min_y/tic_sep);
  ymax = tic_sep * ceil(graph->max_y/tic_sep);
  graph->num_tics_y = 1 + (ymax - ymin) / tic_sep;
  if( ymin < graph->min_y ){
    graph->tic_start_y = ymin + tic_sep;
    graph->num_tics_y--;
    ymin += tic_sep;
  }
  else{
    graph->tic_start_y = ymin;
  }
  if( ymax > graph->max_y ){
    graph->num_tics_y--;
  }
/*  if( ymax > graph->max_y ){
    graph->num_tics_y--;
  }
*/
}


/****************************************************************************
 *
 *                   Procedure count_curves
 *
 * Arguments: instring: pointer to type char
 *          num_curves: pointer to type int
 *
 * Returns: none
 *
 * Action:  Figures out the number of curves in 'instring.
 *
 *  This is just the number of entries in the line minus 1.
 *
 ****************************************************************************/
void count_curves(char *instring,int *num_curves)
{
  char instring_copy[240];
  char *foo;

  strcpy(instring_copy,instring);

  *num_curves = 0;
  foo = (char *)strtok(instring_copy," ");

  if( !foo ) return;

  while( strtok(0," ") ) (*num_curves)++;
}




/****************************************************************************
 *
 *                   Procedure read_graph_data
 *
 * Arguments: infile: pointer to type FILE
 *               graph: a pointer to graph_type
 *            sideways: a char
 *
 * Returns: none
 *
 * Action:  Reads the data out of 'infile.
 *       This reads until it hits a line beginning with a # mark.
 *
 *   If 'sideways is non-zero then it is assumed that the data columns
 *    occur first in the file, followed by the independant variable.
 *
 ****************************************************************************/
void read_graph_data(FILE *infile,graph_type *graph,char sideways,char parse_header)
{
  char instring[MAX_STR_LEN],foostring[MAX_STR_LEN];
  int num_curves,num_names;
  int i,j;
  int done;
  int max_p,num_p;
  float xval,yval;
  float min_x,min_y,max_x,max_y;

  if( !graph )FATAL_BUG("No graph passed to read_graph_data.");

  if( skipcomments(infile,instring) < 0){
    error("That file is empty!");
    display("Too Bad....");
    return;
  }

  /*******

    make sure that we've read past any header information

    write some code to process this stuff!

  ********/
  while(instring[0] == '#'){
    if(skipcomments(infile,instring)<0){
      error("That file is empty!");
      display("Too Bad....");
      return;
    }
  }
  /* figure out how many curves there are in the file */
  count_curves(instring,&num_curves);

  if(num_curves < 1){
    error("That data file has less than one data set in it... check it please.");
    display("Yikes!");
    return;
  }

  /* parse the header information */
  if( parse_header ){
    rewind(infile);
    graph->curve_names = (char *)D_CALLOC(NORMAL_STR_LEN*num_curves,sizeof(char));
    num_names = 1;
    skipcomments(infile,instring);
    while(instring[0] == '#'){
      skipcomments(infile,instring);
      if( strstr(instring,"Curve") ){
        sscanf(instring,"%s %s",foostring,
               &(graph->curve_names[num_names*NORMAL_STR_LEN]));
        num_names++;
      } else if (strstr(instring,"Title")){
        sscanf(instring,"%s %s",foostring,graph->curve_names);
      }
    }
  }

  /* get space to keep track of which curves should be displayed */
  graph->curves_to_display = (char *)D_CALLOC(num_curves,sizeof(char));
  if( !graph->curves_to_display )
    fatal("Can't get space for curves_to_display.");

  /* get space to store linestyles */
  graph->styles = (char *)D_CALLOC(num_curves,sizeof(char));
  if( !graph->styles )fatal("Can't get space for curve styles.");

  /* get space to store fill toggles */
  graph->fills = (char *)D_CALLOC(num_curves,sizeof(char));
  if( !graph->fills )fatal("Can't get space for curve fills.");

  /* initialize the styles */
  for( i=0;i<num_curves; i++){
    graph->styles[i] = i;
  }

  /* default to showing just the first curve */
  graph->curves_to_display[0] = 1;

  graph->do_x_tics = 1;
  graph->do_y_tics = 1;

  /*******
    get initial space for the data (we don't really know how much there
    will be yet, but we'll deal with that.
  *******/
  max_p = 100;
  num_p = 0;
  graph->raw_data = (point_type2D *)D_CALLOC(max_p*num_curves,sizeof(point_type2D));
  if( !graph->raw_data ){
    error("Can't get space for graph data....");
    display("Oh well...");
    return;
  }

  /* read in data until either EOF or a # mark is hit */
  done = 1;
  while(done >= 0 ){

    /* now use strtok to read out space delimited numbers */

    if( !sideways ){
      sscanf((const char *)strtok(instring," "),"%lf",&xval);
      for(i=0;i<num_curves;i++){
        graph->raw_data[num_p*num_curves+i].x = xval;
        sscanf((const char *)strtok(0," "),"%lf",&(graph->raw_data[num_p*num_curves+i].y));
      }
    }
    else{
      /*******
        assume that the x coordinate is at the end of the line for properties
        data....
      *******/
      sscanf((const char *)strtok(instring," "),"%lf",&(graph->raw_data[num_p*num_curves].x));
      for(i=1;i<num_curves;i++){
        sscanf((const char *)strtok(0," "),"%lf",&(graph->raw_data[num_p*num_curves+i].x));
      }
      sscanf((const char *)strtok(0," "),"%lf",&yval);
      for(i=0;i<num_curves;i++){
        graph->raw_data[num_p*num_curves+i].y = yval;
      }
    }
    num_p++;
    /* if we need more space, get it now */
    if( num_p == max_p ){
      max_p += 100;
      graph->raw_data = (point_type2D *)D_REALLOC((void *)graph->raw_data,
                                            (unsigned)max_p*num_curves*
                                            sizeof(point_type2D));
      if( !graph->raw_data ){
        error("Can't D_REALLOCate space for graph data....");
        display("Oh well...");
        return;
      }
    }
    done = skipcomments(infile,instring);

    /* check to see if this line begins with a # */
    if( done >= 0 && instring[0] == '#' ) done = -1;
  }

  graph->num_p = num_p;
  graph->num_curves = num_curves;

  graph->data = (point_type2D *)D_CALLOC((unsigned)num_p*num_curves,sizeof(point_type2D));
  if( !graph->data ) fatal("Can't get space for graph data.");

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

  /* that's it */
}


/****************************************************************************
 *
 *                   Procedure preprocess_graph_data
 *
 * Arguments: graph: a pointer to graph_type
 *
 * Returns: none
 *
 * Action: This does any preprocessing which is needed on the graph data.
 *
 ****************************************************************************/
void preprocess_graph_data(graph_type *graph)
{
  float xscale,yscale;
  int i,j;

  /* find the tic mark spacing */
  find_tic_sep(graph);

  /* scale all the points to fit in a DEF_GRAPH_X x DEF_GRAPH_Y box */
  xscale = DEF_GRAPH_X/(graph->max_x - graph->min_x);
  yscale = DEF_GRAPH_Y/(graph->max_y - graph->min_y);


  /**********

    now scale the points, clipping any which are outside the range

  **********/
  for( i=0; i<graph->num_p; i++ ){
    for(j=0;j<graph->num_curves;j++){
      if( graph->raw_data[i*graph->num_curves+j].x > graph->max_x ){
        graph->data[i*graph->num_curves+j].x = graph->max_x * xscale;
      }
      else    if( graph->raw_data[i*graph->num_curves+j].x < graph->min_x ){
        graph->data[i*graph->num_curves+j].x = graph->min_x * xscale;
      }
      else{
        graph->data[i*graph->num_curves+j].x = graph->raw_data[i*graph->num_curves+j].x *
          xscale;
      }
      if( graph->raw_data[i*graph->num_curves+j].y > graph->max_y ){
        graph->data[i*graph->num_curves+j].y = graph->max_y * yscale;
      }
      else    if( graph->raw_data[i*graph->num_curves+j].y < graph->min_y ){
        graph->data[i*graph->num_curves+j].y = graph->min_y * yscale;
      }
      else{
        graph->data[i*graph->num_curves+j].y = graph->raw_data[i*graph->num_curves+j].y *
          yscale;
      }
    }
  }

  /* adjust the tic mark spacing to fit the new scaling */
  graph->tic_sep_x *= xscale;
  graph->tic_sep_y *= yscale;
  graph->tic_start_x *= xscale;
  graph->tic_start_y *= yscale;
}


/****************************************************************************
 *
 *                   Procedure new_graph
 *
 * Arguments: filename: pointer to type char (optional)
 *
 * Returns: none
 *
 * Action: does everything to get space for and read in a new graph
 *
 ****************************************************************************/
void new_graph(char *filename)
{
  char file_name[80];
  char *theinline;
  char failed;
  FILE *infile;

  failed = 0;
  /* set up a new object to hold the graph */
  makenewobject();
  whichobj = head->obj;

  /* now build the graph primitive */
  whichobj->prim = (prim_type *)D_CALLOC(1,sizeof(prim_type));
  if( !whichobj->prim )fatal("Can't get space for graph primitive.");
  whichobj->prim->which = GRAPH;

  whichobj->prim->graph = (graph_type *)D_CALLOC(1,sizeof(graph_type));
  if( !whichobj->prim->graph )
    fatal("Can't get space for graph.");


  if( !filename ){
    display("Look in the xterm...");
#ifndef USE_READLINE
    printf("Enter the file name containing the graph data: ");
    scanf("%s",file_name);
#else
    theinline= readline("Enter the file name containing the graph data: ");
    add_history(theinline);
    if( theinline ){
      sscanf(theinline,"%s",file_name);
      free(theinline);
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

  if( !failed ){
    read_graph_data(infile,whichobj->prim->graph,0,1);
    strcpy(whichobj->prim->graph->filename,file_name);
  }

  /* check to see if any curves were actually read in.... */
  if(!whichobj->prim->graph->num_curves || failed ){
    /* no... free the memory that we asked for */
    D_FREE(whichobj->prim->graph);
    D_FREE(whichobj->prim);
    D_FREE(whichobj);
    whichobj=0;
    head->obj = 0;
    head = head->next;
  }
  else{
    whichobj->scale.x=whichobj->scale.y=whichobj->scale.z=2.0*GRAPHICS_SCALE;
    whichobj->cent.x=-180;whichobj->cent.y=180;
    whichobj->cent.z=0*GRAPHICS_SCALE;
    whichobj->trans.x=0;whichobj->trans.y=0;
    whichobj->trans.z=0;

    /* make the button window */
#ifdef INTERACTIVE_USE
    build_graph_button_window(&button_wins,whichobj->prim->graph);
    whichobj->prim->but_win = button_wins;
#endif

  }

}


/****************************************************************************
 *
 *                   Procedure draw_graph
 *
 * Arguments: prim: pointer to primitive_type
 *             obj: pointer to object_type
 *
 * Returns: none
 *
 * Action: Draws in a graph.
 *   This takes care of the user specified scaling and translating.
 *
 ****************************************************************************/
void draw_graph(prim_type *prim,object_type *obj)
{
  static XPoint *points=0;
  static point_type2D *fpoints=0;
  static int num_points=0;
  static char integ_last_time=0;
  int i,j;
  graph_type *the_graph;
  point_type2D origin,dim;
  char numstring[20];
  float xloc,yloc;
  float inv_xscale,inv_yscale;
  float scaled_fermi;
  float xscale,yscale;
  int lineshown;
  float xref,yref;
  float xval,yval;
  float *max_x,*max_y,*min_x,*min_y;
  float *old_max_x,*old_max_y,*old_min_x,*old_min_y;

  /* check to see what type of graph this is */
  if(prim->prop_graph){
    the_graph = prim->prop_graph->the_data;
    max_x = &(prim->prop_graph->max_x);
    max_y = &(prim->prop_graph->max_y);
    min_x = &(prim->prop_graph->min_x);
    min_y = &(prim->prop_graph->min_y);
    old_max_x = &(prim->prop_graph->old_max_x);
    old_max_y = &(prim->prop_graph->old_max_y);
    old_min_x = &(prim->prop_graph->old_min_x);
    old_min_y = &(prim->prop_graph->old_min_y);
  } else{
    if( prim->graph )  the_graph = prim->graph;
    else if(prim->band_graph) the_graph = prim->band_graph->the_data;
    else if(prim->walsh_graph) the_graph = prim->walsh_graph->the_data;
    else FATAL_BUG("Bogus primitive passed to draw_graph.");
    max_x = &(the_graph->max_x);
    max_y = &(the_graph->max_y);
    min_x = &(the_graph->min_x);
    min_y = &(the_graph->min_y);
    old_max_x = &(the_graph->old_max_x);
    old_max_y = &(the_graph->old_max_y);
    old_min_x = &(the_graph->old_min_x);
    old_min_y = &(the_graph->old_min_y);
  }

  /* check to see if we need to re-determine the location of tic marks */
  if( *old_max_x != *max_x ||
     *old_max_y != *max_y ||
     *old_min_x != *min_x ||
     *old_min_y != *min_y ||
     (prim->prop_graph && (prim->prop_graph->integs_for_tics && !integ_last_time))){

    /********

      if we're using the integration for the X axis, then we need to
      update the graph max and min values.

      *******/
    if( prim->prop_graph && prim->prop_graph->integs_for_tics ){
      /*******

        check to see if we've already updated them

        ********/
      if(!integ_last_time){
        *max_x = prim->prop_graph->the_integration->max_x;
        *max_y = prim->prop_graph->the_integration->max_y;
        *min_x = prim->prop_graph->the_integration->min_x;
        *min_y = prim->prop_graph->the_integration->min_y;

        /******

          make sure that we aren't gonna get horked when we display COOP
          integrations by displaying COOP's on the same scale as the
          integration

          *******/
        if( prim->prop_graph->type == COOP_PROP ){
          prim->prop_graph->the_data->max_x = *max_x;
          prim->prop_graph->the_data->min_x = *min_x;
        }
        prim->prop_graph->the_data->max_y = *max_y;
        prim->prop_graph->the_data->min_y = *min_y;

        integ_last_time = 1;
      } /* if(!integ_last_time) */
      else{
        /*******

          if we got this far, then it means that the user updated the
          max and/or min values, so we need to change the integration
          max and min values to reflect the change.

          ********/
        prim->prop_graph->the_integration->max_x = *max_x;
        prim->prop_graph->the_integration->min_x = *min_x;
        prim->prop_graph->the_integration->max_y = *max_y;
        prim->prop_graph->the_integration->min_y = *min_y;

        if( prim->prop_graph->type == COOP_PROP ){
          prim->prop_graph->the_data->max_x = *max_x;
          prim->prop_graph->the_data->min_x = *min_x;
        }
        prim->prop_graph->the_data->max_y = *max_y;
        prim->prop_graph->the_data->min_y = *min_y;
      }
    } /*  if( prim->prop_graph && prim->prop_graph->integs_for_tics ) */
    else{
      /******

        reset the max and min values to those of the
        data...

        ******/
      if( prim->prop_graph ){
        if( integ_last_time ){
          *max_x = prim->prop_graph->the_data->max_x;
          *min_x = prim->prop_graph->the_data->min_x;
          *max_y = prim->prop_graph->the_data->max_y;
          *min_y = prim->prop_graph->the_data->min_y;
        }
        if( prim->prop_graph->the_integration){
          if( prim->prop_graph->type == COOP_PROP ){
            prim->prop_graph->the_integration->max_x = *max_x;
            prim->prop_graph->the_integration->min_x = *min_x;
          }
          prim->prop_graph->the_integration->max_y = *max_y;
          prim->prop_graph->the_integration->min_y = *min_y;
        }
      }
      the_graph->max_x = *max_x;
      the_graph->min_x = *min_x;
      the_graph->max_y = *max_y;
      the_graph->min_y = *min_y;

      integ_last_time = 0;
    }
    preprocess_graph_data(the_graph);
    if( prim->prop_graph && prim->prop_graph->the_integration){
      preprocess_graph_data(prim->prop_graph->the_integration);
    }

/*
    the_graph->old_max_x = the_graph->max_x;
    the_graph->old_max_y = the_graph->max_y;
    the_graph->old_min_x = the_graph->min_x;
    the_graph->old_min_y = the_graph->min_y;
*/
    *old_max_x = *max_x;
    *old_max_y = *max_y;
    *old_min_x = *min_x;
    *old_min_y = *min_y;

  }


  /* check to see if we need memory for the Xpoints */
  if( !points || num_points < the_graph->num_p ){
    if( points ) free(points);
    num_points = the_graph->num_p;
    points = (XPoint *)calloc(num_points+2,sizeof(XPoint));
    if( !points ) fatal("Can't allocate memory for point storage in draw_graph.");
    if( fpoints ) free(fpoints);
    fpoints = (point_type2D *)calloc(num_points+2,sizeof(point_type2D));
    if( !fpoints ) fatal("Can't allocate memory for fpoint storage in draw_graph.");
  }

  /* determine the location of the origin on screen */
  origin.x = obj->cent.x + obj->trans.x + g_xmax / 2;
  origin.y = obj->cent.y - obj->trans.y + g_ymax / 2;
  dim.x = DEF_GRAPH_X * obj->scale.x;
  dim.y = DEF_GRAPH_Y * obj->scale.y;

  /* scaling terms */
  inv_xscale = (the_graph->max_x - the_graph->min_x) / DEF_GRAPH_X;
  inv_yscale = (the_graph->max_y - the_graph->min_y) / DEF_GRAPH_Y;
  xscale = 1/inv_xscale;
  yscale = 1/inv_yscale;

  /* find the point in data space which will appear at the origin */
  xref = the_graph->min_x * xscale;
  yref = the_graph->min_y * yscale;

  /* determine the bounding box for this object */
  localmin.x = obj->bmin.x = origin.x;
  localmin.y = obj->bmin.y = origin.y - dim.y;
  localmax.x = obj->bmax.x = origin.x + dim.x;
  localmax.y = obj->bmax.y = origin.y;

  /* put in the tic marks now... */
  if( the_graph->do_x_tics && (!prim->prop_graph || !prim->prop_graph->integs_for_tics)){
    for(i=0;i<(int)rint(the_graph->num_tics_x);i++){
      xloc = origin.x + obj->scale.x * (the_graph->tic_start_x +
                                        i * the_graph->tic_sep_x - xref);
      g_line(xloc,origin.y+obj->scale.y*TIC_DIM,xloc,origin.y);

      localmin.y += obj->scale.y*TIC_DIM;

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

  /* Now do the title */
  if( the_graph->title[0] != 0 && the_graph->do_title){
    xloc = origin.x + dim.x/2;
    yloc = origin.y - dim.y;

    g_title(xloc,yloc,the_graph->title);
  }

  /* now do the plot itself */
  lineshown = 0;
  for(i=0;i<the_graph->num_curves;i++){
    if( the_graph->curves_to_display[i] ){
      /* adjust the line styles to get distinguishable graphs */
      g_change_linestyle(the_graph->styles[i]);

      for(j=0;j<the_graph->num_p;j++){
        points[j].x = fpoints[j].x = origin.x +
          (the_graph->data[j*the_graph->num_curves + i].x - xref)*obj->scale.x;

        points[j].y = fpoints[j].y = origin.y +
          (yref - the_graph->data[j*the_graph->num_curves + i].y)* obj->scale.y;
      }
      /*********

        Deal with the fact that DOS projections may be filled

      **********/
      if( prim->prop_graph && ((prim->prop_graph->type != COOP_PROP &&
         fill_projections && i != 0 ) || prim->prop_graph->the_data->fills[i])){

        /**********

          There are some complications in the filling if the
          end points of the graph are not at the same x level.

          deal with these potential problems by inserting 2 extra
          points into the curve.
        ************/
        if( *min_x == 0.0 || (*max_x >= 0 && *min_x < 0.0) ){
          points[the_graph->num_p].x = fpoints[the_graph->num_p].x = origin.x -
            xref*obj->scale.x;
          points[the_graph->num_p+1].x = fpoints[the_graph->num_p+1].x = origin.x -
            xref*obj->scale.x;
        }else if( *min_x > 0.0 ){
          points[the_graph->num_p].x = fpoints[the_graph->num_p].x = origin.x;
          points[the_graph->num_p+1].x = fpoints[the_graph->num_p+1].x = origin.x;
        }
        else{
          points[the_graph->num_p].x = fpoints[the_graph->num_p].x = origin.x +
            (*max_x - xref)*obj->scale.x;
          points[the_graph->num_p+1].x = fpoints[the_graph->num_p+1].x = origin.x +
            (*max_x - xref)*obj->scale.x;
        }

        points[the_graph->num_p].y = fpoints[the_graph->num_p].y = origin.y +
          (yref - the_graph->data[(the_graph->num_p-1)*the_graph->num_curves + i].y)
            * obj->scale.y;
        points[the_graph->num_p+1].y = fpoints[the_graph->num_p+1].y = origin.y +
          (yref - the_graph->data[i].y)* obj->scale.y;
        g_lines(points,fpoints,the_graph->num_p+2,1);
      }
      else{
        g_lines(points,fpoints,the_graph->num_p,0);
      }

      /* set the line style back to the default value */
      g_change_linestyle(0);
    }
  }

  /******

    put in the fermi level if this is a property graph and the user
    asked for it

  ******/
  if( prim->prop_graph && prim->prop_graph->show_fermi ){
    if( prim->prop_graph->Fermi_E > the_graph->min_y &&
       prim->prop_graph->Fermi_E < the_graph->max_y ){
      scaled_fermi = prim->prop_graph->Fermi_E * yscale;
      yloc = origin.y + (yref - scaled_fermi)*obj->scale.y;
      xloc = origin.x;
      /* use a dashed line */
      g_change_linestyle(1);
      /* draw a line */
      g_line(xloc,yloc,xloc+dim.x,yloc);
      /* set the line style back to the default */
      g_change_linestyle(0);
    }
  }

  /* check to see if we should draw a line through zero */
  if( the_graph->max_x > 0.0 && the_graph->min_x < 0.0 ){
    xloc = origin.x - xref * obj->scale.x;
    yloc = origin.y;
    g_line(xloc,yloc,xloc,yloc-dim.y);
  }


  /* check to see if we need to draw in an integration */
  if( prim->prop_graph && prim->prop_graph->the_integration ){
    the_graph = prim->prop_graph->the_integration;

    inv_xscale = (the_graph->max_x - the_graph->min_x) / DEF_GRAPH_X;
    inv_yscale = (the_graph->max_y - the_graph->min_y) / DEF_GRAPH_Y;
    xscale = 1/inv_xscale;
    yscale = 1/inv_yscale;
    /* find the point in data space which will appear at the origin */
    xref = the_graph->min_x * xscale;
    yref = the_graph->min_y * yscale;

    /*****

      if we are using the integration to provide the tic marks then
      show them now.

    *****/
    if( prim->prop_graph->the_data->do_x_tics && prim->prop_graph->integs_for_tics){
      for(i=0;i<(int)rint(the_graph->num_tics_x);i++){
        xloc = origin.x + obj->scale.x * (the_graph->tic_start_x +
                                          i * the_graph->tic_sep_x - xref);
        g_line(xloc,origin.y+obj->scale.y*TIC_DIM,xloc,origin.y);

        localmin.y += obj->scale.y*TIC_DIM;

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


    lineshown = 0;

    for(i=0;i<the_graph->num_curves;i++){
      if( the_graph->curves_to_display[i] ){
        /* adjust the line styles to get distinguishable graphs */
        g_change_linestyle(the_graph->styles[i]);

        for(j=0;j<the_graph->num_p;j++){
          points[j].x = fpoints[j].x = origin.x +
            (the_graph->data[j*the_graph->num_curves + i].x - xref)*obj->scale.x;
          points[j].y = fpoints[j].y = origin.y +
            (yref - the_graph->data[j*the_graph->num_curves + i].y)*obj->scale.y;
        }
        g_lines(points,fpoints,the_graph->num_p,0);

        g_change_linestyle(0);
      }
    }
    the_graph = prim->prop_graph->the_data;
  }

  /* draw in a box surrounding the plot */
  g_rectangle(origin.x,origin.y,dim.x,dim.y);

}

