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

/********

  this has got the stuff for dealing with band graphs....

  Created by greg Landrum June 1994
*********/

/***
  Recent Edit History:

  03.05.98 gL:
    bounding boxes adjusted.

  14.06.98 gL:
    added support for reading out Fermi levels (if there)
    added some EOF checking in read_band_data

  26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal).  Hopefully this will fix
     a nasty crashing bug though.
  18.10.98 gL:
     support for fatbands added.  draw_band_graph still needs to be modified
     a little bit to be able to deal with multiple fatbands plots in one
     file.  Uh, not that there is currently any software which generates such
     a data file... there will be though, yes there will... you'll see.

***/

#include "viewkel.h"

/****************************************************************************
 *
 *                   Procedure preprocess_band_graph_data
 *
 * Arguments: band_graph: a pointer to band_graph_type
 *
 * Returns: none
 *
 * Action: This does any preprocessing which is needed on the band graph data.
 *
 ****************************************************************************/
void preprocess_band_graph_data( band_graph_type *band_graph )
{
  float xscale,yscale;
  float oldmax_x,oldmin_x;
  int i,j,k;
#ifdef SUPPORT_FATBANDS
  int fatband_offset;
#endif

  /* save the x axis maxima and minima */
  oldmax_x = band_graph->the_data->max_x;
  oldmin_x = band_graph->the_data->min_x;

  /* find the tic mark spacing */
  find_tic_sep(band_graph->the_data);

  band_graph->the_data->max_x = oldmax_x;
  band_graph->the_data->min_x = oldmin_x;

  /* scale all the points to fit in a DEF_GRAPH_X x DEF_GRAPH_Y box */
  xscale = DEF_GRAPH_X/(band_graph->the_data->max_x - band_graph->the_data->min_x);
  yscale = DEF_GRAPH_Y/(band_graph->the_data->max_y - band_graph->the_data->min_y);

#ifdef SUPPORT_FATBANDS
  fatband_offset = band_graph->the_data->num_p*band_graph->the_data->num_curves;
#endif

  for( i=0; i<band_graph->the_data->num_p; i++ ){
    for(j=0;j<band_graph->the_data->num_curves;j++){
      band_graph->the_data->data[i*band_graph->the_data->num_curves+j].x =
        band_graph->the_data->raw_data[i*band_graph->the_data->num_curves+j].x*xscale;

      if( band_graph->the_data->raw_data[i*band_graph->the_data->num_curves+j].y >
         band_graph->the_data->max_y ){
        band_graph->the_data->data[i*band_graph->the_data->num_curves+j].y =
          band_graph->the_data->max_y * yscale;
#ifdef SUPPORT_FATBANDS
        for(k=1;k<=band_graph->num_fatbands;k++){
          band_graph->the_data->data[k*fatband_offset +
                                    i*band_graph->the_data->num_curves+j].y = 0;
        }
#endif
      }
      else if( band_graph->the_data->raw_data[i*band_graph->the_data->num_curves+j].y <
                 band_graph->the_data->min_y ){
        band_graph->the_data->data[i*band_graph->the_data->num_curves+j].y =
          band_graph->the_data->min_y * yscale;
#ifdef SUPPORT_FATBANDS
        for(k=1;k<=band_graph->num_fatbands;k++){
          band_graph->the_data->data[k*fatband_offset +
                                    i*band_graph->the_data->num_curves+j].y = 0;
        }
#endif

      }
      else{
        band_graph->the_data->data[i*band_graph->the_data->num_curves+j].y =
          band_graph->the_data->raw_data[i*band_graph->the_data->num_curves+j].y *
            yscale;
#ifdef SUPPORT_FATBANDS
        for(k=1;k<=band_graph->num_fatbands;k++){
          band_graph->the_data->data[k*fatband_offset +
                                    i*band_graph->the_data->num_curves+j].y =
            band_graph->the_data->raw_data[k*fatband_offset +
                                          i*band_graph->the_data->num_curves+j].y*yscale;
        }
#endif

      }
    }
  }

  /* change the x tic marks */
  band_graph->the_data->tic_sep_x =  band_graph->points_per_line;

  band_graph->the_data->num_tics_x = band_graph->num_special_points;


  /* adjust the tic mark spacing to fit the new scaling */
  band_graph->the_data->tic_sep_x *= xscale;
  band_graph->the_data->tic_sep_y *= yscale;
  band_graph->the_data->tic_start_x *= xscale;
  band_graph->the_data->tic_start_y *= yscale;
}



/****************************************************************************
 *
 *                   Procedure read_band_data
 *
 * Arguments: infile: pointer to FILE
 *        band_graph: pointer to band_graph_type
 *
 * Returns: none
 *
 * Action: Reads in all the data needed to construct a band structure
 *   graph from 'infile and stores it in the proper format in 'band_graph
 *
 ****************************************************************************/
void read_band_data(FILE *infile,band_graph_type *band_graph)
{
  graph_type *the_data;
  char instring[MAX_STR_LEN],foostring[MAX_STR_LEN];
  int eof_hit;
  int i,j,k;
  int points_per_line,num_special_points,num_orbs;
  int tot_num_points;
  float max_x,max_y,min_x,min_y;
#ifdef SUPPORT_FATBANDS
  int num_fatbands=0;
  int points_required;
#endif
  /* read out the number of special points */
  skipcomments(infile,instring);

  if( instring[0] == '#' || instring[1] == '#'){
#ifdef SUPPORT_FATBANDS
    upcase(instring);
    if( strstr(instring,"FATBANDS") ){
      sscanf(instring,"%s %d",foostring,&num_fatbands);
      band_graph->num_fatbands = num_fatbands;
    }
#endif
    skipcomments(infile,instring);
  }
  sscanf(instring,"%d",&num_special_points);
  band_graph->num_special_points = num_special_points;

  /* read out the number of k-points per symmetry line */
  skipcomments(infile,instring);
  sscanf(instring,"%d",&points_per_line);
  band_graph->points_per_line = points_per_line;

  /* read out the number of orbitals */
  skipcomments(infile,instring);
  sscanf(instring,"%d",&num_orbs);

  /* get all the memory that we need */
  band_graph->special_points = (special_point_type *)D_CALLOC(num_special_points,
                                                            sizeof(special_point_type));
  if( !band_graph->special_points ) fatal("Can't get memory for special_points.");

  band_graph->the_data = (graph_type *)D_CALLOC(1,sizeof(graph_type));
  if( !band_graph->the_data )fatal("Can't get memory for the_data.");

  the_data = band_graph->the_data;
  tot_num_points = the_data->num_p = (num_special_points-1)*points_per_line + 1;
  the_data->num_curves = num_orbs;
#ifdef SUPPORT_FATBANDS
  points_required = the_data->num_curves*the_data->num_p +
    num_fatbands*the_data->num_curves*the_data->num_p;
  the_data->raw_data = (point_type2D *)D_CALLOC(points_required,sizeof(point_type2D));

#else
  the_data->raw_data = (point_type2D *)D_CALLOC(the_data->num_curves*the_data->num_p,
                                              sizeof(point_type2D));
#endif
  if(!the_data->raw_data)fatal("Can't allocate the_data->raw_data");

  the_data->curves_to_display = (char *)D_CALLOC(1,sizeof(char));
  if( !the_data->curves_to_display )fatal("Can't get space for curves_to_display.");

  the_data->curves_to_display[0] = 1;

  the_data->styles = (char *)D_CALLOC(1,sizeof(char));
  if( !the_data->styles )fatal("Can't get space for styles.");

  the_data->styles[0] = 0;

  /* now read out the special point information */
  for(i=0;i<num_special_points;i++){
    skipcomments(infile,instring);
    sscanf(instring,"%s %lf %lf %lf",band_graph->special_points[i].name,
           &band_graph->special_points[i].loc.x,&band_graph->special_points[i].loc.y,
           &band_graph->special_points[i].loc.z);
  }


  /* now read out the actual information */
  for(i=0;i<tot_num_points;i++){
    for(j=0;j<num_orbs;j++){
      the_data->raw_data[i*num_orbs+j].x = i;
      eof_hit = skipcomments(infile,instring);
      if( eof_hit < 0 ){
        fprintf(stderr,"EOF hit whilst reading band file\n");
        fprintf(stderr,"We had read %d points of %d,\n",
                i,tot_num_points);
        fprintf(stderr,"\tand died on orbital %d of %d\n",
                j,num_orbs);
        error("Early EOF");
        the_data->num_curves = 0;
        return;
      }
#ifdef SUPPORT_FATBANDS
      sscanf((const char *)strtok(instring," "),"%lf",&the_data->raw_data[i*num_orbs+j].y);
      for(k=1;k<=num_fatbands;k++){
        sscanf((const char *)strtok(0," "),"%lf",
               &the_data->raw_data[k*num_orbs*tot_num_points + i*num_orbs+j].y);
      }
#else
      sscanf(instring,"%lf",&the_data->raw_data[i*num_orbs+j].y);
#endif
    }
  }

  /* check if there's a Fermi Energy present */
  while(eof_hit >= 0 && !strstr(instring,"FERMI_ENERGY")){
    eof_hit = skipcomments(infile,instring);
    upcase(instring);
  }
  if (!eof_hit){
    sscanf(instring,"%s %lf",foostring,&band_graph->Fermi_E);
    band_graph->show_fermi = 1;
  }


#ifdef SUPPORT_FATBANDS
  points_required = the_data->num_curves*the_data->num_p +
    num_fatbands*the_data->num_curves*the_data->num_p;
  the_data->data = (point_type2D *)D_CALLOC((unsigned)points_required,sizeof(point_type2D));
#else
  the_data->data = (point_type2D *)D_CALLOC((unsigned)tot_num_points*num_orbs,
                                          sizeof(point_type2D));
#endif
  if( !the_data->data ) fatal("Can't get space for graph data.");


  /* find the minimum and maximum values (to use in scaling) */
  min_y = min_x = 1e10;
  max_y = max_x = -1e10;
  for( i=0; i<tot_num_points; i++ ){
    if(the_data->raw_data[i*num_orbs].x > max_x)
      max_x = the_data->raw_data[i*num_orbs].x;
    if(the_data->raw_data[i*num_orbs].x < min_x)
      min_x = the_data->raw_data[i*num_orbs].x;

    for(j=0;j<num_orbs;j++){
      if(the_data->raw_data[i*num_orbs+j].y > max_y)
        max_y = the_data->raw_data[i*num_orbs+j].y;
      if(the_data->raw_data[i*num_orbs+j].y < min_y)
        min_y = the_data->raw_data[i*num_orbs+j].y;
    }
  }
  the_data->max_x = max_x;
  the_data->max_y = max_y;
  the_data->min_x = min_x;
  the_data->min_y = min_y;
}



/****************************************************************************
 *
 *                   Procedure new_band_graph
 *
 * Arguments: filename: pointer to char
 *
 * Returns: none
 *
 * Action: does everything to get space for and read in a new band_graph
 *
 ****************************************************************************/
void new_band_graph(char *filename)
{
  char failed = 0;
  char file_name[80];
  char *theinline;
  FILE *infile;
  graph_type *the_graph;

  /* set up a new object to hold the thing */
  makenewobject();
  whichobj = head->obj;

  /* now build the band_graph primitive */
  whichobj->prim = (prim_type *)D_CALLOC(1,sizeof(prim_type));
  if( !whichobj->prim )fatal("Can't get space for prop_graph primitive.");
  whichobj->prim->which = BAND_GRAPH;

  whichobj->prim->band_graph =
    (band_graph_type *)D_CALLOC(1,sizeof(band_graph_type));
  if( !whichobj->prim->band_graph )
    fatal("Can't get space for band_graph.");


#ifndef USING_THE_MAC
  if( !filename ){
    display("Look in the xterm...");
#ifndef USE_READLINE
    printf("Enter the file name containing the band structure data: ");
    scanf("%s",file_name);
#else
    theinline=
      readline("Enter the file name containing the band structure data: ");
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
    strcpy(whichobj->prim->band_graph->filename,file_name);
    read_band_data(infile,whichobj->prim->band_graph);
    preprocess_band_graph_data(whichobj->prim->band_graph);
  }

  /* check to see if any curves were actually read in.... */
  if(!whichobj->prim->band_graph->the_data->num_curves || failed){
    /* no... free the memory that we asked for */
    D_FREE(whichobj->prim->band_graph->the_data);
    D_FREE(whichobj->prim->band_graph);
    D_FREE(whichobj->prim);
    D_FREE(whichobj);
    whichobj=0;
    head->obj=0;
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
    the_graph = whichobj->prim->band_graph->the_data;
    strcpy(the_graph->ylegend,"Energy (eV)");

    the_graph->do_x_tics = 1;
    the_graph->do_y_tics = 1;

#ifdef SUPPORT_FATBANDS
    if( whichobj->prim->band_graph->num_fatbands ){
      whichobj->prim->band_graph->fatband_scale = 1.0;
      whichobj->prim->band_graph->fatband_fill = FATBANDS_LINE;
      whichobj->prim->band_graph->fatbands_on = 1;
    }
#endif

    build_band_button_window(&button_wins,whichobj->prim->band_graph);
    whichobj->prim->but_win = button_wins;

  }
}

/****************************************************************************
 *
 *                   Procedure draw_band_graph
 *
 * Arguments: prim: pointer to primitive_type
 *             obj: pointer to object_type
 *
 * Returns: none
 *
 * Action: Draws in a band structure graph.
 *   This takes care of the user specified scaling and translating.
 *
 ****************************************************************************/
void draw_band_graph(prim_type *prim,object_type *obj)
{
  static XPoint *points=0;
  static point_type2D *fpoints=0;
  static int num_points=0;
  int i,j;
  graph_type *the_graph;
  band_graph_type *band_graph;
  point_type2D origin,dim;
  char numstring[20];
  float xloc,yloc;
  float scaled_fermi;
  float xscale,yscale;
  float inv_xscale,inv_yscale;
  int max_str_len;
  float xref,yref;
  float yval;
#ifdef SUPPORT_FATBANDS
  int tot_num_p;
  float fat_width1,fat_width2;
  point_type2D fb_p1,fb_p2;
#endif

  the_graph = prim->band_graph->the_data;
  band_graph = prim->band_graph;

  /* check to see if we need to re-determine the location of tic marks */
  if( the_graph->old_max_x != the_graph->max_x || the_graph->old_max_y != the_graph->max_y ||
     the_graph->old_min_x != the_graph->min_x || the_graph->old_min_y != the_graph->min_y ){

    preprocess_band_graph_data(band_graph);
    the_graph->old_max_x = the_graph->max_x;
    the_graph->old_max_y = the_graph->max_y;
    the_graph->old_min_x = the_graph->min_x;
    the_graph->old_min_y = the_graph->min_y;
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

  xscale = 1/inv_xscale;
  yscale = 1/inv_yscale;

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

  /******

    Instead of tic marks on the X axis, put in vertical lines delineating
    the symmetry lines

    ******/
  if(the_graph->do_x_tics){
    for(i=0;i<the_graph->num_tics_x;i++){
      xloc = origin.x + obj->scale.x * i * the_graph->tic_sep_x,
      g_line(xloc,origin.y+obj->scale.y*TIC_DIM,xloc,origin.y-dim.y);

      /*****

        do the label

        *****/
      strcpy(numstring,band_graph->special_points[i].name);
      yloc = origin.y+obj->scale.y*TIC_DIM;

      g_center_text(xloc,yloc,numstring);
    }
  }

  /* Y tics */
  if( the_graph->do_y_tics ){
    max_str_len = 0;
    for(i=0;i<(int)rint(the_graph->num_tics_y);i++){
      yloc = origin.y + obj->scale.y * (yref - the_graph->tic_start_y -
                                        i * the_graph->tic_sep_y);
      g_line(origin.x-obj->scale.x*TIC_DIM,yloc,origin.x,yloc);

      yval = (the_graph->tic_start_y + i*the_graph->tic_sep_y)*inv_yscale;
      if( fabs(yval) < 1e-12 ) yval = 0.0;
      sprintf(numstring,"%lg",yval);
      xloc = origin.x-obj->scale.x*TIC_DIM-3;

      g_right_text(xloc,yloc,numstring);
    }
    obj->bmin.x -= obj->scale.x*TIC_DIM;
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


  /******

    put in the fermi level if the user asked for it

    ******/
  if(  prim->band_graph->show_fermi ){
    if( prim->band_graph->Fermi_E > the_graph->min_y &&
       prim->band_graph->Fermi_E < the_graph->max_y ){
      scaled_fermi = prim->band_graph->Fermi_E * yscale;
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

  /* now do the plot itself */
  if( the_graph->curves_to_display[0] ){
#ifdef SUPPORT_FATBANDS
    tot_num_p = the_graph->num_curves*the_graph->num_p;
#endif
    for(i=0;i<the_graph->num_curves;i++){
      g_change_linestyle(the_graph->styles[0]);
#ifdef SUPPORT_FATBANDS
      if( prim->band_graph->num_fatbands && prim->band_graph->fatbands_on ){
        g_change_linewidth(0.5);
        for(j=0;j<the_graph->num_p-1;j++){
          fat_width1 = the_graph->data[tot_num_p + j*the_graph->num_curves + i].y;
          fat_width2 = the_graph->data[tot_num_p + (j+1)*the_graph->num_curves + i].y;
          fat_width1 *= prim->band_graph->fatband_scale/2.0;
          fat_width2 *= prim->band_graph->fatband_scale/2.0;
          fb_p1.x = the_graph->data[j*the_graph->num_curves + i].x - xref;
          fb_p1.y = yref - the_graph->data[j*the_graph->num_curves + i].y;
          fb_p2.x = the_graph->data[(j+1)*the_graph->num_curves + i].x - xref;
          fb_p2.y = yref - the_graph->data[(j+1)*the_graph->num_curves + i].y;
          points[0].x = (int)(origin.x + fb_p1.x*obj->scale.x);
          points[0].y = (int)(origin.y + (fb_p1.y+fat_width1)*obj->scale.y);
          points[1].x = (int)(origin.x + fb_p1.x*obj->scale.x);
          points[1].y = (int)(origin.y + (fb_p1.y - fat_width1)*obj->scale.y);
          points[2].x = (int)(origin.x + fb_p2.x*obj->scale.x);
          points[2].y = (int)(origin.y + (fb_p2.y - fat_width2)*obj->scale.y);
          points[3].x = (int)(origin.x + fb_p2.x*obj->scale.x);
          points[3].y = (int)(origin.y + (fb_p2.y + fat_width2)*obj->scale.y);
          points[4].x = points[0].x;
          points[4].y = points[0].y;

          switch(prim->band_graph->fatband_fill){
          case FATBANDS_SHADE:
            g_change_color(2);
            g_filled_polygon(points,4);
            break;
          case FATBANDS_LINE:
            g_change_color(0);
            g_open_polygon(points,4);
            break;
          }
        }
        g_change_linewidth(1.0);
      }
      g_change_color(0);
      for(j=0;j<the_graph->num_p;j++){
        points[j].x = fpoints[j].x = origin.x +
          (the_graph->data[j*the_graph->num_curves + i].x - xref)*obj->scale.x;
        points[j].y = fpoints[j].y = origin.y +
          (yref - the_graph->data[j*the_graph->num_curves + i].y)*obj->scale.y;
      }
      g_lines(points,fpoints,the_graph->num_p,0);
      g_change_linewidth(1.0);
#else
      for(j=0;j<the_graph->num_p;j++){
        points[j].x = fpoints[j].x = origin.x +
          (the_graph->data[j*the_graph->num_curves + i].x - xref)*obj->scale.x;
        points[j].y = fpoints[j].y = origin.y +
          (yref - the_graph->data[j*the_graph->num_curves + i].y)*obj->scale.y;
      }
      g_lines(points,fpoints,the_graph->num_p,0);
#endif
      /* set the line style back to the default value */
      g_change_linestyle(0);
    }
  }
}

