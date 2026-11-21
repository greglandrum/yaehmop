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
 *                   Procedure add_spin_down_props
 *
 *   Spin-polarized run: reads the curves of the spin-down channel from
 *   'partner (same energy grid) and appends them to the graph.  DOS curves
 *   of spin down are drawn on the negative side (mirrored), COOP curves of
 *   spin down are drawn dashed.  Curve i+n is the spin-down partner of
 *   curve i; both are shown or hidden together at the start.
 *
 ****************************************************************************/
static void add_spin_down_props(prop_graph_type *prop_graph,char *partner)
{
  graph_type *u,*d;
  FILE *infile;
  char instring[MAX_STR_LEN];
  point_type2D *raw;
  char *disp,*styles,*fills,*names;
  int nu,nd,nt,np,i,j,ok;
  float sgn;

  infile = fopen(partner,"r");
  if( !infile ) return;
  d = (graph_type *)D_CALLOC(1,sizeof(graph_type));
  if( !d ) fatal("Can't get space for the spin-down curves.");
  skipcomments(infile,instring);
  read_graph_data(infile,d,1,1);
  fclose(infile);
  u = prop_graph->the_data;
  ok = d->num_curves > 0 && d->num_p == u->num_p;
  for(i=0;ok && i<u->num_p;i++){
    if( fabs(d->raw_data[i*d->num_curves].y - u->raw_data[i*u->num_curves].y) > 1e-4 ) ok = 0;
  }
  if( !ok ){
    printf("%s doesn't have the energy grid of %s: showing one spin channel only.\n",
           partner,prop_graph->filename);
  } else {
    nu = u->num_curves; nd = d->num_curves; nt = nu+nd; np = u->num_p;
    sgn = prop_graph->type == DOS_PROP ? -1.0 : 1.0;
    raw = (point_type2D *)D_CALLOC(np*nt,sizeof(point_type2D));
    disp = (char *)D_CALLOC(nt,sizeof(char));
    styles = (char *)D_CALLOC(nt,sizeof(char));
    fills = (char *)D_CALLOC(nt,sizeof(char));
    names = (char *)D_CALLOC((nt+1)*NORMAL_STR_LEN,sizeof(char));
    if( !raw || !disp || !styles || !fills || !names )
      fatal("Can't get space for the spin-polarized curves.");
    for(i=0;i<np;i++){
      for(j=0;j<nu;j++) raw[i*nt+j] = u->raw_data[i*nu+j];
      for(j=0;j<nd;j++){
        raw[i*nt+nu+j].x = sgn*d->raw_data[i*nd+j].x;
        raw[i*nt+nu+j].y = d->raw_data[i*nd+j].y;
      }
    }
    for(j=0;j<nu;j++){
      disp[j] = u->curves_to_display[j];
      styles[j] = u->styles[j];
      fills[j] = u->fills ? u->fills[j] : 0;
    }
    for(j=0;j<nd;j++){
      disp[nu+j] = j < nu ? u->curves_to_display[j] : 0;
      styles[nu+j] = prop_graph->type == DOS_PROP ? (j < nu ? u->styles[j] : j) : 2;
      fills[nu+j] = j < nu ? fills[j] : 0;
    }
    if( u->curve_names ){
      /* the title is kept at 0, the name of curve i at i+1 (as read_graph_data) */
      strncpy(names,u->curve_names,NORMAL_STR_LEN-1);
      for(j=0;j<nt;j++){
        int k = j < nu ? j : j-nu;
        char *src = (k+1 < nu) ? &u->curve_names[(k+1)*NORMAL_STR_LEN] : "";
        snprintf(&names[(j+1)*NORMAL_STR_LEN],NORMAL_STR_LEN,"%.*s%s%s",NORMAL_STR_LEN-8,
                 src,src[0] ? "_" : "",j < nu ? "up" : "down");
      }
      D_FREE(u->curve_names);
    }
    u->curve_names = names;
    D_FREE(u->raw_data); u->raw_data = raw;
    D_FREE(u->data);
    u->data = (point_type2D *)D_CALLOC(np*nt,sizeof(point_type2D));
    if( !u->data ) fatal("Can't get space for the spin-polarized curves.");
    D_FREE(u->curves_to_display); u->curves_to_display = disp;
    D_FREE(u->styles); u->styles = styles;
    if( u->fills ) D_FREE(u->fills);
    u->fills = fills;
    u->num_curves = nt;
    for(i=0;i<np*nt;i++){
      if( raw[i].x > u->max_x ) u->max_x = raw[i].x;
      if( raw[i].x < u->min_x ) u->min_x = raw[i].x;
    }
    prop_graph->spin_down_start = nu;
    printf("Spin-polarized curves: %s (curves 1-%d) and %s (curves %d-%d%s).\n",
           prop_graph->filename,nu,partner,nu+1,nt,
           prop_graph->type == DOS_PROP ? ", drawn as negative DOS" : ", dashed");
  }
  if( d->raw_data ) D_FREE(d->raw_data);
  if( d->data ) D_FREE(d->data);
  if( d->curves_to_display ) D_FREE(d->curves_to_display);
  if( d->styles ) D_FREE(d->styles);
  if( d->fills ) D_FREE(d->fills);
  if( d->curve_names ) D_FREE(d->curve_names);
  D_FREE(d);
}

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
    /* spin-polarized run: add the spin-down channel */
    {
      char partner[MAX_STR_LEN];
      if( whichobj->prim->prop_graph->the_data->num_curves &&
          spin_partner_file(file_name,partner) )
        add_spin_down_props(whichobj->prim->prop_graph,partner);
    }
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
    if( !whichobj->prim->prop_graph->spin_down_start &&
        skipcomments(infile,instring) >= 0 ){
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
    /* spin-polarized run: Fermi level of the .spin summary */
    if( whichobj->prim->prop_graph->spin_down_start &&
        spin_summary_fermi(file_name,&(whichobj->prim->prop_graph->Fermi_E)) )
      whichobj->prim->prop_graph->show_fermi = 1;

    /* make the button window */
#ifdef INTERACTIVE_USE
    build_prop_graph_button_window(&button_wins,whichobj->prim->prop_graph);
    whichobj->prim->but_win = button_wins;
#endif
  }
}
