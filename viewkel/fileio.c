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
   06.04.98 gL:
    modified code to create rayshade file so that it
     (a) just returns if it's not handed a 3D object.
     (b) uses readline to get the filename.
   09.07.98 gL:
    added code to generate VRML files (based on rayshade code)
   05.08.98 gL:
    modifications to VRML generation code.
   26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)
   26.01.99 gL:
     fixed up a bunch of stuff to make file saving/reading work better.
   28.01.99 gL:
     added general_read, which allows only a single file read button
     to appear in the main options window.  yay.  (secretly, this is
     all part of the plan to make the Mac port work better and to
     facilitate the appearance of a Windoze port)
   30.01.99 gL:
     added printing of "#VFIG_FILE" in save_all to allow general_read
     to be able to also deal with save files.
   12.03.99 gL:
     For some reason the tic mark settings weren't being saved
     with property graphs.  that's fixed now.
   19.06.1999 gL:
     trying to save a CONT_PLOT in save_all generated a fatal_bug.
     That was wrong.  Now it just says that it doesn't know how to
     save contour plots.
   30.06.1999 gL:
     added support for binary saves and reads of molecules.
     added dumping of printing options controlling molecule drawing
       to VFIG files.
   15.07.1999 gL:
     saving fatbands plots now saves the fatbands information as well.
     read_bands was changed to reflect this
***/


#include "viewkel.h"
#include "VRML_stuff.h"
#include <fcntl.h>

#ifndef HPUX
extern double atan2(double,double);
#else
extern double atan2();
#endif

#define RECOG_STRING "Viewkel Molecule\n"

/****************************************************************************
 *
 *                   Procedure general_read
 *
 * Arguments: none
 * Returns: none
 *
 * Action: gets a filename from the user and tries to figure out
 *   what kind of file it is.  If successful, the appropriate
 *   read routine is called.  If not, we either give up or allow
 *   the user to pick from a list.
 *
 ****************************************************************************/
void general_read(char *filename)
{
  FILE *infile,*tempfile1,*tempfile2;
  char *theinline,file_name[240],instring[MAX_STR_LEN];
  char tempname1[240],tempname2[240];
  char failed = 0;
  int where,choice,num_choices;

  if(!filename){
#ifndef USING_THE_MAC
    display("Look in the xterm...");
#ifndef USE_READLINE
    printf("Enter the file name to be read: ");
    scanf("%s\n",file_name);
#else /* USE_READLINE */
    theinline= readline("Enter the file name to be read: ");
    if( theinline && *theinline ){
      add_history(theinline);
      sscanf(theinline,"%s",file_name);
      free(theinline);
    } else {
      error("Bad file name");
      file_name[0] = 0;
      return;
    }
  } else {
    strncpy(file_name,filename,240);
  }
#endif /* USE_READLINE */
  if( file_name[0] ) infile = fopen(file_name,"r");
#else /* USING_THE_MAC */
  infile = choose_mac_file(file_name,MAC_FOPEN_OPEN_CD);
#endif

  if( !infile ){
    printf("Problems opening file: %s\n",file_name);
    display("oooooops!");
    return;
  }


  /*****

        okay, grab the first non-blank, non-comment line
        and work with that to see what kind of file we've
    got.  This is relying on matching to strings which
    may not be in files produced by older versions of
    bind or the utilities, so there is an extra bit at
    the bottom to let the user decide what the heck
    kind of file something is.

  ******/
  skipcomments(infile,instring);
  upcase(instring);
  if( strstr(instring,"#BIND_OUTPUT") ||
             strstr(instring,"MOVIE-FILE")) {
    printf("%s looks like a molecule file.\n",file_name);
    fclose(infile);
    new_molecule(file_name);
#ifdef INCLUDE_ADF_PLOTS
  } else if( strstr(instring,"#ADF_MO_DATA") ){
    printf("%s looks like an ADF MO output file.\n",file_name);
    fclose(infile);
    new_ADF_MO_surface(file_name);
  } else if(strstr(instring,"#ADF_VIBRATION_DATA")){
    printf("%s looks like an ADF vibration output file.\n",file_name);
    fclose(infile);
    new_vibration(file_name);
#endif
  } else if(strstr(instring,"#GRAPH_DATA")){
    printf("%s looks like a graph file.\n",file_name);
    fclose(infile);
    new_graph(file_name);
  } else if(strstr(instring,"#FCO_DATA")){
    printf("%s looks like an FCO file.\n",file_name);
    fclose(infile);
    new_cont_plot(file_name);
  } else if(strstr(instring,"#MO_DATA")){
    printf("%s looks like an MO grid file.\n",file_name);
    fclose(infile);
    new_cont_plot(file_name);
  } else if(strstr(instring,"#DOS DATA") || strstr(instring,"# DOS DATA") ||
            strstr(instring,"#COOP DATA")){
    printf("%s looks like a property file.\n",file_name);
    fclose(infile);
    new_prop_graph(file_name);
  } else if(strstr(instring,"#BAND_DATA") ||strstr(instring,"#BAND DATA") ||
            strstr(instring,"#FATBANDS")){
    printf("%s looks like a band file.\n",file_name);
    fclose(infile);
    new_band_graph(file_name);
  } else if(strstr(instring,"#WALSH_DATA")){
    printf("%s looks like a walsh file.\n",file_name);
    fclose(infile);
    new_walsh_graph(file_name);
  } else if(strstr(instring,"#FMO RESULTS") ||
            strstr(instring,"# FMO RESULTS") ){
    printf("%s looks like an FMO file.\n",file_name);
    fclose(infile);
    new_FMO_diagram(file_name);
  } else if(strstr(instring,"#VFIG_FILE")){
    printf("%s looks like a vfig file.\n",file_name);
    fclose(infile);
    read_all(file_name);
  } else if(strstr(instring,"#BEGIN_PARMS")){
    printf("%s looks like an MO file.\n",file_name);
    /*******
      now things get a bit wacky...
      we've got to cut the .MO off the name of the file
      in order for the MO reading stuff (which wants
      the name of an input file) to work.

      NOTE: the way this is wired, the file name better
      end with ".MO" or it's not going to work at all.
      This is a silly assumption, but it is "justified" by
      the way read_MO_data works.
    ********/
    fclose(infile);
    where = strlen(file_name)-3;
    file_name[where] = 0;
    strcpy(tempname1,file_name);
    strcat(tempname1,".out");
    infile = fopen(tempname1,"r");
    if (!infile ) {
      printf("\tbut I can't find: %s.\n",tempname1);
      display("Dogged!");
      return;
    }
    fclose(infile);
    new_MO_surface(file_name);
  } else{
    printf("I don't know how to deal with file: %s\n",file_name);
    printf("What should I do with it?\n");
    choice = 0;
    while(!choice){
      printf("\t(-1) Forget it (Cancel).\n");
      printf("\t(1) Read Molecule\n");
      printf("\t(2) Read MO\n");
      printf("\t(3) Read Graph\n");
      printf("\t(4) Read Contours\n");
      printf("\t(5) Read Properties\n");
      printf("\t(6) Read Bands\n");
      printf("\t(7) Read Walsh\n");
      printf("\t(8) Read FMO\n");
      printf("\t(9) Read All (save file)\n");
      num_choices = 9;
#ifdef INCLUDE_ADF_PLOTS
      printf("\t(10) Read ADF MO\n");
      printf("\t(11) Read ADF Vibration\n");
      num_choices = 11;
#endif
      scanf("%d",&choice);
      if( choice > num_choices ){
        printf("Invalid choice, try again.\n");
        choice = 0;
      }
    }
    switch(choice){
    case -1:
      printf("Giving up.\n");
      display("Cancelled!");
      break;
    case 1:
      new_molecule(file_name);break;
    case 2:
      new_MO_surface(file_name);break;
    case 3:
      new_graph(file_name);break;
    case 4:
      new_cont_plot(file_name);break;
    case 5:
      new_prop_graph(file_name);break;
    case 6:
      new_band_graph(file_name);break;
    case 7:
      new_walsh_graph(file_name);break;
    case 8:
      new_FMO_diagram(file_name);break;
    case 9:
      read_all(file_name);break;
#ifdef INCLUDE_ADF_PLOTS
    case 10:
      new_ADF_MO_surface(file_name);break;
    case 11:
      new_vibration(file_name);break;
#endif
    default:
      FATAL_BUG("bad choice parameter\n");
    }
  }
}

/****************************************************************************
 *
 *                   Procedure save_bands
 *
 * Arguments: vfig_file: pointer to type FILE
 *            prim: a pointer to a primitive type
 *             obj: a pointer to object_type
 * Returns: none
 *
 * Action: writes the band graph in 'prim to 'vfig_file.
 *
 ****************************************************************************/
void save_bands( FILE *vfig_file,prim_type *prim,object_type *obj )
{
  int i;
  band_graph_type *the_graph;
  graph_type *the_data;

  the_graph = prim->band_graph;
  the_data = the_graph->the_data;
  fprintf(vfig_file,"\n\nTYPE:  BAND_GRAPH\n");
  fprintf(vfig_file,"FILENAME: %s\n",the_graph->filename);
#ifdef SUPPORT_FATBANDS
  if( the_graph->num_fatbands > 0 ){
    fprintf(vfig_file,"FATBAND_FILL: %d",the_graph->fatband_fill);
    fprintf(vfig_file,"FATBANDS_ON: %d",the_graph->fatbands_on);
    fprintf(vfig_file,"FATBAND_SCALE: %lf",the_graph->fatband_scale);
  }
#endif
  fprintf(vfig_file,"XLEGEND: ");
  fputs(the_data->xlegend,vfig_file);
  fprintf(vfig_file,"\n");
  fprintf(vfig_file,"YLEGEND: ");
  fputs(the_data->ylegend,vfig_file);
  fprintf(vfig_file,"\n");
  fprintf(vfig_file,"TITLE:\n");
  for(i=0;i<NUM_TITLE_LINES;i++){
    fputs(the_data->title[i],vfig_file);
    fprintf(vfig_file,"\n");
  }
  fprintf(vfig_file,"DO_X_TICS: %d\nDO_Y_TICS: %d\nDO_TITLE: %d\n",
          the_data->do_x_tics,the_data->do_y_tics,the_data->do_title);
  fprintf(vfig_file,"STYLES:\n");
  fprintf(vfig_file,"\t%d\n",the_data->styles[0]);

  fprintf(vfig_file,"MIN_X: %lf\nMAX_X: %lf\nMAX_Y: %lf\nMIN_Y: %lf\n",
          the_data->min_x,the_data->max_x,the_data->max_y,the_data->min_y);
  fprintf(vfig_file,"TIC_SEP_X: %lf\nNUM_TICS_X: %d\nTIC_START_X: %lf\n",
          the_data->tic_sep_x,(int)the_data->num_tics_x,the_data->tic_start_x);
  fprintf(vfig_file,"TIC_SEP_Y: %lf\nNUM_TICS_Y: %d\nTIC_START_Y: %lf\n",
          the_data->tic_sep_y,(int)the_data->num_tics_y,the_data->tic_start_y);
  fprintf(vfig_file,"SHOW_FERMI: %d\nFERMI_E: %lf\n",the_graph->show_fermi,
          the_graph->Fermi_E);

}

/****************************************************************************
 *
 *                   Procedure save_bin_molecule
 *
 * Arguments: molec: pointer to molec_type
 *
 * Returns: none
 *
 * Action: writes a molecule in binary format
 *
 ****************************************************************************/
int save_bin_molecule(molec_type *molec)
{
  int outfile;
  char filename[240];
  int size,num_written;
  int i;

  sprintf(filename,"%s.BIN",molec->filename);
  outfile = creat(filename,0644);
  if(!outfile){
    error("can't open binary molecule save file.");
    return(-1);
  }
  /* deselect all the atoms... otherwise we'll have big problems */
  for(i=0;i<molec->num_frames*molec->num_atoms;i++){
    molec->atoms[i].is_selected = 0;
  }

#define TRY_WRITE(_a_,_b_,_c_) {num_written=write(_a_,_b_,_c_);\
                if(num_written!=_c_){error("wrong size written");\
                                     unlink(filename);return(-2);}}
  TRY_WRITE(outfile,RECOG_STRING,strlen(RECOG_STRING));
  TRY_WRITE(outfile,molec,sizeof(molec_type));
  size=molec->num_atoms*molec->num_frames*sizeof(atom_type);
  TRY_WRITE(outfile,molec->atoms,size);
  for(i=0;i<molec->num_atoms;i++){
    size = molec->atoms[i].num_lines_out*sizeof(int);
    TRY_WRITE(outfile,molec->atoms[i].linesto,size);
#ifdef INCLUDE_ADF_PLOTS
    if( molec->num_vibrations ){
      size = molec->num_vibrations*sizeof(point_type);
      TRY_WRITE(outfile,molec->atoms[i].displacements,size);
    }
#endif
  }
  size = molec->num_frames * sizeof(int);
  TRY_WRITE(outfile,molec->num_lines,size);
  size = 0;
  for(i=0;i<molec->num_frames;i++) size += molec->num_lines[i];
  size *= sizeof(line_type);
  TRY_WRITE(outfile,molec->lines,size);
  close(outfile);
  return(0);
}

/****************************************************************************
 *
 *                   Procedure read_bin_molecule
 *
 * Arguments: molec: pointer to molec_type
 *
 * Returns: none
 *
 * Action: reads in a binary format save of a molecule
 *
 ****************************************************************************/
int read_bin_molecule(molec_type *molec)
{
  int infile;
  char filename[240],instring[240];
  int size,num_written;
  int i;

  sprintf(filename,"%s.BIN",molec->filename);
  infile = open(filename,O_RDONLY,0644);
  if(!infile){
    error("can't open binary molecule save file.");
    return(-1);
  }

#define TRY_READ(_a_,_b_,_c_) {num_written=read(_a_,_b_,_c_);\
                if(num_written!=_c_){error("wrong size written");\
                                     unlink(filename);return(-2);}}
  TRY_READ(infile,instring,strlen(RECOG_STRING));
  if(strcmp(instring,RECOG_STRING)){
    error("Recognition string does not match.  Mabye the file is corrupted.");
    return(-3);
  }
  TRY_READ(infile,molec,sizeof(molec_type));
  size=molec->num_atoms*molec->num_frames*sizeof(atom_type);
  molec->atoms = (atom_type *)D_MALLOC(size);
  if(!molec->atoms) fatal("Can't allocate space for atoms");
  TRY_READ(infile,molec->atoms,size);
  for(i=0;i<molec->num_atoms;i++){
    size = molec->atoms[i].num_lines_out*sizeof(int);
    if( size > 0 ){
      molec->atoms[i].linesto = (int *)D_MALLOC(size);
      if(!molec->atoms[i].linesto) fatal("Can't allocate space for linesto");
      TRY_READ(infile,molec->atoms[i].linesto,size);
    } else {
      molec->atoms[i].linesto = 0;
    }
#ifdef INCLUDE_ADF_PLOTS
    if( molec->num_vibrations ){
      size = molec->num_vibrations*sizeof(point_type);
      molec->atoms[i].displacements = (point_type *)D_MALLOC(size);
      if(!molec->atoms[i].displacements) fatal("Can't allocate space for displacements");
      TRY_READ(infile,molec->atoms[i].displacements,size);
    }
#endif
  }
  size = molec->num_frames * sizeof(int);
  molec->num_lines = (int *)malloc(size);
  if(!molec->num_lines)fatal("Can't allocate molec->num_lines");
  TRY_READ(infile,molec->num_lines,size);
  size = 0;
  for(i=0;i<molec->num_frames;i++) size += molec->num_lines[i];
  size *= sizeof(line_type);
  if( size > 0 ){
    molec->lines = (line_type *)malloc(size);
    if(!molec->lines)fatal("Can't allocate molec->lines");
    TRY_READ(infile,molec->lines,size);
  } else{
    molec->lines = 0;
  }

  close(infile);
  return(0);
}


/****************************************************************************
 *
 *                   Procedure save_molecule
 *
 * Arguments: vfig_file: pointer to type FILE
 *            prim: a pointer to a primitive type
 *             obj: a pointer to object_type
 * Returns: none
 *
 * Action: writes the molecule in 'prim to 'vfig_file.
 *
 ****************************************************************************/
void save_molecule( FILE *vfig_file,prim_type *prim,object_type *obj )
{
  molec_type *the_molec;
  int bin_failed;

  the_molec = prim->molec;
  fprintf(vfig_file,"\n\nTYPE:  MOLECULE\n");
  fprintf(vfig_file,"FILENAME: %s\n",the_molec->filename);
  bin_failed = save_bin_molecule(the_molec);
  if( !bin_failed ){
    fprintf(vfig_file,"BINFILENAME: %s.BIN\n",the_molec->filename);
  } else{
    fprintf(vfig_file,"NUMBERS_ON: %d\n",the_molec->numbers_on);
    fprintf(vfig_file,"DRAW_CONNECTORS: %d\n",the_molec->draw_connectors);
    fprintf(vfig_file,"OUTLINES_ON: %d\n",the_molec->outlines_on);
    fprintf(vfig_file,"SHADING_ON: %d\n",the_molec->shading_on);
    fprintf(vfig_file,"HYDROGENS_ON: %d\n",the_molec->hydrogens_on);
    fprintf(vfig_file,"FANCY_LINES: %d\n",the_molec->fancy_lines);
    fprintf(vfig_file,"SYMBOLS_ON: %d\n",the_molec->symbols_on);
    fprintf(vfig_file,"AXES_ON: %d\n",the_molec->axes_on);
    fprintf(vfig_file,"DUMMIES_ON: %d\n",the_molec->dummies_on);
    fprintf(vfig_file,"CROSSES_ON: %d\n",the_molec->crosses_on);
    fprintf(vfig_file,"TUBES_ON: %d\n",the_molec->tubes_on);
    fprintf(vfig_file,"BREAKING_LINES: %d\n",the_molec->breaking_lines);
    fprintf(vfig_file,"BOND_TOL: %lf\n",the_molec->bond_tol);
    fprintf(vfig_file,"RAD_MULT: %lf\n",the_molec->rad_mult);
  }
}

/****************************************************************************
 *
 *                   Procedure save_MO_surf
 *
 * Arguments: vfig_file: pointer to type FILE
 *            prim: a pointer to a primitive type
 *             obj: a pointer to object_type
 * Returns: none
 *
 * Action: writes the MO surf in 'prim to 'vfig_file.
 *
 ****************************************************************************/
void save_MO_surf( FILE *vfig_file,prim_type *prim,object_type *obj )
{
  MO_surface_type *the_surf;
  molec_type *the_molec;

  the_surf = prim->MO_surf;
  the_molec = the_surf->molec;
  fprintf(vfig_file,"\n\nTYPE:  MO_SURF\n");
  fprintf(vfig_file,"FILENAME: %s\n",the_surf->filename);
  fprintf(vfig_file,"NUMBERS_ON: %d\n",the_molec->numbers_on);
  fprintf(vfig_file,"DRAW_CONNECTORS: %d\n",the_molec->draw_connectors);
  fprintf(vfig_file,"OUTLINES_ON: %d\n",the_molec->outlines_on);
  fprintf(vfig_file,"SHADING_ON: %d\n",the_molec->shading_on);
  fprintf(vfig_file,"HYDROGENS_ON: %d\n",the_molec->hydrogens_on);
  fprintf(vfig_file,"FANCY_LINES: %d\n",the_molec->fancy_lines);
  fprintf(vfig_file,"SYMBOLS_ON: %d\n",the_molec->symbols_on);
  fprintf(vfig_file,"AXES_ON: %d\n",the_molec->axes_on);
  fprintf(vfig_file,"DUMMIES_ON: %d\n",the_molec->dummies_on);
  fprintf(vfig_file,"CROSSES_ON: %d\n",the_molec->crosses_on);
  fprintf(vfig_file,"TUBES_ON: %d\n",the_molec->tubes_on);
  fprintf(vfig_file,"BREAKING_LINES: %d\n",the_molec->breaking_lines);
  fprintf(vfig_file,"BOND_TOL: %lf\n",the_molec->bond_tol);
  fprintf(vfig_file,"RAD_MULT: %lf\n",the_molec->rad_mult);
  fprintf(vfig_file,"ACTIVE_MO: %d\n",the_surf->active_MO);
  fprintf(vfig_file,"DISPLAY_SURF: %d\n",the_surf->display_surf);
  fprintf(vfig_file,"DISPLAY_MOLEC: %d\n",the_surf->display_molec);
  fprintf(vfig_file,"DO_SHADING: %d\n",the_surf->do_shading);
  fprintf(vfig_file,"DO_LINES: %d\n",the_surf->do_lines);
  fprintf(vfig_file,"NUM_LOOKUP_ENTRIES: %d\n",the_surf->num_lookup_entries);
  fprintf(vfig_file,"LOOKUP_MIN: %lf\n",the_surf->lookup_min);
  fprintf(vfig_file,"LOOKUP_MAX: %lf\n",the_surf->lookup_max);
  fprintf(vfig_file,"CONTRACTION: %lf\n",the_surf->contraction);
  fprintf(vfig_file,"SURFACE_VALUE: %lf\n",the_surf->surface_value);
  fprintf(vfig_file,"SURFACE_TOLERANCE: %lf\n",the_surf->surface_tolerance);
  fprintf(vfig_file,"SLOP: %lf\n",the_surf->slop);
  fprintf(vfig_file,"SEARCH_RADIUS: %lf\n",the_surf->search_radius);
  fprintf(vfig_file,"VOXEL_SIZE: %lf\n",the_surf->voxel_size);
}

/****************************************************************************
 *
 *                   Procedure save_graph
 *
 * Arguments: vfig_file: pointer to type FILE
 *            prim: a pointer to a primitive type
 *             obj: a pointer to object_type
 * Returns: none
 *
 * Action: writes the graph in 'prim to 'vfig_file.
 *
 ****************************************************************************/
void save_graph( FILE *vfig_file,prim_type *prim,object_type *obj )
{
  int i;
  graph_type *the_graph;

  the_graph = prim->graph;
  fprintf(vfig_file,"\n\nTYPE:  GRAPH\n");
  fprintf(vfig_file,"FILENAME: %s\n",the_graph->filename);
  fprintf(vfig_file,"XLEGEND: ");
  fputs(the_graph->xlegend,vfig_file);
  fprintf(vfig_file,"\n");
  fprintf(vfig_file,"YLEGEND: ");
  fputs(the_graph->ylegend,vfig_file);
  fprintf(vfig_file,"\n");
  fprintf(vfig_file,"TITLE:\n");
  for(i=0;i<NUM_TITLE_LINES;i++){
    fputs(the_graph->title[i],vfig_file);
    fprintf(vfig_file,"\n");
  }
  fprintf(vfig_file,"DO_X_TICS: %d\nDO_Y_TICS: %d\nDO_TITLE: %d\n",
          the_graph->do_x_tics,the_graph->do_y_tics,the_graph->do_title);
  fprintf(vfig_file,"CURVES:\n");
  for(i=0;i<the_graph->num_curves;i++){
    fprintf(vfig_file,"\t%d\n",the_graph->curves_to_display[i]);
  }
  fprintf(vfig_file,"STYLES:\n");
  for(i=0;i<the_graph->num_curves;i++){
    fprintf(vfig_file,"\t%d\n",the_graph->styles[i]);
  }
  fprintf(vfig_file,"MIN_X: %lf\nMAX_X: %lf\nMAX_Y: %lf\nMIN_Y: %lf\n",
          the_graph->min_x,the_graph->max_x,the_graph->max_y,the_graph->min_y);
  fprintf(vfig_file,"TIC_SEP_X: %lf\nNUM_TICS_X: %d\nTIC_START_X: %lf\n",
          the_graph->tic_sep_x,(int)the_graph->num_tics_x,the_graph->tic_start_x);
  fprintf(vfig_file,"TIC_SEP_Y: %lf\nNUM_TICS_Y: %d\nTIC_START_Y: %lf\n",
          the_graph->tic_sep_y,(int)the_graph->num_tics_y,the_graph->tic_start_y);
}


/****************************************************************************
 *
 *                   Procedure save_FMO
 *
 * Arguments: vfig_file: pointer to type FILE
 *            prim: a pointer to a primitive type
 *             obj: a pointer to object_type
 * Returns: none
 *
 * Action: saves the FMO diagram in 'prim to 'vfig_file.
 *
 ****************************************************************************/
void save_FMO( FILE *vfig_file,prim_type *prim,object_type *obj )
{
  int i;
  FMO_diagram_type *the_graph;

  the_graph = prim->FMO_diagram;
  fprintf(vfig_file,"\n\nTYPE:  FMO\n");
  fprintf(vfig_file,"FILENAME: %s\n",the_graph->filename);
  fprintf(vfig_file,"DO_TITLE: %d\n",the_graph->do_title);
  fprintf(vfig_file,"SHOW_BOX: %d\n",the_graph->show_box);
  fprintf(vfig_file,"SHOW_DATA: %d\n",the_graph->show_data);
  fprintf(vfig_file,"SHOW_CONNECTS: %d\n",the_graph->show_connects);
  fprintf(vfig_file,"FILLING_MODE: %d\n",the_graph->electron_filling_mode);
  fprintf(vfig_file,"DO_Y_TICS: %d\n",the_graph->do_y_tics);
  fprintf(vfig_file,"LEFT_FRAGMENT: %d\n",the_graph->left_fragment);
  fprintf(vfig_file,"RIGHT_FRAGMENT: %d\n",the_graph->right_fragment);
  fprintf(vfig_file,"LEVEL_WIDTH: %d\n",the_graph->level_width);
  fprintf(vfig_file,"LEVEL_THICKNESS: %d\n",the_graph->level_thickness);
  fprintf(vfig_file,"ELECTRON_LENGTH: %d\n",the_graph->electron_length);
  fprintf(vfig_file,"MIN_Y: %lf\n",the_graph->min_y);
  fprintf(vfig_file,"MAX_Y: %lf\n",the_graph->max_y);
  fprintf(vfig_file,"TIC_SEP_Y: %lf\n",the_graph->tic_sep_y);
  fprintf(vfig_file,"TIC_START_Y: %lf\n",the_graph->tic_start_y);
  fprintf(vfig_file,"NUM_TICS_Y: %d\n",(int)the_graph->num_tics_y);
  fprintf(vfig_file,"TITLE\n");
  for(i=0;i<NUM_TITLE_LINES;i++){
    fprintf(vfig_file,"%s\n",the_graph->title[i]);
  }
  fprintf(vfig_file,"LABELS\n");
  for(i=0;i<the_graph->num_frags;i++){
    fprintf(vfig_file,"%s\n",the_graph->frags[i].label);
  }
  fprintf(vfig_file,"\nLABEL\n");
  fprintf(vfig_file,"%s\n",the_graph->label);

}


/****************************************************************************
 *
 *                   Procedure save_prop_graph
 *
 * Arguments: vfig_file: pointer to type FILE
 *            prim: a pointer to a primitive type
 *             obj: a pointer to object_type
 * Returns: none
 *
 * Action: writes the prop_graph in 'prim to 'vfig_file.
 *
 ****************************************************************************/
void save_prop_graph( FILE *vfig_file,prim_type *prim,object_type *obj )
{
  int i;
  prop_graph_type *prop_graph;
  graph_type *the_graph,*the_integration;

  prop_graph = prim->prop_graph;
  the_graph = prim->prop_graph->the_data;
  the_integration = prim->prop_graph->the_integration;
  fprintf(vfig_file,"TYPE:  PROPERTY_GRAPH\n");
  fprintf(vfig_file,"FILENAME: %s\n",prop_graph->filename);
  fprintf(vfig_file,"XLEGEND: ");
  fputs(the_graph->xlegend,vfig_file);
  fprintf(vfig_file,"\n");
  fprintf(vfig_file,"YLEGEND: ");
  fputs(the_graph->ylegend,vfig_file);
  fprintf(vfig_file,"\n");
  fprintf(vfig_file,"TITLE:\n");
  for(i=0;i<NUM_TITLE_LINES;i++){
    fputs(the_graph->title[i],vfig_file);
    fprintf(vfig_file,"\n");
  }
  fprintf(vfig_file,"TITLE:\n");
  for(i=0;i<NUM_TITLE_LINES;i++) fprintf(vfig_file,"%s\n",the_graph->title[i]);
  fprintf(vfig_file,"CURVES:\n");
  for(i=0;i<the_graph->num_curves;i++){
    fprintf(vfig_file,"\t%d\n",the_graph->curves_to_display[i]);
  }
  fprintf(vfig_file,"STYLES:\n");
  for(i=0;i<the_graph->num_curves;i++){
    fprintf(vfig_file,"\t%d\n",the_graph->styles[i]);
  }
  fprintf(vfig_file,"FILLS:\n");
  for(i=0;i<the_graph->num_curves;i++){
    fprintf(vfig_file,"\t%d\n",the_graph->fills[i]);
  }
  fprintf(vfig_file,"INTEG_CURVES:\n");
  if( the_integration){
    for(i=0;i<the_integration->num_curves;i++){
      fprintf(vfig_file,"\t%d\n",the_integration->curves_to_display[i]);
    }
    fprintf(vfig_file,"INTEG_STYLES:\n");
    for(i=0;i<the_integration->num_curves;i++){
      fprintf(vfig_file,"\t%d\n",the_integration->styles[i]);
    }
  }
  fprintf(vfig_file,"INTEGS_FOR_TICS: %d\n",prop_graph->integs_for_tics);
  fprintf(vfig_file,"MIN_X: %lf\nMAX_X: %lf\nMAX_Y: %lf\nMIN_Y: %lf\n",
          prop_graph->min_x,prop_graph->max_x,prop_graph->max_y,prop_graph->min_y);
  fprintf(vfig_file,"DO_X_TICS: %d\nDO_Y_TICS: %d\nDO_TITLE: %d\n",
          the_graph->do_x_tics,the_graph->do_y_tics,the_graph->do_title);

  fprintf(vfig_file,"TIC_SEP_X: %lf\nNUM_TICS_X: %d\nTIC_START_X: %lf\n",
          the_graph->tic_sep_x,(int)the_graph->num_tics_x,the_graph->tic_start_x);
  fprintf(vfig_file,"TIC_SEP_Y: %lf\nNUM_TICS_Y: %d\nTIC_START_Y: %lf\n",
          the_graph->tic_sep_y,(int)the_graph->num_tics_y,the_graph->tic_start_y);

  fprintf(vfig_file,"PROP_TYPE: %d\n",prop_graph->type);
  fprintf(vfig_file,"SHOW_FERMI: %d\nFERMI_E: %lf\n",prop_graph->show_fermi,
          prop_graph->Fermi_E);
}



/****************************************************************************
 *
 *                   Procedure save_prim
 *
 * Arguments: vfig_file: pointer to type FILE
 *            prim: a pointer to a primitive type
 *             obj: a pointer to object_type
 * Returns: none
 *
 * Action: writes the 'prim to 'vfig_file.  This currently doesn't do
 *   anything with molecules or MO surfaces
 *
 ****************************************************************************/
void save_prim( FILE *vfig_file,prim_type *prim,object_type *obj )
{
  fprintf(vfig_file,"\nBEGIN_OBJECT\n");
  switch(prim->which){
  case GRAPH:
    save_graph(vfig_file,prim,obj);
    break;
  case PROP_GRAPH:
    save_prop_graph(vfig_file,prim,obj);
    break;
  case MOLECULE:
    save_molecule(vfig_file,prim,obj);
    break;
  case MO_SURF:
    save_MO_surf(vfig_file,prim,obj);
    break;
  case BAND_GRAPH:
    save_bands(vfig_file,prim,obj);
    break;
  case WALSH_GRAPH:
/*
    save_walsh_graph(vfig_file,prim,obj);
*/
    error("Can't save walsh graphs yet.");
    break;
  case FMO_DIAGRAM:
    save_FMO(vfig_file,prim,obj);
    break;
  case CONT_PLOT:
    error("Can't save contour plots yet.");
    break;
  case LABEL:
    break;
  default:
    FATAL_BUG("Bogus primitive in save_prim.");
    break;
  }

  fprintf(vfig_file,"CENT: %lf %lf %lf\n",
          obj->cent.x,obj->cent.y,obj->cent.z);
  fprintf(vfig_file,"TRANS: %lf %lf %lf\n",
          obj->trans.x,obj->trans.y,obj->trans.z);
  fprintf(vfig_file,"SCALE: %lf %lf %lf\n",
          obj->scale.x,obj->scale.y,obj->scale.z);
  fprintf(vfig_file,"ROT: %lf %lf %lf\n",
          obj->rot.x,obj->rot.y,obj->rot.z);

  fprintf(vfig_file,"END_OBJECT\n");
}



/****************************************************************************
 *
 *                   Procedure save_all
 *
 * Arguments: none
 * Returns: none
 *
 * Action: writes everything to an output file
 *
 ****************************************************************************/
void save_all(void)
{
  FILE *vfig_file;
  char filename[240];
  char *theinline;
  head_type *temphead;
  object_type *obj;
  int num_written;

  display("Look in xterm");
#ifndef USE_READLINE
  printf("Enter save file name: ");
  scanf("%s\n",filename);
#else
  theinline= readline("Enter save file name: ");
  add_history(theinline);
  if( theinline ){
    sscanf(theinline,"%s",filename);
    free(theinline);
  } else {
    error("Bad file name");
    filename[0] = 0;
  }
#endif

  vfig_file = fopen(filename,"w+");
  if(!vfig_file){
    error("Can't open save file.");
    return;
  }

  /* write the global variables */
  fprintf(vfig_file,"#VFIG_FILE\n");
  fprintf(vfig_file,"FILL_PROJECTIONS: %d\n",fill_projections);
  fprintf(vfig_file,"PS_WHERE_TO_PRINT: %d\n",PS_options.where_to_print);
  fprintf(vfig_file,"PS_FONTNAME: %s\n",PS_options.fontname);
  fprintf(vfig_file,"PS_FONTSIZE: %lf\n",PS_options.fontsize);
  fprintf(vfig_file,"PS_PRINTSCALE: %lf\n",PS_options.printscale);
  fprintf(vfig_file,"ATOM_TYPE: %d\n",PS_options.atom_sphere_type);
  fprintf(vfig_file,"BOND_TYPE: %d\n",PS_options.bond_type);
  num_written = 0;
  temphead = head;
  while(temphead){
    obj = temphead->obj;
    if(obj->prim){
      save_prim(vfig_file,obj->prim,obj);
      num_written++;
    }
    temphead = temphead->next;
  }
  fclose(vfig_file);
  printf("%d objects saved\n",num_written);
}




/****************************************************************************
 *
 *                   Function read_molec_object
 *
 * Arguments: vfig_file: pointer to type FILE
 *             filename: pointer to type char
 * Returns: int
 *
 * Action: reads a molec object out of 'filename, then sets its options
 *    according to the keywords in 'vfig_file.
 *
 *    returns a negative number on failure.
 *
 ****************************************************************************/
int read_molec_object(FILE *vfig_file,char *filename)
{
  int temp;
  molec_type *the_molec;
  char instring[MAX_STR_LEN],keyword[80];
  int read_error = 0;

  /* start off by reading in the data */
  new_molecule(filename);
  if( !whichobj ) return -1;

  the_molec = whichobj->prim->molec;

  /* now read the information out of 'vfig_file */
  read_error = skipcomments(vfig_file,instring);
  while(read_error >= 0 && !strstr(instring,"END_OBJECT")){
    sscanf(instring,"%s",keyword);
    if( strstr(keyword,"BINFILENAME") ){
      printf("There's a binary file for this molecule.  I'll use that.\n");
    }else if(strstr(keyword,"RAD_MULT")){
      sscanf(instring,"%s %lf",keyword,&the_molec->rad_mult);
    }else if(strstr(keyword,"BOND_TOL")){
      sscanf(instring,"%s %lf",keyword,&the_molec->bond_tol);
    }else if(strstr(keyword,"BREAKING_LINES")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->breaking_lines = temp;
    }else if(strstr(keyword,"DUMMIES_ON")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->dummies_on = temp;
    }else if(strstr(keyword,"CROSSES_ON")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->crosses_on = temp;
    }else if(strstr(keyword,"TUBES_ON")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->tubes_on = temp;
    }else if(strstr(keyword,"AXES_ON")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->axes_on = temp;
    }else if(strstr(keyword,"SYMBOLS_ON")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->symbols_on = temp;
    }else if(strstr(keyword,"FANCY_LINES")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->fancy_lines = temp;
    }else if(strstr(keyword,"HYDROGENS_ON")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->hydrogens_on = temp;
    }else if(strstr(keyword,"SHADING_ON")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->shading_on = temp;
    }else if(strstr(keyword,"OUTLINES_ON")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->outlines_on = temp;
    }else if(strstr(keyword,"DRAW_CONNECTORS")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->draw_connectors = temp;
    }else if(strstr(keyword,"NUMBERS_ON")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->numbers_on = temp;
    } else if(strstr(keyword,"CENT")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->cent.x,
             &whichobj->cent.y,&whichobj->cent.z);
    } else if(strstr(keyword,"TRANS")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->trans.x,
             &whichobj->trans.y,&whichobj->trans.z);
    } else if(strstr(keyword,"SCALE")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->scale.x,
             &whichobj->scale.y,&whichobj->scale.z);
    } else if(strstr(keyword,"ROT")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->rot.x,
             &whichobj->rot.y,&whichobj->rot.z);
    } else{
      fprintf(stderr,"Unrecognized keyword: %s\n",keyword);
    }
    read_error = skipcomments(vfig_file,instring);
  }
  if( read_error >= 0 ) return(1);
}


/****************************************************************************
 *
 *                   Function read_MO_surf_object
 *
 * Arguments: vfig_file: pointer to type FILE
 *             filename: pointer to type char
 * Returns: int
 *
 * Action: reads an MO surface object out of 'filename, then sets its options
 *    according to the keywords in 'vfig_file.
 *
 *    returns a negative number on failure.
 *
 ****************************************************************************/
int read_MO_surf_object(FILE *vfig_file,char *filename)
{
  int temp;
  MO_surface_type *the_surf;
  molec_type *the_molec;
  char instring[MAX_STR_LEN],keyword[80];
  int read_error = 0;

  /* start off by reading in the data */
  new_MO_surface(filename);
  if( !whichobj ) return -1;

  the_surf = whichobj->prim->MO_surf;
  the_molec = the_surf->molec;

  /* now read the information out of 'vfig_file */
  read_error = skipcomments(vfig_file,instring);
  while(read_error >= 0 && !strstr(instring,"END_OBJECT")){
    sscanf(instring,"%s",keyword);
    if(strstr(keyword,"VOXEL_SIZE")){
      sscanf(instring,"%s %lf",keyword,&the_surf->voxel_size);
    }else if(strstr(keyword,"SEARCH_RADIUS")){
      sscanf(instring,"%s %lf",keyword,&the_surf->search_radius);
    } else if(strstr(keyword,"SLOP")){
      sscanf(instring,"%s %lf",keyword,&the_surf->slop);
    } else if(strstr(keyword,"SURFACE_TOLERANCE")){
      sscanf(instring,"%s %lf",keyword,&the_surf->surface_tolerance);
    } else if(strstr(keyword,"SURFACE_VALUE")){
      sscanf(instring,"%s %lf",keyword,&the_surf->surface_value);
    } else if(strstr(keyword,"CONTRACTION")){
      sscanf(instring,"%s %lf",keyword,&the_surf->contraction);
    } else if(strstr(keyword,"LOOKUP_MAX")){
      sscanf(instring,"%s %lf",keyword,&the_surf->lookup_max);
    } else if(strstr(keyword,"LOOKUP_MIN")){
      sscanf(instring,"%s %lf",keyword,&the_surf->lookup_min);
    } else if(strstr(keyword,"NUM_LOOKUP_ENTRIES")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_surf->num_lookup_entries = temp;
    }else if(strstr(keyword,"DO_LINES")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_surf->do_lines = temp;
    } else if(strstr(keyword,"DO_SHADING")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_surf->do_shading = temp;
    } else if(strstr(keyword,"DISPLAY_MOLEC")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_surf->display_molec = temp;
    } else if(strstr(keyword,"DISPLAY_SURF")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_surf->display_surf = temp;
    } else if(strstr(keyword,"ACTIVE_MO")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_surf->active_MO = temp;
    } else if(strstr(keyword,"RAD_MULT")){
      sscanf(instring,"%s %lf",keyword,&the_molec->rad_mult);
    } else if(strstr(keyword,"BOND_TOL")){
      sscanf(instring,"%s %lf",keyword,&the_molec->bond_tol);
    }else if(strstr(keyword,"BREAKING_LINES")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->breaking_lines = temp;
    }else if(strstr(keyword,"CROSSES_ON")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->crosses_on = temp;
    }else if(strstr(keyword,"TUBES_ON")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->tubes_on = temp;
    }else if(strstr(keyword,"DUMMIES_ON")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->dummies_on = temp;
    }else if(strstr(keyword,"AXES_ON")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->axes_on = temp;
    }else if(strstr(keyword,"SYMBOLS_ON")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->symbols_on = temp;
    }else if(strstr(keyword,"FANCY_LINES")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->fancy_lines = temp;
    }else if(strstr(keyword,"HYDROGENS_ON")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->hydrogens_on = temp;
    }else if(strstr(keyword,"SHADING_ON")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->shading_on = temp;
    }else if(strstr(keyword,"OUTLINES_ON")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->outlines_on = temp;
    }else if(strstr(keyword,"DRAW_CONNECTORS")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->draw_connectors = temp;
    }else if(strstr(keyword,"NUMBERS_ON")){
      sscanf(instring,"%s %d",keyword,&temp);
      the_molec->numbers_on = temp;
    } else if(strstr(keyword,"CENT")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->cent.x,
             &whichobj->cent.y,&whichobj->cent.z);
    } else if(strstr(keyword,"TRANS")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->trans.x,
             &whichobj->trans.y,&whichobj->trans.z);
    } else if(strstr(keyword,"SCALE")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->scale.x,
             &whichobj->scale.y,&whichobj->scale.z);
    } else if(strstr(keyword,"ROT")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->rot.x,
             &whichobj->rot.y,&whichobj->rot.z);
    } else{
      fprintf(stderr,"Unrecognized keyword: %s\n",keyword);
    }
    read_error = skipcomments(vfig_file,instring);
  }
  if( read_error >= 0 ) return(1);
}

/****************************************************************************
 *
 *                   Function read_graph_object
 *
 * Arguments: vfig_file: pointer to type FILE
 *             filename: pointer to type char
 * Returns: int
 *
 * Action: reads a graph object out of 'filename, then sets its options
 *    according to the keywords in 'vfig_file.
 *
 *    returns a negative number on failure.
 *
 ****************************************************************************/
int read_graph_object(FILE *vfig_file,char *filename)
{
  int i;
  int temp;
  graph_type *the_graph;
  char instring[MAX_STR_LEN],keyword[80];
  char *start_point;
  int read_error = 0;


  /* start off by reading in the data */
  new_graph(filename);
  if( !whichobj ) return -1;

  the_graph = whichobj->prim->graph;
  /* now read the information out of 'vfig_file */
  read_error = skipcomments(vfig_file,instring);
  while(read_error >= 0 && !strstr(instring,"END_OBJECT")){
    sscanf(instring,"%s",keyword);
    if( strstr(keyword,"XLEGEND") ){
      start_point = &(instring[strspn(instring,"XLEGEND:")]);
      start_point[strlen(start_point)-1] = 0;
      strcpy(the_graph->xlegend,start_point);
    } else if( strstr(keyword,"YLEGEND") ){
      start_point = &(instring[strspn(instring,"YLEGEND:")]);
      start_point[strlen(start_point)-1] = 0;
      strcpy(the_graph->ylegend,start_point);
    } else if( strstr(keyword,"DO_TITLE") ){
      sscanf(instring,"%s %d",keyword,&temp);
      the_graph->do_title = temp;
    } else if( strstr(keyword,"DO_X_TICS") ){
      sscanf(instring,"%s %d",keyword,&temp);
      the_graph->do_x_tics = temp;
    } else if( strstr(keyword,"DO_Y_TICS") ){
      sscanf(instring,"%s %d",keyword,&temp);
      the_graph->do_y_tics = temp;
    } else if( strstr(keyword,"MIN_X") ){
      sscanf(instring,"%s %lf",keyword,&the_graph->min_x);
    } else if( strstr(keyword,"MIN_Y") ){
      sscanf(instring,"%s %lf",keyword,&the_graph->min_y);
    } else if( strstr(keyword,"MAX_X") ){
      sscanf(instring,"%s %lf",keyword,&the_graph->max_x);
    } else if( strstr(keyword,"MAX_Y") ){
      sscanf(instring,"%s %lf",keyword,&the_graph->max_y);
    } else if( strstr(keyword,"TIC_SEP_X") ){
      sscanf(instring,"%s %lf",keyword,&the_graph->tic_sep_x);
    } else if( strstr(keyword,"TIC_SEP_Y") ){
      sscanf(instring,"%s %lf",keyword,&the_graph->tic_sep_y);
    } else if( strstr(keyword,"TIC_START_X") ){
      sscanf(instring,"%s %lf",keyword,&the_graph->tic_start_x);
    } else if( strstr(keyword,"TIC_START_Y") ){
      sscanf(instring,"%s %lf",keyword,&the_graph->tic_start_y);
    } else if( strstr(keyword,"NUM_TICS_X") ){
      sscanf(instring,"%s %lf",keyword,&the_graph->num_tics_x);
    } else if( strstr(keyword,"NUM_TICS_Y") ){
      sscanf(instring,"%s %lf",keyword,&the_graph->num_tics_y);
    } else if( strstr(keyword,"TITLE") ){
      for(i=0;i<NUM_TITLE_LINES;i++){
        fgets(instring,MAX_STR_LEN,vfig_file);
        instring[strlen(instring)-1] = 0;
        strcpy(the_graph->title[i],instring);
      }
    } else if( strstr(keyword,"CURVES") ){
      for(i=0;i<the_graph->num_curves;i++){
        skipcomments(vfig_file,instring);
        sscanf(instring,"%d",&temp);
        the_graph->curves_to_display[i] = temp;
      }
    } else if( strstr(keyword,"STYLES") ){
      for(i=0;i<the_graph->num_curves;i++){
        skipcomments(vfig_file,instring);
        sscanf(instring,"%d",&temp);
        the_graph->styles[i] = temp;
      }
    } else if(strstr(keyword,"CENT")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->cent.x,
             &whichobj->cent.y,&whichobj->cent.z);
    } else if(strstr(keyword,"TRANS")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->trans.x,
             &whichobj->trans.y,&whichobj->trans.z);
    } else if(strstr(keyword,"SCALE")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->scale.x,
             &whichobj->scale.y,&whichobj->scale.z);
    } else if(strstr(keyword,"ROT")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->rot.x,
             &whichobj->rot.y,&whichobj->rot.z);
    } else{
      fprintf(stderr,"Unrecognized keyword: %s\n",keyword);
    }
    read_error = skipcomments(vfig_file,instring);
  }
  if( read_error >= 0 ) return(1);

}


/****************************************************************************
 *
 *                   Function read_band_object
 *
 * Arguments: vfig_file: pointer to type FILE
 *             filename: pointer to type char
 * Returns: int
 *
 * Action: reads a band graph object out of 'filename, then sets its options
 *    according to the keywords in 'vfig_file.
 *
 *    returns a negative number on failure.
 *
 ****************************************************************************/
int read_band_object(FILE *vfig_file,char *filename)
{
  int i;
  int temp;
  band_graph_type *the_graph;
  graph_type *the_data;
  char instring[MAX_STR_LEN],keyword[80];
  char *start_point;
  int read_error = 0;


  /* start off by reading in the data */
  new_band_graph(filename);
  if( !whichobj ) return -1;

  the_graph = whichobj->prim->band_graph;
  the_data = the_graph->the_data;

  /* now read the information out of 'vfig_file */
  read_error = skipcomments(vfig_file,instring);
  while(read_error >= 0 && !strstr(instring,"END_OBJECT")){
    sscanf(instring,"%s",keyword);
    if( strstr(keyword,"XLEGEND") ){
      start_point = &(instring[strspn(instring,"XLEGEND:")]);
      start_point[strlen(start_point)-1] = 0;
      strcpy(the_data->xlegend,start_point);
    } else if( strstr(keyword,"YLEGEND") ){
      start_point = &(instring[strspn(instring,"YLEGEND:")]);
      start_point[strlen(start_point)-1] = 0;
      strcpy(the_data->ylegend,start_point);
    } else if( strstr(keyword,"DO_TITLE") ){
      sscanf(instring,"%s %d",keyword,&temp);
      the_data->do_title = temp;
#ifdef SUPPORT_FLATBANDS
    } else if( strstr(keyword,"FATBAND_FILL") ){
      sscanf(instring,"%s %d",keyword,&the_data->fatband_fill);
    } else if( strstr(keyword,"FATBANDS_ON") ){
      sscanf(instring,"%s %d",keyword,&temp);
      the_data->fatbands_on = temp;
    } else if( strstr(keyword,"FATBAND_SCALE") ){
      sscanf(instring,"%s %lf",keyword,&the_data->fatband_scale);
#endif
    } else if( strstr(keyword,"DO_X_TICS") ){
      sscanf(instring,"%s %d",keyword,&temp);
      the_data->do_x_tics = temp;
    } else if( strstr(keyword,"DO_Y_TICS") ){
      sscanf(instring,"%s %d",keyword,&temp);
      the_data->do_y_tics = temp;
    } else if( strstr(keyword,"MIN_X") ){
      sscanf(instring,"%s %lf",keyword,&the_data->min_x);
    } else if( strstr(keyword,"MIN_Y") ){
      sscanf(instring,"%s %lf",keyword,&the_data->min_y);
    } else if( strstr(keyword,"MAX_X") ){
      sscanf(instring,"%s %lf",keyword,&the_data->max_x);
    } else if( strstr(keyword,"MAX_Y") ){
      sscanf(instring,"%s %lf",keyword,&the_data->max_y);
    } else if( strstr(keyword,"TIC_SEP_X") ){
      sscanf(instring,"%s %lf",keyword,&the_data->tic_sep_x);
    } else if( strstr(keyword,"TIC_SEP_Y") ){
      sscanf(instring,"%s %lf",keyword,&the_data->tic_sep_y);
    } else if( strstr(keyword,"TIC_START_X") ){
      sscanf(instring,"%s %lf",keyword,&the_data->tic_start_x);
    } else if( strstr(keyword,"TIC_START_Y") ){
      sscanf(instring,"%s %lf",keyword,&the_data->tic_start_y);
    } else if( strstr(keyword,"NUM_TICS_X") ){
      sscanf(instring,"%s %lf",keyword,&the_data->num_tics_x);
    } else if( strstr(keyword,"NUM_TICS_Y") ){
      sscanf(instring,"%s %lf",keyword,&the_data->num_tics_y);
    } else if( strstr(keyword,"SHOW_FERMI") ){
      sscanf(instring,"%s %d",keyword,&temp);
      the_graph->show_fermi = temp;
    } else if( strstr(keyword,"FERMI_E") ){
      sscanf(instring,"%s %lf",keyword,&the_graph->Fermi_E);
    } else if( strstr(keyword,"TITLE") ){
      for(i=0;i<NUM_TITLE_LINES;i++){
        fgets(instring,MAX_STR_LEN,vfig_file);
        instring[strlen(instring)-1] = 0;
        strcpy(the_data->title[i],instring);
      }
    } else if( strstr(keyword,"STYLES") ){
        skipcomments(vfig_file,instring);
        sscanf(instring,"%d",&temp);
        the_data->styles[0] = temp;
    } else if(strstr(keyword,"CENT")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->cent.x,
             &whichobj->cent.y,&whichobj->cent.z);
    } else if(strstr(keyword,"TRANS")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->trans.x,
             &whichobj->trans.y,&whichobj->trans.z);
    } else if(strstr(keyword,"SCALE")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->scale.x,
             &whichobj->scale.y,&whichobj->scale.z);
    } else if(strstr(keyword,"ROT")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->rot.x,
             &whichobj->rot.y,&whichobj->rot.z);
    } else{
      fprintf(stderr,"Unrecognized keyword: %s\n",keyword);
    }
    read_error = skipcomments(vfig_file,instring);
  }

  if( read_error >= 0 ) return(1);
}


/****************************************************************************
 *
 *                   Function read_FMO_object
 *
 * Arguments: vfig_file: pointer to type FILE
 *             filename: pointer to type char
 * Returns: int
 *
 * Action: reads an FMO object out of 'filename, then sets its options
 *    according to the keywords in 'vfig_file.
 *
 *    returns a negative number on failure.
 *
 ****************************************************************************/
int read_FMO_object(FILE *vfig_file,char *filename)
{
  int i;
  int temp;
  FMO_diagram_type *the_graph;
  char instring[MAX_STR_LEN],keyword[80];
  char *start_point;
  int read_error = 0;


  /* start off by reading in the data */
  new_FMO_diagram(filename);
  if( !whichobj ) return -1;

  the_graph = whichobj->prim->FMO_diagram;
  /* now read the information out of 'vfig_file */
  read_error = skipcomments(vfig_file,instring);
  while(read_error >= 0 && !strstr(instring,"END_OBJECT")){
    sscanf(instring,"%s",keyword);
    if( strstr(keyword,"XLEGEND") ){
      start_point = &(instring[strspn(instring,"XLEGEND:")]);
      start_point[strlen(start_point)-1] = 0;
      strcpy(the_graph->xlegend,start_point);
    } else if( strstr(keyword,"YLEGEND") ){
      start_point = &(instring[strspn(instring,"YLEGEND:")]);
      start_point[strlen(start_point)-1] = 0;
      strcpy(the_graph->ylegend,start_point);
    } else if( strstr(keyword,"DO_TITLE") ){
      sscanf(instring,"%s %d",keyword,&temp);
      the_graph->do_title = temp;
    } else if( strstr(keyword,"SHOW_BOX") ){
      sscanf(instring,"%s %d",keyword,&temp);
      the_graph->show_box = temp;
    } else if( strstr(keyword,"SHOW_CONNECTS") ){
      sscanf(instring,"%s %d",keyword,&temp);
      the_graph->show_connects = temp;
    } else if( strstr(keyword,"SHOW_DATA") ){
      sscanf(instring,"%s %d",keyword,&temp);
      the_graph->show_data = temp;
    } else if( strstr(keyword,"FILLING_MODE") ){
      sscanf(instring,"%s %d",keyword,&temp);
      the_graph->electron_filling_mode = temp;
    } else if( strstr(keyword,"DO_Y_TICS") ){
      sscanf(instring,"%s %d",keyword,&temp);
      the_graph->do_y_tics = temp;
    } else if( strstr(keyword,"LEFT_FRAGMENT") ){
      sscanf(instring,"%s %d",keyword,&the_graph->left_fragment);
    } else if( strstr(keyword,"RIGHT_FRAGMENT") ){
      sscanf(instring,"%s %d",keyword,&the_graph->right_fragment);
    } else if( strstr(keyword,"LEVEL_WIDTH") ){
      sscanf(instring,"%s %d",keyword,&the_graph->level_width);
    } else if( strstr(keyword,"LEVEL_THICKNESS") ){
      sscanf(instring,"%s %d",keyword,&the_graph->level_thickness);
    } else if( strstr(keyword,"ELECTRON_LENGTH") ){
      sscanf(instring,"%s %d",keyword,&the_graph->electron_length);
    } else if( strstr(keyword,"MIN_Y") ){
      sscanf(instring,"%s %lf",keyword,&the_graph->min_y);
    } else if( strstr(keyword,"MAX_Y") ){
      sscanf(instring,"%s %lf",keyword,&the_graph->max_y);
    } else if( strstr(keyword,"TIC_SEP_Y") ){
      sscanf(instring,"%s %lf",keyword,&the_graph->tic_sep_y);
    } else if( strstr(keyword,"TIC_START_Y") ){
      sscanf(instring,"%s %lf",keyword,&the_graph->tic_start_y);
    } else if( strstr(keyword,"NUM_TICS_Y") ){
      sscanf(instring,"%s %lf",keyword,&the_graph->num_tics_y);
    } else if( strstr(keyword,"TITLE") ){
      for(i=0;i<NUM_TITLE_LINES;i++){
        fgets(instring,MAX_STR_LEN,vfig_file);
        instring[strlen(instring)-1] = 0;
        strcpy(the_graph->title[i],instring);
      }
    } else if( strstr(keyword,"LABELS") ){
      for(i=0;i<the_graph->num_frags;i++){
        fgets(instring,MAX_STR_LEN,vfig_file);
        instring[strlen(instring)-1] = 0;
        strcpy(the_graph->frags[i].label,instring);
      }
    } else if( strstr(keyword,"LABEL") ){
      fgets(instring,MAX_STR_LEN,vfig_file);
      instring[strlen(instring)-1] = 0;
      strcpy(the_graph->label,instring);
    } else if(strstr(keyword,"CENT")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->cent.x,
             &whichobj->cent.y,&whichobj->cent.z);
    } else if(strstr(keyword,"TRANS")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->trans.x,
             &whichobj->trans.y,&whichobj->trans.z);
    } else if(strstr(keyword,"SCALE")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->scale.x,
             &whichobj->scale.y,&whichobj->scale.z);
    } else if(strstr(keyword,"ROT")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->rot.x,
             &whichobj->rot.y,&whichobj->rot.z);
    } else{
      fprintf(stderr,"Unrecognized keyword: %s\n",keyword);
    }
    read_error = skipcomments(vfig_file,instring);
  }

  if( read_error >= 0 ) return(1);
}


/****************************************************************************
 *
 *                   Function read_prop_graph_object
 *
 * Arguments: vfig_file: pointer to type FILE
 *             filename: pointer to type char
 * Returns: int
 *
 * Action: reads a property graph object out of 'filename, then sets it's options
 *    according to the keywords in 'vfig_file.
 *
 *    returns a negative number on failure.
 *
 ****************************************************************************/
int read_prop_graph_object(FILE *vfig_file,char *filename)
{
  int i;
  int temp;
  prop_graph_type *the_graph;
  graph_type *the_data,*the_integration;
  char instring[MAX_STR_LEN],keyword[80];
  char *start_point;
  int read_error = 0;

  /* start off by reading in the data */
  new_prop_graph(filename);
  if( !whichobj ) return -1;

  the_graph = whichobj->prim->prop_graph;
  the_data = the_graph->the_data;
  the_integration = the_graph->the_integration;

  /* now read the information out of 'vfig_file */
  read_error = skipcomments(vfig_file,instring);
  while(read_error >= 0 && !strstr(instring,"END_OBJECT")){
    sscanf(instring,"%s",keyword);
    if( strstr(keyword,"XLEGEND") ){
      start_point = &(instring[strspn(instring,"XLEGEND:")]);
      start_point[strlen(start_point)-1] = 0;
      strcpy(the_data->xlegend,start_point);
    } else if( strstr(keyword,"YLEGEND") ){
      start_point = &(instring[strspn(instring,"YLEGEND:")]);
      start_point[strlen(start_point)-1] = 0;
      strcpy(the_data->ylegend,start_point);
    } else if( strstr(keyword,"DO_TITLE") ){
      sscanf(instring,"%s %d",keyword,&temp);
      the_data->do_title = temp;
    } else if( strstr(keyword,"DO_X_TICS") ){
      sscanf(instring,"%s %d",keyword,&temp);
      the_data->do_x_tics = temp;
    } else if( strstr(keyword,"DO_Y_TICS") ){
      sscanf(instring,"%s %d",keyword,&temp);
      the_data->do_y_tics = temp;
    } else if( strstr(keyword,"MIN_X") ){
      sscanf(instring,"%s %lf",keyword,&the_graph->min_x);
    } else if( strstr(keyword,"MIN_Y") ){
      sscanf(instring,"%s %lf",keyword,&the_graph->min_y);
    } else if( strstr(keyword,"MAX_X") ){
      sscanf(instring,"%s %lf",keyword,&the_graph->max_x);
    } else if( strstr(keyword,"MAX_Y") ){
      sscanf(instring,"%s %lf",keyword,&the_graph->max_y);
    } else if( strstr(keyword,"TIC_SEP_X") ){
      sscanf(instring,"%s %lf",keyword,&the_data->tic_sep_x);
    } else if( strstr(keyword,"TIC_SEP_Y") ){
      sscanf(instring,"%s %lf",keyword,&the_data->tic_sep_y);
    } else if( strstr(keyword,"TIC_START_X") ){
      sscanf(instring,"%s %lf",keyword,&the_data->tic_start_x);
    } else if( strstr(keyword,"TIC_START_Y") ){
      sscanf(instring,"%s %lf",keyword,&the_data->tic_start_y);
    } else if( strstr(keyword,"NUM_TICS_X") ){
      sscanf(instring,"%s %lf",keyword,&the_data->num_tics_x);
    } else if( strstr(keyword,"NUM_TICS_Y") ){
      sscanf(instring,"%s %lf",keyword,&the_data->num_tics_y);
    } else if( strstr(keyword,"INTEGS_FOR_TICS") ){
      sscanf(instring,"%s %d",keyword,&temp);
      the_graph->integs_for_tics = temp;
    } else if( strstr(keyword,"PROP_TYPE") ){
      sscanf(instring,"%s %d",keyword,&temp);
      the_graph->type = temp;
    } else if( strstr(keyword,"SHOW_FERMI") ){
      sscanf(instring,"%s %d",keyword,&temp);
      the_graph->show_fermi = temp;
    } else if( strstr(keyword,"FERMI_E") ){
      sscanf(instring,"%s %lf",keyword,&the_graph->Fermi_E);
    } else if( strstr(keyword,"TITLE") ){
      for(i=0;i<NUM_TITLE_LINES;i++){
        fgets(instring,MAX_STR_LEN,vfig_file);
        instring[strlen(instring)-1] = 0;
        strcpy(the_data->title[i],instring);
      }
    } else if( strstr(keyword,"INTEG_CURVES") ){
      if(the_integration){
        for(i=0;i<the_data->num_curves;i++){
          skipcomments(vfig_file,instring);
          sscanf(instring,"%d",&temp);
          the_integration->curves_to_display[i] = temp;
        }
      }
    } else if( strstr(keyword,"INTEG_STYLES") ){
      for(i=0;i<the_data->num_curves;i++){
        if( the_integration ){
          skipcomments(vfig_file,instring);
          sscanf(instring,"%d",&temp);
          the_integration->styles[i] = temp;
        }
      }
    }else if( strstr(keyword,"CURVES") ){
      for(i=0;i<the_data->num_curves;i++){
        skipcomments(vfig_file,instring);
        sscanf(instring,"%d",&temp);
        the_data->curves_to_display[i] = temp;
      }
    } else if( strstr(keyword,"STYLES") ){
      for(i=0;i<the_data->num_curves;i++){
        skipcomments(vfig_file,instring);
        sscanf(instring,"%d",&temp);
        the_data->styles[i] = temp;
      }
    } else if( strstr(keyword,"FILLS") ){
      for(i=0;i<the_data->num_curves;i++){
        skipcomments(vfig_file,instring);
        sscanf(instring,"%d",&temp);
        the_data->fills[i] = temp;
      }
    } else if(strstr(keyword,"CENT")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->cent.x,
             &whichobj->cent.y,&whichobj->cent.z);
    } else if(strstr(keyword,"TRANS")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->trans.x,
             &whichobj->trans.y,&whichobj->trans.z);
    } else if(strstr(keyword,"SCALE")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->scale.x,
             &whichobj->scale.y,&whichobj->scale.z);
    } else if(strstr(keyword,"ROT")){
      sscanf(instring,"%s %lf %lf %lf",keyword,&whichobj->rot.x,
             &whichobj->rot.y,&whichobj->rot.z);
    } else{
      fprintf(stderr,"Unrecognized keyword: %s\n",keyword);
    }

    read_error = skipcomments(vfig_file,instring);
  }
  return read_error;

  if( read_error >= 0 ) return(1);
}



/****************************************************************************
 *
 *                   Function read_object
 *
 * Arguments: vfig_file: pointer to type FILE
 * Returns: int
 *
 * Action: gets space for an object and reads it out of 'vfig_file
 *    returns a negative number on failure.
 *
 ****************************************************************************/
int read_object(FILE *vfig_file)
{
  char keyword[80];
  char instring[MAX_STR_LEN],foostring[80];
  char filename[80];
  int success = 0;
  int read_error;
  int object_type;

  /* start by finding out what kind of object this is */
  read_error = skipcomments(vfig_file,instring);
  if( read_error >= 0 ){
    sscanf(instring,"%s %s\n",keyword,foostring);
    if( strstr(keyword,"TYPE:") ){
      if(strstr(foostring,"MOLECULE"))
        object_type = MOLECULE;
      else if(strstr(foostring,"MO_SURF"))
        object_type = MO_SURF;
      else if(strstr(foostring,"PROPERTY_GRAPH"))
        object_type = PROP_GRAPH;
      else if(strstr(foostring,"BAND_GRAPH"))
        object_type = BAND_GRAPH;
      /* this case must be near the bottom */
      else if(strstr(foostring,"GRAPH"))
        object_type = GRAPH;
      else if(strstr(foostring,"FMO"))
        object_type = FMO_DIAGRAM;
      else{
        error("Unrecognized object type hit in read_object.");
        success = -1;
      }
    }else{
      error("The first field in an object must be the type");
      success = -1;
    }

    if(success >= 0 && read_error >= 0 ){
      if( skipcomments(vfig_file,instring) >= 0){
        if( strstr(instring,"FILENAME:") ){
          sscanf(instring,"%s %s",keyword,filename);
        }else{
          error("The second field in an object must be the filename");
          success = -1;
        }
      } else{
        read_error = -1;
      }
    }
  }

  if( read_error >= 0 && success >= 0 ){
    switch(object_type){
    case MOLECULE:
      success = read_molec_object(vfig_file,filename);
      break;
    case MO_SURF:
      success = read_MO_surf_object(vfig_file,filename);
      break;
    case PROP_GRAPH:
      success = read_prop_graph_object(vfig_file,filename);
      break;
    case BAND_GRAPH:
      success = read_band_object(vfig_file,filename);
      break;
    case GRAPH:
      success = read_graph_object(vfig_file,filename);
      break;
    case FMO_DIAGRAM:
      success = read_FMO_object(vfig_file,filename);
      break;
    }
  }
  if( read_error < 0 ) success = -1;

  return success;
}




/****************************************************************************
 *
 *                   Procedure read_all
 *
 * Arguments: none
 * Returns: none
 *
 * Action: reads in the contents of an input file
 *
 ****************************************************************************/
void read_all(char *file_name)
{
  FILE *vfig_file;
  char *theinline;
  char filename[240];
  char keyword[80];
  char instring[MAX_STR_LEN],foostring[80];
  int num_read,temp;

  if( !file_name ){
#ifndef USE_READLINE
    printf("Enter input file name: ");
    scanf("%s\n",filename);
#else
    theinline= readline("Enter input file name: ");
    add_history(theinline);
    if( theinline ){
      sscanf(theinline,"%s",filename);
      free(theinline);
    } else {
      error("Bad file name");
      filename[0] = 0;
    }
#endif
  } else {
    strcpy(filename,file_name);
  }

  vfig_file = fopen(filename,"r");
  if(!vfig_file){
    error("Can't open input file.");
    return;
  }

  num_read = 0;
  while(skipcomments(vfig_file,instring) >= 0){
    sscanf(instring,"%s",keyword);
    if(strstr(keyword,"FILL_PROJECTIONS")){
      sscanf(instring,"%s %d",foostring,&temp);
      fill_projections = temp;
    }
    else if(strstr(keyword,"PS_WHERE_TO_PRINT")){
      sscanf(instring,"%s %d",foostring,&temp);
      PS_options.where_to_print = temp;
    }
    else if(strstr(keyword,"PS_FONTNAME"))
      sscanf(instring,"%s %s",foostring,PS_options.fontname);
    else if(strstr(keyword,"PS_FONTSIZE"))
      sscanf(instring,"%s %lf",foostring,&PS_options.fontsize);
    else if(strstr(keyword,"PS_PRINTSCALE"))
      sscanf(instring,"%s %lf",foostring,&PS_options.printscale);
    else if(strstr(keyword,"ATOM_TYPE"))
      sscanf(instring,"%s %d",foostring,&PS_options.atom_sphere_type);
    else if(strstr(keyword,"BOND_TYPE"))
      sscanf(instring,"%s %d",foostring,&PS_options.bond_type);
    else if(strstr(keyword,"BEGIN_OBJECT")){
      if( read_object(vfig_file) >= 0 ){
        num_read++;
      }
    }
  }
  fclose(vfig_file);
  printf("%d objects read\n",num_read);
}


/****************************************************************************
 *
 *                   Procedure dump_3D_objects
 *
 * Arguments: prim: pointer to prim_type
 *             obj: pointer to object_type
 *
 * Returns: none
 *
 * Action: Dumps the 3D objects contained in 'prim into an input file
 *         for rayshade.
 *
 *    At the moment, this is not particularly sophisticated in terms of
 *    how it does it's work.  All atoms are the same color, etc...
 *    in the future this should probably be remedied.
 *
 ****************************************************************************/
void dump_3D_objects(prim_type *prim,object_type *obj)
{
  char filename[80],*theinline;
  int i,j;
  int objects_so_far;
  molec_type *molec;
  MO_surface_type *surf;
  int num_atoms,num_triangles;
  atom_type *atom,*atom2,*atoms;
  int num_unique;
  char found;
  atom_type *unique_atoms;
  triangle_type *triangles,*triangle;
  vertex_type *vertices,*v1,*v2,*v3;


  /* just return if there's no 3D object here */
  if(!prim || !obj ||
     !(prim->which==MOLECULE || prim->which == MO_SURF)){
    return;
  }


  display("Look in xterm!");
#ifndef USE_READLINE
  printf("Enter rayshade input file name: ");
  scanf("%s\n",filename);
#else
  theinline= readline("Enter rayshade input file name: ");
  add_history(theinline);
  if( theinline ){
    sscanf(theinline,"%s",filename);
    free(theinline);
  } else {
    error("Bad file name");
    filename[0] = 0;
  }
#endif

  rayfile = fopen(filename,"w+");
  if(!rayfile){
    error("Can't open output file.");
    return;
  }


  fprintf(rayfile,"eyep  0 0 -100\n");
  fprintf(rayfile,"lookp  0 0 0\n");
  fprintf(rayfile,"up  0 1 0\n");
  fprintf(rayfile,"fov  2\n");
  fprintf(rayfile,"screen 256 256 \n");
  fprintf(rayfile,"light 1 point 50 50 -100\n");
  fprintf(rayfile,"background 1 1 1\n");
  fprintf(rayfile,"surface s1\n");
  fprintf(rayfile,"        ambient 0.25 0.25 0.25\n");
  fprintf(rayfile,"        diffuse 0.55 0.55 0.55 \n");
  fprintf(rayfile,"        specular 0.2 0.2 0.2 \n");
  fprintf(rayfile,"        specpow 3 \n");
  fprintf(rayfile,"        reflect 0.1 \n");
  fprintf(rayfile,"        transp 0.3 \n");
  fprintf(rayfile,"surface s2\n");
  fprintf(rayfile,"          ambient 0.15 0.15 0.15 \n");
  fprintf(rayfile,"        diffuse 0.25 0.25 0.25 \n");
  fprintf(rayfile,"        specular 0.05 0.05 0.05  \n");
  fprintf(rayfile,"        specpow 3  \n");
  fprintf(rayfile,"        reflect 0.1 \n");
  fprintf(rayfile,"        transp 0.3 \n");
  fprintf(rayfile,"surface atomsurf\n");
  fprintf(rayfile,"        ambient 0.35 0.15 0.15\n");
  fprintf(rayfile,"        diffuse 0.65 0.15 0.15\n");
  fprintf(rayfile,"        specular 0.6 0.6 0.6 \n");
  fprintf(rayfile,"        specpow 3 \n");
  fprintf(rayfile,"        reflect 0.3 \n");

  /* write the atomic surface values */
  molec = 0;
  switch(prim->which){
  case MO_SURF:
    if(prim->MO_surf->display_molec) molec = prim->MO_surf->molec;
    break;
  case MOLECULE:
    molec = prim->molec;
    break;
  }
  if(molec){
    num_atoms = molec->num_atoms;
    num_unique = 0;
    unique_atoms = (atom_type *)D_CALLOC(num_atoms,sizeof(atom_type));
    if( !unique_atoms ) fatal("Can't get space for unique_atoms");
    for(i=0;i<num_atoms;i++){
      found = 0;
      for(j=0;j<num_unique;j++){
        if( !strcmp(unique_atoms[j].type,molec->atoms[i].type) ){
          found = 1;
        }
      }
      if( ! found ){
        safe_strcpy(unique_atoms[num_unique].type,molec->atoms[i].type);
        if( !strstr(unique_atoms[num_unique].type,"&") ){
          fprintf(rayfile,"surface %ssurf\n",unique_atoms[num_unique].type);
        }
        else{
          fprintf(rayfile,"surface Dummysurf\n");
        }
        fprintf(rayfile,"        ambient 0.35 0.15 0.15\n");
        fprintf(rayfile,"        diffuse 0.65 0.15 0.15\n");
        fprintf(rayfile,"        specular 0.6 0.6 0.6 \n");
        fprintf(rayfile,"        specpow 3 \n");
        fprintf(rayfile,"        reflect 0.3 \n");
        num_unique++;
      }
    }
  }


  fprintf(rayfile,"surface cylsurf\n");
  fprintf(rayfile,"        ambient 0.25 0.25 0.25\n");
  fprintf(rayfile,"        diffuse 0.45 0.45 0.45\n");
  fprintf(rayfile,"        specular 0.1 0.1 0.1 \n");
  fprintf(rayfile,"        specpow 3 \n");
  fprintf(rayfile,"        reflect 0.3 \n");


  printf("Dumping objects to  rayshade file\n");

  molec = 0;
  switch(prim->which){
  case MO_SURF:
    surf = prim->MO_surf;
    if( surf->display_molec ){
      molec = surf->molec;
    }
    if( surf->display_surf ){
      fprintf(rayfile,"name MOsurf grid 3 3 3\n");
      triangles = surf->triangles;
      vertices = surf->triangle_vertices;
      num_triangles = surf->num_triangles;
      for(i=0;i<num_triangles;i++,objects_so_far++){
        triangle = &(triangles[i]);
        v1 = &vertices[triangle->vertices[0]];
        v2 = &vertices[triangle->vertices[1]];
        v3 = &vertices[triangle->vertices[2]];

        if( triangle->color ){
          fprintf(rayfile,"triangle s1 %4.2lf %4.2lf %4.2lf %4.2lf %4.2lf %4.2lf %4.2lf %4.2lf %4.2lf\n",
                  v1->position.x,v1->position.y,v1->position.z,
                  v2->position.x,v2->position.y,v2->position.z,
                  v3->position.x,v3->position.y,v3->position.z);
        }else{
          fprintf(rayfile,"triangle s2 %4.2lf %4.2lf %4.2lf %4.2lf %4.2lf %4.2lf %4.2lf %4.2lf %4.2lf\n",
                  v1->position.x,v1->position.y,v1->position.z,
                  v2->position.x,v2->position.y,v2->position.z,
                  v3->position.x,v3->position.y,v3->position.z);
        }

      }

      fprintf(rayfile,"end\n");
      fprintf(rayfile,"object MOsurf\n");
      fprintf(rayfile,"rotate 0 0 1 %4.2lf\n",180.0*obj->rot.z/PI);
      fprintf(rayfile,"rotate 0 1 0 %4.2lf\n",180.0*obj->rot.y/PI);
      fprintf(rayfile,"rotate 1 0 0 %4.2lf\n",180.0*obj->rot.x/PI);

      fprintf(rayfile,"scale %4.2lf %4.2lf %4.2lf\n",obj->scale.x,obj->scale.y,obj->scale.z);
      fprintf(rayfile,"translate %4.2lf %4.2lf %4.2lf\n",obj->trans.x,obj->trans.y,obj->trans.z);
    }
    break;
  case MOLECULE:
    molec = prim->molec;
    break;
  }

  if( molec ){
    num_atoms = molec->num_atoms;
    fprintf(rayfile,"name molec grid 10 10 10\n");
    if(molec->num_frames > 1){
      atoms = &(molec->atoms[(molec->current_frame%molec->num_frames)*molec->num_atoms]);
    }
    else{
      atoms = molec->atoms;
    }
    for(i=0;i<num_atoms;i++){
      atom = &(atoms[i]);
      if( !atom->exclude && (molec->hydrogens_on ||
                             atom->type[0] != 'H' || atom->type[1] != 0) &&
         (molec->dummies_on || atom->type[0] != '&') ){
        if( !strstr(atom->type,"&") ){
          fprintf(rayfile,"sphere %ssurf %4.2lf %4.2lf %4.2lf %4.2lf\n",
                  atom->type,
                  atom->rad*.5*molec->rad_mult,atom->loc.x,atom->loc.y,atom->loc.z);
        }
        else{
          fprintf(rayfile,"sphere Dummysurf %4.2lf %4.2lf %4.2lf %4.2lf\n",
                  atom->rad*.5*molec->rad_mult,atom->loc.x,atom->loc.y,atom->loc.z);
        }
        if( (molec->draw_connectors && atom->num_lines_out) ){
          for(j=0;j<atom->num_lines_out;j++){
            if( atom->linesto[j] < i ){
              atom2 = &(atoms[atom->linesto[j]]);
              fprintf(rayfile,"/* %s - %s */\n",atom->type,atom2->type);
              fprintf(rayfile,"cylinder cylsurf 0.075 %4.2lf %4.2lf %4.2lf %4.2lf %4.2lf %4.2lf\n",
                      atom->loc.x,atom->loc.y,atom->loc.z,
                      atom2->loc.x,atom2->loc.y,atom2->loc.z);
            }
          }
        }
      }
    }
    fprintf(rayfile,"end\n");
    fprintf(rayfile,"object molec\n");
    fprintf(rayfile,"rotate 0 0 1 %4.2lf\n",180.0*obj->rot.z/PI);
    fprintf(rayfile,"rotate 0 1 0 %4.2lf\n",180.0*obj->rot.y/PI);
    fprintf(rayfile,"rotate 1 0 0 %4.2lf\n",180.0*obj->rot.x/PI);

    fprintf(rayfile,"scale %4.2lf %4.2lf %4.2lf\n",obj->scale.x,obj->scale.y,obj->scale.z);
    fprintf(rayfile,"translate %4.2lf %4.2lf %4.2lf\n",obj->trans.x/20.0,obj->trans.y/20.0,
            obj->trans.z/20.0);
  }

  fclose(rayfile);
  printf("done\n");
}


/****************************************************************************
 *
 *                   Procedure dump_VRML
 *
 * Arguments: prim: pointer to prim_type
 *             obj: pointer to object_type
 *
 * Returns: none
 *
 * Action: Dumps the 3D objects contained in 'prim into an input file
 *         for VRML.
 *
 *    At the moment, this is not particularly sophisticated in terms of
 *    how it does it's work.  All atoms are the same color, etc...
 *    in the future this should probably be remedied.
 *
 ****************************************************************************/
void dump_VRML(prim_type *prim,object_type *obj)
{
  FILE *outfile;
  char *stringarr[3];
  char instring[80],label[80],outstring[80];
  char filename[80],*theinline;
  float theta_y,theta_z,len,innerlen;
  int i,j;
  molec_type *molec;
  MO_surface_type *surf;
  int num_atoms;
  atom_type *atom,*atom2,*atoms;
  int num_unique;
  int new_file;
  point_type translation,color;
  char found;
  atom_type *unique_atoms;
  float LOD_range=10;
  char use_LOD=1;
  char do_axes=1;
  line_type *the_line;

  /* just return if there's no 3D object here */
  if(!prim || !obj ||
     !(prim->which==MOLECULE || prim->which == MO_SURF)){
    return;
  }


  display("Look in xterm!");
#ifndef USE_READLINE
  printf("Enter VRML input file name: ");
  scanf("%s\n",filename);
#else
  theinline= readline("Enter VRML input file name: ");
  add_history(theinline);
  if( theinline ){
    sscanf(theinline,"%s",filename);
    free(theinline);
  } else {
    error("Bad file name");
    filename[0] = 0;
  }
#endif

  outfile = fopen(filename,"r");
  if(!outfile){
    new_file = 1;
  } else {
    fclose(outfile);
    printf("file: %s exits already.\n",filename);
    strcpy(instring,"n");
    readcharparm("overwrite file (otherwise we append)",instring);
    if( instring[0] == 'y' || instring[0] == 'Y'){
      new_file = 1;
    } else{
      new_file = 0;
    }
  }

  strcpy(instring,"n");
  readcharparm("Use LOD with this object? (speeds drawing)",instring);
  if( instring[0] == 'y' || instring[0] == 'Y'){
    use_LOD = 1;
    LOD_range = 10.0;
    readfloatparm("LOD range",&LOD_range);
  } else{
    use_LOD = 0;
  }


  strcpy(instring,"y");
  readcharparm("show an axis system?",instring);
  if( instring[0] == 'n' || instring[0] == 'N'){
    do_axes = 0;
  } else{
    do_axes = 1;
  }


  if( new_file ) outfile = fopen(filename,"w+");
  else outfile = fopen(filename,"a");
  if(!outfile){
    error("Can't open output file.");
    return;
  }

  stringarr[0] = instring;

  if( !new_file ){
    strcpy(instring,"10 10 -10");
    readstringparm("translation of this object",stringarr);
    sscanf(instring,"%lf %lf %lf",&translation.x,&translation.y,
           &translation.z);
  } else{
    i=0;
    while(VRML_header[i] != 0){
      fprintf(outfile,"%s\n",VRML_header[i]);
      i++;
    }
  }

  /* write the atomic surface values */
  molec = 0;
  switch(prim->which){
  case MO_SURF:
    if(prim->MO_surf->display_molec) molec = prim->MO_surf->molec;
    break;
  case MOLECULE:
    molec = prim->molec;
    break;
  }
  if(molec){
    num_atoms = molec->num_atoms;
    num_unique = 0;
    unique_atoms = (atom_type *)D_CALLOC(num_atoms,sizeof(atom_type));
    if( !unique_atoms ) fatal("Can't get space for unique_atoms");
    for(i=0;i<num_atoms;i++){
      found = 0;
      for(j=0;j<num_unique;j++){
        if( !strcmp(unique_atoms[j].type,molec->atoms[i].type) ){
          found = 1;
        }
      }
      if( ! found ){
        safe_strcpy(unique_atoms[num_unique].type,molec->atoms[i].type);
        if( !strstr(unique_atoms[num_unique].type,"&") ){
          fprintf(outfile,"PROTO %sball [ field SFFloat MyRadius %5.4lf\n",
                  unique_atoms[num_unique].type,
                  0.5*molec->atoms[i].rad*molec->rad_mult);
          sprintf(outstring,"color of atom type: %s",
                  unique_atoms[num_unique].type);
          strcpy(instring,(char *)"0.6 0.6 0.6");
          readstringparm(outstring,stringarr);
          sscanf(instring,"%lf %lf %lf",&color.x,&color.y,&color.z);
          if( color.x < 0 ) color.x = 0;
          if( color.x > 1 ) color.x = 1;
          if( color.y < 0 ) color.y = 0;
          if( color.y > 1 ) color.y = 1;
          if( color.z < 0 ) color.z = 0;
          if( color.z > 1 ) color.z = 1;
          fprintf(outfile,
                  "       field SFColor MyColor %4.2lf %4.2lf %4.2lf]\n",
                  color.x,color.y,color.z);
        }
        else{
          fprintf(outfile,"PROTO Dummyball [ field SFFloat MyRadius 0.1\n");
          fprintf(outfile,"       field SFColor MyColor 0.6 0.6 0.6]\n");
        }
        fprintf(outfile,"{Shape {\n");
        fprintf(outfile,"geometry Sphere { radius IS MyRadius }\n");
        fprintf(outfile,"appearance Appearance {\n");
        fprintf(outfile,"material Material {diffuseColor IS MyColor}}}}\n");
        num_unique++;
      }
    }
  }

  switch(prim->which){
  case MO_SURF:
    strcpy(label,prim->MO_surf->filename); break;
  case MOLECULE:
    strcpy(label,prim->molec->filename); break;
  }

  stringarr[0] = label;
  readstringparm("molecule label",stringarr);


  if( new_file ){
    printf("Dumping objects to VRML file\n");
  } else{
    printf("Adding objects to VRML file\n");
    fprintf(outfile,
            "Transform {translation %4.2lf %4.2lf %4.2lf children [\n",
            translation.x,translation.y,translation.z);
  }

  if( use_LOD ){
    fprintf(outfile,"LOD{ center 0 0 0 range [%4.2lf]\n",LOD_range);
    fprintf(outfile,"level [ \n");
    fprintf(outfile,"Group{ children [\n");
  }


  molec = 0;
  switch(prim->which){
  case MO_SURF:
    surf = prim->MO_surf;
    if( surf->display_molec ){
      molec = surf->molec;
    }
    if( surf->display_surf ){
    }
    break;
  case MOLECULE:
    molec = prim->molec;
    break;
  }

  if( molec ){
    num_atoms = molec->num_atoms;
    if(molec->num_frames > 1){
      atoms = &(molec->atoms[(molec->current_frame%molec->num_frames)*molec->num_atoms]);
    }
    else{
      atoms = molec->atoms;
    }
    for(i=0;i<num_atoms;i++){
      atom = &(atoms[i]);
      if( !atom->exclude && (molec->hydrogens_on ||
                             atom->type[0] != 'H' || atom->type[1] != 0) &&
         (molec->dummies_on || atom->type[0] != '&') ){
        fprintf(outfile,
                "Transform { translation %4.2lf %4.2lf %4.2lf children[\n",
                atom->loc.x,atom->loc.y,atom->loc.z);
        if( !strstr(atom->type,"&") ){
          fprintf(outfile,"%sball {}",atom->type);
        }
        else{
          fprintf(outfile,"Dummyball{}\n");
        }
        fprintf(outfile,"]}\n");
        if( (molec->draw_connectors && atom->num_lines_out) ){
          for(j=0;j<atom->num_lines_out;j++){
            if( atom->linesto[j] < i ){
              atom2 = &(atoms[atom->linesto[j]]);
              the_line = find_the_line(atom->num,atom2->num,molec);
              if(!the_line) FATAL_BUG("whoa!  cant't find the line.");
              if( !the_line->custom || the_line->drawn ){
                angles_from_bond_vect(&(atom->loc),&(atom2->loc),
                                      &theta_y,&theta_z,&len);
                fprintf(outfile,"# /* %s - %s */\n",atom->type,atom2->type);
                fprintf(outfile,
                        "Transform { translation %5.4lf %5.4lf %5.4lf \n",
                        atom->loc.x,atom->loc.y,atom->loc.z);
                fprintf(outfile,
                        "\trotation 0 1 0 %4.2lf children[\n",
                        -theta_y);
                fprintf(outfile,
                        "Transform { rotation 0 0 1 %4.2lf children[\n",
                        -theta_z);
                innerlen = len - 0.5*molec->rad_mult*(atom->rad+atom2->rad);
                fprintf(outfile,
                        "Transform { translation 0 %5.4lf 0 children[\n",
                        innerlen / 2.0 + 0.5*molec->rad_mult*atom->rad);
                fprintf(outfile,
                        "diffcyl {MyHeight %5.4lf} ]}  ]}  ]}\n",
                        innerlen);
              }
            }
          }
        }
      }
    }
  }

  fprintf(outfile,"Transform { translation 5. 0. 0. children[\n");
  fprintf(outfile,
          " Billboard { children[Shape{geometry Text {string \"%s\"}}]}]}\n",
          label);

  if( do_axes ){
    fprintf(outfile,"Transform {translation -5 0 0 children[axes{}]}\n");
  }

  if( use_LOD ){
    fprintf(outfile,"]}\n");
    fprintf(outfile,"Transform{ translation %4.2lf %4.2lf %4.2lf children [\n",
            translation.x,translation.y,translation.z);
    fprintf(outfile,"Shape{ geometry Text {string \"%s\"}}]}\n",
            label);
    fprintf(outfile,"]}\n");
  }
  if( !new_file ){
    fprintf(outfile,"]}\n");
  }

  i=0;
  while(VRML_trailer[i] != 0){
    fprintf(outfile,"%s\n",VRML_trailer[i]);
    i++;
  }

  fclose(outfile);
  printf("done\n");
}


