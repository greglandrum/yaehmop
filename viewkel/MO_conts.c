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

#define CLOSE_TO_ZERO 1e-3

/********

  this has got the stuff for dealing with contour plots of MOs
   (the 3D kind)
*********/

/***
  Recent Edit History:
   10.05.98 gL:
     added support for arbitrary MO planes
   26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
      under gcc (yeah yeah... it's anal)
     extensive modifications to use symmetry (mirror planes at least)
      in MO plots
   22.01.99 gL:
     dump_grids functionality added.
   28.06.1999 gL:
     Fixed EVIL overflow in drawing of hidden lines plots. (Note: it was
     of course just a silly mistake, but it was a cast iron bitch to track
     down.  I'd be proud of myself if the initial mistake weren't so stupid.)
     Fixed some core leaks in contour plot stuff.
  29.06.1999 gL:
     added support for dumping MO volumes (binary files which contain MO data)

***/

#include "viewkel.h"

#define MY_MAX(a,b)  ((a)>(b)?(a):(b))
#define MY_MIN(a,b)  ((a)<(b)?(a):(b))

#ifdef DEBUG_HIDDEN_LINE
extern void g_clines PROTO((XPoint *,point_type2D *,int,char));
#endif


/****************************************************************************
 *
 *                   Function determine_plane_bounds
 *
 * Arguments:  cont_plot: pointer to MO_contour_plot_type
 *                 atoms: pointer to atom_type
 *             num_atoms: int
 *   width,height,offset: floats
 *
 * Returns: int
 *
 * Action:  Determines the bounding box for the arbitrary plane
 *   determined by the three selected atoms.  This routine should
 *   not be called unless there are three atoms selected.
 *
 *  The return value is 0 for success, 1 for failure
 *   (failure due to linear points, etc)
 *
 *  The plane is calculated as follows.
 *   The three selected atoms are P1,P2,P3 (selection order)
 *   V1 = P1 - P2;  V2 = P2 - P3;
 *   the plane normal is Normal = V1 x V2
 *   The (orthonormal) basis vectors for the plane are then:
 *    Bas1 = V1/|V1|; Bas2 = Bas1 x Normal/|Normal|
 *
 ****************************************************************************/
int determine_plane_bounds( MO_contour_plot_type *cont_plot,
                            atom_type *atoms, int num_atoms,
                            float width, float height, float offset)
{
  int i;
  atom_type *selected_atoms[3];
  int num_found;
  point_type V1,V2,Normal;
  point_type *Bas1,*Bas2,*Origin;

  if( num_selected != 3 )
    FATAL_BUG("determine_plane_bounds called with num_selected != 3");

  num_found = 0;
  /* find the three selected atoms */
  for( i=0;i<num_atoms;i++){
    if( atoms[i].is_selected ){
      selected_atoms[atoms[i].is_selected-1] = &(atoms[i]);
      num_found++;
    }
  }
  if( num_found != 3 )
    FATAL_BUG("can't find all selected atoms in determine_plane_bounds");

  Bas1 = &cont_plot->Bas1;
  Bas2 = &cont_plot->Bas2;
  Origin = &cont_plot->Origin;

  /* determine the vectors connecting the atoms */
  V3Sub(&(selected_atoms[0]->loc),&(selected_atoms[1]->loc),
        &V1);
  V3Sub(&(selected_atoms[2]->loc),&(selected_atoms[1]->loc),
        &V2);
  /* normalize them */
  V3Normalize(&V1);
  V3Normalize(&V2);

  /* make sure they aren't colinear */
  if( (fabs(V3Dot(&V1,&V2))-1.0) > -1e-4 ){
    fprintf(stderr,"The atoms are colinear... this is bad.\n");
    return(1);
  }

  /****
    calculate the plane normal, this is
     normalized because V1 and V2 are
  ****/
  V3Cross(&V1,&V2,&Normal);

  /* and determine the "basis vectors" for the plane */
  Bas1->x=V1.x;Bas1->y=V1.y;Bas1->z=V1.z;
  V3Cross(Bas1,&Normal,Bas2);

  /****
    okay, that was pretty easy, now that we have the definition
     of the plane, we can generate the bounding box.  Put
     the middle atom in the center of the plane.
   ****/
  Origin->x=selected_atoms[1]->loc.x;
  Origin->y=selected_atoms[1]->loc.y;
  Origin->z=selected_atoms[1]->loc.z;
  cont_plot->left_top_corner.x = Origin->x -
    (width * Bas1->x + height * Bas2->x)/2.0 +
    offset * Normal.x;
  cont_plot->left_top_corner.y = Origin->y -
    (width * Bas1->y + height * Bas2->y)/2.0 +
    offset * Normal.y;
  cont_plot->left_top_corner.z = Origin->z -
    (width * Bas1->z + height * Bas2->z)/2.0 +
    offset * Normal.z;

  cont_plot->right_bottom_corner.x = Origin->x +
    (width * Bas1->x + height * Bas2->x)/2.0 +
    offset * Normal.x;
  cont_plot->right_bottom_corner.y = Origin->y +
    (width * Bas1->y + height * Bas2->y)/2.0 +
    offset * Normal.y;
  cont_plot->right_bottom_corner.z = Origin->z +
    (width * Bas1->z + height * Bas2->z)/2.0 +
    offset * Normal.z;

  /*******

    okay, now the final trick... divide the two basis
    vectors by the number of steps in each direction so
    that we can determine a location by simple multiplication.

    crafty eh?

  ********/
  V3Scale(Bas1,width/(float)cont_plot->num_a);
  V3Scale(Bas2,height/(float)cont_plot->num_b);


  /* I believe that is that */
  return(0);
}

#ifdef SUPPORT_VOLUMES
#define CONFIRM_STRING "VIEWKEL_MO_DATA"
#define BIN_WRITE(_a_,_b_,_c_) {numwritten=write(_a_,_b_,_c_);\
  if(numwritten!=_c_)fatal("numwritten!=size");}


/****************************************************************************
 *
 *                   Procedure eval_MO_volume
 *
 * Arguments:  surf: pointer to MO_surface_type
 *
 *
 * Returns: none
 *
 * Action:  Calculates the value of an MO in a user specified volume
 *   and dumps it to a binary file on disk.
 *
 *  File format:
 *   Header:
 *     Confirmation string: "VIEWKEL_MO_DATA"
 *     doubles: xrange,yrange,zrange
 *     ints: numx,numy,numz
 *   Data:
 *     numx*numy*numz doubles in order:
 *      x[0]y[0]z[0], x[0]y[0]z[1]... x[0]y[0]z[numz]
 *      x[0]y[1]z[0], x[0]y[1]z[1]... x[0]y[1]z[numz]
 *      ...
 *      x[numx]y[numy]z[0], x[numx]y[numy]z[1]... x[numx]y[numy]z[numz]
 *
 * Notes: the volume is centered on the origin
 *
 ****************************************************************************/
void eval_MO_volume(MO_surface_type *surf)
{
  double xr,yr,zr,*data;
  float stepx,stepy,stepz,tempf;
  MO_info_type MO_info;
  int numx,numy,numz;
  char tempfilename[240],*stringarr[3];
  char confirm_string[80];
  int outfile,confirm;
  int numwritten,size;
  int i,j,k;

  confirm = 0;
  readintparm("Dump a volume? ",&confirm);
  if( !confirm ){
    display("cancelled!");
    return;
  }
  sprintf(tempfilename,"%s.VOL",surf->filename);
  stringarr[0] = tempfilename;
  readstringparm("Volume file name:",stringarr);

  /* open the file */
#ifndef USING_THE_MAC
#include <sys/mode.h>
  outfile = open(tempfilename,
                 O_RDWR|O_TRUNC|O_CREAT,S_IRUSR|S_IWUSR);
#else
  outfile = open(tempfilename,O_RDWR|O_TRUNC|O_CREAT);
#endif

  strcpy(confirm_string,CONFIRM_STRING);

  size = strlen(CONFIRM_STRING)*sizeof(char);
  BIN_WRITE(outfile,confirm_string,size);

  tempf = 10;
  readfloatparm("x range",&tempf);
  xr = tempf;
  numx = 100;
  readintparm("numx",&numx);
  readfloatparm("y range",&tempf);
  yr = tempf;
  numy = 100;
  readintparm("numy",&numx);
  readfloatparm("z range",&tempf);
  zr = tempf;
  numz = 100;
  readintparm("numz",&numz);
  size = sizeof(double);
  BIN_WRITE(outfile,&xr,size);
  BIN_WRITE(outfile,&yr,size);
  BIN_WRITE(outfile,&zr,size);
  size = sizeof(int);
  BIN_WRITE(outfile,&numx,size);
  BIN_WRITE(outfile,&numy,size);
  BIN_WRITE(outfile,&numz,size);


  /* the header is written... allocate the data array and go */
  data = (double *)D_CALLOC(numz,sizeof(double));
  if(!data)fatal("can't allocate data array");

  stepx = (float)xr /(float)numx;
  stepy = (float)yr /(float)numy;
  stepz = (float)zr /(float)numz;

  size = numz * sizeof(double);

  MO_info.loc.x = -xr / 2.0;
  for(i=0;i<numx;i++,MO_info.loc.x += stepx){
    printf("X=%lf\n",MO_info.loc.x);
    MO_info.loc.y = -yr / 2.0;
    for(j=0;j<numy;j++,MO_info.loc.y += stepy){

      MO_info.loc.z = -zr / 2.0;
      for(k=0;k<numy;k++,MO_info.loc.z += stepz){
        calc_MO_value(surf->active_MO,&MO_info,surf->MO_centers,
                      surf->num_centers,surf->adf_plot);
        data[k] = MO_info.val;
      }
      BIN_WRITE(outfile,data,size);
    }
  }
  close(outfile);
  D_FREE(data);
}
#endif

/****************************************************************************
 *
 *                   Procedure eval_MO_plane
 *
 * Arguments:  surf: pointer to MO_surface_type
 *
 *
 * Returns: none
 *
 * Action:  Evaluates the value of the MO in the plane
 *          defined in 'surf.
 *
 ****************************************************************************/
void eval_MO_plane(MO_surface_type *surf)
{
  int i,j;
  int symm_i,symm_j;
  int use_symm_a,use_symm_b;
  int num_a_to_do,num_b_to_do;
  MO_contour_plot_type *cont_plot;
  iso_curve_type *data,*forward_ptr,*back_ptr;
  int num_isocurves;
  MO_info_type MO_info;
  point_type step;
  point_type *origin;

  cont_plot = surf->MO_contours;

  /* clean out any old data if it's there and is too small */
  data = cont_plot->data;
  num_isocurves=0;
  while(data){
    if( !(data->points) ||  data->p_max < cont_plot->num_b ){
      if( data->points ) D_FREE(data->points);
      data->p_max = cont_plot->num_b;
      data->points = (point_type *)D_CALLOC(data->p_max,sizeof(point_type));
      if(!data->points)fatal("can't get memory for data->points");
    }
    data = data->next;
    num_isocurves++;
  }
  /* get extra isocurves if we need them */
  printf("num_isocurves: %d, cont_plot->num_a: %d\n");
  while( num_isocurves < cont_plot->num_a ){
    printf("\t looping: %d\n",num_isocurves);
    data = (iso_curve_type *)D_CALLOC(1,sizeof(iso_curve_type));
    if( !data ) fatal("can't get memory for new iso_curve");
    data->p_max = cont_plot->num_b;
    data->points = (point_type *)D_CALLOC(data->p_max,sizeof(point_type));
    if(!data->points)fatal("can't get memory for data->points");
    if( cont_plot->data ) cont_plot->data->prev = data;
    data->next = cont_plot->data;
    cont_plot->data = data;
    num_isocurves++;
  }

  /**************

    hokay, we've got space to store the isocurves now... so
    go head and generate the data

  **************/
  origin = &cont_plot->left_top_corner;
  switch(cont_plot->orientation){
  case ORIENT_Z:
    step.x = (cont_plot->right_bottom_corner.x -
              cont_plot->left_top_corner.x) / (float)(cont_plot->num_a-1);
    step.y = (cont_plot->right_bottom_corner.y -
              cont_plot->left_top_corner.y) / (float)(cont_plot->num_b-1);
    step.z = 0;

    /* check the character, it will be zero if we're not using symmetry */
    use_symm_a = surf->characters[surf->active_MO].planes[X_AX];
    if( use_symm_a ){
#ifdef REALLY_CHATTY_SYMMETRY
      fprintf(stderr,"Oh good, there's an X mirror with character: %d\n",
              use_symm_a);
#endif
      /* find the end of the list */
      back_ptr = cont_plot->data;
      for(i=0;i<cont_plot->num_a-1;i++){
        back_ptr = back_ptr->next;
        if( !back_ptr )FATAL_BUG("back_ptr ran off the end of the list!");
      }
      num_a_to_do = cont_plot->num_a / 2;
    } else{
      num_a_to_do = cont_plot->num_a;
    }
    use_symm_b = surf->characters[surf->active_MO].planes[Y_AX];
    if( use_symm_b ){
#ifdef REALLY_CHATTY_SYMMETRY
      fprintf(stderr,"Oh good, there's a Y mirror with character: %d\n",
              use_symm_b);
#endif
      num_b_to_do = cont_plot->num_b / 2;
    } else{
      num_b_to_do = cont_plot->num_b;
    }


    /*****
      Memo to self:
       don't *ever* do something like this again:
        for(i=0;i<cont_plot->num_a;i++,data=data->next){
       you oughtn't to be working with the linked list
       in that for statement you dingus.
    *****/

    data = cont_plot->data;
    MO_info.loc.z = origin->z;
    for(i=0;i<num_a_to_do;i++){
      MO_info.loc.x = origin->x + (float)i*step.x;
      for(j=0;j<num_b_to_do;j++){
        MO_info.loc.y = origin->y + (float)j*step.y;
        calc_MO_value(surf->active_MO,&MO_info,surf->MO_centers,
                      surf->num_centers,surf->adf_plot);
        data->points[j].x = MO_info.loc.x;
        data->points[j].y = MO_info.loc.y;
        data->points[j].z = MO_info.val;
        if( use_symm_b ){
          symm_j = cont_plot->num_b-(j+1);
          data->points[symm_j].x = MO_info.loc.x;
          data->points[symm_j].y = origin->y + (float)symm_j*step.y;
          data->points[symm_j].z =
            MO_info.val*(float)use_symm_b;
        }
      }
      data->p_count = cont_plot->num_b;
      if( use_symm_a ){
        symm_i = cont_plot->num_a-(i+1);
        for(j=0;j<cont_plot->num_b;j++){
          back_ptr->points[j].x = origin->x + (float)symm_i*step.x;
          back_ptr->points[j].y = data->points[j].y;
          back_ptr->points[j].z = data->points[j].z*(float)use_symm_a;
        }
        back_ptr->p_count = cont_plot->num_b;
        back_ptr = back_ptr->prev;
      }
      data = data->next;
    }
    break;

  case ORIENT_X:
    step.x = 0;
    step.y = (cont_plot->right_bottom_corner.y - cont_plot->left_top_corner.y) /
      (float)(cont_plot->num_a-1);
    step.z = (cont_plot->right_bottom_corner.z - cont_plot->left_top_corner.z) /
      (float)(cont_plot->num_b-1);

    /* check the character, it will be zero if we're not using symmetry */
    use_symm_a = surf->characters[surf->active_MO].planes[Y_AX];
    if( use_symm_a ){
#ifdef REALLY_CHATTY_SYMMETRY
      fprintf(stderr,"Oh good, there's a Y mirror with character: %d\n",
              use_symm_a);
#endif
      /* find the end of the list */
      back_ptr = cont_plot->data;
      for(i=0;i<cont_plot->num_a-1;i++){
        back_ptr = back_ptr->next;
        if( !back_ptr )FATAL_BUG("back_ptr ran off the end of the list!");
      }
      num_a_to_do = cont_plot->num_a / 2;
    } else{
      num_a_to_do = cont_plot->num_a;
    }
    use_symm_b = surf->characters[surf->active_MO].planes[Z_AX];
    if( use_symm_b ){
#ifdef REALLY_CHATTY_SYMMETRY
      fprintf(stderr,"Oh good, there's a Z mirror with character: %d\n",
              use_symm_b);
#endif
      num_b_to_do = cont_plot->num_b / 2;
    } else{
      num_b_to_do = cont_plot->num_b;
    }

    data = cont_plot->data;
    MO_info.loc.x = origin->x;
    for(i=0;i<num_a_to_do;i++){
      MO_info.loc.y = origin->y + (float)i*step.y;
      for(j=0;j<num_b_to_do;j++){
        MO_info.loc.z = origin->z + (float)j*step.z;
        calc_MO_value(surf->active_MO,&MO_info,surf->MO_centers,
                      surf->num_centers,surf->adf_plot);
        data->points[j].x = MO_info.loc.y;
        data->points[j].y = MO_info.loc.z;
        data->points[j].z = MO_info.val;
        if( use_symm_b ){
          symm_j = cont_plot->num_b-(j+1);
          data->points[symm_j].x = MO_info.loc.y;
          data->points[symm_j].y = origin->z + (float)symm_j*step.z;
          data->points[symm_j].z =
            MO_info.val*(float)use_symm_b;
        }
      }
      data->p_count = cont_plot->num_b;
      if( use_symm_a ){
        symm_i = cont_plot->num_a-(i+1);
        for(j=0;j<cont_plot->num_b;j++){
          back_ptr->points[j].x = origin->y + (float)symm_i*step.y;
          back_ptr->points[j].y = data->points[j].y;
          back_ptr->points[j].z = data->points[j].z*(float)use_symm_a;
        }
        back_ptr->p_count = cont_plot->num_b;
        back_ptr = back_ptr->prev;
      }
      data = data->next;
    }
    break;

  case ORIENT_Y:
    step.y = 0;
    step.x = (cont_plot->right_bottom_corner.x - cont_plot->left_top_corner.x) /
      (float)(cont_plot->num_a-1);
    step.z = (cont_plot->right_bottom_corner.z - cont_plot->left_top_corner.z) /
      (float)(cont_plot->num_b-1);

    /* check the character, it will be zero if we're not using symmetry */
    use_symm_a = surf->characters[surf->active_MO].planes[X_AX];
    if( use_symm_a ){
#ifdef REALLY_CHATTY_SYMMETRY
      fprintf(stderr,"Oh good, there's an X mirror with character: %d\n",
              use_symm_a);
#endif
      /* find the end of the list */
      back_ptr = cont_plot->data;
      for(i=0;i<cont_plot->num_a-1;i++){
        back_ptr = back_ptr->next;
        if( !back_ptr )FATAL_BUG("back_ptr ran off the end of the list!");
      }
      num_a_to_do = cont_plot->num_a / 2;
    } else{
      num_a_to_do = cont_plot->num_a;
    }
    use_symm_b = surf->characters[surf->active_MO].planes[Z_AX];
    if( use_symm_b ){
#ifdef REALLY_CHATTY_SYMMETRY
      fprintf(stderr,"Oh good, there's a Z mirror with character: %d\n",
              use_symm_b);
#endif
      num_b_to_do = cont_plot->num_b / 2;
    } else{
      num_b_to_do = cont_plot->num_b;
    }

    data = cont_plot->data;
    MO_info.loc.y = origin->y;
    for(i=0;i<num_a_to_do;i++){
      MO_info.loc.x = origin->x + (float)i*step.x;
      for(j=0;j<num_b_to_do;j++){
        MO_info.loc.z = origin->z + (float)j*step.z;
        calc_MO_value(surf->active_MO,&MO_info,surf->MO_centers,
                      surf->num_centers,surf->adf_plot);
        data->points[j].x = MO_info.loc.x;
        data->points[j].y = MO_info.loc.z;
        data->points[j].z = MO_info.val;
        if( use_symm_b ){
          symm_j = cont_plot->num_b-(j+1);
          data->points[symm_j].x = MO_info.loc.x;
          data->points[symm_j].y = origin->z + (float)symm_j*step.z;
          data->points[symm_j].z =
            MO_info.val*(float)use_symm_b;
        }
      }
      data->p_count = cont_plot->num_b;
      if( use_symm_a ){
        symm_i = cont_plot->num_a-(i+1);
        for(j=0;j<cont_plot->num_b;j++){
          back_ptr->points[j].x = origin->x + (float)symm_i*step.x;
          back_ptr->points[j].y = data->points[j].y;
          back_ptr->points[j].z = data->points[j].z*(float)use_symm_a;
        }
        back_ptr->p_count = cont_plot->num_b;
        back_ptr = back_ptr->prev;
      }
      data = data->next;
    }
    break;
  case ORIENT_ARBITRARY:
    data = cont_plot->data;
    for(i=0;i<cont_plot->num_a;i++){
      for(j=0;j<cont_plot->num_b;j++){
        MO_info.loc.x = origin->x + (float)i*(cont_plot->Bas1.x)+
          (float)j*(cont_plot->Bas2.x);
        MO_info.loc.y = origin->y + (float)i*(cont_plot->Bas1.y)+
          (float)j*(cont_plot->Bas2.y);
        MO_info.loc.z = origin->z + (float)i*(cont_plot->Bas1.z)+
          (float)j*(cont_plot->Bas2.z);


        calc_MO_value(surf->active_MO,&MO_info,surf->MO_centers,
                      surf->num_centers,surf->adf_plot);
        data->points[j].x = (float)i;
        data->points[j].y = (float)j;
        data->points[j].z = MO_info.val;
      }
      data->p_count = cont_plot->num_b;
      data=data->next;
    }
    break;

  default:
    FATAL_BUG("bogus plane orientation in eval_MO_plane");
  }
}




/****************************************************************************
 *
 *                   Procedure construct_MO_contours
 *
 * Arguments:  surf: pointer to MO_surface_type
 *
 *
 * Returns: none
 *
 * Action:  Sets up and runs everything needed to slap a set of contours
 *        on the current MO.
 *
 ****************************************************************************/
void construct_MO_contours(int num_args,char **MO_surf_ptr)
{
  static char contour_string[240]="-.1 -.09 -.08 -.07 -.06 -.05 -.04 -.03 -.02 -.01 .01 .02 .03 .04 .05 .06 .07 .08 .09 .1\n";
  static float width=10.0,height=10.0;
  static int num_levels=DEFAULT_NUM_OF_ZLEVELS;
  static float offset=0.0;
  static int first_call = 1;
  static int stack_num = 10;
  static int num_a_store = 100;
  static int num_b_store = 100;
  static float  stack_step = 0.5;
  static float lower_val = -0.045;
  static float val_step = 0.01;
  static float start_stack;
  char filename[240],*theinline;
  FILE *grid_file;
  iso_curve_type *data;
  char instring[240],*str_ptr;
  point_type *origin;
  MO_surface_type *surf;
  MO_contour_plot_type *cont_plot;
  MO_contours_type *MO_cont;
  gnuplot_contour_type *gnu_cont;
  int real_num_levels;
  int i,j,k;
  int num_read;
  char do_stack;
  float startx,starty,endx,endy,stepx,stepy;
  float last_val;

  surf = (MO_surface_type *)MO_surf_ptr[0];

#ifdef DEBUG
  for(i=0;i<surf->num_centers;i++){
    printf("Center: %d\n",i);
    for(j=0;j<surf->MO_centers[i].num_AOs;j++){
      printf("% 4d % 4d % 4d % 4d % -6.4lf % -6.4lf\n",
             surf->MO_centers[i].AO_list[j].kx,surf->MO_centers[i].AO_list[j].ky,
             surf->MO_centers[i].AO_list[j].kz,surf->MO_centers[i].AO_list[j].kr,
             surf->MO_centers[i].AO_list[j].zeta1,surf->MO_centers[i].AO_list[j].norm_fact);
    }
    printf("----\n");
    for(j=0;j<surf->MO_centers[i].num_AOs;j++){
      printf("% -6.4lf ",surf->MO_centers[i].AO_list[j].coeff[surf->active_MO]);
      if( !((j+1)%10) ) printf("\n");
    }
    printf("\n");
  }
#endif


  /* get memory for the contour plot if we need it */
  if(!surf->MO_contours){
    cont_plot = (MO_contour_plot_type *)D_CALLOC(1,sizeof(MO_contour_plot_type));
    if( !cont_plot ) fatal("Can't get space for cont_plot");
    surf->MO_contours = cont_plot;
  }

  cont_plot = surf->MO_contours;

  /*******

    just read everything in from the user..... mmmmmmm
    quick hack.

  ********/
  if( first_call ){
    cont_plot->orientation = ORIENT_Z;
  }

  /* allow arbitrary planes if three atoms are selected */
  if(num_selected == 3 ){
    readintparm("Orientation (0:X 1:Y 2:Z 3:Other)",&cont_plot->orientation);
    if(cont_plot->orientation < ORIENT_X ||
       cont_plot->orientation > ORIENT_ARBITRARY){
      printf("Bogus value... setting to Z\n");
      cont_plot->orientation = ORIENT_Z;
    }
  } else{
    readintparm("Orientation (0:X 1:Y 2:Z)",&cont_plot->orientation);
    if(cont_plot->orientation < ORIENT_X ||
       cont_plot->orientation > ORIENT_Z ){
      printf("Bogus value... setting to Z\n");
      cont_plot->orientation = ORIENT_Z;
    }
  }
  readfloatparm("Plane Width",&width);
  readfloatparm("Plane Height",&height);
  readfloatparm("Offset",&offset);

  readintparm("num_a",&num_a_store);
  if( num_a_store % 2 ){
    printf("I'm making that an even number, that's happier.\n");
    num_a_store++;
  }
  readintparm("num_b",&num_b_store);
  if( num_b_store % 2 ){
    printf("I'm making that an even number, that's happier.\n");
    num_b_store++;
  }
  cont_plot->num_a = num_a_store;
  cont_plot->num_b = num_b_store;


  readintparm("number of contours",&num_levels);
  cont_plot->num_levels = num_levels;
  if( cont_plot->levels_list ) D_FREE(cont_plot->levels_list);
  cont_plot->levels_list = (double *)D_CALLOC(cont_plot->num_levels,sizeof(double));
  if( !cont_plot->levels_list ) fatal("can't allocate cont_plot->levels_list");

  cont_plot->num_approx_pts = DEFAULT_NUM_APPROX_PTS;
  cont_plot->interp_kind = INTERP_NOTHING;
  cont_plot->order = DEFAULT_BSPLINE_ORDER;

  cont_plot->levels_kind = surf->levels_kind;
  if( cont_plot->levels_kind == LEVELS_INCREMENTAL ){
    display("Look in the xterm");
    if( !cont_plot->levels_list ){
      cont_plot->levels_list = (float *)D_CALLOC(2,sizeof(float));
      if( !(cont_plot->levels_list) ) fatal("Can't allocate levels list");

      /* set some defaults */
      cont_plot->levels_list[0] = 0.01;
      cont_plot->levels_list[1] = 0.01;
    }
    readfloatparm("lower limit on contours",&lower_val);
    cont_plot->levels_list[0] = lower_val;
    readfloatparm("contour increment",&val_step);
    cont_plot->levels_list[1] = val_step;
  } else if (cont_plot->levels_kind == LEVELS_DISCRETE){
    printf("Enter the %d contour values:\n",num_levels);
    if( strlen(contour_string) > 1 ){
      printf("<return> for last set: %s",contour_string);
      fgets(instring,240,stdin);
      if(instring[0] == '\n' ) strcpy(instring,contour_string);
      else strcpy(contour_string,instring);
    } else{
      instring[0] = 0;
      while (instring[0] == '\n' || instring[0] == 0)
        fgets(instring,240,stdin);
      strcpy(contour_string,instring);
    }
    str_ptr = strtok(instring," ");
    num_read = 0;
    sscanf(str_ptr,"%lf",&(cont_plot->levels_list[num_read++]));
    while(num_read < num_levels){
      if(!str_ptr){
        while(!str_ptr){
          fgets(instring,240,stdin);
          if(instring[0] != '\n' && instring[0] != 0)
            str_ptr = strtok(instring," ");
        }
      }else str_ptr = strtok(NULL," ");
      sscanf(str_ptr,"%lf",&(cont_plot->levels_list[num_read]));
      num_read++;
    }
  }
  cont_plot->num_levels = num_levels;

  /* blow out any old contour data */
  while(cont_plot->contours){
    MO_cont = cont_plot->contours;
    if(MO_cont->hidden_points)D_FREE(MO_cont->hidden_points);
    if(MO_cont->inv_slope)D_FREE(MO_cont->inv_slope);
    if(MO_cont->intercept)D_FREE(MO_cont->intercept);
    if(MO_cont->coords)D_FREE(MO_cont->coords);
    MO_cont = cont_plot->contours->next;
    D_FREE(cont_plot->contours);
    cont_plot->contours = MO_cont;
  }
  cont_plot->gnu_contours = 0;
  cont_plot->contours = 0;

  do_stack = 'n';
  readcharparm("do a stack",&do_stack);

  if( do_stack == 'y' || do_stack == 'Y' ){
    do_stack = 1;
    start_stack = offset;
    readintparm("Stack num",&stack_num);
    readfloatparm("Start Stack",&start_stack);
    readfloatparm("Stack Step",&stack_step);
  }
  else{
    do_stack = 0;
    start_stack = offset;
    stack_step = 0.0;
    stack_num = 1;
  }

  cont_plot->num_conts = 0;
  for(i=0;i<stack_num;i++){
    if( i != 0 ){
      cont_plot->levels_kind = LEVELS_DISCRETE;
    }
    offset = start_stack + i*stack_step;
    switch(cont_plot->orientation){
    case ORIENT_X:
      cont_plot->left_top_corner.x = offset;
      cont_plot->left_top_corner.y = -width/2.0;
      cont_plot->left_top_corner.z = height/2.0;
      cont_plot->right_bottom_corner.x = offset;
      cont_plot->right_bottom_corner.y = width/2.0;
      cont_plot->right_bottom_corner.z = -height/2.0;
      break;
    case ORIENT_Y:
      cont_plot->left_top_corner.x = -width/2.0;
      cont_plot->left_top_corner.y = offset;
      cont_plot->left_top_corner.z = height/2.0;
      cont_plot->right_bottom_corner.x = width/2.0;
      cont_plot->right_bottom_corner.y = offset;
      cont_plot->right_bottom_corner.z = -height/2.0;
      break;
    case ORIENT_ARBITRARY:
      if( determine_plane_bounds(cont_plot,surf->molec->atoms,
                                 surf->molec->num_atoms,width,height,
                                 offset) == 0 ){
        break;
      }else{
        fprintf(stderr,"Problems determining arbitrary plane.\n");
        fprintf(stderr,"Z orientation used instead.\n");
        cont_plot->orientation = ORIENT_Z;
        /* we'll just slide gently into the ORIENT_Z case now */
      }
    case ORIENT_Z:
      cont_plot->left_top_corner.x = -width/2.0;
      cont_plot->left_top_corner.y = height/2.0;
      cont_plot->left_top_corner.z = offset;
      cont_plot->right_bottom_corner.x = width/2.0;
      cont_plot->right_bottom_corner.y = -height/2.0;
      cont_plot->right_bottom_corner.z = offset;
      break;

    }

    grid_file = 0;

    /********

      okay... we've got everything set up.  First evaluate the
      MO in the plane, then contour

      ********/
    printf("Evaluating plane\n");
    eval_MO_plane(surf);

    if( dump_grids_on ){
#ifndef USE_READLINE
      printf("\nEnter grid file name: ");
      scanf("%s\n",filename);
#else
      theinline= readline("Enter the grid file name: ");
      add_history(theinline);
      if( theinline ){
        sscanf(theinline,"%s",filename);
        free(theinline);
      } else {
        error("Bad file name");
        filename[0] = 0;
      }
#endif
      if( filename[0] ){
        grid_file = fopen(filename,"w+");
        if( grid_file ){
          fprintf(grid_file,"#MO_DATA\n");
          fprintf(grid_file,"#NUM_CURVES: 1\n");
          fprintf(grid_file,"#LEFT_TOP: %lf %lf %lf\n",cont_plot->left_top_corner.x,
                  cont_plot->left_top_corner.y,cont_plot->left_top_corner.z);
          fprintf(grid_file,"#RIGHT_BOTTOM: %lf %lf %lf\n",cont_plot->right_bottom_corner.x,
                  cont_plot->right_bottom_corner.y,cont_plot->right_bottom_corner.z);
          fprintf(grid_file,"#STEPS: %d %d\n",cont_plot->num_a,cont_plot->num_b);
          fprintf(grid_file,"#NUM_X: %d\n",cont_plot->num_a);
          fprintf(grid_file,"#NUM_Y: %d\n",cont_plot->num_b);
          switch(cont_plot->orientation){
          case ORIENT_X:
            startx = cont_plot->left_top_corner.y;
            endx = cont_plot->right_bottom_corner.y;
            starty = cont_plot->left_top_corner.z;
            endy = cont_plot->right_bottom_corner.z;
            break;
          case ORIENT_Y:
            startx = cont_plot->left_top_corner.x;
            endx = cont_plot->right_bottom_corner.x;
            starty = cont_plot->left_top_corner.z;
            endy = cont_plot->right_bottom_corner.z;
            break;
          case ORIENT_Z:
            startx = cont_plot->left_top_corner.x;
            endx = cont_plot->right_bottom_corner.x;
            starty = cont_plot->left_top_corner.y;
            endy = cont_plot->right_bottom_corner.y;
            break;
          }
          stepx = (endx - startx)/(float)cont_plot->num_a;
          stepy = (endy - starty)/(float)cont_plot->num_b;
          fprintf(grid_file,"#MIN_X: %lf\n",startx);
          fprintf(grid_file,"#MAX_X: %lf\n",endx);
          fprintf(grid_file,"#STEP_X: %lf\n",stepx);
          fprintf(grid_file,"#MIN_Y: %lf\n",starty);
          fprintf(grid_file,"#MAX_Y: %lf\n",endy);
          fprintf(grid_file,"#STEP_Y: %lf\n",stepy);

          fprintf(grid_file,"\n#BEGIN_DATA\n");
          data = cont_plot->data;
          for(j=0;j<cont_plot->num_a;j++){
            for(k=0;k<cont_plot->num_b;k++){
              fprintf(grid_file,"% -12.8lg\n",data->points[k].z);
            }
            fprintf(grid_file,"\n");
            data = data->next;
          }
          fprintf(grid_file,"#that's all folks\n");
          fclose(grid_file);
        }else{
          error("Can't open output file.");
        }
      }
    }



    printf("Contouring\n");
    cont_plot->gnu_contours = contour_data(cont_plot->num_a,0,cont_plot->data,
                                           cont_plot->num_levels,
                                           cont_plot->num_approx_pts,
                                           cont_plot->interp_kind,
                                           cont_plot->order,
                                           cont_plot->levels_kind,
                                           cont_plot->levels_list);

    /**********

      now postprocess the data so that we can actually plot it

      **********/
    printf("Postprocessing\n");
    gnu_cont = cont_plot->gnu_contours;

    origin = &cont_plot->left_top_corner;

    last_val = -10;

    if(!i) real_num_levels = 0;
    while(gnu_cont){
      MO_cont = (MO_contours_type *)D_CALLOC(1,sizeof(MO_contours_type));
      if(!MO_cont)fatal("Can't allocate MO_cont");
      MO_cont->num_pts = gnu_cont->num_pts;
      MO_cont->value = gnu_cont->coords[0].z;
      if(!i){
        if( MO_cont->value != last_val ){
          cont_plot->levels_list[real_num_levels] = MO_cont->value;
          last_val = MO_cont->value;
          real_num_levels++;
        }
      }
      MO_cont->coords = (point_type *)
        D_CALLOC(MO_cont->num_pts,sizeof(point_type));
      if(!MO_cont->coords) fatal("can't allocate MO_cont->coords");

      for(j=0;j<gnu_cont->num_pts;j++){
        switch(cont_plot->orientation){
        case ORIENT_X:
          MO_cont->coords[j].x = cont_plot->left_top_corner.x;
          MO_cont->coords[j].y = gnu_cont->coords[j].x;
          MO_cont->coords[j].z = gnu_cont->coords[j].y;
          break;
        case ORIENT_Y:
          MO_cont->coords[j].x = gnu_cont->coords[j].x;
          MO_cont->coords[j].y = cont_plot->left_top_corner.y;
          MO_cont->coords[j].z = gnu_cont->coords[j].y;
          break;
        case ORIENT_Z:
          MO_cont->coords[j].x = gnu_cont->coords[j].x;
          MO_cont->coords[j].y = gnu_cont->coords[j].y;
          MO_cont->coords[j].z = cont_plot->left_top_corner.z;
          break;
        case ORIENT_ARBITRARY:
          MO_cont->coords[j].x = origin->x +
            gnu_cont->coords[j].x*(cont_plot->Bas1.x) +
            gnu_cont->coords[j].y*(cont_plot->Bas2.x);
          MO_cont->coords[j].y = origin->y +
            gnu_cont->coords[j].x*(cont_plot->Bas1.y) +
            gnu_cont->coords[j].y*(cont_plot->Bas2.y);
          MO_cont->coords[j].z = origin->z +
            gnu_cont->coords[j].x*(cont_plot->Bas1.z) +
            gnu_cont->coords[j].y*(cont_plot->Bas2.z);

          break;
        }
      }
      MO_cont->next = surf->MO_contours->contours;
      surf->MO_contours->contours = MO_cont;

      /* since we don't need it anymore, free up the gnu contour */
      D_FREE(gnu_cont->coords);
      cont_plot->gnu_contours = gnu_cont;
      gnu_cont = gnu_cont->next;
      D_FREE(cont_plot->gnu_contours);

      cont_plot->num_conts++;
    }
    cont_plot->gnu_contours = 0;
  }
  printf("Done!\n");

  /* figure out what the contour levels are */
  cont_plot->num_levels = real_num_levels;
  printf("There are %d contours at values:\n",cont_plot->num_levels);
  for(i=0;i<cont_plot->num_levels;i++)
    printf("%lg ",cont_plot->levels_list[i]);
  printf("\n");

  first_call = 0;
}


/****************************************************************************
 *
 *                   Procedure MO_contour_surf
 *
 * Arguments:  surf: pointer to MO_surface_type
 *
 *
 * Returns: none
 *
 * Action:  Sets up and runs everything needed to slap a set of contours
 *        on the current MO.
 *
 ****************************************************************************/
void MO_contour_surf(int num_args,char **MO_surf_ptr)
{
  static float width=10.0,height=10.0;
  static int num_levels=DEFAULT_NUM_OF_ZLEVELS;
  static float offset=-1.0;
  static int stack_num = 10;
  static float start_stack=-3;
  static float stack_step=0.7;
  MO_surface_type *surf;
  MO_contour_plot_type *cont_plot;
  MO_contours_type *MO_cont;
  iso_curve_type *iso_curve;
  gnuplot_contour_type *gnu_cont;
  int i,j;
  int num_stacks_to_do,use_symm;
  surf = (MO_surface_type *)MO_surf_ptr[0];

  /* get memory for the contour plot if we need it */
  if(!surf->MO_contours){
    cont_plot = (MO_contour_plot_type *)D_CALLOC(1,sizeof(MO_contour_plot_type));
    if( !cont_plot ) fatal("Can't get space for cont_plot");
    surf->MO_contours = cont_plot;
  }

  cont_plot = surf->MO_contours;

  /* if there's already data stored, blow that out too */
  if( cont_plot->data ){
    iso_curve = cont_plot->data;
    while(iso_curve){
      iso_curve = iso_curve->next;
      D_FREE(cont_plot->data->points);
      D_FREE(cont_plot->data);
      cont_plot->data = iso_curve;
    }
    cont_plot->data = 0;
  }


  /* if we already have some contours, blow them out now */
  MO_cont = cont_plot->contours;
  while(MO_cont){
    if(MO_cont->coords) D_FREE(MO_cont->coords);
    if(MO_cont->hidden_points) free(MO_cont->hidden_points);
    if(MO_cont->inv_slope) free(MO_cont->inv_slope);
    if(MO_cont->intercept) free(MO_cont->intercept);
    cont_plot->contours = MO_cont->next;
    D_FREE(MO_cont);
    MO_cont = cont_plot->contours;
  }
  cont_plot->contours = 0;

  num_levels = cont_plot->num_levels = 2;

  /*******

    just read everything in from the user..... mmmmmmm
    quick hack.

   ********/


  readfloatparm("Plane Width",&width);
  readfloatparm("Plane Height",&height);

  if( cont_plot->num_a == 0 )cont_plot->num_a = 100;
  if( cont_plot->num_b == 0 )cont_plot->num_b = 100;
  readintparm("num_a",&cont_plot->num_a);
  if( cont_plot->num_a % 2 ){
    printf("I'm making that an even number, that's happier.\n");
    cont_plot->num_a++;
  }
  readintparm("num_b",&cont_plot->num_b);
  if( cont_plot->num_b % 2 ){
    printf("I'm making that an even number, that's happier.\n");
    cont_plot->num_b++;
  }

  cont_plot->levels_kind = LEVELS_DISCRETE;
  cont_plot->interp_kind = INTERP_NOTHING;
  cont_plot->num_approx_pts = 0;
  if(cont_plot->levels_list)D_FREE(cont_plot->levels_list);
  cont_plot->levels_list = (double *)D_CALLOC(2,sizeof(double));
  if(!cont_plot->levels_list) fatal("Can't get space for cont_plot->levels_list");
  cont_plot->levels_list[0] = -surf->surface_value;
  cont_plot->levels_list[1] = surf->surface_value;

  cont_plot->gnu_contours = 0;

  readintparm("Stack num",&stack_num);
  if( stack_num % 2 ){
    printf("I'm making that an even number, that's happier.\n");
    stack_num++;
  }
  readfloatparm("Start Stack",&start_stack);
  stack_step = -(2*start_stack/(float)(stack_num-1));
  readfloatparm("Stack Step",&stack_step);

  cont_plot->num_conts = 0;

  if( strstr(surf->plane_dirs,"Z") ){
    cont_plot->orientation = ORIENT_Z;
    use_symm = surf->characters[surf->active_MO].planes[Z_AX];
    if( use_symm ){
#ifdef CHATTY_SYMMETRY
      fprintf(stderr,"That Z mirror cuts our work by half!\n");
#endif
      num_stacks_to_do = stack_num / 2;
    } else {
      num_stacks_to_do = stack_num;
    }

    for(i=0;i<num_stacks_to_do;i++){
      offset = start_stack + i*stack_step;
      cont_plot->left_top_corner.x = -width/2.0;
      cont_plot->left_top_corner.y = height/2.0;
      cont_plot->left_top_corner.z = offset;
      cont_plot->right_bottom_corner.x = width/2.0;
      cont_plot->right_bottom_corner.y = -height/2.0;
      cont_plot->right_bottom_corner.z = offset;

      /********

        okay... we've got everything set up.  First evaluate the
        MO in the plane, then contour

        ********/
      printf("Evaluating plane\n");
      eval_MO_plane(surf);

      printf("Contouring\n");
      cont_plot->gnu_contours = contour_data(cont_plot->num_a,0,cont_plot->data,
                                             cont_plot->num_levels,
                                             cont_plot->num_approx_pts,
                                             cont_plot->interp_kind,
                                             cont_plot->order,
                                             cont_plot->levels_kind,
                                             cont_plot->levels_list);

      /**********

        now postprocess the data so that we can actually plot it

        **********/
      printf("Postprocessing\n");
      gnu_cont = cont_plot->gnu_contours;

      while(gnu_cont){
        MO_cont = (MO_contours_type *)D_CALLOC(1,sizeof(MO_contours_type));
        if(!MO_cont)fatal("Can't allocate MO_cont");
        MO_cont->num_pts = gnu_cont->num_pts;
        MO_cont->value = gnu_cont->coords[0].z;
        MO_cont->orientation = cont_plot->orientation;
        MO_cont->coords = (point_type *)D_CALLOC(MO_cont->num_pts,sizeof(point_type));
        if(!MO_cont->coords) fatal("can't allocate MO_cont->coords");

        for(j=0;j<gnu_cont->num_pts;j++){
          MO_cont->coords[j].x = gnu_cont->coords[j].x;
          MO_cont->coords[j].y = gnu_cont->coords[j].y;
          MO_cont->coords[j].z = cont_plot->left_top_corner.z;
        }
        MO_cont->next = surf->MO_contours->contours;
        surf->MO_contours->contours = MO_cont;
        cont_plot->num_conts++;

        if( use_symm ){
          MO_cont = (MO_contours_type *)D_CALLOC(1,sizeof(MO_contours_type));
          if(!MO_cont)fatal("Can't allocate MO_cont");
          MO_cont->num_pts = gnu_cont->num_pts;
          MO_cont->value = gnu_cont->coords[0].z*use_symm;
          MO_cont->orientation = cont_plot->orientation;
          MO_cont->coords = (point_type *)D_CALLOC(MO_cont->num_pts,sizeof(point_type));
          if(!MO_cont->coords) fatal("can't allocate MO_cont->coords");

          for(j=0;j<gnu_cont->num_pts;j++){
            MO_cont->coords[j].x = gnu_cont->coords[j].x;
            MO_cont->coords[j].y = gnu_cont->coords[j].y;
            MO_cont->coords[j].z = start_stack +
              (float)(stack_num - (i+1))*stack_step;
          }
          MO_cont->next = surf->MO_contours->contours;
          surf->MO_contours->contours = MO_cont;
          cont_plot->num_conts++;
        }
        /* since we don't need it anymore, free up the gnu contour */
        D_FREE(gnu_cont->coords);
        cont_plot->gnu_contours = gnu_cont;
        gnu_cont = gnu_cont->next;
        D_FREE(cont_plot->gnu_contours);
      }
      cont_plot->gnu_contours = 0;
    }
  }
  if( strstr(surf->plane_dirs,"Y") ){
    cont_plot->orientation = ORIENT_Y;
    use_symm = surf->characters[surf->active_MO].planes[Y_AX];
    if( use_symm ){
#ifdef CHATTY_SYMMETRY
      fprintf(stderr,"That Y mirror cuts our work by half!\n");
#endif
      num_stacks_to_do = stack_num / 2;
    } else {
      num_stacks_to_do = stack_num;
    }
    for(i=0;i<num_stacks_to_do;i++){
      offset = start_stack + i*stack_step;
      cont_plot->left_top_corner.x = -width/2.0;
      cont_plot->left_top_corner.y = offset;
      cont_plot->left_top_corner.z = height/2.0;
      cont_plot->right_bottom_corner.x = width/2.0;
      cont_plot->right_bottom_corner.y = offset;
      cont_plot->right_bottom_corner.z = -height/2.0;

      printf("Evaluating plane\n");
      eval_MO_plane(surf);

      printf("Contouring\n");
      cont_plot->gnu_contours = contour_data(cont_plot->num_a,0,cont_plot->data,
                                             cont_plot->num_levels,
                                             cont_plot->num_approx_pts,
                                             cont_plot->interp_kind,
                                             cont_plot->order,
                                             cont_plot->levels_kind,
                                             cont_plot->levels_list);
      printf("Postprocessing\n");
      gnu_cont = cont_plot->gnu_contours;

      while(gnu_cont){
        MO_cont = (MO_contours_type *)D_CALLOC(1,sizeof(MO_contours_type));
        if(!MO_cont)fatal("Can't allocate MO_cont");
        MO_cont->num_pts = gnu_cont->num_pts;
        MO_cont->value = gnu_cont->coords[0].z;
        MO_cont->orientation = cont_plot->orientation;
        MO_cont->coords = (point_type *)D_CALLOC(MO_cont->num_pts,sizeof(point_type));
        if(!MO_cont->coords) fatal("can't allocate MO_cont->coords");

        for(j=0;j<gnu_cont->num_pts;j++){
          MO_cont->coords[j].x = gnu_cont->coords[j].x;
          MO_cont->coords[j].y = cont_plot->left_top_corner.y;
          MO_cont->coords[j].z = gnu_cont->coords[j].y;
        }
        MO_cont->next = surf->MO_contours->contours;
        surf->MO_contours->contours = MO_cont;
        cont_plot->num_conts++;

        if( use_symm ){
          MO_cont = (MO_contours_type *)D_CALLOC(1,sizeof(MO_contours_type));
          if(!MO_cont)fatal("Can't allocate MO_cont");
          MO_cont->num_pts = gnu_cont->num_pts;
          MO_cont->value = gnu_cont->coords[0].z*use_symm;
          MO_cont->orientation = cont_plot->orientation;
          MO_cont->coords = (point_type *)D_CALLOC(MO_cont->num_pts,sizeof(point_type));
          if(!MO_cont->coords) fatal("can't allocate MO_cont->coords");

          for(j=0;j<gnu_cont->num_pts;j++){
            MO_cont->coords[j].x = gnu_cont->coords[j].x;
            MO_cont->coords[j].y = start_stack +
              (float)(stack_num - (i+1))*stack_step;
            MO_cont->coords[j].z = gnu_cont->coords[j].y;
          }
          MO_cont->next = surf->MO_contours->contours;
          surf->MO_contours->contours = MO_cont;
          cont_plot->num_conts++;
        }

        D_FREE(gnu_cont->coords);
        cont_plot->gnu_contours = gnu_cont;
        gnu_cont = gnu_cont->next;
        D_FREE(cont_plot->gnu_contours);
      }
      cont_plot->gnu_contours = 0;
    }
  }

  if( strstr(surf->plane_dirs,"X") ){
    cont_plot->orientation = ORIENT_X;
    use_symm = surf->characters[surf->active_MO].planes[X_AX];
    if( use_symm ){
#ifdef CHATTY_SYMMETRY
      fprintf(stderr,"That X mirror cuts our work by half!\n");
#endif
      num_stacks_to_do = stack_num / 2;
    } else {
      num_stacks_to_do = stack_num;
    }
    for(i=0;i<num_stacks_to_do;i++){
      offset = start_stack + i*stack_step;
      cont_plot->left_top_corner.x = offset;
      cont_plot->left_top_corner.y = -width/2.0;
      cont_plot->left_top_corner.z = height/2.0;
      cont_plot->right_bottom_corner.x = offset;
      cont_plot->right_bottom_corner.y = width/2.0;
      cont_plot->right_bottom_corner.z = -height/2.0;

      printf("Evaluating plane\n");
      eval_MO_plane(surf);

      printf("Contouring\n");
      cont_plot->gnu_contours = contour_data(cont_plot->num_a,0,cont_plot->data,
                                             cont_plot->num_levels,
                                             cont_plot->num_approx_pts,
                                             cont_plot->interp_kind,
                                             cont_plot->order,
                                             cont_plot->levels_kind,
                                             cont_plot->levels_list);
      printf("Postprocessing\n");
      gnu_cont = cont_plot->gnu_contours;

      while(gnu_cont){
        MO_cont = (MO_contours_type *)D_CALLOC(1,sizeof(MO_contours_type));
        if(!MO_cont)fatal("Can't allocate MO_cont");
        MO_cont->num_pts = gnu_cont->num_pts;
        MO_cont->value = gnu_cont->coords[0].z;
        MO_cont->orientation = cont_plot->orientation;
        MO_cont->coords = (point_type *)D_CALLOC(MO_cont->num_pts,sizeof(point_type));
        if(!MO_cont->coords) fatal("can't allocate MO_cont->coords");

        for(j=0;j<gnu_cont->num_pts;j++){
          MO_cont->coords[j].x = cont_plot->left_top_corner.x;
          MO_cont->coords[j].y = gnu_cont->coords[j].x;
          MO_cont->coords[j].z = gnu_cont->coords[j].y;
        }
        MO_cont->next = surf->MO_contours->contours;
        surf->MO_contours->contours = MO_cont;
        cont_plot->num_conts++;
        if( use_symm ){
          MO_cont = (MO_contours_type *)D_CALLOC(1,sizeof(MO_contours_type));
          if(!MO_cont)fatal("Can't allocate MO_cont");
          MO_cont->num_pts = gnu_cont->num_pts;
          MO_cont->value = gnu_cont->coords[0].z*use_symm;
          MO_cont->orientation = cont_plot->orientation;
          MO_cont->coords = (point_type *)D_CALLOC(MO_cont->num_pts,sizeof(point_type));
          if(!MO_cont->coords) fatal("can't allocate MO_cont->coords");

          for(j=0;j<gnu_cont->num_pts;j++){
            MO_cont->coords[j].x = start_stack +
              (float)(stack_num - (i+1))*stack_step;
            MO_cont->coords[j].y = gnu_cont->coords[j].x;
            MO_cont->coords[j].z = gnu_cont->coords[j].y;
          }
          MO_cont->next = surf->MO_contours->contours;
          surf->MO_contours->contours = MO_cont;
          cont_plot->num_conts++;
        }

        /* since we don't need it anymore, free up the gnu contour */
        D_FREE(gnu_cont->coords);
        cont_plot->gnu_contours = gnu_cont;
        gnu_cont = gnu_cont->next;
        D_FREE(cont_plot->gnu_contours);


      }
      cont_plot->gnu_contours = 0;
    }
  }


  printf("Done!\n");

  /* figure out what the contour levels are */
  cont_plot->num_levels = num_levels;
  printf("There are %d contours at values:\n",cont_plot->num_levels);
  for(i=0;i<cont_plot->num_levels;i++)
    printf("%lg ",cont_plot->levels_list[i]);
  printf("\n");
}


/****************************************************************************
 *
 *                   Procedure draw_contours
 *
 * Arguments: surf: pointer to MO_surface_type
 *        contours: pointer to MO_contours_type
 *
 * Returns: none
 *
 * Action:
 *   Draws in the contours in 'surf.
 *
 *   If hidden line removal is being used, uses Jorgensen's algorithm
 *      (from psi88) to do the removal.
 *
 *   NOTE:  The contours should already have been transformed into
 *     the screen frame by now....
 *
 *
 *
 ****************************************************************************/
void draw_contours(MO_surface_type  *surf,MO_contours_type *contours)
{
  static num_Xpoints_allocated=0;
  static XPoint *Xpoints=0;
  static point_type2D *points_to_draw=0;
  int i,j,k,l;
  int num_conts;
  MO_contours_type *contour,*contour2;
  float z1,z2,plane_z,denom;
  point_type *loc;
  char bool1,bool2;
  char hidden,even;
  float xp;
  int points_drawn;
  int points_culled;
#ifdef DEBUG_HIDDEN_LINE
  static XPoint *Xhidpoints=0;
  int hidpoints_drawn;
#endif
  num_conts = surf->MO_contours->num_conts;

  /* this is the non-hidden line removal (quick) way to do it */
  if( !surf->hidden_surf ){

#ifdef CONT_DEBUG
    printf("Points per contour (non-hidden):\n");
    for(i=0;i<num_conts;i++){
      printf("(%d: %d) ",i,contours[i].num_pts);
    }
    printf("\n");
#endif

    for(i=0;i<num_conts;i++){
      contour = &contours[i];
      if( contour->num_pts > num_Xpoints_allocated ){
        if(Xpoints) free(Xpoints);
        if(points_to_draw) free(points_to_draw);
        Xpoints = (XPoint *)calloc(contour->num_pts,sizeof(XPoint));
        points_to_draw = (point_type2D *)calloc(contour->num_pts,
                                                sizeof(point_type2D));
        if(!Xpoints || !points_to_draw) fatal("can't get memory for Xpoints\n");
        num_Xpoints_allocated = contour->num_pts;
#ifdef DEBUG_HIDDEN_LINE
        if(Xhidpoints) free(Xhidpoints);
        Xhidpoints = (XPoint *)calloc(contour->num_pts,sizeof(XPoint));
        if(!Xhidpoints) fatal("can't get memory for Xpoints\n");
#endif

      }
      for(j=0;j<contour->num_pts;j++){
        Xpoints[j].x = points_to_draw[j].x = contour->coords[j].x;
        Xpoints[j].y = points_to_draw[j].y = contour->coords[j].y;
      }
      g_change_linewidth(1);
      if( contour->value > 0 ){
        if( !surf->invert_phase )
          g_change_linestyle(0);
        else
          g_change_linestyle(2);
      }else if( contour->value < 0 ){
        if( !surf->invert_phase )
          g_change_linestyle(2);
        else
          g_change_linestyle(0);
      } else{
        g_change_linestyle(1);
      }
      g_lines(Xpoints,points_to_draw,contour->num_pts,0);
    }
  } else{
    /* okay, we're doing the REAL-STUFF (TM) */

    /*******

      start out by looping over the contours and filling in the
      equations for their planes. and their max and min vals.

      The planes are defined as:
      z = -( A_C X + B_C Y + D_C )  where
      A_C = A/C, B_C = B/C, D_C = D/C

      ********/
#ifdef CONT_DEBUG
    printf("Points per contour (hidden):\n");
    for(i=0;i<num_conts;i++){
      printf("(%d: %d) ",i,contours[i].num_pts);
    }
    printf("\n");
#endif
    for(i=0;i<num_conts;i++){

      contour = &contours[i];
      find_contour_plane(contour);

      contour->min_vals.x = contour->min_vals.y = contour->min_vals.z = 1000.0;
      contour->max_vals.x = contour->max_vals.y = contour->max_vals.z = -1000.0;

      if( !contour->inv_slope ){
#ifdef CONT_DEBUG
        printf("allocating inv_slope\n");
#endif
        contour->inv_slope = (float *)calloc(contour->num_pts,sizeof(float));
        contour->intercept = (float *)calloc(contour->num_pts,sizeof(float));
        if(!contour->intercept) fatal("can't get memory for contour->intercept");
      }


      for(j=0;j<contour->num_pts;j++){
        contour->max_vals.x = MY_MAX(contour->max_vals.x,contour->coords[j].x);
        contour->max_vals.y = MY_MAX(contour->max_vals.y,contour->coords[j].y);
        contour->max_vals.z = MY_MAX(contour->max_vals.z,contour->coords[j].z);
        contour->min_vals.x = MY_MIN(contour->min_vals.x,contour->coords[j].x);
        contour->min_vals.y = MY_MIN(contour->min_vals.y,contour->coords[j].y);
        contour->min_vals.z = MY_MIN(contour->min_vals.z,contour->coords[j].z);

        /* while we're at it, set the intercepts and inverse slopes too (used later) */
        if( j != 0 ){
          denom = contour->coords[j].y - contour->coords[j-1].y;
          if( fabs(denom) > 1e-4 ){
            contour->inv_slope[j] = (contour->coords[j].x - contour->coords[j-1].x)/denom;
          }else{
            contour->inv_slope[j] = 1e4;
          }
          if(fabs(contour->inv_slope[j]) > 1e-4){
            contour->intercept[j] = contour->coords[j].y -
              contour->coords[j].x/contour->inv_slope[j];
          }else{
            /*            contour->intercept[j] = -1e4;*/
            contour->intercept[j] = contour->coords[j].y;
          }
        }
      }
      /*  set the intercept and inverse slope of first point */
      j = contour->num_pts-1;
      denom = contour->coords[0].y - contour->coords[j].y;
      if( fabs(denom) > 1e-4 ){
        contour->inv_slope[0] = (contour->coords[0].x - contour->coords[j].x)/denom;
      }else{
        contour->inv_slope[0] = 1e4;
      }
      if(fabs(contour->inv_slope[0]) > 1e-4){
        contour->intercept[0] = contour->coords[0].y -
          contour->coords[0].x/contour->inv_slope[0];
      }else{
        /*            contour->intercept[0] = -1e4;*/
        contour->intercept[0] = contour->coords[0].y;
      }

    }

    /***************

      okay... everything is set to do the hidden line removal stuff
      now.  Here's how it works...
      -loop over each contour i
      -loop over each contour j  (i != j)
      -if contour j is parallel to i and has a higher z value
      at (1,1), or if the min and max of j are outside those of i
      we know that it can't be hidden by i.
      -else  contour j may have parts hidden by i.
      swing through the points in j and check to
      see which are hidden (using the equation of i's plane).
      mark any points which are hidden.

      *****************/
    points_culled = 0;
    for(i=0;i<num_conts;i++){
/*      printf("culling contour %d\n",i);*/
      contour = &contours[i];
      if( !contour->hidden_points ){
#ifdef CONT_DEBUG
        printf("allocating hidden_points\n");
#endif
        contour->hidden_points = (char *)calloc(contour->num_pts,sizeof(char));
        if(!contour->hidden_points) fatal("Can't get memory for contour->hidden_points");
      } else{
        bzero(contour->hidden_points,contour->num_pts*sizeof(char));
      }

      /* figure out the z value at (1,1) */
      z1 = (contour->A_C + contour->B_C + contour->D_C);

      for(j=0;j<num_conts;j++){
        if( j != i ){
          contour2 = &contours[j];
          /* figure out the z value at (1,1) */
          z2 = (contour2->A_C + contour2->B_C + contour2->D_C);

          /*******

            check the z, max and min values

            the max and min check is for totally disjoint planes.

            *******/
          if( (z2 < z1 || contour->orientation!=contour2->orientation) &&
             !( contour2->min_vals.x-contour->max_vals.x > CLOSE_TO_ZERO ||
               contour2->min_vals.y-contour->max_vals.y > CLOSE_TO_ZERO ) &&
             !(contour2->max_vals.x - contour->min_vals.x < -CLOSE_TO_ZERO||
               contour2->max_vals.y - contour->min_vals.y < -CLOSE_TO_ZERO)){

            /* okay... if we got this far, it means that there's a chance of overlap */
            for(k=0;k<contour->num_pts;k++){
              loc = &contour->coords[k];

              /******

                before we evaluate the equation of the plane, or, god forbid,
                have to go in and do the nasty (compare point-wise),
                check max and min values on this guy

                *******/
              if( !contour->hidden_points[k] &&
                 loc->x <= contour2->max_vals.x && loc->x >= contour2->min_vals.x &&
                 loc->y <= contour2->max_vals.y && loc->y >= contour2->min_vals.y ){

                /********

                  use the equation of the plane to try and avoid culling
                  this point...

                  ********/
                plane_z = (contour2->A_C * loc->x + contour2->B_C * loc->y+
                            contour2->D_C);

                if( plane_z < loc->z ){
                  even = 1;
                  hidden = 0;
                  bool1 = (loc->y - contour2->coords[0].y) < -CLOSE_TO_ZERO;
                  for(l=1;l<contour2->num_pts;l++){
                    bool2 = bool1;
                    bool1 = (loc->y - contour2->coords[l].y) < -CLOSE_TO_ZERO;

                    /* if we're between the y values of the points */
                    if( (bool1 || bool2) && !(bool1 && bool2) ){
                      even = !even;
                      if( fabs(contour2->inv_slope[l]) > 1e-4 ){
                        xp = (loc->y - contour2->intercept[l])*contour2->inv_slope[l];
                        if( loc->x - xp >= -CLOSE_TO_ZERO/1000.0 ){
                          hidden = !hidden;
                        }
                      } else{
                        hidden = !hidden;
                      }
                    }
                  }
                  hidden = hidden && even;
                  if( hidden ){
                    contour->hidden_points[k]=1;
/*printf("Curve %d, point %d clipped by curve %d\n",i,k,j);                    */
                    points_culled++;
                  }
                }
              }
            }
          }
        }
      }
    } /* end of hidden surface spoo */

/*    printf("%d points were culled\n",points_culled);*/
    /* okay.. now we get to draw */
    for(i=0;i<num_conts;i++){
      contour = &contours[i];
      if( contour->num_pts > num_Xpoints_allocated ){
        if(Xpoints) free(Xpoints);
        if(points_to_draw) free(points_to_draw);
        Xpoints = (XPoint *)calloc(contour->num_pts,sizeof(XPoint));
        points_to_draw = (point_type2D *)calloc(contour->num_pts,
                                                sizeof(point_type2D));
        if(!Xpoints || !points_to_draw) fatal("can't get memory for Xpoints\n");
        num_Xpoints_allocated = contour->num_pts;
#ifdef DEBUG_HIDDEN_LINE
        if(Xhidpoints) free(Xhidpoints);
        Xhidpoints = (XPoint *)calloc(contour->num_pts,sizeof(XPoint));
        if(!Xhidpoints) fatal("can't get memory for Xpoints\n");
#endif
      }

      g_change_linewidth(1);
      points_drawn = 0;
#ifdef DEBUG_HIDDEN_LINE
      hidpoints_drawn = 0;
#endif
      for(j=0;j<contour->num_pts;j++){
        if( !contour->hidden_points[j] ){
#ifdef DEBUG_HIDDEN_LINE
          if( hidpoints_drawn ){
            if( contour->value > 0 ){
              if( !surf->invert_phase )
                g_change_linestyle(0);
              else
                g_change_linestyle(2);
            }else if( contour->value < 0 ){
              if( !surf->invert_phase )
                g_change_linestyle(2);
              else
                g_change_linestyle(0);
            } else{
              g_change_linestyle(1);
            }
            g_clines(Xhidpoints,points_to_draw,hidpoints_drawn,0);
            hidpoints_drawn = 0;
          }
#endif
          Xpoints[points_drawn].x = points_to_draw[points_drawn].x =
            contour->coords[j].x;
          Xpoints[points_drawn].y = points_to_draw[points_drawn].y =
            contour->coords[j].y;
#ifdef DEBUG_HIDDEN_LINE
          XFillArc(disp,gpix,graphgc,Xpoints[points_drawn].x-3,Xpoints[points_drawn].y-3,
                   6,6,0,23040);
#endif
          points_drawn++;
        } else{
#ifdef DEBUG_HIDDEN_LINE
          Xhidpoints[hidpoints_drawn].x = contour->coords[j].x;
          Xhidpoints[hidpoints_drawn].y = contour->coords[j].y;
          hidpoints_drawn++;
#endif
          if( points_drawn ){
            if( contour->value > 0 ){
              if( !surf->invert_phase )
                g_change_linestyle(0);
              else
                g_change_linestyle(2);
            }else if( contour->value < 0 ){
              if( !surf->invert_phase )
                g_change_linestyle(2);
              else
                g_change_linestyle(0);
            } else{
              g_change_linestyle(1);
            }
            g_lines(Xpoints,points_to_draw,points_drawn,0);
            points_drawn = 0;
          }
        }
      }
      if( points_drawn ){
        if( contour->value > 0 ){
          if( !surf->invert_phase )
            g_change_linestyle(0);
          else
            g_change_linestyle(2);
        }else if( contour->value < 0 ){
          if( !surf->invert_phase )
            g_change_linestyle(2);
          else
            g_change_linestyle(0);
        } else{
          g_change_linestyle(1);
        }
        g_lines(Xpoints,points_to_draw,points_drawn,0);
      }
#ifdef DEBUG_HIDDEN_LINE
      if( hidpoints_drawn ){
        if( contour->value > 0 ){
          if( !surf->invert_phase )
            g_change_linestyle(0);
          else
            g_change_linestyle(2);
        }else if( contour->value < 0 ){
          if( !surf->invert_phase )
            g_change_linestyle(2);
          else
            g_change_linestyle(0);
        } else{
          g_change_linestyle(1);
        }
        g_clines(Xhidpoints,points_to_draw,hidpoints_drawn,0);
        hidpoints_drawn = 0;
      }
#endif

    }
  }
}

void find_contour_plane(MO_contours_type *contour)
{
  point_type V1,V2,V3;
  float D;
  point_type *a,*b,*c;

  a = &(contour->coords[0]);
  b = &(contour->coords[contour->num_pts/3+1]);
  c = &(contour->coords[2*contour->num_pts/3+1]);

  V1.x = b->x-a->x;V1.y = b->y-a->y;V1.z = b->z-a->z;
  V2.x = c->x-a->x;V2.y = c->y-a->y;V2.z = c->z-a->z;


  V3Cross(&V1,&V2,&V3);
  D = (V3.x*a->x + V3.y*a->y + V3.z*a->z);

  if( fabs(V3.z) < 1e-8 ){
    contour->A_C = 0;
    contour->B_C = 0;
    contour->D_C = 1e6;
  }else{
    contour->A_C = -V3.x/V3.z;
    contour->B_C = -V3.y/V3.z;
    contour->D_C = D/V3.z;
  }

}


