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
   22.02.98 gL:
     General cleanup of commented out code.
     increased verbosity of error message in find_numbered_atom_in_objects
   30.03.98 gL:
     Fixed broken stuff with find_the_line and animations.
   26.04.98 gL:
     added near plane clipping
       this has been needed forever
     tried to fix problems with lines being drawn
      to bogus positions on deep atoms where things are almost
      directly pointed at the camera.  this has been around for
      most of forever.  It seems to be mostly taken care of, though
      there are still some problems that might be tracked down.
     The next thing that really should be done is stopping the white
      parts of breaking lines and tubes from penetrating the atom
      they are drawn from.
      This should be fairly straightforward, but I'm a little bit
      tired to worry about it right now.
   02.05.98 gL:
     Yeah, well, the bogus lines to deep atoms was fixed, but
       it looked bad.  I cleared things up a bit and made
       everything simpler by using the "right" answer with a
       fudge factor.  So far things look okay.
   03.05.98 gL:
     The aforementioned bit about stopping breaking lines from
       penetrating the deep atoms is taken care of.  Again, another
       fudge factor was introduced.
   18.05.98 gL:
     When doing PS code, polyhedra are outlined now.  This should
     be made a toggle somehow.

     triangle normal z coordinates are now divided by the z scaling
      of the object to allow proper shading with scaling
   30.05.98 gL:
     Polyhedra outlines are toggle-able (grossly, using
     a global, but it works).
   04.09.98 gL:
     don't transform atoms which aren't drawn.
   08.09.98 gL:
     support for colored and/or shaded atoms
   24.09.98 gL:
     support for colored and/or shaded atoms updated
   26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)
   29.01.99 gL:
     changes to lattice drawing code, added drawing of cell_bounds.
   30.01.99 gL:
     cell_bounds and lattice drawing cleanup (draw atoms *before* cell/lattice)
     This is still wrong sometimes... atoms end up being drawn before lines when
     they really shouldn't be.  The fix for this puppy is somewhat nontrivial.
   26.06.1999 gL:
     Turning on either draw_axes or draw_box for solids with surfaces would
     cause either seg faults (linux) or calls to FATAL_BUG (AIX).  This was
     the result of a stupid oversight on my part... it has been fixed.
   28.06.1999 gL:
     Fixed EVIL overflow in drawing of hidden lines plots. (Note: it was
     of course just a silly mistake, but it was a cast iron bitch to track
     down.)
***/

/********

  this has got the stuff for dealing with drawing 3D objects.

  created by gL June 1995

*********/
#include "viewkel.h"

point_type *displacements=0;

#define ATOM_RAD 20
#define CLIP_VAL 10
#define LOC_SHOULD_BE_CLIPPED(_p_) (_p_.z < CLIP_VAL ? 1 : 0)
#define V3POINT_ASSIGN(_a_,_b_) {_a_.x=_b_.x;_a_.y=_b_.y;_a_.z=_b_.z;}
line_type *find_the_line(int num1,int num2,molec_type *molec)
{
  int i;
  int num_lines;

  if( !molec ) FATAL_BUG("find_the_line called with null molec");
#if 0
  if( molec->num_frames <= 1 )
    num_lines = molec->num_lines[0];
  else
    num_lines = molec->num_lines[molec->current_frame % molec->num_frames];
#endif
  num_lines = molec->num_lines[0];

  for(i=0;i<num_lines;i++){
    if((molec->lines[i].end1 == num1 || molec->lines[i].end2 == num1) &&
       (molec->lines[i].end1 == num2 || molec->lines[i].end2 == num2))
      return(&molec->lines[i]);
  }

  return(0);
}



int find_numbered_atom_in_objects(generic_3D_object *array,int length,
                                  int number)
{
  int i;

  for(i=0;i<length;i++){
    if( array[i].type == ATOM_3D && array[i].object.atom->num == number )
      return i;
  }
  error("find_number_atom_in_objects can't find an atom for a connector.");
  fprintf(stderr,"offending atom is: %d, length is %d\n",number,length);
}


/****************************************************************************
 *
 *                   Procedure draw_axes
 *
 * Arguments: axes: pointer to axis_type
 *
 * Returns: none
 *
 * Action: Draws in the set of axes at the position indicated in 'axes
 *
 ****************************************************************************/
void draw_axes( axis_type *axes )
{
  float dx,dy;

  g_change_linewidth(2);
  g_change_color(0);

  g_line(axes->points[0].x,axes->points[0].y,axes->points[1].x,axes->points[1].y);
  dx = axes->points[1].x - axes->points[0].x;
  dy = axes->points[1].y - axes->points[0].y;
  g_center_text(axes->points[1].x + .1*dx,axes->points[1].y + .1*dy,"X");

  g_line(axes->points[0].x,axes->points[0].y,axes->points[2].x,axes->points[2].y);
  dx = axes->points[2].x - axes->points[0].x;
  dy = axes->points[2].y - axes->points[0].y;
  g_center_text(axes->points[2].x + .1*dx,axes->points[2].y + .1*dy,"Y");

  g_line(axes->points[0].x,axes->points[0].y,axes->points[3].x,axes->points[3].y);
  dx = axes->points[3].x - axes->points[0].x;
  dy = axes->points[3].y - axes->points[0].y;
  g_center_text(axes->points[3].x + .1*dx,axes->points[3].y + .1*dy,"Z");

}


/****************************************************************************
 *
 *                   Procedure draw_atom
 *
 * Arguments: atom: pointer to atom_type
 *        atom_num: int
 *         objects: pointer to generic_3D_object type
 *     num_objects: int
 *           molec: pointer to molec_type
 *
 * Returns: none
 *
 * Action: Draws in the atom pointed to by 'atom
 *
 ****************************************************************************/
void draw_atom(atom_type *atom,int atom_num,generic_3D_object *objects,
               int num_objects,molec_type *molec)
{
  static char first_call=1, pix_copy_needed=1;
  static float base_radius_scale;
  atom_type *atom2;

  int j;
  char num_string[80];
  point_type origin;
  float radius_scale;
  float xcoord,ycoord,zcoord;
  float dx,dy,dxy,dz,newx,newy,dist;
  float xs,ys;
  line_type *the_line;
  int the_linestyle;

  int radius,radius2,sqr_rad;
  int num_atoms;
  int tab2;
  char shading_on,crosses_on,outlines_on;
  char tube_on,breaking_on,thickness,dashed_on;
  num_atoms = molec->num_atoms;

  /********

    if this is the first call, transform a point at the origin
    to use for radius scaling

  ********/
  if( first_call ){
    origin.x = 0.0;origin.y=0.0;origin.z=0.0;
    transform(&origin);
    base_radius_scale = ATOM_RAD*origin.z;
    first_call = 0;
  }
  radius_scale = molec->rad_mult*base_radius_scale;

  if(near_plane_clipping_on && LOC_SHOULD_BE_CLIPPED(atom->loc)){
    return;
  }
  if( !atom->exclude && (molec->hydrogens_on ||
       atom->type[0] != 'H' || atom->type[1] != 0) &&
     (molec->dummies_on || atom->type[0] != '&') ){

    zcoord = atom->loc.z;
    radius = (int)ceil(atom->rad*radius_scale/zcoord);


    /********
      find the projected coordinates and translate them
      so that the circles appear in the right place.

      (X takes the upper left corner of the rectangle around an
       arc as the argument)
    **********/
    xcoord =  atom->loc.x - radius;
    ycoord =  atom->loc.y - radius;


    if( atom->custom ){
      shading_on = atom->shading_on;
      crosses_on = atom->crosses_on;
      outlines_on = atom->outlines_on;
    } else {
      shading_on = molec->shading_on;
      crosses_on = molec->crosses_on;
      outlines_on = molec->outlines_on;
    }


    /* choose a color */
    g_change_color(atom->color);

    /* draw a filled circle */
    if( shading_on ){
      g_filled_circle(atom->loc.x,atom->loc.y,(float)radius,
                      atom->atom_shade,atom->atom_color,
                      atom->Gpixel_val,atom->Cpixel_val);
#ifdef X_GRAPHICS
#ifdef SUPPORT_COLOR_X
      if( refresh_all_colormaps || pix_copy_needed ){
        memcpy(molec->atoms[atom->num].Gpixel_val,atom->Gpixel_val,
               NUM_X_SHADES*sizeof(long int));
        memcpy(molec->atoms[atom->num].Cpixel_val,atom->Cpixel_val,
               NUM_X_SHADES*sizeof(long int));
        pix_copy_needed = 1;
      }
#endif
#endif
    }else if( crosses_on ){
      g_white_circle(atom->loc.x,atom->loc.y,(float)radius);
    }

    /* draw in the outer circle */
    if( outlines_on ){
      g_change_linewidth(3);
      g_change_color(0);
      if( atom->is_selected && mainmode == CHOOSE){
        g_change_linestyle(2);
      }
      if( crosses_on ){
        g_crossed_circle(atom->loc.x,atom->loc.y,(float)radius);
      } else{
        g_open_circle(atom->loc.x,atom->loc.y,(float)radius);
      }
      g_change_linewidth(1);
      g_change_linestyle(0);
    }

    /* update the atom's screen location (for selecting) */
    molec->atoms[atom->num].screen_loc.x = atom->loc.x;
    molec->atoms[atom->num].screen_loc.y = atom->loc.y;
    molec->atoms[atom->num].screen_rad = radius;

    /* draw in a line (if we need to) */
    if( (molec->draw_connectors && atom->num_lines_out) ){

      /* find the end point */
      for(j=0;j<atom->num_lines_out;j++){
        tab2 = find_numbered_atom_in_objects(objects,num_objects,
                                             atom->linesto[j]);

        /* if the other end point is closer to the camera, then draw it */
        if(tab2 < atom_num && !objects[tab2].object.atom->exclude &&
           ((molec->hydrogens_on ||
             objects[tab2].object.atom->type[0] != 'H' ||
             objects[tab2].object.atom->type[1] != 0) &&
            (molec->dummies_on ||
             objects[tab2].object.atom->type[0] != '&')) ){

          atom2 = objects[tab2].object.atom;
          radius2 = (int)ceil(atom2->rad*radius_scale/atom2->loc.z);
          sqr_rad = radius+radius2;
          sqr_rad *= sqr_rad;

          if( !near_plane_clipping_on ||
              !LOC_SHOULD_BE_CLIPPED(atom2->loc)){
            the_linestyle =0;

          /****

            okay, if we need to determine which line this is
             to determine how to draw it... do so now

          ****/
            the_line = find_the_line(atom->num,objects[tab2].object.atom->num,
                                     molec);
            if(the_line){
              if( the_line->custom ){
                tube_on = the_line->tube;
                breaking_on = the_line->breaking;
                thickness = the_line->thickness;
                dashed_on = the_line->dashed;
              } else{
                tube_on = molec->tubes_on;
                breaking_on = molec->breaking_lines;
                thickness = molec->line_width;
                dashed_on = 0;
#ifdef INCLUDE_BOND_VALENCE
                if( molec->valence_for_bonds ) dashed_on = 1;
#endif
              }

              if(dashed_on){
                the_linestyle = the_line->type;
              }

              if( !the_line->custom || the_line->drawn ){
                if( molec->fancy_lines ){

                  /****
                    figure out how to break the line at the edge
                    of the bottom atom
                    *****/

#if 0
                  dx = atom2->loc.x-atom->loc.x;
                  dy = atom2->loc.y-atom->loc.y;
                  dxy = dx*dx+dy*dy;
                  if( dxy > 1e-2 ){
                    dz = objects[tab2].object.atom->loc.z-atom->loc.z;
                    dtot = dxy + dz*dz;
                    if( fabs(dz) < .1 ){
                      dist = 1.05*radius;
                    }else{
                      dist = .8*radius*dxy/dtot;
                    }
                    if( dx != 0.0 ){
                      slope = dy / dx;
                      newx = dist/sqrt(1.0+slope*slope);
                      if(dx < 0.0) newx *= -1.0;
                      newy = slope*newx;
                    }
                    else{
                      newx = 0.0;
                      newy = dist;
                      if( dy < 0.0 ) newy *= -1;
                    }
                  }
#else
                  dx = atom2->loc.x - atom->loc.x;
                  dy = atom2->loc.y - atom->loc.y;
                  dxy = dx*dx+dy*dy;
                  /****

                    here we have a fudge factor to take into
                    account the fact that the z values we are working with
                    have not been altered by the homogenous transformation.
                    this is *wrong*, but it looks okay most of the time, so
                    we'll keep it.

                  ****/
                  dz = fabs(atom->loc.z)*0.5*(atom2->loc.z - atom->loc.z);
                  dist = sqrt(dx*dx+dy*dy+dz*dz);
                  dxy = sqrt(dxy);
                  newx = radius * dx / dist;
                  newy = radius * dy / dist;
                  /* calculate where the edge of the atom is */
                  if( dxy < radius + radius2 ){
                    xs = dx;
                    ys = dy;
                  }else{
                    xs = 1.02*radius * dx / dxy;
                    ys = 1.02*radius * dy / dxy;
                  }

#endif
                }
                else{
                  newx = 0.0;
                  newy = 0.0;
                  dxy = 1.0;
                }
                /******

                  draw if there is sufficient xy change.
                  the test for this is fairly simple, there's
                  a lower limit and then a test to see if the front
                  atom is inside the circle of the back one.
                  in either of these cases, don't bother drawing.

                ****/
#if 0
                if( dxy > 1e-2 && dxy > sqr_rad ){
#else
                if( dxy > 1e-2 ){
#endif
                  if( breaking_on ){
                    if( !the_linestyle )
                      g_draw_stop_line(atom->loc.x+newx,
                                       atom->loc.y+newy,
                                       atom2->loc.x,
                                       atom2->loc.y,
                                       thickness,
                                       atom->loc.x+xs,
                                       atom->loc.y+ys);
                    else
                      g_draw_dashed_stop_line(atom->loc.x+newx,
                                              atom->loc.y+newy,
                                              atom2->loc.x,
                                              atom2->loc.y,
                                              thickness,the_linestyle,
                                              atom->loc.x+xs,
                                              atom->loc.y+ys);

                  } else if( tube_on ){
                    if(  !the_linestyle )
                      g_draw_stop_tube(atom->loc.x+newx,
                                       atom->loc.y+newy,
                                       atom2->loc.x,
                                       atom2->loc.y,
                                       thickness,
                                       atom->loc.x+xs,
                                       atom->loc.y+ys);

                    else{
                      if( !the_line->custom )
                        g_draw_dashed_stop_line(atom->loc.x+newx,
                                                atom->loc.y+newy,
                                                atom2->loc.x,
                                                atom2->loc.y,
                                                thickness,the_linestyle,
                                                atom->loc.x+xs,
                                                atom->loc.y+ys);
                      else
                        g_draw_dashed_stop_tube(atom->loc.x+newx,
                                                atom->loc.y+newy,
                                                atom2->loc.x,
                                                atom2->loc.y,
                                                thickness,the_linestyle,
                                                atom->loc.x+xs,
                                                atom->loc.y+ys);
                    }
                  }else{
                    if( !the_linestyle )
                      g_line(atom->loc.x+newx,
                             atom->loc.y+newy,
                             atom2->loc.x,
                             atom2->loc.y);
                  }
                }
              }
            }
          }
        }
      }
    }
    g_change_color(0);
    /* draw in a number for the atom if it is needed */
    if(molec->numbers_on){
      sprintf(num_string,"%d",atom->num+1);
      g_right_text(xcoord,ycoord,num_string);
    }
    if(molec->symbols_on){
      sprintf(num_string,"%s",atom->type);
      g_center_text(xcoord+radius,ycoord+.75*radius,num_string);
    }

#ifdef INCLUDE_ADF_PLOTS
    /* draw in the vibration if that is needed */
    if( molec->num_vibrations ){
      g_change_linewidth(2);
      g_line(atom->loc.x,atom->loc.y,displacements[atom->num].x,displacements[atom->num].y);
    }
#endif
  }
}

/****************************************************************************
 *
 *                   Procedure draw_triangle
 *
 * Arguments: triangle: pointer to triangle_type
 *            vertices: pointer to vertex_type
 *                surf: pointer to MO_surf_type
 * Returns: none
 *
 * Action: Draws in the triangle pointed to by 'triangle
 *
 ****************************************************************************/
void draw_triangle(triangle_type *triangle,vertex_type *vertices,
                   MO_surface_type *surf)
{
  static XPoint xpoints[4];

  vertex_type *v1,*v2,*v3;
  int num_vis_vertices;

  v1 = &vertices[triangle->vertices[0]];
  v2 = &vertices[triangle->vertices[1]];
  v3 = &vertices[triangle->vertices[2]];
  num_vis_vertices = 2;

  if( num_vis_vertices > 0 ){
    /* set up the vertices */
    xpoints[0].x = (int)v1->position.x;
    xpoints[0].y = (int)v1->position.y;
    xpoints[1].x = (int)v2->position.x;
    xpoints[1].y = (int)v2->position.y;
    xpoints[2].x = (int)v3->position.x;
    xpoints[2].y = (int)v3->position.y;
    xpoints[3].x = xpoints[0].x;
    xpoints[3].y = xpoints[0].y;

    if( surf ){
      if( surf->do_shading ){
        /* draw a filled triangle of the appropriate color */
        if( triangle->color ){
          g_change_color(3);
        }
        else{
          g_change_color(6);
        }
        g_filled_polygon(xpoints,3);
      }

      if( surf->do_lines ){
        g_change_color(0);
        g_open_polygon(xpoints,3);
      }
    } else{
#ifdef CULL_POLYGONS
      if( triangle->normal.z < 0 ){
#endif
        g_shaded_polygon(xpoints,3,&triangle->normal);

        if( outline_polyhed_on ){
          g_change_linewidth(2);
          g_change_color(0);
          g_open_polygon(xpoints,3);
        }
#ifdef CULL_POLYGONS
      }
#endif
    }
  }

}



/* this function is used by qsort to depth sort the objects */
int object_zcompare(const void *obj1p,const void *obj2p)
{
  generic_3D_object *obj1,*obj2;
  float z1,z2,diff;

  obj1 = (generic_3D_object *)obj1p;
  obj2 = (generic_3D_object *)obj2p;


  switch(obj1->type){
  case ATOM_3D: z1 = obj1->object.atom->loc.z;break;
  case TRIANGLE_3D:z1 = obj1->object.triangle->center.z;break;
  case CONT_POINT_3D:z1=obj1->object.contour_point.loc->z;break;
  case LINE_3D:
    z1 = (obj1->object.line.end1->z > obj1->object.line.end2->z ?
          obj1->object.line.end1->z : obj1->object.line.end2->z);
    break;
  default: FATAL_BUG("Invalid type in 3D_object_zcompare");
  }
  switch(obj2->type){
  case ATOM_3D: z2 = obj2->object.atom->loc.z;break;
  case TRIANGLE_3D:z2 = obj2->object.triangle->center.z;break;
  case CONT_POINT_3D:z2=obj2->object.contour_point.loc->z;break;
  case LINE_3D:
    z2 = (obj2->object.line.end1->z > obj2->object.line.end2->z ?
          obj2->object.line.end1->z : obj2->object.line.end2->z);
    break;
  default: FATAL_BUG("Invalid type in 3D_object_zcompare");
  }

  /* atoms are always closer than lines at the same z value */
  diff = z1-z2;
  if( fabs(diff) >= 0.01){
    return (int)(100.0*(z1 - z2));
  } else if( obj1->type == LINE_3D && obj2->type == ATOM_3D ){
    return(1);
  } else if(obj2->type == LINE_3D && obj1->type == ATOM_3D ){
    return(-1);
  }
}

/****************************************************************************
 *
 *                   Procedure draw_3D_objects
 *
 * Arguments: prim: pointer to prim_type
 *             obj: pointer to object_type
 *
 * Returns: none
 *
 * Action: Draws in the 3D objects contained in 'prim.
 *
 *   The algorithm used is:
 *    1) Do the perspective transformation on all the objects
 *    2) Depth sort the atoms by distance from camera
 *    3) Draw from back to front so that the front objects cover the
 *       back ones
 *        (primitive hidden surface removal).
 *
 *    <this is the Painter's algorithm>
 *
 ****************************************************************************/
void draw_3D_objects(prim_type *prim,object_type *obj)
{
  static generic_3D_object *objects=0;
  static int num_objects_allocated=0;
  static atom_type *atom_store=0;
  static triangle_type *triangle_store=0;
  static vertex_type *vertex_store=0;

  static int num_atoms_allocated=0;
  static int num_triangles_allocated=0;
  static int num_vertices_allocated=0;
  static MO_contours_type *contour_store=0;
  static int num_conts_allocated=0;
  static axis_type axes;
  static point_type lattice[8];
  static point_type cell_box[24];
  int i,j;
  int num_objects,objects_so_far,frame;
  int num_lattice_p,num_box_p;
  molec_type *molec;
  MO_surface_type *surf;
  int num_atoms,num_triangles,num_vertices;
  atom_type *atoms;
  triangle_type *triangles;
  float xcoord,ycoord;
  MO_contours_type *contour;
  int num_conts, num_cont_p;
  char do_axes=0;

  /* initialize this object's bounding box */
  obj->bmin.x = obj->bmin.y = 10000;
  obj->bmax.x = obj->bmax.y = 0;


  /**********

    figure out the total number of objects so we can get the correct
    amount of memory.

  ***********/
  num_objects = 0;
  num_atoms = 0;
  num_triangles = 0;
  switch(prim->which){
  case MOLECULE:
    molec = prim->molec;
    num_atoms = molec->num_atoms;
    num_objects += num_atoms;
    if( molec->num_frames <= 1){
      atoms = molec->atoms;
    }else{
      frame = molec->current_frame % molec->num_frames;
      atoms = &(molec->atoms[frame*molec->num_atoms]);
    }
    num_vertices = 0;
    num_conts = 0;
    if( molec->draw_polyhed ){
      num_triangles = molec->num_triangles;
      num_objects += num_triangles;
      if( num_triangles ) num_vertices = molec->num_polyhed_verts;
    } else{
      num_triangles = 0;
      num_vertices = 0;
    }

    if( molec->axes_on ) do_axes = 1;

    /* do we need memory? */
#ifdef INCLUDE_ADF_PLOTS
    if( num_atoms > num_atoms_allocated || !atom_store ||
        (molec->num_vibrations && !displacements) ){
#else
    if( num_atoms > num_atoms_allocated || !atom_store ){
#endif
      if(atom_store) free(atom_store);
      atom_store = (atom_type *)calloc(num_atoms,sizeof(atom_type));
      if( !atom_store ) fatal("Can't get memory for atom_store.");
#ifdef INCLUDE_ADF_PLOTS
      if( molec->num_vibrations ){
        if( displacements ) free(displacements);
        displacements = (point_type *)calloc(num_atoms,sizeof(point_type));
        if(!displacements) fatal("can't allocated displacements");
      }
#endif

      num_atoms_allocated = num_atoms;
    }

    /* copy the atoms over */
    memcpy((char *)atom_store,(char *)atoms,num_atoms*sizeof(atom_type));

#ifdef INCLUDE_ADF_PLOTS
    /* copy in the displacement array of the vibrations if that is needed */
    if( molec->num_vibrations ){
      for(i=0;i<num_atoms;i++){
        displacements[i].x = molec->atoms[i].loc.x +
          molec->atoms[i].displacements[molec->active_vibn-1].x *
            molec->vibration_scale;
        displacements[i].y = molec->atoms[i].loc.y +
          molec->atoms[i].displacements[molec->active_vibn-1].y *
            molec->vibration_scale;
        displacements[i].z = molec->atoms[i].loc.z +
          molec->atoms[i].displacements[molec->active_vibn-1].z *
            molec->vibration_scale;
      }
    }
#endif
    if( num_triangles > num_triangles_allocated ){
      if(triangle_store) free(triangle_store);
      triangle_store = (triangle_type *)
        calloc(num_triangles,sizeof(triangle_type));
      if( !triangle_store ) fatal("Can't get memory for triangle_store.");
      num_triangles_allocated = num_triangles;
    }

    if( num_vertices > num_vertices_allocated ){
      if(vertex_store) free(vertex_store);
      vertex_store = (vertex_type *)
        calloc(num_vertices,sizeof(vertex_type));
      if( !vertex_store ) fatal("Can't get memory for vertex_store.");
      num_vertices_allocated = num_vertices;
    }

    /* copy the triangles & vertices over */
    if(num_triangles){
      memcpy((void *)triangle_store,(void *)molec->triangles,
             num_triangles*sizeof(triangle_type));
      memcpy((void *)vertex_store,(void *)molec->polyhed_verts,
              num_vertices*sizeof(vertex_type));
    }

    surf = 0;
    break;

  case MO_SURF:
    surf = prim->MO_surf;
    if( surf->display_molec ){
      molec = surf->molec;
      num_atoms = molec->num_atoms;
      num_objects += num_atoms;
      atoms = molec->atoms;

      if( molec->axes_on ) do_axes = 1;

      /* do we need memory? */
      if( num_atoms > num_atoms_allocated ){
        if(atom_store) free(atom_store);
        atom_store = (atom_type *)calloc(num_atoms,sizeof(atom_type));
        if( !atom_store ) fatal("Can't get memory for atom_store.");
        num_atoms_allocated = num_atoms;
      }

      /* copy the atoms over */
      memcpy((char *)atom_store,(char *)atoms,num_atoms*sizeof(atom_type));

    } else{
      molec = 0;
      num_atoms = 0;
    }
    if( surf->display_surf ){
      triangles = surf->triangles;
      num_triangles = surf->num_triangles;
      num_vertices = surf->num_vertices;
      num_objects += num_triangles;

      /* do we need memory? */
      if( num_triangles > num_triangles_allocated ){
        if(triangle_store) free(triangle_store);
        triangle_store = (triangle_type *)
          calloc(num_triangles,sizeof(triangle_type));
        if( !triangle_store ) fatal("Can't get memory for triangle_store.");
        num_triangles_allocated = num_triangles;
      }
      if( num_vertices > num_vertices_allocated ){
        if(vertex_store) free(vertex_store);
        vertex_store = (vertex_type *)
          calloc(num_vertices,sizeof(vertex_type));
        if( !vertex_store ) fatal("Can't get memory for vertex_store.");
        num_vertices_allocated = num_vertices;
      }

      /* copy the triangles & vertices over */
      memcpy((char *)triangle_store,(char *)surf->triangles,
             num_triangles*sizeof(triangle_type));

      memcpy((char *)vertex_store,(char *)surf->triangle_vertices,
             num_vertices*sizeof(vertex_type));

    } else{
      triangles = 0;
      num_triangles = 0;
    }
    if( surf->display_conts && surf->MO_contours){
      num_conts = surf->MO_contours->num_conts;
#ifdef CONT_DEBUG
      fprintf(stderr,"surf from file: %s has %d contours\n",
              surf->filename,num_conts);
#endif
      if( num_conts_allocated < num_conts ){
        if( contour_store ) free(contour_store);
        num_conts_allocated = num_conts;
        contour_store = (MO_contours_type *)calloc(num_conts,
                                                    sizeof(MO_contours_type));
        if(!contour_store) fatal("can't get memory for contour_store");
        for(i=0;i<num_conts-1;i++){
          contour_store[i].next = &contour_store[i+1];
        }
      }
      /* copy the data over now */
      i = 0;
      num_cont_p = 0;
      contour = surf->MO_contours->contours;
      while(contour){
        if( contour_store[i].num_pts != contour->num_pts ){
          if( contour_store[i].coords) free(contour_store[i].coords);
          contour_store[i].coords =
            (point_type *)calloc(contour->num_pts,sizeof(point_type));
          if( !contour_store[i].coords )
            fatal("can't get memory for contour_store[i].coords");
          contour_store[i].num_pts = contour->num_pts;
          if( contour_store[i].inv_slope ){
            free(contour_store[i].inv_slope);
            contour_store[i].inv_slope = 0;
          }
          if( contour_store[i].intercept ){
            free(contour_store[i].intercept);
            contour_store[i].intercept = 0;
          }
          if( contour_store[i].hidden_points ){
            free(contour_store[i].hidden_points);
            contour_store[i].hidden_points = 0;
          }
        }
        contour_store[i].orientation = contour->orientation;
        memcpy((char *)contour_store[i].coords,(char *)contour->coords,
               contour->num_pts*sizeof(point_type));
        contour_store[i].value = contour->value;
        num_cont_p += contour->num_pts;
        contour = contour->next;
        i++;
      }
      if( surf->sort_conts )
        num_objects += num_cont_p;
    } else{
      num_conts = 0;
    }
    break;

  default:
    FATAL_BUG("Invalid primitive passed to draw_3D_objects");
  }


  if( molec && molec->draw_lattice ){
    for(i=0;i<=molec->num_dim;i++){
      V3POINT_ASSIGN(lattice[i],molec->lattice_vect[i]);
    }
    num_lattice_p = molec->num_dim+1;
    num_objects += molec->num_dim;
  }

  if( molec && molec->draw_box ){
    V3POINT_ASSIGN(cell_box[0],molec->cell_box[0]);
    V3POINT_ASSIGN(cell_box[1],molec->cell_box[1]);
    num_box_p = 2;
    if( molec->num_dim > 1 ){
      V3POINT_ASSIGN(cell_box[2],cell_box[1]);
      V3POINT_ASSIGN(cell_box[3],molec->cell_box[2]);
      V3POINT_ASSIGN(cell_box[4],cell_box[3]);
      V3POINT_ASSIGN(cell_box[5],molec->cell_box[3]);
      V3POINT_ASSIGN(cell_box[6],cell_box[5]);
      V3POINT_ASSIGN(cell_box[7],molec->cell_box[0]);
      num_box_p = 8;
      if( molec->num_dim > 2 ){
        V3POINT_ASSIGN(cell_box[8],cell_box[7]);
        V3POINT_ASSIGN(cell_box[9],molec->cell_box[4]);
        V3POINT_ASSIGN(cell_box[10],cell_box[9]);
        V3POINT_ASSIGN(cell_box[11],molec->cell_box[5]);
        V3POINT_ASSIGN(cell_box[12],cell_box[11]);
        V3POINT_ASSIGN(cell_box[13],molec->cell_box[1]);
        V3POINT_ASSIGN(cell_box[14],cell_box[12]);
        V3POINT_ASSIGN(cell_box[15],molec->cell_box[6]);
        V3POINT_ASSIGN(cell_box[16],cell_box[15]);
        V3POINT_ASSIGN(cell_box[17],molec->cell_box[2]);
        V3POINT_ASSIGN(cell_box[18],cell_box[16]);
        V3POINT_ASSIGN(cell_box[19],molec->cell_box[7]);
        V3POINT_ASSIGN(cell_box[20],cell_box[19]);
        V3POINT_ASSIGN(cell_box[21],molec->cell_box[3]);
        V3POINT_ASSIGN(cell_box[22],cell_box[20]);
        V3POINT_ASSIGN(cell_box[23],molec->cell_box[4]);
        num_box_p = 24;
      }
    }
    num_objects += num_box_p/2;
  }
  /******

    okay, check to see if we need to get more memory
    for the 3D objects.

  *******/
  if( num_objects > num_objects_allocated ){
    if( objects ) free(objects);
    objects = (generic_3D_object *)
      calloc(num_objects,sizeof(generic_3D_object));
    num_objects_allocated = num_objects;
    if( !objects ) fatal("Can't get memory for objects.");
  }

  /********

    do we need to re-determine bond locations?

  *********/
#ifdef INCLUDE_BOND_VALENCE
  if(molec && (molec->bond_tol2 != molec->old_bond_tol2 ||
               molec->bond_tol != molec->old_bond_tol)){
    determine_connections(molec);
    molec->old_bond_tol = molec->bond_tol;
    molec->old_bond_tol2 = molec->bond_tol2;
  }
#else
  if(molec && molec->bond_tol != molec->old_bond_tol){
    determine_connections(molec);
    molec->old_bond_tol = molec->bond_tol;
  }
#endif

  /**********

    we're set, transform all the objects and set the
    pointers

  **********/
  objects_so_far = 0;
  if( atoms ){
    for(i=0;i<num_atoms;i++,objects_so_far++){
      if( !atom_store[i].exclude && (molec->hydrogens_on ||
                                     atom_store[i].type[0] != 'H' ||
                                     atom_store[i].type[1] != 0) &&
          (molec->dummies_on || atom_store[i].type[0] != '&') ){
        transform(&(atom_store[i].loc));
        atom_store[i].loc.x += obj->cent.x;
        atom_store[i].loc.y += obj->cent.y;
        xcoord = atom_store[i].loc.x;
        ycoord = atom_store[i].loc.y;

#ifdef INCLUDE_ADF_PLOTS
        /* transform the displacement if that is needed */
        if( molec->num_vibrations ){
          transform(&(displacements[i]));
          displacements[i].x += obj->cent.x;
          displacements[i].y += obj->cent.y;
        }
#endif

        /* update the bounding box */
        if( xcoord < obj->bmin.x ) obj->bmin.x = xcoord;
        if( ycoord < obj->bmin.y ) obj->bmin.y = ycoord;
        if( xcoord > obj->bmax.x ) obj->bmax.x = xcoord;
        if( ycoord > obj->bmax.y ) obj->bmax.y = ycoord;
      }
      objects[objects_so_far].type = ATOM_3D;
      objects[objects_so_far].object.atom =
        &(atom_store[i]);
    }
  }

  if(molec && molec->draw_lattice){
    transform(&lattice[0]);
    lattice[0].x += obj->cent.x;
    lattice[0].y += obj->cent.y;
    for(i=1;i<num_lattice_p;i++){
      transform(&lattice[i]);
      lattice[i].x += obj->cent.x;
      lattice[i].y += obj->cent.y;
      objects[objects_so_far].type = LINE_3D;
      objects[objects_so_far].object.line.end1 = &lattice[0];
      objects[objects_so_far].object.line.end2 = &lattice[i];
      objects_so_far++;
    }
  }

  if(molec && molec->draw_box){
    for(i=0;i<num_box_p;i+=2){
      transform(&cell_box[i]);
      transform(&cell_box[i+1]);
      cell_box[i].x += obj->cent.x;
      cell_box[i].y += obj->cent.y;
      cell_box[i+1].x += obj->cent.x;
      cell_box[i+1].y += obj->cent.y;
      objects[objects_so_far].type = LINE_3D;
      objects[objects_so_far].object.line.end1 = &cell_box[i];
      objects[objects_so_far].object.line.end2 = &cell_box[i+1];
      objects_so_far++;
    }
  }


  if( triangles ){
    for(i=0;i<num_triangles;i++,objects_so_far++){
      triangle_store[i].normal.x += triangle_store[i].center.x;
      triangle_store[i].normal.y += triangle_store[i].center.y;
      triangle_store[i].normal.z += triangle_store[i].center.z;
      transform(&(triangle_store[i].center));
      transform(&(triangle_store[i].normal));
      triangle_store[i].normal.x += obj->cent.x;
      triangle_store[i].normal.y += obj->cent.y;
      triangle_store[i].center.x += obj->cent.x;
      triangle_store[i].center.y += obj->cent.y;

      triangle_store[i].normal.z -= triangle_store[i].center.z;
      triangle_store[i].normal.z /= obj->scale.z;

      objects[objects_so_far].type = TRIANGLE_3D;

      objects[objects_so_far].object.triangle =
        &(triangle_store[i]);
    }

    /* the vertices and normals need to be transformed too */
    for(i=0;i<num_vertices;i++){
      transform(&(vertex_store[i].position));
      vertex_store[i].position.x += obj->cent.x;
      vertex_store[i].position.y += obj->cent.y;
    }
  }


  if( num_conts && surf->display_conts ){
    for(i=0;i<num_conts;i++){
      contour = &contour_store[i];
      for(j=0;j<contour->num_pts;j++){
        transform(&(contour->coords[j]));
        contour->coords[j].x += obj->cent.x;
        contour->coords[j].y += obj->cent.y;
        xcoord = contour->coords[j].x;
        ycoord = contour->coords[j].y;

        /* update the bounding box */
        if( xcoord < obj->bmin.x ) obj->bmin.x = xcoord;
        if( ycoord < obj->bmin.y ) obj->bmin.y = ycoord;
        if( xcoord > obj->bmax.x ) obj->bmax.x = xcoord;
        if( ycoord > obj->bmax.y ) obj->bmax.y = ycoord;


        if( surf->sort_conts ){
          objects[objects_so_far].type = CONT_POINT_3D;
          objects[objects_so_far].object.contour_point.loc =
            &(contour->coords[j]);
          if( j < contour->num_pts-1)
            objects[objects_so_far].object.contour_point.next =
              &(contour->coords[j+1]);
          else
            objects[objects_so_far].object.contour_point.next = 0;

          objects[objects_so_far].object.contour_point.val = contour->value;
          objects_so_far++;
        }
      }
    }
  }



  /* internal consistency check */
  if( objects_so_far != num_objects )
    FATAL_BUG("objects_so_far != num_objects in draw_3D_objects");

  /* depth sort the objects */
  qsort(objects,num_objects,sizeof(generic_3D_object),
        object_zcompare);

  /* okay, everything is ready to draw.... do so */
  for(i=num_objects-1;i>=0;i--){
    switch(objects[i].type){
    case ATOM_3D:
      draw_atom(objects[i].object.atom,i,objects,num_objects,molec);

      if( prim->molec )
        prim->molec->atoms[objects[i].object.atom->num].draw_order = i;
      else if (prim->MO_surf )
        prim->MO_surf->molec->atoms[objects[i].object.atom->num].draw_order = i;
      break;
    case TRIANGLE_3D:
      draw_triangle(objects[i].object.triangle,vertex_store,surf);
      break;
    case CONT_POINT_3D:
      g_change_linewidth(1);
      if( objects[i].object.contour_point.val > 0 ){
        g_change_linestyle(0);
      }else if( objects[i].object.contour_point.val < 0 ){
        g_change_linestyle(2);
      } else{
        g_change_linestyle(1);
      }
      if( objects[i].object.contour_point.next ){
        g_line(objects[i].object.contour_point.loc->x,
               objects[i].object.contour_point.loc->y,
               objects[i].object.contour_point.next->x,
               objects[i].object.contour_point.next->y);
      }
      break;
    case LINE_3D:
      g_change_linestyle(0);
      g_draw_breaking_line(objects[i].object.line.end1->x,
                           objects[i].object.line.end1->y,
                           objects[i].object.line.end2->x,
                           objects[i].object.line.end2->y,2);
      break;

    }

  }



  /* draw in the contours now if they were not sorted in */
  if( surf && surf->MO_contours && surf->display_conts && !surf->sort_conts )
    draw_contours(surf,contour_store);



  /* draw in the axes if we need to */
  if( do_axes ){
    /* set up and transform the points */
    axes.points[0].x = 0.5;axes.points[0].y=0.5;axes.points[0].z=0.5;
    axes.points[1].x = 1.5;axes.points[1].y=0.5;axes.points[1].z=0.5;
    axes.points[2].x = 0.5;axes.points[2].y=1.5;axes.points[2].z=0.5;
    axes.points[3].x = 0.5;axes.points[3].y=0.5;axes.points[3].z=1.5;
    for(i=0;i<4;i++){
      transform(&(axes.points[i]));
      axes.points[i].x = 100 - 10 + axes.points[i].x;
      axes.points[i].y = 100 - 10 + axes.points[i].y;
    }
    draw_axes(&axes);
  }

#if 0
  if( molec && molec->draw_lattice ){
    g_change_linestyle(0);
    g_change_linewidth(2);
    for(i=1;i<num_lattice_p;i++){
      g_line(lattice[0].x,lattice[0].y,lattice[i].x,lattice[i].y);
    }
  }
#endif
  if( molec && animating_molecule ){
    if( animating_molecule ) molec->current_frame++;
  }
}

