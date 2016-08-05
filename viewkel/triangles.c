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

#include "viewkel.h"


/*******************

  this has got the stuff for dealing with surfaces made up
   of triangles.

  re-created by greg Landrum  Jun 1995

*******************/

/***
  Recent Edit History
   26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)
***/


/****************************************************************************
 *
 *                   Procedure convert_explicit_triangles
 *
 * Arguments: explicit: pointer to explicit_triangle_type
 *        num_explicit: int
 *           triangles: pointer to pointer to triangle_type
 *             num_tri: pointer to int
 *              points: pointer to pointer to point_type
 *          num_points: pointer to int
 *
 * Returns: none
 *
 * Action: This takes the triangles in 'explicit, which have all of
 *    their vertices specified, and convert them to triangles with
 *    pointers to vertices stored in 'points.  This is to allow more
 *    efficient drawing.
 *
 ****************************************************************************/
void convert_explicit_triangles(explicit_triangle_type *explicit,
                                int num_explicit,triangle_type **triangles,
                                int *num_tri,point_type **points,
                                int *num_points)
{
  int i,j,k;
  int done;
  int max_tri,max_points;
  int which_vert;

  /* start out by allocating some memory to store the points */
  max_points = num_explicit;
  *points = (point_type *)D_CALLOC(max_points,sizeof(point_type));
  if( !(*points) ) fatal("Can't get space for points in convert_explicit_triangles");

  /*******

    get memory for the triangles

    for now there will be as many triangles as there are explicit
    triangles.  Later, however, it may be desirable to add some culling
    of small or strange triangles.

  *******/
  max_tri = num_explicit;
  *triangles = (triangle_type *)D_CALLOC(max_tri,sizeof(triangle_type));
  if( !(*triangles) )
    fatal("Can't get space for triangles in convert_explicit_triangles");

  /******

    now loop through the explicit triangles:
      add the vertices to the list stored in 'points and
      set up the corresponding triangles

  *******/
  *num_tri = 0;
  *num_points = 0;
  for(i=0;i<num_explicit;i++){
    /* loop over vertices */
    for(j=0;j<3;j++){

      /* check to see if this vertex is already in the list */
      done = 0;
      for(k=0;k<(*num_points)&&!done;k++){
        if( (*points)[k].x == explicit[i].vertices[j].x &&
           (*points)[k].y == explicit[i].vertices[j].y &&
           (*points)[k].z == explicit[i].vertices[j].z ){

          /* okay, this is a match, set the vertex number */
          done = 1;
          which_vert = k;
        }
      }

      if( !done ){
        /* add this vertex to the points list */
        (*points)[*num_points].x = explicit[i].vertices[j].x;
        (*points)[*num_points].y = explicit[i].vertices[j].y;
        (*points)[*num_points].z = explicit[i].vertices[j].z;

        /* set the vertex number */
        which_vert = (*num_points);
        (*num_points)++;

        /* do we need more memory? */
        if( (*num_points) == max_points ){
          max_points += 512;
          *points = (point_type *)D_REALLOC((void *)(*points),
                                          max_points*sizeof(point_type));
          if( !(*points) )
            fatal("Can't D_REALLOC points in convert_explicit_triangles");
        }
      }

      /* now set up the triangle vertex */
      (*triangles)[(*num_tri)].vertices[j] = which_vert;

      /* that's it for now */
    }


    /* set the center of each triangle */
    (*triangles)[(*num_tri)].center.x =
      ((*points)[(*triangles)[(*num_tri)].vertices[0]].x +
       (*points)[(*triangles)[(*num_tri)].vertices[1]].x +
       (*points)[(*triangles)[(*num_tri)].vertices[2]].x)/3.0;
    (*triangles)[(*num_tri)].center.y =
      ((*points)[(*triangles)[(*num_tri)].vertices[0]].y +
       (*points)[(*triangles)[(*num_tri)].vertices[1]].y +
       (*points)[(*triangles)[(*num_tri)].vertices[2]].y)/3.0;
    (*triangles)[(*num_tri)].center.z =
      ((*points)[(*triangles)[(*num_tri)].vertices[0]].z +
       (*points)[(*triangles)[(*num_tri)].vertices[1]].z +
       (*points)[(*triangles)[(*num_tri)].vertices[2]].z)/3.0;


    /**** TEMPORARY HACK ********/
    (*triangles)[(*num_tri)].color = 0;

    /*****

      increment the number of triangles and see if we need more
      memory

    ******/
    (*num_tri)++;
    if( (*num_tri) == max_tri ){
      max_tri += 512;
      (*triangles) = (triangle_type *)D_REALLOC((void *)(*triangles),
                                              max_tri*sizeof(triangle_type));
      if( !(*triangles) )
        fatal("Can't D_REALLOC triangles in convert_explicit_triangles.");
    }
  }

  /* that's it */
}


/****************************************************************************
 *
 *                   Procedure remove_degen_triangles
 *
 * Arguments: MO_surf: pointer to MO_surface_type
 *
 * Returns: none
 *
 * Action: removes all the triangles in 'part_surf that are degenerate
 *   (identical or have 2 parallel sides).
 *
 ****************************************************************************/
void remove_degen_triangles(MO_surface_type *MO_surf)
{
  int i,j;
  int *ditched;
  triangle_type *triangle_store,*curr_tri;
  vertex_type *vertices;
  point_type *c1,*c2,temp_v;
  point_type edge1,edge2;
  int num_kept;

  return;

  /* get memory for the temporary triangle and ditched arrays */
  triangle_store = (triangle_type *)
    D_MALLOC(MO_surf->num_triangles*sizeof(triangle_type));
  ditched = (int *)D_CALLOC(MO_surf->num_triangles,sizeof(triangle_type));
  if( !triangle_store || !ditched ) fatal("Can't allocate temporary storage.\n");


  vertices = MO_surf->triangle_vertices;

  num_kept = 0;
  /* loop over the original triangle list */
  for( i=0; i<MO_surf->num_triangles; i++ ){
    if( !ditched[i] ){
      curr_tri = &(MO_surf->triangles[i]);
      c1 = &curr_tri->center;

      for( j=i+1; j<MO_surf->num_triangles; j++ ){
        if( !ditched[j] ){
          /* check the distance between the centers */
          c2 = &(MO_surf->triangles[j].center);
          if( V3SquaredLength(V3Sub(c1,c2,&temp_v)) <= .00001 ){
            ditched[j] = 1;
          }
        }
      }

      /* check the sides of triangle i to see if we should ditch them */
      V3Sub(&(vertices[curr_tri->vertices[1]].position),
            &(vertices[curr_tri->vertices[0]].position),
            &edge1);
      V3Sub(&(vertices[curr_tri->vertices[2]].position),
            &(vertices[curr_tri->vertices[0]].position),
            &edge2);

      /******

        do the cross product between the two edges and see if it's
        zero.

      *******/
      if( V3SquaredLength(V3Cross(&edge1,&edge2,&temp_v)) > .000001 ){
        bcopy(&(MO_surf->triangles[i]), &(triangle_store[num_kept]),
              sizeof(triangle_type));
        num_kept++;
      }else{
        ditched[i] = 1;
      }
    }
  }
  /* free the original array */
  D_FREE(MO_surf->triangles);
  D_FREE(ditched);

  /* D_REALLOC the new data and set a pointer*/
  MO_surf->triangles = (triangle_type *)
    D_REALLOC(triangle_store,num_kept*sizeof(triangle_type));
  if( !MO_surf->triangles ) fatal("Can't D_REALLOC triangle array.  This is weird.");

  fprintf(stderr,"%d of %d triangles were kept (%4.2lf %%).\n",num_kept,
          MO_surf->num_triangles,
          100.0*(float)num_kept/(float)MO_surf->num_triangles);
  MO_surf->num_triangles = num_kept;


}




/****************************************************************************
 *
 *                   Procedure calc_triangle_centers
 *
 * Arguments: MO_surf: pointer to MO_surface_type
 *
 * Returns: none
 *
 * Action: finds the centers of all the triangles in 'part_surf.
 *
 ****************************************************************************/
void calc_triangle_centers(MO_surface_type *MO_surf)
{
  int i;
  point_type *verts[3];
  triangle_type *temp_tri;




  for( i=0; i<MO_surf->num_triangles; i++ ){
    /* get pointers to the triangle and its vertices */
    temp_tri = &MO_surf->triangles[i];

    verts[0] =
      &(MO_surf->triangle_vertices[temp_tri->vertices[0]].position);
    verts[1] =
      &(MO_surf->triangle_vertices[temp_tri->vertices[1]].position);
    verts[2] =
      &(MO_surf->triangle_vertices[temp_tri->vertices[2]].position);


    if( temp_tri->color != -1 ){
      temp_tri->center.x =
        (verts[0]->x + verts[1]->x + verts[2]->x)/3.0;
      temp_tri->center.y =
        (verts[0]->y + verts[1]->y + verts[2]->y)/3.0;
      temp_tri->center.z =
        (verts[0]->z + verts[1]->z + verts[2]->z)/3.0;
    }
  }
}


/****************************************************************************
 *
 *                   Procedure save_triangle_locs
 *
 * Arguments:   num_args:  an integer
 *             MO_surf_p: array of pointers to char
 *
 * Returns: none
 *
 * Action: writes the locations of the triangles in 'MO_surf
 *    to an output file.
 *
 *  'num_args is just used because this function is intended to be
 *    called from a function button
 *
 ****************************************************************************/
void save_triangle_locs(int num_args,char *MO_surf_p[MAX_ARGS])
{
  int outfile;
  char filename[240];
  MO_surface_type *MO_surf;

  MO_surf = (MO_surface_type *)MO_surf_p[0];

/*  strcpy(filename,MO_surf->filename);
  strcat(filename,".tri");
*/
  sprintf(filename,"%s.num%d.tri",MO_surf->filename,MO_surf->active_MO);

#ifdef X_GRAPHICS
  outfile = open(filename,O_RDWR|O_APPEND|O_CREAT,"rw");
  if( outfile == -1 ){
    error("Can't open triangle file for binary I/O.");
    return;
  }

  /* write the number of triangles and number of vertices */
  write(outfile,&(MO_surf->num_vertices),sizeof(int));
  write(outfile,&(MO_surf->num_triangles),sizeof(int));

  /* write the vertices */
  write(outfile,MO_surf->triangle_vertices,
        MO_surf->num_vertices*sizeof(vertex_type));

  /* write the triangles */
  write(outfile,MO_surf->triangles,
        MO_surf->num_triangles*sizeof(triangle_type));

  printf("Wrote %d vertices and %d triangles\n",MO_surf->num_vertices,
         MO_surf->num_triangles);

  close(outfile);
#else
   fprintf(stderr,"Triangle Saves not yet implemented on this system.  Sorry\n");
#endif
}





/****************************************************************************
 *
 *                   Procedure read_triangle_locs
 *
 * Arguments:   num_args:  an integer
 *             MO_surf_p: array of pointers to char
 *
 * Returns: none
 *
 * Action: reads in the locations of the triangles in the output file
 *         which corresponds to 'MO_surf
 *
 *  'num_args is just used because this function is intended to be
 *    called from a function button
 *
 ****************************************************************************/
void read_triangle_locs(int num_args, char *MO_surf_p[MAX_ARGS])
{
  int infile;
  char filename[240];
  MO_surface_type *MO_surf;

  MO_surf = (MO_surface_type *)MO_surf_p[0];

/*  strcpy(filename,MO_surf->filename);
  strcat(filename,".tri");
*/
  sprintf(filename,"%s.num%d.tri",MO_surf->filename,MO_surf->active_MO);

#ifdef X_GRAPHICS
  infile = open(filename,O_RDONLY,"r");
  if( infile == -1 ){
    error("Can't open triangle file for binary I/O.");
    return;
  }

  /* if there's any residual crap allocated, get rid of it now */
  if( MO_surf->num_triangles ) D_FREE(MO_surf->triangles);
  if( MO_surf->num_vertices ) D_FREE(MO_surf->triangle_vertices);


  /* read out the number of triangles and number of vertices */
  read(infile,&(MO_surf->num_vertices),sizeof(int));
  read(infile,&(MO_surf->num_triangles),sizeof(int));

  /* get space for the vertices and triangles */
  MO_surf->triangle_vertices =
    (vertex_type *)D_CALLOC(MO_surf->num_vertices,sizeof(vertex_type));
  if( !(MO_surf->triangle_vertices) ) fatal("Can't get memory for vertices");
  MO_surf->triangles =
    (triangle_type *)D_CALLOC(MO_surf->num_triangles,sizeof(triangle_type));
  if( !(MO_surf->triangles) ) fatal("Can't get memory for triangles");


  /* read the vertices */
  read(infile,MO_surf->triangle_vertices,
        MO_surf->num_vertices*sizeof(vertex_type));

  /* write the triangles */
  read(infile,MO_surf->triangles,
        MO_surf->num_triangles*sizeof(triangle_type));

  printf("Read in %d vertices and %d triangles\n",MO_surf->num_vertices,
         MO_surf->num_triangles);

  close(infile);
#else
   fprintf(stderr,"Triangle Saves not yet implemented on this system.  Sorry\n");
#endif

}




