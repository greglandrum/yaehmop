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

  this has got the stuff for dealing with coordination
  polyhedra

*********/
#include "viewkel.h"
#include "polyhed.h"

/***
 Recent Edit History
  18.05.98 gL: added saving of polyhedra centers and radii to
    help with the fileio process, which must eventually be finished.
  26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)

***/

void set_tri_parms(triangle_list_type *temp_tri,atom_type *active_atoms){
  int i,j,k;
  point_type v1,v2;
  float proj;

  i = temp_tri->tri.vertices[0];
  j = temp_tri->tri.vertices[1];
  k = temp_tri->tri.vertices[2];

  /* set center */
  temp_tri->tri.center.x = (active_atoms[i].loc.x +
                            active_atoms[j].loc.x +
                            active_atoms[k].loc.x ) / 3;
  temp_tri->tri.center.y = (active_atoms[i].loc.y +
                            active_atoms[j].loc.y +
                            active_atoms[k].loc.y ) / 3;
  temp_tri->tri.center.z = (active_atoms[i].loc.z +
                            active_atoms[j].loc.z +
                            active_atoms[k].loc.z ) / 3;

  /* set normal */
  V3Sub(&active_atoms[i].loc,&active_atoms[j].loc,&v1);
  V3Sub(&active_atoms[i].loc,&active_atoms[k].loc,&v2);

  /* check to make sure the points aren't co-linear */
  V3Cross(&v1,&v2,&temp_tri->tri.normal);
  if(V3Length(&temp_tri->tri.normal) > DOT_TOL){
    /* they aren't.  normalize the normal */
    V3Normalize(&temp_tri->tri.normal);

    /* find the bonding box */
    temp_tri->bmax.x = -10000;
    temp_tri->bmax.y = -10000;
    temp_tri->bmax.z = -10000;
    temp_tri->bmin.x = 10000;
    temp_tri->bmin.y = 10000;
    temp_tri->bmin.z = 10000;
    if(active_atoms[i].loc.x > temp_tri->bmax.x )
      temp_tri->bmax.x = active_atoms[i].loc.x;
    if(active_atoms[j].loc.x > temp_tri->bmax.x )
      temp_tri->bmax.x = active_atoms[j].loc.x;
    if(active_atoms[k].loc.x > temp_tri->bmax.x )
      temp_tri->bmax.x = active_atoms[k].loc.x;
    if(active_atoms[i].loc.y > temp_tri->bmax.y )
      temp_tri->bmax.y = active_atoms[i].loc.y;
    if(active_atoms[j].loc.y > temp_tri->bmax.y )
      temp_tri->bmax.y = active_atoms[j].loc.y;
    if(active_atoms[k].loc.y > temp_tri->bmax.y )
      temp_tri->bmax.y = active_atoms[k].loc.y;
    if(active_atoms[i].loc.z > temp_tri->bmax.z )
      temp_tri->bmax.z = active_atoms[i].loc.z;
    if(active_atoms[j].loc.z > temp_tri->bmax.z )
      temp_tri->bmax.z = active_atoms[j].loc.z;
    if(active_atoms[k].loc.z > temp_tri->bmax.z )
      temp_tri->bmax.z = active_atoms[k].loc.z;

    if(active_atoms[i].loc.x < temp_tri->bmin.x )
      temp_tri->bmin.x = active_atoms[i].loc.x;
    if(active_atoms[j].loc.x < temp_tri->bmin.x )
      temp_tri->bmin.x = active_atoms[j].loc.x;
    if(active_atoms[k].loc.x < temp_tri->bmin.x )
      temp_tri->bmin.x = active_atoms[k].loc.x;
    if(active_atoms[i].loc.y < temp_tri->bmin.y )
      temp_tri->bmin.y = active_atoms[i].loc.y;
    if(active_atoms[j].loc.y < temp_tri->bmin.y )
      temp_tri->bmin.y = active_atoms[j].loc.y;
    if(active_atoms[k].loc.y < temp_tri->bmin.y )
      temp_tri->bmin.y = active_atoms[k].loc.y;
    if(active_atoms[i].loc.z < temp_tri->bmin.z )
      temp_tri->bmin.z = active_atoms[i].loc.z;
    if(active_atoms[j].loc.z < temp_tri->bmin.z )
      temp_tri->bmin.z = active_atoms[j].loc.z;
    if(active_atoms[k].loc.z < temp_tri->bmin.z )
      temp_tri->bmin.z = active_atoms[k].loc.z;


    /* make sure the normal points away from the central atom */
    V3Sub(&temp_tri->tri.center,&active_atoms[0].loc,&v1);
    proj = V3Dot(&temp_tri->tri.normal,&v1);
    if(proj < 0 ){
      V3Negate(&temp_tri->tri.normal);
    }
  } else{
    error("somehow a bad triangle slipped through");
  }
}


/****************************************************************************
 *
 *                   Procedure gen_coord_polyhed
 *
 * Arguments: molec: pointer to molec_type
 *
 * Returns: none
 *
 * Action:  generates coordination polyhedra for all selected atoms.
 *
 ****************************************************************************/
void gen_coord_polyhed(molec_type *molec)
{
  int i,num_polyhed;
  atom_type *atoms;
  static float cut_off=3.0;
  int first = 1;
  if( !molec ) FATAL_BUG("gen_coord_polyhed called with bogus molec");

  if(molec->num_frames > 1)
    atoms = &(molec->atoms[molec->current_frame*molec->num_atoms]);
  else atoms = molec->atoms;

  readfloatparm("coordination cut off",&cut_off);

  num_polyhed=0;
  for(i=0;i<molec->num_atoms;i++){
    if( atoms[i].is_selected ){
      num_polyhed++;
      find_coord_polyhed(molec,i,cut_off,first);
      first = 0;
    }
  }
  if( num_polyhed ){
    molec->polyhed_centers = (int *)D_CALLOC(num_polyhed,sizeof(int));
    molec->polyhed_rads = (float *)D_CALLOC(num_polyhed,sizeof(float));
    if(!molec->polyhed_centers || !molec->polyhed_rads){
      error("can't allocate polyhed centers or rads\n");
      return;
    } else{
      num_polyhed =0;
      for(i=0;i<molec->num_atoms;i++){
        if( atoms[i].is_selected ){
          molec->polyhed_centers[num_polyhed] = i;
          molec->polyhed_rads[num_polyhed] = cut_off;
          num_polyhed++;
        }
      }
    }
  }

}
/****************************************************************************
 *
 *                   Procedure find_coord_polyhed
 *
 * Arguments: molec: pointer to molec_type
 *       which_atom: int
 *           cutoff: float
 *       initialize: int
 *
 *
 * Returns: none
 *
 * Action:  determines the faces of the coordination
 *  polyhedron around 'which_atom.  'cutoff determines
 *  the maximum length for things to be considered
 *  coordinated.
 *
 *  if 'initialize is nonzero, existing triangles will be
 *   blown out.
 *
 ****************************************************************************/
void find_coord_polyhed(molec_type *molec,int which_atom,float cut_off,
                        int initialize)
{
  int i,j;
  atom_type *atoms,*active_atoms,*temp_atoms;
  int num_active_atoms,max_active_atoms;
  triangle_list_type *tri_list;
  triangle_list_type *tri1;
  point_type *p1,*p2;
  int num_tris;
  float dist;

  active_atoms = 0;
  tri_list =0;

  if( !molec ) FATAL_BUG("find_coord_polyhed called with bogus molec");

  if(molec->num_frames > 1)
    atoms = &(molec->atoms[molec->current_frame*molec->num_atoms]);
  else atoms = molec->atoms;

  /******

    find the atoms which make up the coordination polyhedron

  ******/

  /* first allocate space to store them */
  max_active_atoms = 16;
  active_atoms = (atom_type *)D_CALLOC(max_active_atoms,sizeof(atom_type));
  if(!active_atoms) fatal("can't allocate active atom array.");



  memcpy((void *)&(active_atoms[0]),
         (void *)&(atoms[which_atom]),
         sizeof(atom_type));

  /* now find the other atoms */
  p1 = &atoms[which_atom].loc;
  num_active_atoms = 1;
  for(i=0;i<molec->num_atoms;i++){
    if( i != which_atom && atoms[i].type[0] != '&'){
      p2 = &atoms[i].loc;
      dist = V3DistanceBetween2Points(p1,p2);
      if( dist < cut_off ){
        /* yep, this one is good, copy it in */
        memcpy((void *)&(active_atoms[num_active_atoms]),
               (void *)&(atoms[i]),
               sizeof(atom_type));
        num_active_atoms++;
        /* do we need more memory?  I hope not */
        if(num_active_atoms == max_active_atoms){
          max_active_atoms += 16;
          temp_atoms = active_atoms;
          active_atoms = (atom_type *)D_CALLOC(max_active_atoms,
                                             sizeof(atom_type));
          if(!active_atoms){
            fatal("can't D_REALLOCate active atom array.");
          }
          memcpy((void *)temp_atoms,(void *)active_atoms,
                 num_active_atoms*sizeof(atom_type));
        }
      }
    }
  }

  /***********

    okay, we have the group of atoms which are within
    the cut_off radius.  The next step is to go through
    and set up the polyhedron using the chull code

  **************/
  gen_chull(&(active_atoms[1]),num_active_atoms-1,&tri_list);


  tri1 = tri_list;
  num_tris = 0;
  while(tri1){
    num_tris++;
    tri1 = tri1->next;
  }
  if(initialize ){
    if( molec->triangles) D_FREE(molec->triangles);

    molec->triangles = (triangle_type *)D_CALLOC(num_tris,
                                               sizeof(triangle_type));

    if(!molec->triangles) fatal("Can't get space to store triangles.");

    if(molec->polyhed_verts)D_FREE(molec->polyhed_verts);
    molec->polyhed_verts = (vertex_type *)D_CALLOC(num_active_atoms,
                                               sizeof(vertex_type));
    if(!molec->polyhed_verts)
      fatal("Can't get space to store polyedron vertices.");
    molec->num_triangles = 0;
    molec->num_polyhed_verts = 0;
  } else{

    molec->triangles = (triangle_type *)
      D_REALLOC((void *)molec->triangles,
              (num_tris+molec->num_triangles)*sizeof(triangle_type));
    if(!molec->triangles)fatal("can't D_REALLOCated molec->triangles.");
    molec->polyhed_verts = (vertex_type *)
      D_REALLOC((void *)molec->polyhed_verts,
              (num_active_atoms+molec->num_polyhed_verts)*sizeof(vertex_type));
    if(!molec->polyhed_verts)fatal("can't D_REALLOCated molec->polyhed_verts.");

  }
  for(i=0;i<num_active_atoms;i++){
    memcpy((void *)&molec->polyhed_verts[i+molec->num_polyhed_verts].position,
           (void *)&active_atoms[i].loc,sizeof(point_type));
  }



  i=molec->num_triangles;

  tri1 = tri_list;
  while(tri1){
    for(j=0;j<3;j++){
      tri1->tri.vertices[j]++;
    }
    set_tri_parms(tri1,active_atoms);
    for(j=0;j<3;j++){
      tri1->tri.vertices[j] += molec->num_polyhed_verts;
    }
    memcpy((void *)&molec->triangles[i],(void *)&(tri1->tri),
           sizeof(triangle_type));
    molec->triangles[i].color = atoms[which_atom].color;
    i++;
    tri1 = tri1->next;
  }

  molec->num_triangles += num_tris;
  molec->num_polyhed_verts += num_active_atoms;



fprintf(stderr,"Find polyhedron\n");
for(i=0;i<num_active_atoms;i++){
  fprintf(stderr,"%d %s\n",i,active_atoms[i].type);
}
tri1 = tri_list;
while(tri1){
  fprintf(stderr,"tri: %d %d %d",
          tri1->tri.vertices[0],tri1->tri.vertices[1],tri1->tri.vertices[2]);
  fprintf(stderr,"\tcenter: (%f %f %f)\n",
          tri1->tri.center.x,tri1->tri.center.y,tri1->tri.center.z);
  fprintf(stderr,"\tnormal: (%f %f %f)\n",
          tri1->tri.normal.x,tri1->tri.normal.y,tri1->tri.normal.z);
  tri1 = tri1->next;
}


  /* clean up the memory we allocated */
  tri1 = tri_list;
  while(tri1){
    tri_list = tri1->next;
    D_FREE(tri1);
    tri1 = tri_list;
  }
  if( active_atoms ) D_FREE(active_atoms);
}
