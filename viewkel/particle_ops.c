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

  this has got the stuff for dealing with the particles used to
    construct particle_surfaces as well as the triangles used
    to build those surfaces

  created by greg Landrum  March 1995

*******************/




/****************************************************************************
 *
 *                   Procedure read_triangle_locs
 *
 * Arguments:   num_args:  an integer
 *             part_surf_p: array of pointers to char
 *
 * Returns: none
 *
 * Action: reads the locations of the triangles in an input file
 *     into 'part_surf
 *
 *  'num_args is just used because this function is intended to be
 *    called from a function button
 *
 ****************************************************************************/
void read_triangle_locs(num_args,part_surf_p)
  int num_args;
  char *part_surf_p[MAX_ARGS];
{
  int i;
  FILE *infile;
  char filename[240];
  triangle_type *temp_tri;
  particle_surface_type *part_surf;
  int phase;

  part_surf = (particle_surface_type *)part_surf_p[0];

  strcpy(filename,part_surf->filename);
  strcat(filename,".tri");

  infile = fopen(filename,"r");
  if( !infile ){
    error("Can't open output file.");
    return;
  }

  fscanf(infile,"%d",&part_surf->num_triangles);

  part_surf->triangles = (triangle_type *)calloc(part_surf->num_triangles,
                                                 sizeof(triangle_type));
  if( !part_surf->triangles ) fatal("can't get triangle storage");

  fprintf(stderr,"Reading %d triangles\n",part_surf->num_triangles);
  for(i=0;i<part_surf->num_triangles;i++){
    temp_tri = &(part_surf->triangles[i]);
    fscanf(infile,"%d %d %d %d",
            &temp_tri->vertices[0],&temp_tri->vertices[1],
            &temp_tri->vertices[2],&phase);
    temp_tri->color = phase;
  }
}


/****************************************************************************
 *
 *                   Procedure save_triangle_locs
 *
 * Arguments:   num_args:  an integer
 *             part_surf_p: array of pointers to char
 *
 * Returns: none
 *
 * Action: writes the locations of the triangles in 'part_surf
 *    to an output file.
 *
 *  'num_args is just used because this function is intended to be
 *    called from a function button
 *
 ****************************************************************************/
void save_triangle_locs(num_args,part_surf_p)
  int num_args;
  char *part_surf_p[MAX_ARGS];
{
  int i;
  FILE *outfile;
  char filename[240];
  triangle_type *temp_tri;
  particle_surface_type *part_surf;

  part_surf = (particle_surface_type *)part_surf_p[0];

  strcpy(filename,part_surf->filename);
  strcat(filename,".tri");

  outfile = fopen(filename,"w+");
  if( !outfile ){
    error("Can't open output file.");
    return;
  }

  fprintf(outfile,"%d\n",part_surf->num_triangles);

  for(i=0;i<part_surf->num_triangles;i++){
    temp_tri = &(part_surf->triangles[i]);
    fprintf(outfile,"%d %d %d %d\n",
            temp_tri->vertices[0],temp_tri->vertices[1],
            temp_tri->vertices[2],temp_tri->color);
  }
  fclose(outfile);
}





/****************************************************************************
 *
 *                   Procedure save_particle_locs
 *
 * Arguments:   num_args:  an integer
 *             part_surf_p: array of pointers to char
 *
 * Returns: none
 *
 * Action: writes the locations of the particles in 'part_surf
 *    to an output file.
 *
 *  'num_args is just used because this function is intended to be
 *    called from a function button
 *
 ****************************************************************************/
void save_particle_locs(num_args,part_surf_p)
  int num_args;
  char *part_surf_p[MAX_ARGS];
{
  int i;
  FILE *outfile;
  char filename[240];
  particle_type *temp_part;
  particle_surface_type *part_surf;

  part_surf = (particle_surface_type *)part_surf_p[0];

  strcpy(filename,part_surf->filename);
  strcat(filename,".part");

  outfile = fopen(filename,"w+");
  if( !outfile ){
    error("Can't open output file.");
    return;
  }

  fprintf(outfile,"%d\n",part_surf->num_particles);
  fprintf(outfile,"%lf %lf\n",part_surf->surface_value,
          part_surf->surface_tolerance);

  temp_part = part_surf->bound_particles;
  for(i=0;i<part_surf->num_particles;i++){
    if( !temp_part ) FATAL_BUG("ran out of particles in save_particle_locs");

    fprintf(outfile,"%6.4lf %6.4lf %6.4lf %d\n",
            temp_part->loc.x,temp_part->loc.y,temp_part->loc.z,
            temp_part->phase);
    temp_part = temp_part->next;
  }
  fclose(outfile);
}




/****************************************************************************
 *
 *                   Procedure read_particle_locs
 *
 * Arguments:   num_args:  an integer
 *             part_surf_p: array of pointers to char
 *
 * Returns: none
 *
 * Action: reads particles into 'part_surf from an output file.
 *
 *  'num_args is just used because this function is intended to be
 *    called from a function button
 *
 ****************************************************************************/
void read_particle_locs(num_args,part_surf_p)
  int num_args;
  char *part_surf_p[MAX_ARGS];
{
  int i;
  FILE *infile;
  char filename[240];
  particle_type *temp_part,*next_part;
  float val;
  particle_surface_type *part_surf;
  int phase;

  part_surf = (particle_surface_type *)part_surf_p[0];

  strcpy(filename,part_surf->filename);
  strcat(filename,".part");

  infile = fopen(filename,"r");
  if( !infile ){
    error("Can't open input file.");
    return;
  }

  /* free up any memory that's currently being used */
  if( part_surf->bound_particles ){
    temp_part = part_surf->bound_particles;
    while(temp_part){
      next_part = temp_part->next;
      free(temp_part);
      temp_part = next_part;
    }
  }

  if( part_surf->triangles ){
    part_surf->triangles = 0;
    part_surf->num_triangles = 0;
  }

  fscanf(infile,"%d",&(part_surf->num_particles));
  fscanf(infile,"%lf %lf",&(part_surf->surface_value),
         &(part_surf->surface_tolerance));

  printf("Reading in %d particles.\n",part_surf->num_particles);

  for(i=0; i<part_surf->num_particles;i++){
    temp_part = (particle_type *)calloc(1,sizeof(particle_type));
    if(!temp_part) fatal("can't get memory to read in a particle");

    temp_part->next = part_surf->bound_particles;
    part_surf->bound_particles = temp_part;

    fscanf(infile,"%lf %lf %lf %d",
            &(temp_part->loc.x),&(temp_part->loc.y),
            &(temp_part->loc.z),&phase);
    temp_part->phase = phase;
    temp_part->num = i;
  }
  display("Done Reading");
}



/****************************************************************************
 *
 *                   function compare_triangles
 *
 * Arguments: tri1,tri2: pointers to triangle_type
 *
 * Returns: int
 *
 * Action: returns non-zero if 'tri1 and 'tri2 are the same
 *
 ****************************************************************************/
int compare_triangles(tri1,tri2)
  triangle_type *tri1,*tri2;
{
  if( tri1->vertices[0] == tri2->vertices[0] ||
     tri1->vertices[0] == tri2->vertices[1]  ||
     tri1->vertices[0] == tri2->vertices[2] ){
    if( tri1->vertices[1] == tri2->vertices[0] ||
       tri1->vertices[1] == tri2->vertices[1]  ||
       tri1->vertices[1] == tri2->vertices[2] ){
      if( tri1->vertices[2] == tri2->vertices[0] ||
         tri1->vertices[2] == tri2->vertices[1]  ||
         tri1->vertices[2] == tri2->vertices[2] ){
        return(1);
      }
    }
  }
  return(0);
}





/****************************************************************************
 *
 *                   Procedure build_part_lookup_list
 *
 * Arguments: particles: pointer to particle_type
 *        num_particles: int
 *          lookup_list:  pointer to int
 *
 * Returns: none
 *
 * Action:  builds 'lookup_list, an array which allows a particular
 *   numbered particle to be found in 'particles instantly.  (it's a
 *   hash table)
 *
 ****************************************************************************/
void build_part_lookup_list(particles,num_particles,lookup_list)
  particle_type *particles;
  int num_particles, *lookup_list;
{
  int i;

  for(i=0;i<num_particles;i++) lookup_list[particles[i].num] = i;
}



/****************************************************************************
 *
 *                   Procedure find_and_connect_nearest_neighbors
 *
 * Arguments: particles: pointer to particle_type
 *         num_particles: int
 *
 * Returns: none
 *
 * Action:  goes through the list of particles 'particles and connects
 *  up nearest neighbors for each one.
 *
 *   At the moment this is only intended to work within a slice of
 *    data and will only connect to the two nearest neighbors in
 *    the current slice.
 *
 ****************************************************************************/
void find_and_connect_nearest_neighbors(particles,num_particles)
  particle_type *particles;
  int num_particles;
{
  int i,j;
  particle_type *i_part,*j_part,*closest1,*closest2;
  float dist;
  float dist1,dist2;

  i_part = particles;
  for(i=0;i<num_particles;i++){
    dist1 = dist2 = 1e10;
    closest1 = closest2 = 0;
    if( !i_part ) FATAL_BUG("ran out of i_parts in connection routine.");
    if( i_part->num_neighbors < 2 ){
      j_part = i_part->next;
      for(j=i+1;j<num_particles;j++){
        if( !j_part ) FATAL_BUG("ran out of j_parts in connection routine.");

        if( j_part->num_neighbors < 2 && j_part->phase == i_part->phase ){
          dist = (i_part->loc.x-j_part->loc.x)*(i_part->loc.x-j_part->loc.x) +
            (i_part->loc.y-j_part->loc.y)*(i_part->loc.y-j_part->loc.y) +
              (i_part->loc.z-j_part->loc.z)*(i_part->loc.z-j_part->loc.z);
          if( dist < dist1 && dist < .10){
            dist2 = dist1;
            closest2 = closest1;
            dist1 = dist;
            closest1 = j_part;
          } else if(dist < dist2 && dist < .10){
            dist2 = dist;
            closest2 = j_part;
          }
        }
        j_part = j_part->next;
      }
      if( closest1 ){
        i_part->neighbors[i_part->num_neighbors++] = closest1;
        closest1->neighbors[closest1->num_neighbors++] = i_part;
      }
      if( i_part->num_neighbors < 2 && closest2 ){
        i_part->neighbors[i_part->num_neighbors++] = closest2;
        closest2->neighbors[closest2->num_neighbors++] = i_part;
      }
    }
    i_part = i_part->next;
  }
}


/****************************************************************************
 *
 *                   Procedure split_a_particle
 *
 * Arguments: particle_surf: pointer to particle_surface_type.
 *                 particle: pointer to particle_type
 *
 * Returns: none
 *
 * Action:  divides 'particle into 2 pieces
 *
 ****************************************************************************/
void split_a_particle(particle_surf,particle)
  particle_surface_type *particle_surf;
  particle_type *particle;
{
  particle_type *new_part;
  if( !particle_surf || !particle )
    FATAL_BUG("Bad arguments to split_a_particle");

  particle_surf->num_particles++;

  new_part = (particle_type *)calloc(1,sizeof(particle_type));
  if( !new_part )
    fatal("Can't allocate a particle.\n");

  new_part->next = particle_surf->bound_particles;
  particle_surf->bound_particles = new_part;

  new_part->loc.x = particle->loc.x + (float)my_drand(.05);
  new_part->loc.y = particle->loc.y + (float)my_drand(.05);
  new_part->loc.z = particle->loc.z + (float)my_drand(.05);

  particle->rad *= .65;
  new_part->rad = particle->rad;
}




/****************************************************************************
 *
 *                   Procedure split_bound_particles
 *
 * Arguments: particle_surf: pointer to particle_surface_type.
 *
 * Returns: none
 *
 * Action:  divides all the particles on the surface into 2 particles.
 *   in doing this, the repulsion radius of all the particles is lowered.
 *
 ****************************************************************************/
void split_bound_particles(particle_surf)
  particle_surface_type *particle_surf;
{
  particle_type *curr_part,*new_part;

  if( !particle_surf )
    FATAL_BUG("Bad particle_surf in split_bound_particles.");

  curr_part = particle_surf->bound_particles;
  particle_surf->num_particles *= 2;

  while(curr_part){
    new_part = (particle_type *)calloc(1,sizeof(particle_type));
    if( !new_part )
      fatal("Can't allocate a particle.\n");

    new_part->next = particle_surf->bound_particles;
    particle_surf->bound_particles = new_part;

    new_part->loc.x = curr_part->loc.x + (float)my_drand((float).05);
    new_part->loc.y = curr_part->loc.y + (float)my_drand((float).05);
    new_part->loc.z = curr_part->loc.z + (float)my_drand((float).05);

    curr_part->rad *= .25;
    new_part->rad = curr_part->rad;
    curr_part = curr_part->next;
  }
}



/****************************************************************************
 *
 *                   Procedure remove_outlying_particles
 *
 * Arguments: particle_surf: pointer to particle_surface_type.
 *
 * Returns: none
 *
 * Action:  frees all the bound particles in 'particle_surf that have
 *   managed to crawl off of the surface.
 *
 ****************************************************************************/
void remove_outlying_particles(particle_surf)
  particle_surface_type *particle_surf;
{
  particle_type *curr_part,*next_part,*last_part;

  if( !particle_surf )
    FATAL_BUG("Bad particle_surf in remove_outlying_particles.");

  curr_part = particle_surf->bound_particles;

  while(curr_part){
    next_part = curr_part->next;
    if( fabs(fabs(curr_part->val) - particle_surf->surface_value) >
       10*particle_surf->surface_tolerance ){

      if( curr_part == particle_surf->bound_particles ){
        particle_surf->bound_particles = next_part;
        last_part = particle_surf->bound_particles;
      }else{
        last_part->next = next_part;
        last_part = curr_part;
      }
      free(curr_part);
      particle_surf->num_particles--;
      curr_part = next_part;
    }
    else{
      last_part = curr_part;
      curr_part = next_part;
    }
  }
}



/****************************************************************************
 *
 *                   Procedure remove_free_particles
 *
 * Arguments: particle_surf: pointer to particle_surface_type.
 *
 * Returns: none
 *
 * Action:  frees all the non surface bound particles in 'particle_surf
 *
 ****************************************************************************/
void remove_free_particles(particle_surf)
  particle_surface_type *particle_surf;
{
  particle_type *curr_part,*next_part;

  if( !particle_surf )
    FATAL_BUG("Bad particle_surf in remove_free_particles.");

  curr_part = particle_surf->free_particles;
  particle_surf->free_particles = 0;

  while(curr_part){
    next_part = curr_part->next;
    particle_surf->num_particles--;
    free(curr_part);
    curr_part = next_part;
  }
}


/****************************************************************************
 *
 *                   Procedure add_free_particles
 *
 * Arguments: particle_surf: pointer to particle_surface_type.
 *               num_to_add: integer
 *
 * Returns: none
 *
 * Action:  adds 'num_to_add free particles to 'particle_surf
 *
 ****************************************************************************/
void add_free_particles(particle_surf,num_to_add)
  particle_surface_type *particle_surf;
  int num_to_add;
{
  int i;
  particle_type *temp_part;
  MO_center_list_type *center;
  int num_at_this_center;

  if( !particle_surf )
    FATAL_BUG("Bad particle_surf in add_free_particles.");
  particle_surf->num_particles += num_to_add;

  /* start at random locations surrounding the centers */
  num_at_this_center = 0;
  center = particle_surf->MO_centers;
  for( i=0;i<num_to_add;i++ ) {
    temp_part = (particle_type *)calloc(1,sizeof(particle_type));
    if( !temp_part )
      fatal("Can't allocate a particle.\n");

    temp_part->next = particle_surf->free_particles;
    particle_surf->free_particles = temp_part;

    temp_part->rad = 0.3;

    temp_part->loc.x = center->loc->x + my_drand(1.0);
    temp_part->loc.y = center->loc->y + my_drand(1.0);
    temp_part->loc.z = center->loc->z + my_drand(1.0);


    num_at_this_center++;
    if( num_at_this_center == num_to_add / particle_surf->num_centers ){
      num_at_this_center = 0;
#if 0
      if( center->next )center = center->next;
#endif
    }
  }
}



/****************************************************************************
 *
 *                   Procedure calc_bound_particle_forces
 *
 * Arguments: particle_surf: pointer to particle_surface_type.
 *
 * Returns: none
 *
 * Action:  Calculates the forces on bound particles.
 *  This is done by figuring out the net force on each particle
 *   arising from the repulsion of other particles, then projecting
 *   out the component of that force which lies along the gradient
 *   of the iso-surface at that point (the gradient of the iso-surface
 *   is a normal to the surface).  At the end a small spring force
 *   is added to drag the particle back towards the surface (this is to
 *   correct for numerical error in the projection).
 *
 ****************************************************************************/
void calc_bound_particle_forces(particle_surf)
  particle_surface_type *particle_surf;
{
  static point_type *part_force_accum=0;
  static int num_particles=0;

  float dist_from_isosurf;
  float dist,dist_squared,cut_off,cut_off2;
  MO_info_type MO_info;
  point_type dist_vect,grad;
  particle_type *part1,*part2;
  int which_part1,which_part2;
  float repulse_strength,rad_squared;
  float rad_squared2,force_mag;

  /* primitive error checking */
  if( !particle_surf || !particle_surf->MO_centers ){
    FATAL_BUG("Bogus particle_surf passed to calc_bound_particle_forces.");
  }

  /* if there are no bound particles, just return now */
  if( !particle_surf->bound_particles )  return;


  /* see if we need more space for the part_force_accum array */
  if( !part_force_accum || num_particles != particle_surf->num_particles ){
    if( part_force_accum ) free(part_force_accum);

    num_particles = particle_surf->num_particles;
    part_force_accum =
      (point_type *)calloc(num_particles,sizeof(point_type));
    if(!part_force_accum) fatal("Can't get memory for part_force_accum.");
  }

  /* zero out the force accumulation array */
  bzero(part_force_accum,num_particles*sizeof(point_type));

  /*******

    first accumulate the forces, this involves two trips
    through the list of particles

  *******/
  part1 = particle_surf->bound_particles;
  which_part1 = 0;
  while(part1){
    /********

      first determine the distance cut-off for the repulsion
      we'll make this 3*'particle_surf->particle_rad, because by
      that point the repulsion is down way low

    ********/
    rad_squared = part1->rad*part1->rad;
    cut_off = 9.0*rad_squared;

    part2 = part1->next;
    which_part2 = which_part1+1;
    while(part2){
/*      if( part1->phase == part2->phase ){*/
        rad_squared2 = part2->rad*part2->rad;
        cut_off2 = 9.0 * rad_squared2;

        dist_vect.x = part1->loc.x - part2->loc.x;
        dist_vect.y = part1->loc.y - part2->loc.y;
        dist_vect.z = part1->loc.z - part2->loc.z;

        /********

          calculate the distance between the two particles.

          *********/
        dist_squared = dist_vect.x*dist_vect.x + dist_vect.y*dist_vect.y +
          dist_vect.z*dist_vect.z;

        /* normalize the inter-particle vector */
        dist = sqrt(dist_squared);
        dist_vect.x /= dist;
        dist_vect.y /= dist;
        dist_vect.z /= dist;

        if(dist_squared < cut_off ){
          /* add the forces into the accum array */
          repulse_strength = exp(-dist_squared/(2.0*rad_squared));
          part_force_accum[which_part1].x += dist_vect.x*repulse_strength;
          part_force_accum[which_part1].y += dist_vect.y*repulse_strength;
          part_force_accum[which_part1].z += dist_vect.z*repulse_strength;
        }

        /* now do the same thing for the other particle */
        if( dist_squared < cut_off2 ){
          repulse_strength = exp(-dist_squared/(2.0*rad_squared2));
          part_force_accum[which_part2].x -= dist_vect.x*repulse_strength;
          part_force_accum[which_part2].y -= dist_vect.y*repulse_strength;
          part_force_accum[which_part2].z -= dist_vect.z*repulse_strength;
        }
/*      }         */
      part2 = part2->next;
      which_part2++;
    }
    part1 = part1->next;
    which_part1++;
  }

  /*****

    Now step through the list and project out the components of the
    forces which point off of the surface.

  ******/
  part1 = particle_surf->bound_particles;
  which_part1 = 0;

  while(part1){

    /* find the value and gradient of the wave func */
    bzero((char *)&MO_info,sizeof(MO_info_type));
    MO_info.loc.x = part1->loc.x;
    MO_info.loc.y = part1->loc.y;
    MO_info.loc.z = part1->loc.z;
    calc_MO_value(&MO_info,particle_surf->MO_centers);

    /* normalize the gradient */
    grad.x = MO_info.grad.x;
    grad.y = MO_info.grad.y;
    grad.z = MO_info.grad.z;
    dist_squared = grad.x*grad.x + grad.y*grad.y + grad.z*grad.z;
    dist = sqrt(dist_squared);

    grad.x /= dist;
    grad.y /= dist;
    grad.z /= dist;


    /* okay, do the projection and subtraction */
    part_force_accum[which_part1].x -=
      part_force_accum[which_part1].x * grad.x;
    part_force_accum[which_part1].y -=
      part_force_accum[which_part1].y * grad.y;
    part_force_accum[which_part1].z -=
      part_force_accum[which_part1].z * grad.z;

    /* update the radius of the particle */
    force_mag =
      part_force_accum[which_part1].x*part_force_accum[which_part1].x +
        part_force_accum[which_part1].y*part_force_accum[which_part1].y +
          part_force_accum[which_part1].z*part_force_accum[which_part1].z;

/*printf("%lf\n",force_mag);*/
    if( force_mag > 1.00 ){
      part1->rad *= .85;
    } else if( force_mag < .5 ){
      part1->rad *= 1.20;
      if( part1->rad > .35 && particle_surf->okay_to_split ){
#if 0
         && fabs(fabs(part1->val) - particle_surf->surface_value) <
         100.0*particle_surf->surface_tolerance){
#endif
        printf(".\n");
        split_a_particle(particle_surf,part1);
      }
    }



    /* add in the projected force */
    part1->loc.x +=
      particle_surf->bound_step_size*part_force_accum[which_part1].x;
    part1->loc.y +=
      particle_surf->bound_step_size*part_force_accum[which_part1].y;
    part1->loc.z +=
      particle_surf->bound_step_size*part_force_accum[which_part1].z;


    /* find the value and gradient of the wave func */
    MO_info.loc.x = part1->loc.x;
    MO_info.loc.y = part1->loc.y;
    MO_info.loc.z = part1->loc.z;
    calc_MO_value(&MO_info,particle_surf->MO_centers);

    /* check to see if we need to be pulled back onto the surface */
    dist_from_isosurf = fabs(MO_info.val) - particle_surf->surface_value;
    if( fabs(dist_from_isosurf) > particle_surf->surface_tolerance ){
      if( MO_info.val > 0.0 ){
        part1->loc.x += particle_surf->free_step_size*dist_from_isosurf*
          MO_info.grad.x;
        part1->loc.y += particle_surf->free_step_size*dist_from_isosurf*
          MO_info.grad.y;
        part1->loc.z += particle_surf->free_step_size*dist_from_isosurf*
          MO_info.grad.z;
      } else{
        part1->loc.x -= particle_surf->free_step_size*dist_from_isosurf*
          MO_info.grad.x;
        part1->loc.y -= particle_surf->free_step_size*dist_from_isosurf*
          MO_info.grad.y;
        part1->loc.z -= particle_surf->free_step_size*dist_from_isosurf*
          MO_info.grad.z;
      }
    }
    part1->val = MO_info.val;
    if( MO_info.val >= 0.0 ) part1->phase = 1;
    else part1->phase = 0;

    part1 = part1->next;
    which_part1++;
  } /* end of loop over particles */
/*printf("\n\n",force_mag);*/
}



/****************************************************************************
 *
 *                   Procedure calc_free_particle_forces
 *
 * Arguments: particle_surf: pointer to particle_surface_type.
 *
 * Returns: none
 *
 * Action:  Calculates the forces on free particles.  If the particle
 *   has not yet made it onto the desired isosurface, there there will be
 *   a force pulling it towards the surface.
 *
 ****************************************************************************/
void calc_free_particle_forces(particle_surf)
  particle_surface_type *particle_surf;
{
  int i;

  float dist_from_isosurf;
  particle_type *temp_part,*temp_part2;
  particle_surface_type *particle_surface;
  MO_center_list_type *center;
  MO_info_type MO_info;
  AO_list_type *AO;
  point_type loc;
  point_type ang_grad,rad_grad;
  point_type ang_grad_accum,rad_grad_accum;


  /* primitive error checking */
  if( !particle_surf || !particle_surf->MO_centers ){
    FATAL_BUG("Bogus particle_surf passed to calc_free_article_forces.");
  }


  /* step through the list of free particles */
  temp_part = particle_surf->free_particles;
  while(temp_part){

    /* find the value and gradient of the wave func */
    MO_info.loc.x = temp_part->loc.x;
    MO_info.loc.y = temp_part->loc.y;
    MO_info.loc.z = temp_part->loc.z;
    calc_MO_value(&MO_info,particle_surf->MO_centers);

    /******

      if we're not on the desired iso-surface yet, then take a step
      along the gradient direction

    ******/
    dist_from_isosurf = fabs(MO_info.val) - particle_surf->surface_value;

    /*****

      before stepping, check to see if the last step put us on the
      surface

    *****/
    if( fabs(dist_from_isosurf) <= particle_surf->surface_tolerance){
      /******

        we're there, move the particle to the bound list and
        set its phase

      *******/
      temp_part->on_surf = 1;
      if( MO_info.val > 0.0 ) temp_part->phase = 1;
      else temp_part->phase = 0;

      /* find the particle in the free_particle list */
      if( particle_surf->free_particles == temp_part ){
        /* this is the easiest case */
        particle_surf->free_particles = temp_part->next;
        /* now move the particle */
        temp_part->next = particle_surf->bound_particles;
        particle_surf->bound_particles = temp_part;
        temp_part = particle_surf->free_particles;
      } else{
        temp_part2 = particle_surf->free_particles;
        while(temp_part2 && temp_part2->next != temp_part )
          temp_part2 = temp_part2->next;

        if( !temp_part2 )
          FATAL_BUG("Can't find a particle to move it.");

        temp_part2->next = temp_part->next;
        /* now move the particle */
        temp_part->next = particle_surf->bound_particles;
        particle_surf->bound_particles = temp_part;
        temp_part = temp_part2->next;
      }

    } else{
      /* we're not there yet, take a step */
      if( MO_info.val > 0.0 ){
        temp_part->loc.x += particle_surf->free_step_size*dist_from_isosurf*
          MO_info.grad.x;
        temp_part->loc.y += particle_surf->free_step_size*dist_from_isosurf*
          MO_info.grad.y;
        temp_part->loc.z += particle_surf->free_step_size*dist_from_isosurf*
          MO_info.grad.z;
      } else{
        temp_part->loc.x -= particle_surf->free_step_size*dist_from_isosurf*
          MO_info.grad.x;
        temp_part->loc.y -= particle_surf->free_step_size*dist_from_isosurf*
          MO_info.grad.y;
        temp_part->loc.z -= particle_surf->free_step_size*dist_from_isosurf*
          MO_info.grad.z;
      }

      temp_part = temp_part->next;
    }
  } /* end of loop over particles */
}

