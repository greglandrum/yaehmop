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

/****************************************************************************
 *
 *     this file contains everything needed to deal with systems
 *      containing geom_frags
 *
 *
 *  created:  greg landrum  August 1996
 *
 *****************************************************************************/
#include "bind.h"




/****************************************************************************
 *
 *                   Procedure insert_geom_frags
 *
 * Arguments: cell: pointer to cell_type
 *
 *
 * Returns: none
 *
 * Action:  Loops through the geom_frags that are part of 'cell and
 *        inserts them into the geometry.
 *
 *    Each geom_frag is independant and it keeps track of which atom
 *     it replaces.  When the insertion is done, the atom to be
 *     replaced (which logically should be a dummy) is used to
 *     determine the translation of the geom_frag.  No orientation
 *     changes of the geom_frags are handled.  All changes in
 *     orientation must be handled *within* the definition of each fragment.
 *
 *  NOTE:  everything must be in CARTESIAN COORDINATES by the time
 *     this function is called.
 *  NOTE 2:  If cell->atoms has not been reallocated to hold
 *     the fragment atoms, this will be done here.
 *  NOTE 3:  The lattice vectors (if any) are moved to the end of
 *     the array of atoms within this function.
 *
 ****************************************************************************/
void insert_geom_frags(cell_type *cell)
{
  geom_frag_type *geom_frag;
  point_type loc;
  int i,j;
  int next_atom,num_atoms;

  if( !cell->geom_frags ) FATAL_BUG("insert_geom_frags called with no frags.");

  /* count the total number of atoms */
  num_atoms = cell->num_raw_atoms;
  geom_frag = cell->geom_frags;
  while(geom_frag){
    num_atoms += geom_frag->num_atoms;
    geom_frag = geom_frag->next;
  }

  /* check to see if we need to reallocate cell->atoms */
  if( cell->num_atoms != num_atoms ){
    cell->atoms = (atom_type *)my_realloc((int *)cell->atoms,(num_atoms+cell->dim)*sizeof(atom_type));
    if( !cell->atoms ) fatal("Can't realloc cell->atoms in insert_geom_frags.");


    /* move the lattice vects (if any) to the end */
    for(i=0;i<cell->dim;i++){
      bcopy(&cell->atoms[cell->num_atoms+i],
            &cell->atoms[num_atoms+i], sizeof(atom_type));
      bzero(&cell->atoms[cell->num_atoms+i],sizeof(atom_type));
    }
    cell->num_atoms = num_atoms;
  }

  geom_frag = cell->geom_frags;
  next_atom = cell->num_raw_atoms;
  while(geom_frag){
    bcopy(&cell->atoms[geom_frag->which].loc,&loc,sizeof(point_type));

fprintf(stderr,"Inserting a geom_frag at atom %d\n",geom_frag->which+1);
    /* copy the atoms in, then translate them */
    bcopy(geom_frag->atoms,&(cell->atoms[next_atom]),
          geom_frag->num_atoms*sizeof(atom_type));
    for(i=0;i<geom_frag->num_atoms;i++){
      cell->atoms[next_atom].loc.x += loc.x;
      cell->atoms[next_atom].loc.y += loc.y;
      cell->atoms[next_atom].loc.z += loc.z;
      next_atom++;
      /* some error checking */
      if( next_atom > cell->num_atoms )
        FATAL_BUG("insert_geom_frags ran out of atom space");
    }
    geom_frag = geom_frag->next;
  }
}






/****************************************************************************
 *
 *                   Procedure process_geom_frags
 *
 * Arguments: cell: pointer to cell_type
 *
 *
 * Returns: none
 *
 * Action:  Loops through the geom_frags that are part of 'cell and
 *        makes certain that their geometry is in cartesian coordinates.
 *        The fragments are then inserted into the list of atoms.
 *
 ****************************************************************************/
void process_geom_frags(cell_type *cell)
{
  geom_frag_type *geom_frag;

  if( !cell->geom_frags )
    FATAL_BUG("process_geom_frags called without geom_frags.");

  geom_frag = cell->geom_frags;
  while(geom_frag){
    if( geom_frag->using_Z_mat )
      eval_Zmat_locs(geom_frag->atoms,geom_frag->num_atoms,0,0);

    geom_frag = geom_frag->next;
  }

  /* now insert the geom frags */
  insert_geom_frags(cell);
}
