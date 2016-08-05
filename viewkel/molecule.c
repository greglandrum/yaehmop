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
     modifications to determine_connections:
      cleaned up memory allocation, extra bond insertion at beginning
      added persistence of custom bonds
   30.03.98 gL:
     modifications to determine_connections:
      patched up at least some of the horrible bugs the aforementioned
      changes wrought with animations.  what *was* I thinking?
   06.04.98 gL:
     more mods to the bug farm that is determine_connections:
      sooprise sooprise sooprise, that modification was still completely
      screwed up.  This time it worked for animations, but not for crystals.
      things seem to be okay now.   It still gets broken sometimes
      when crystals are grown with custom bonds, but I haven't tracked
      that little taste of evil down yet... at least it doesn't crash.
      I suspect it has something to do with atoms being renumbered.
      This probably needs to be changed.
   02.05.98 gL:
     crashing bug due to new_molecule returning too soon on failure
      fixed
   19.05.98 gL:
     YA boundary condition fixed in determine_connections.
   30.05.98 gL:
     determine_connections is still screwed up when custom lines
     are turned on.  Work around is simple: don't customize lines until
     after you're done finding them.
     Will the bugs in this function never end?  A smart
     boy would just give up and rewrite the damn thing entirely, then
     mabye it would be easier to write the code to allow individual
     lines to be added.
   18.08.98 gL:
     modified center_molecule so that it translates coordination polyhedra
     as well as atoms.... this was a dumb oversight.
   07.09.98 gL:
     determine_mol_bounds adds a bit to the bounding box for
     linear molecules, so now it is (hopefully) possible to select them.
   08.09.98 gL:
     support for colored and/or shaded atoms
   24.09.98 gL:
     changed definition of compare_atoms to remove a compiler warning
   26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)
   29.01.99 gL:
     modifications to deal with drawing the lattice more properly
   21.05.1999 gL:
     some of the problems in determine_connectors with respect to
      custom bonds seem to have been due to a stupid oversight on
      my part.  This appears to have been remedied.  In theory
      that stuff should all work now.  Also, in theory, I've now
      grokked enough of the code (again) to be able to add bonds.
      We'll see if that actually works.

      Frighteningly, that actually seemed to work.  So, now, at last,
      the oft-requested add-a-bond(TM) feature seems to function.
      Like I won't find a crashing bug two days after the release,
      but for now I get to WOO HOO!

      Also: the revision history is now Y2K compliant.

   30.06.1999 gL:
     added support for binary saves and reads of molecules.

***/

/********

  this has got the stuff for dealing with molecules

*********/
#include "viewkel.h"

#define ATOM_RAD 20

#ifndef PARM_FILE
#define PARM_FILE "new_atomic_parms.dat"
#endif

/* only needed locally */
typedef struct{
  int num;
  float dist;
} atom_dist_type;


/* to compare atom_dist_type structures in qsort */
int comp_atom_dists(const void *d1,const void *d2)
{
  if(((atom_dist_type *)d1)->dist > ((atom_dist_type *)d2)->dist){
    return(1);
  } else if(((atom_dist_type *)d1)->dist < ((atom_dist_type *)d2)->dist){
    return(-1);
  } else{
    return(0);
  }
}



/****************************************************************************
 *
 *                   Procedure show_coord_env
 *
 * Arguments: num_selected: int
 *                     obj: pointer to object_type
 *
 * Returns: none
 *
 * Action: Hides all atoms other than those within the coordination
 *    environment of the atom selected
 *
 ***************************************************************************/
void show_coord_env(int num_selected,object_type *obj)
{
  int i;
  molec_type *molec;
  atom_type *atoms;
  atom_dist_type *atom_dists;
  point_type *loc1,*loc2;
  int active_atom;
  float curr_dist,target_dist;
  int num_found;

  if(!num_selected || num_selected > 1 )return;

  if(obj->prim->molec){
    molec = obj->prim->molec;
    if(molec->num_frames > 1)
      atoms = &(molec->atoms[molec->current_frame*molec->num_atoms]);
    else atoms = molec->atoms;
  }else if( obj->prim->MO_surf && obj->prim->MO_surf->molec ){
    molec = obj->prim->MO_surf->molec;
    atoms = molec->atoms;
  } else{
    return;
  }

  /* find the selected atom */
  num_found = 0;
  for(i=0;i<molec->num_atoms && num_found < num_selected;i++){
    if(atoms[i].is_selected){
      active_atom = i;
      num_found++;
    }
  }

  /* get the target distance from the user */
  target_dist = 3.5;
  readfloatparm("coordination cut off",&target_dist);

  atom_dists = (atom_dist_type *)D_CALLOC(molec->num_atoms,sizeof(atom_dist_type));
  if(!atom_dists)fatal("Can't get memory for atom_dists");

  /* loop over all the atoms, hiding the ones outside the cutoff */
  num_found = 0;
  loc1 = &atoms[active_atom].loc;
  for(i=0;i<molec->num_atoms;i++){
    if( i != active_atom && molec->atoms[i].type[0] != '&' ){
      loc2 = &atoms[i].loc;
      curr_dist = V3DistanceBetween2Points(loc1,loc2);
      if( curr_dist > target_dist ){
        atoms[i].exclude = 1;
      } else{
        atom_dists[num_found].num = i;
        atom_dists[num_found].dist = curr_dist;
        num_found++;
      }
    }
  }
  /* sort the list of atom distances */
  qsort((void *)atom_dists,num_found,sizeof(atom_dist_type),comp_atom_dists);

  printf("Atoms within %6.4f Angstroms of: %s(%d)\n",target_dist,
        atoms[active_atom].type,atoms[active_atom].num);

  for(i=0;i<num_found;i++){
    printf("\t%s(%d): %6.4f\n",atoms[atom_dists[i].num].type,
           atoms[atom_dists[i].num].num, atom_dists[i].dist);
  }
}



void add_a_bond(molec_type *molec,int end1,int end2,float length){
  line_type *templines;
  int *templinesto;
  int foo;

  templines = molec->lines;
  foo = molec->num_lines[0];
  molec->lines = (line_type *)D_CALLOC(foo+1,sizeof(line_type));
  if(!molec->lines) fatal("can't allocated molec->lines");
  if( templines ){
    memcpy(molec->lines,templines,foo*sizeof(line_type));
    D_FREE(templines);
  }

  molec->lines[foo].end1 = end1;
  molec->lines[foo].end2 = end2;
  molec->lines[foo].type = 1;
  molec->lines[foo].length = length;

  foo = molec->atoms[end1].num_lines_out;
  templinesto = molec->atoms[end1].linesto;
  molec->atoms[end1].linesto=(int *)D_CALLOC(foo+1,sizeof(int));
  if( !molec->atoms[end1].linesto) fatal("can't allocate linesout 1");
  if(templinesto){
    memcpy(molec->atoms[end1].linesto,templinesto,foo*sizeof(int));
    D_FREE(templinesto);
  }
  molec->atoms[end1].linesto[foo] = end2;
  molec->atoms[end1].num_lines_out++;

  foo = molec->atoms[end2].num_lines_out;
  templinesto = molec->atoms[end2].linesto;
  molec->atoms[end2].linesto=(int *)D_CALLOC(foo+1,sizeof(int));
  if( !molec->atoms[end2].linesto) fatal("can't allocate linesout 1");
  if(templinesto){
    memcpy(molec->atoms[end2].linesto,templinesto,foo*sizeof(int));
    D_FREE(templinesto);
  }
  molec->atoms[end2].linesto[foo] = end1;
  molec->atoms[end2].num_lines_out++;

  molec->num_lines[0] += 1;
}

void add_bonds(molec_type *molec,char *type1,char *type2,float target){
  int i,j;
  float length;

  for(i=0;i<molec->num_atoms;i++){
    for(j=0;j<molec->num_atoms;j++){
      if( j != i ){
        if(!strcmp(molec->atoms[i].type,type1) &&
           !strcmp(molec->atoms[j].type,type2)){
          length = V3DistanceBetween2Points(&molec->atoms[i].loc,
                                            &molec->atoms[j].loc);
          if(fabs(length-target)<0.0001){
            add_a_bond(molec,i,j,length);
          }
        }
      }
    }
  }
  return;
}


/****************************************************************************
 *
 *                   Procedure adjust_color
 *
 * Arguments: num_selected: int
 *                     obj: pointer to object_type
 *          shade_or_color: int
 *
 *
 * Returns: none
 *
 * Action:  allows the user to adjust shading and/or color parameters for
 *   the molecule in 'obj.  If there is no molecule there, we
 *   just report a bug and return.
 *
 *   only a single atom should be selected
 *
 ****************************************************************************/
void adjust_color(int num_selected,object_type *obj,int shade_or_color)
{
  char *stringarr[3];
  char instring[80];
  atom_type *selected_atom;
  molec_type *molec;
  atom_type *atoms;
  int num_found;
  float oldshade,oldcolor[3];
  int i,j;
  char changed;

  if(!num_selected || num_selected > 1 )return;
  if(obj->prim->molec){
    molec = obj->prim->molec;
    if(molec->num_frames > 1)
      atoms = &(molec->atoms[molec->current_frame*molec->num_atoms]);
    else atoms = molec->atoms;
  }else if( obj->prim->MO_surf && obj->prim->MO_surf->molec ){
    molec = obj->prim->MO_surf->molec;
    atoms = molec->atoms;
  } else{
    return;
  }

  num_found = 0;
  for(i=0;i<molec->num_atoms && num_found < num_selected;i++){
    if(atoms[i].is_selected){
      /******

        insert a pointer to the atom into the selected_atoms
        array.  Put it in the slot corresponding to its selection
        order.

        ******/
      selected_atom = &(atoms[i]);
      num_found++;
    }
  }
  if( num_found != num_selected )
    FATAL_BUG("Can't find the selected atom!");

  printf("Altering atom: %s(%d)\n",selected_atom->type,
         selected_atom->num);
  changed = 0;

  switch(shade_or_color){
  case ATOM_SHADE_FILL:
    oldshade = selected_atom->atom_shade;
    readfloatparm("atom shading (between 0 and 1)",&selected_atom->atom_shade);
    if( selected_atom->atom_shade < 0 ) selected_atom->atom_shade = 0;
    if( selected_atom->atom_shade > 1 ) selected_atom->atom_shade = 1;
    if( selected_atom->atom_shade != oldshade ) changed = 1;
    break;
  case ATOM_COLOR_FILL:
    for(i=0;i<3;i++) oldcolor[i] = selected_atom->atom_color[i];
    sprintf(instring,"%3.3lf %3.3lf %3.3lf",oldcolor[0],oldcolor[1],
            oldcolor[2]);
    stringarr[0] = instring;
    readstringparm("atom color (RGB triplet)",stringarr);
    sscanf(instring,"%lf %lf %lf",&selected_atom->atom_color[0],
           &selected_atom->atom_color[1],&selected_atom->atom_color[2]);
    if( selected_atom->atom_color[0] < 0 ) selected_atom->atom_color[0] = 0;
    if( selected_atom->atom_color[0] > 1 ) selected_atom->atom_color[0] = 1;
    if( selected_atom->atom_color[1] < 0 ) selected_atom->atom_color[1] = 0;
    if( selected_atom->atom_color[1] > 1 ) selected_atom->atom_color[1] = 1;
    if( selected_atom->atom_color[2] < 0 ) selected_atom->atom_color[2] = 0;
    if( selected_atom->atom_color[2] > 1 ) selected_atom->atom_color[2] = 1;

    changed = 1;
    break;
  }

  if(changed){
#ifdef X_GRAPHICS
    refresh_all_colormaps = 1;
#endif
    if(!selected_atom->custom){
      selected_atom->outlines_on = molec->outlines_on;
      selected_atom->shading_on = molec->shading_on;
      selected_atom->crosses_on = molec->crosses_on;
    }

    selected_atom->custom = 1;
    strcpy(instring,"no");
    readcharparm("propagate change to all similar atoms?",instring);
    if(instring[0] == 'Y' || instring[0] == 'y' ){
      for(i=0;i<molec->num_atoms*molec->num_frames;i++){
        if(&molec->atoms[i] != selected_atom &&
           !strcmp(molec->atoms[i].type,selected_atom->type)){
          molec->atoms[i].custom = 1;
          molec->atoms[i].outlines_on = selected_atom->outlines_on;
          molec->atoms[i].shading_on =  selected_atom->shading_on;
          molec->atoms[i].crosses_on =  selected_atom->crosses_on;

          switch(shade_or_color){
          case ATOM_SHADE_FILL:
            molec->atoms[i].atom_shade = selected_atom->atom_shade;
            break;
          case ATOM_COLOR_FILL:
            for(j=0;j<3;j++)
              molec->atoms[i].atom_color[j] =
                selected_atom->atom_color[j];
            break;
          }
        }
      }
    }
  } else{
    printf("Okay, you didn't change anything...\n");
  }

}

/****************************************************************************
 *
 *                   Procedure adjust_style
 *
 * Arguments: num_selected: int
 *                     obj: pointer to object_type
 *
 *
 * Returns: none
 *
 * Action:  allows the user to adjust drawing parameters for
 *   the molecule in 'obj.  If there is no molecule there, we
 *   just report a bug and return.
 *
 *   if a single atom is selected, the drawing parms for that atom
 *   will be altered, if two atoms are selected, the parms for
 *   a bond will be altered.  otherwise we just return.
 *
 ****************************************************************************/
void adjust_style(int num_selected,object_type *obj)
{
  char instring[80];
  atom_type **selected_atoms;
  molec_type *molec;
  atom_type *atoms;
  float newfloat,length;
  line_type *the_line;
  int newint,num_found;
  char add_similar;
  int num_lines;
  int i;
  char changed,match_distance;

  if(!num_selected || num_selected > 2 )return;
  if(obj->prim->molec){
    molec = obj->prim->molec;
    if(molec->num_frames > 1)
      atoms = &(molec->atoms[molec->current_frame*molec->num_atoms]);
    else atoms = molec->atoms;
  }else if( obj->prim->MO_surf && obj->prim->MO_surf->molec ){
    molec = obj->prim->MO_surf->molec;
    atoms = molec->atoms;
  } else{
    return;
  }

  selected_atoms = (atom_type **)D_CALLOC(num_selected,
                                        sizeof(atom_type *));
  if(!selected_atoms)fatal("Can't get memory for selected atoms array");

  num_found = 0;
  for(i=0;i<molec->num_atoms && num_found < num_selected;i++){
    if(atoms[i].is_selected){
      /******

        insert a pointer to the atom into the selected_atoms
        array.  Put it in the slot corresponding to its selection
        order.

        ******/
      selected_atoms[atoms[i].is_selected-1] = &(atoms[i]);
      num_found++;
    }
  }
  if( num_found != num_selected )
    FATAL_BUG("Can't find all the selected atoms!");

  switch(num_selected){

  case 1:
    printf("Altering atom: %s(%d)\n",selected_atoms[0]->type,
           selected_atoms[0]->num);
    changed = 0;
    if(!selected_atoms[0]->custom){
      selected_atoms[0]->outlines_on = molec->outlines_on;
      selected_atoms[0]->shading_on = molec->shading_on;
      selected_atoms[0]->crosses_on = molec->crosses_on;
    }
    newfloat = selected_atoms[0]->rad;
    readfloatparm("radius",&(newfloat));
    if(newfloat!=selected_atoms[0]->rad){
      changed = 1;
      selected_atoms[0]->rad = newfloat;
    }
    newint = selected_atoms[0]->color;
    readintparm("color",&newint);
    if(newint!=selected_atoms[0]->color){
      changed=1;
      selected_atoms[0]->color = (char)newint;
    }
    newint = selected_atoms[0]->shading_on;
    readintparm("shading on",&newint);
    if(newint!=selected_atoms[0]->shading_on){
      changed=1;
      selected_atoms[0]->shading_on = (char)newint;
    }
    newint = selected_atoms[0]->outlines_on;
    readintparm("outlines on",&newint);
    if(newint!=selected_atoms[0]->outlines_on){
      changed=1;
      selected_atoms[0]->outlines_on = (char)newint;
    }
    newint = selected_atoms[0]->crosses_on;
    readintparm("crosses on",&newint);
    if(newint!=selected_atoms[0]->crosses_on){
      changed=1;
      selected_atoms[0]->crosses_on = (char)newint;
    }
    if(changed){
      selected_atoms[0]->custom = 1;
      strcpy(instring,"no");
      readcharparm("propagate change to all similar atoms?",instring);
      if(instring[0] == 'Y' || instring[0] == 'y' ){
        for(i=0;i<molec->num_atoms*molec->num_frames;i++){
          if(&molec->atoms[i] != selected_atoms[0] &&
             !strcmp(molec->atoms[i].type,selected_atoms[0]->type)){
            molec->atoms[i].custom = 1;
            molec->atoms[i].outlines_on =
              selected_atoms[0]->outlines_on;
            molec->atoms[i].shading_on =
              selected_atoms[0]->shading_on;
            molec->atoms[i].crosses_on =
              selected_atoms[0]->crosses_on;
            molec->atoms[i].rad =
              selected_atoms[0]->rad;
            molec->atoms[i].color =
              selected_atoms[0]->color;
          }
        }
      }
    } else{
      printf("Okay, you didn't change anything...\n");
    }
    break;

  case 2:
    /* here we have to check to see if the line already exists */
    the_line = find_the_line(selected_atoms[0]->num,
                             selected_atoms[1]->num,
                             molec);
    if(!the_line){
      length = V3DistanceBetween2Points(&selected_atoms[0]->loc,
                                        &selected_atoms[1]->loc);
      printf("There is no bond between: %s(%d) and %s(%d) (%f A long)\n",
             selected_atoms[0]->type,selected_atoms[0]->num,
             selected_atoms[1]->type,selected_atoms[1]->num,
             length);
      printf("Should I add one?\n");
      strcpy(instring,"no");
      readcharparm("add the bond?",instring);
      if(instring[0] == 'Y' || instring[0] == 'y'){
        strcpy(instring,"no");
        readcharparm("add all similar bonds?",instring);
        if(instring[0] == 'Y' || instring[0] == 'y' ){
          add_similar = 1;
        } else{
          add_similar = 0;
        }
        num_lines = molec->num_lines[molec->current_frame%molec->num_frames];
        if( add_similar ){
          add_bonds(molec,selected_atoms[0]->type,selected_atoms[1]->type,
                    length);
          the_line = find_the_line(selected_atoms[0]->num,selected_atoms[1]->num,
                                   molec);
        }else{
          add_a_bond(molec,selected_atoms[0]->num,selected_atoms[1]->num,
                     length);
          the_line = find_the_line(selected_atoms[0]->num,selected_atoms[1]->num,
                                   molec);
        }

      }
    }

    if( the_line ){
      printf("Altering the bond between: %s(%d) and %s(%d) (%f A long)\n",
             selected_atoms[0]->type,selected_atoms[0]->num,
             selected_atoms[1]->type,selected_atoms[1]->num,
             the_line->length);
      changed = 0;
      if(!the_line->custom){
        the_line->tube = molec->tubes_on;
        the_line->breaking = molec->breaking_lines;
        the_line->drawn = 1;
        the_line->thickness = molec->line_width;
      }

      newint = the_line->drawn;
      readintparm("show bond",&newint);
      if(newint!=the_line->drawn){
        changed=1;
        the_line->drawn = (char)newint;
      }
      newint = the_line->tube;
      readintparm("tube bond",&newint);
      if(newint!=the_line->tube){
        changed=1;
        the_line->tube = (char)newint;
      }
      newint = the_line->breaking;
      readintparm("breaking line",&newint);
      if(newint!=the_line->breaking){
        changed=1;
        the_line->breaking = (char)newint;
      }
      newint = the_line->dashed;
      readintparm("dashed line",&newint);
      if(newint!=the_line->dashed){
        changed=1;
        the_line->dashed = (char)newint;
      }

#ifdef INCLUDE_BOND_VALENCE
      if( the_line->dashed || molec->valence_for_bonds ){
#else  /* } */
      if( the_line->dashed ){
#endif
        newint = the_line->type;
        readintparm("line style",&newint);
        if(newint!=the_line->type){
          changed=1;
          the_line->type = newint;
        }
      }

      newint = the_line->thickness;
      readintparm("line width",&newint);
      if(newint!=the_line->thickness){
        changed=1;
        the_line->thickness = newint;
      }

    } else{
      printf("gotta add it first\n");
      changed = 0;
    }
    if(changed){
      the_line->custom = 1;

      strcpy(instring,"no");
      readcharparm("propagate change to all similar bonds?",instring);
      if(instring[0] == 'Y' || instring[0] == 'y' ){
        strcpy(instring,"yes");
        readcharparm("use distance to match?",instring);
        if(instring[0] == 'Y' || instring[0] == 'y' ) match_distance = 1;
        else match_distance = 0;
        num_lines = molec->num_lines[molec->current_frame%molec->num_frames];

        for(i=0;i<num_lines;i++){
          if((!match_distance ||
              fabs(molec->lines[i].length-the_line->length)<0.0001 &&
             &(molec->lines[i]) != the_line) &&
             ((!strcmp(molec->atoms[molec->lines[i].end1].type,
                        selected_atoms[0]->type) &&
               !strcmp(molec->atoms[molec->lines[i].end2].type,
                        selected_atoms[1]->type)) ||
              (!strcmp(molec->atoms[molec->lines[i].end1].type,
                        selected_atoms[1]->type) &&
               !strcmp(molec->atoms[molec->lines[i].end2].type,
                        selected_atoms[0]->type)))){
            molec->lines[i].drawn = the_line->drawn;
            molec->lines[i].tube = the_line->tube;
            molec->lines[i].breaking = the_line->breaking;
            molec->lines[i].dashed = the_line->dashed;
            molec->lines[i].type = the_line->type;
            molec->lines[i].thickness = the_line->thickness;
            molec->lines[i].custom = 1;
          }
        }
      }
    } else{
      printf("Okay, you didn't change anything...\n");
    }

    break;
  }

  if(selected_atoms)D_FREE(selected_atoms);
}
/* this function is used by qsort to depth sort the atoms */
int zcompare(atom_type *e1,atom_type *e2)
{
  return((int)(100.0*(e1->loc.z - e2->loc.z)));
}



int find_numbered_atom(atom_type *array,int length,int number)
{
  int i;

  for(i=0;i<length;i++){
    if(array[i].num == number) return(i);
  }

  error("Can't find an atom for a connector.");
  return(-1);
}


/****************************************************************************
 *
 *                   Procedure center_molecule
 *
 * Arguments:    num_args: int
 *              molec_ptr: pointer to pointer to char
 *
 * Returns: none
 *
 * Action: centers the molecule pointed to by molec_ptr about the origin.
 *    this is intended to be called from a button window
 *
 ****************************************************************************/
void center_molecule(int num_args,char **molec_ptr)
{
  int i;
  molec_type *molec;
  atom_type *atom;
  point_type center;
  int num_shown;
  molec = (molec_type *)molec_ptr[0];

  center.x = center.y = center.z = 0;

  /* find the center of "mass" of atoms being displayed. */
  num_shown = 0;
  for(i=0;i<molec->num_atoms;i++){
    atom = &molec->atoms[i];
    if( !atom->exclude && (molec->hydrogens_on ||
        atom->type[0] != 'H' || atom->type[1] != 0) &&
        (molec->dummies_on || atom->type[0] != '&') ){
      center.x += atom->loc.x;
      center.y += atom->loc.y;
      center.z += atom->loc.z;
      num_shown++;
    }
  }

  if( num_shown > 0 ){
    center.x /= num_shown;
    center.y /= num_shown;
    center.z /= num_shown;
  }else{
    return;
  }

  /* move the atoms */
  for(i=0;i<molec->num_atoms*molec->num_frames;i++){
    molec->atoms[i].loc.x -= center.x;
    molec->atoms[i].loc.y -= center.y;
    molec->atoms[i].loc.z -= center.z;
  }

  /* if we have coordination polyhedra on, move those too */
  for(i=0;i<molec->num_polyhed_verts;i++){
    molec->polyhed_verts[i].position.x -= center.x;
    molec->polyhed_verts[i].position.y -= center.y;
    molec->polyhed_verts[i].position.z -= center.z;
  }
  for(i=0;i<molec->num_triangles;i++){
    molec->triangles[i].center.x -= center.x;
    molec->triangles[i].center.y -= center.y;
    molec->triangles[i].center.z -= center.z;
  }


  /* if there's a lattice, please move that too */
  if(molec->num_dim){
    for(i=0;i<=molec->num_dim;i++){
      molec->lattice_vect[i].x -= center.x;
      molec->lattice_vect[i].y -= center.y;
      molec->lattice_vect[i].z -= center.z;
    }
    for(i=0;i<8;i++){
      molec->cell_box[i].x -= center.x;
      molec->cell_box[i].y -= center.y;
      molec->cell_box[i].z -= center.z;
    }
  }
  display("Centered the atoms!");
}


/****************************************************************************
 *
 *                   Procedure determine_mol_bounds
 *
 * Arguments: molec: pointer to molecule
 *          bmin,bmax: pointers to point_type
 *
 * Returns: none
 *
 * Action: determines a bounding box for the molecule 'molec
 *
 ****************************************************************************/
void determine_mol_bounds(molec_type *molec,point_type *bmin,
                          point_type *bmax)
{
  int i;

  bmin->x = 1e10;
  bmin->y = 1e10;
  bmin->z = 1e10;
  bmax->x = -1e10;
  bmax->y = -1e10;
  bmax->z = -1e10;

  for(i=0;i<molec->num_atoms;i++){
    if( molec->atoms[i].loc.x < bmin->x )
      bmin->x = molec->atoms[i].loc.x;
    if( molec->atoms[i].loc.y < bmin->y )
      bmin->y = molec->atoms[i].loc.y;
    if( molec->atoms[i].loc.z < bmin->z )
      bmin->z = molec->atoms[i].loc.z;

    if( molec->atoms[i].loc.x > bmax->x )
      bmax->x = molec->atoms[i].loc.x;
    if( molec->atoms[i].loc.y > bmax->y )
      bmax->y = molec->atoms[i].loc.y;
    if( molec->atoms[i].loc.z > bmax->z )
      bmax->z = molec->atoms[i].loc.z;
  }

  /* check for linear molecules */
  if( bmax->x - bmin->x  < 1 ){
    bmin->x -= 2;
    bmax->x += 2;
  }
  if( bmax->y - bmin->y  < 1 ){
    bmin->y -= 2;
    bmax->y += 2;
  }
}


/****************************************************************************
 *
 *                   Function hide_atoms
 *
 * Arguments:
 *
 * Returns: none
 *
 * Action:  allows the user to exclude particular atoms from the
 *  display of the molecule
 *
 ****************************************************************************/
void hide_atoms(int num_args,char **molec_ptr)
{
  char instring[MAX_STR_LEN];
  int *atoms_to_exclude;
  int num_to_exclude;
  int i;

  molec_type *molec;

  display("Look in Xterm");
  molec = (molec_type *)molec_ptr[0];

  printf("Please enter the numbers of the atoms you wish to hide:\n");
  fgets(instring,MAX_STR_LEN,stdin);

  parse_integer_string(instring,&atoms_to_exclude,&num_to_exclude);

  /*******

    now set the exclusion...

    we set it in the raw_MO_centers array in case this is an extended
    system and the user has grown the crystal

  ********/
  for(i=0;i<num_to_exclude;i++){
    molec->atoms[atoms_to_exclude[i]-1].exclude = 1;
  }

  printf("%d Atoms excluded.\n",num_to_exclude);
  if( molec->num_dim >= 1 ){
    printf("If you have grown this crystal, you will need to grow it again\n");
    printf("  to make these changes take effect.  Sorry\n");
  }
}


/****************************************************************************
 *
 *                   Function show_atoms
 *
 * Arguments:
 *
 * Returns: none
 *
 * Action:  allows the user to include particular atoms in the
 *  display of the molecule.  Note that there is no point
 *  in doing this unless they have explicitly hidden some atoms
 *
 ****************************************************************************/
void show_atoms(int num_args,char **molec_ptr)
{
  char instring[MAX_STR_LEN];
  int *atoms_to_exclude;
  int num_to_exclude;
  int i;

  molec_type *molec;

  display("Look in Xterm");
  molec = (molec_type *)molec_ptr[0];

  printf("Please enter the numbers of the atoms you wish to show:\n");
  fgets(instring,MAX_STR_LEN,stdin);

  parse_integer_string(instring,&atoms_to_exclude,&num_to_exclude);

  /*******

    now set the inclusion...

    we set it in the raw_MO_centers array in case this is an extended
    system and the user has grown the crystal

  ********/
  for(i=0;i<num_to_exclude;i++){
    molec->atoms[atoms_to_exclude[i]-1].exclude = 0;
  }

  printf("%d Atoms shown.\n",num_to_exclude);
  if( molec->num_dim >= 1 ){
    printf("If you have grown this crystal, you will need to grow it again\n");
    printf("  to make these changes take effect.  Sorry\n");
  }
}

/****************************************************************************
 *
 * Procedure determine_connections
 *
 * Arguments: molec: pointer to molecule
 *
 * Returns: none
 *
 * Action: This uses the covalent radii of the atoms in the molecule
 *    to determine which ones should be connected by lines ("bonds")
 *    when the molecule is drawn.
 *
 ****************************************************************************/
void determine_connections(molec_type *molec)
{
  int i,j,k,frame;
  atom_type *atoms;
  int num_connectors,max_connectors,connectors_so_far;
  line_type *templines,*linesave;
  int num_custom;
  float dist,rad_sum;
  int end1,end2;
  int num_atoms;
  float R0_val,valence;
  int add_a_bond;
  int *connectors_per_atom;

  frame = 0;
  num_atoms = molec->num_atoms*molec->num_frames;

  connectors_so_far=0;
  max_connectors = 0;
  num_connectors = 0;

  /* clear out all non-custom connectors already present
     (in case we've been here before) */
  atoms = molec->atoms;
  num_custom = 0;
  if( molec->lines ){
    /* self-consistency check */
    if( !molec->num_lines ) FATAL_BUG("molec->num_lines should not be null");

    templines = (line_type *)D_CALLOC(molec->num_lines[0],sizeof(line_type));
    if(!templines){
      fatal("can't allocate templines.");
    }
    num_custom=0;
    for(i=0;i<molec->num_lines[0];i++){
      if( molec->lines[i].custom ){
        memcpy(&templines[num_custom],&molec->lines[i],sizeof(line_type));
        num_custom++;
      }
    }
    D_FREE(molec->lines);
    if( num_custom ){
      molec->lines = (line_type *)D_CALLOC(num_custom,sizeof(line_type));
      if(!molec->lines) fatal("can't allocate storage for molec->lines");
      memcpy(molec->lines,templines,num_custom*sizeof(line_type));
      max_connectors = num_custom;
      connectors_so_far = num_custom;
      num_connectors = connectors_so_far;
    } else{
      molec->lines = 0;
    }
    D_FREE(templines);
  }
  for(i=0;i<num_atoms;i++){
    if( atoms[i].linesto ){
      D_FREE(atoms[i].linesto);
    }
    atoms[i].linesto = 0;
    atoms[i].num_lines_out = 0;
  }
  molec->num_lines[0] = 0;

  /* allocate memory to store the number of connectors needed
     at each atom */
  connectors_per_atom = (int *)D_CALLOC(num_atoms,sizeof(int));
  if(!connectors_per_atom) fatal("can't allocate connectors_per_atom");

  /* just loop over all of the atoms */
  frame = 0;
  for( i=0;i<num_atoms;i++){
    if( i % molec->num_atoms == 0) frame++;
    for(j=i+1;j<frame*molec->num_atoms;j++){
      dist = V3DistanceBetween2Points(&atoms[i].loc,&atoms[j].loc);
      rad_sum = atoms[i].rad + atoms[j].rad;

      /* now check for a bond */
      add_a_bond = 0;
#ifdef INCLUDE_BOND_VALENCE
      if( molec->valence_for_bonds ){
        R0_val = 0.0;
        bond_length_to_bond_valence(&atoms[i],&atoms[j],
                                    dist,&R0_val,&valence);
        if( valence > molec->bond_tol ){
          add_a_bond = 2;
        }
        else if( valence > molec->bond_tol2 ) {
          add_a_bond = 1;
        }
      }
      else{
#endif
        if( dist - rad_sum <= molec->bond_tol*rad_sum ){
          add_a_bond=1;
        }
#ifdef INCLUDE_BOND_VALENCE
      }
#endif
      if(add_a_bond){
        /* make sure we don't already have a bond between these
           atoms */
        for(k=0;k<num_connectors && add_a_bond;k++){
          if( (molec->lines[k].end1 == i && molec->lines[k].end2 == j) ||
              (molec->lines[k].end1 == j && molec->lines[k].end2 == k) ){
            add_a_bond = 0;
          }
        }
      }
      if(add_a_bond){
        num_connectors++;

        /* check to see if we need more memory */
        if( connectors_so_far >= max_connectors || !molec->lines){
          max_connectors += 10;
          templines = (line_type *)D_CALLOC(max_connectors,sizeof(line_type));
          if( !templines ) fatal("Can't get more memory for connectors");

          /* copy in the lines from the current array, if there is one */
          if( molec->lines ){
            memcpy((char *)templines,(char *)molec->lines,
                   connectors_so_far*sizeof(line_type));
            linesave = molec->lines;
            molec->lines = templines;
            /* free up the old array */
            D_FREE(linesave);
          }
          else{
            molec->lines = templines;
          }
        }
        if(!molec->lines) molec->lines = templines;
        molec->lines[connectors_so_far].end1 = i;
        molec->lines[connectors_so_far].end2 = j;
        molec->lines[connectors_so_far].type = add_a_bond;
        molec->lines[connectors_so_far].length = dist;
        connectors_so_far++;
      }
    }
  }

    /* this is a boundary condition we have to catch */
  /*    if(num_connectors != 0) num_connectors++;*/

  /******
      set the flag in each of the endpoints indicating that they
      are involved in a line (and telling which atom they connect to!).
    *******/
  if( num_connectors ){
    /* start by counting the number of connectors per atom */
    for(i=0;i<num_connectors;i++){
      connectors_per_atom[molec->lines[i].end1]++;
      connectors_per_atom[molec->lines[i].end2]++;
    }

    for(i=0;i<num_connectors;i++){
      end1 = molec->lines[i].end1;
      end2 = molec->lines[i].end2;

      /* check to see if we need memory for the linesto array */
      if( !atoms[end1].linesto ){
        atoms[end1].linesto = (int *)D_CALLOC(connectors_per_atom[end1],
                                              sizeof(int));
        if( !atoms[end1].linesto )
          fatal("Can't get atomic connector memory.");
      }
      if( !atoms[end2].linesto ){
        atoms[end2].linesto = (int *)D_CALLOC(connectors_per_atom[end2],
                                              sizeof(int));
        if( !atoms[end2].linesto )
          fatal("Can't get atomic connector memory.");
      }

      /* set the endpoints */
      atoms[end1].linesto[atoms[end1].num_lines_out++] = end2;
      atoms[end2].linesto[atoms[end2].num_lines_out++] = end1;
      if( atoms[end1].num_lines_out > connectors_per_atom[end1]){
        FATAL_BUG("num_lines_out exceeds connectors_per_atom");
      }
      if( atoms[end2].num_lines_out > connectors_per_atom[end2]){
        FATAL_BUG("num_lines_out exceeds connectors_per_atom");
      }
    }
  }

  if( num_connectors > 0 ) molec->num_lines[0] = num_connectors;
  if(connectors_per_atom) D_FREE(connectors_per_atom);

}


/****************************************************************************
 *
 *                   Procedure fill_atomic_info
 *
 * Arguments: num_atoms: integer
 *            atom: a pointer to atom_type
 *
 * Returns: none
 *
 * Action: This routine reads in the information about the atoms in 'atoms from
 *  the parameter file.
 *
 ****************************************************************************/
void fill_atomic_info(int num_atoms,atom_type *atoms)
{
  int i;
  char found;
  char instring[MAX_STR_LEN],foostring[80];
  char *file_name;
  FILE *infile;
  int eof_hit;

  /* open the parameter file */
  file_name = getenv("VIEWKEL_PARM_FILE");
  if( !file_name ){
    file_name = calloc(240,sizeof(char));
    strcpy(file_name,PARM_FILE);
  }

  infile = fopen(file_name,"r");
  if(!infile){
    fprintf(stderr,"Can't open atomic parameter file: %s\n",file_name);
    fatal("Can't open parameter file.");
  }

  /* look for each atom */
  for(i=0;i<num_atoms;i++){
    found = 0;
    strcpy(foostring,atoms[i].type);
    upcase(foostring);
    rewind(infile);
    eof_hit = skipcomments(infile,instring);
    while( !found && eof_hit >= 0 ){
      upcase(instring);
      if(strstr(instring,foostring)){
#ifdef INCLUDE_BOND_VALENCE
        sscanf(instring,"%s %lf %c %lf %lf",foostring,&(atoms[i].rad),
               &(atoms[i].color),
               &atoms[i].ci,&atoms[i].ri);
#else
        sscanf(instring,"%s %lf %c",foostring,&(atoms[i].rad),
               &(atoms[i].color));
#endif
        found = 1;
        atoms[i].color -= '0';
      }
      else{
        eof_hit = skipcomments(infile,instring);
      }
    }
    /* check to see if we actually found the atom */
    if( eof_hit < 0){
      /* nope... */
      fprintf(stderr,"Can't find atom type: %s in parameter file.",
              atoms[i].type);
      atoms[i].color = 0;
      atoms[i].rad = 1.0;
    }
    if( atoms[i].color != 0 ){
      atoms[i].atom_shade = .9 - (float)(NUM_COLORS-atoms[i].color)/(float)NUM_COLORS;
      atoms[i].atom_color[0] = atoms[i].atom_shade;
      atoms[i].atom_color[1] = atoms[i].atom_shade;
      atoms[i].atom_color[2] = atoms[i].atom_shade;
    } else{
      atoms[i].atom_shade = 0.0;
      atoms[i].atom_color[0] = atoms[i].atom_shade;
      atoms[i].atom_color[1] = atoms[i].atom_shade;
      atoms[i].atom_color[2] = atoms[i].atom_shade;
    }
    atoms[i].Gpixel_val[0] = -1;
    atoms[i].Cpixel_val[0] = -1;
  }

#ifdef X_GRAPHICS
    refresh_all_colormaps = 1;
#endif

}


/****************************************************************************
 *
 *                   Procedure read_molecule_data
 *
 * Arguments:   infile: pointer to type FILE
 *              molec: a pointer to molecule_type
 *
 * Returns: none
 *
 * Action: This routine reads in the information for viewing a molecule from
 *    'infile.  The data is stored in 'molec
 *
 *  The file format assumed is that used in the output files from bind.
 *
 ****************************************************************************/
void read_molecule_data(FILE *infile,molec_type *molec)
{
  FILE *test;
  char foostring[80];
  int bin_read_failed;
  char instring[MAX_STR_LEN];

  int i,j;
  char extended_system;
  int eof_hit;
  int num_atoms,num_steps,which_step,num_frames;
  int reading_movie;
  int num_to_skip;
  int num_skipped;
  int lines_read;

  if(!infile)FATAL_BUG("No file passed to read_molecule_data!");
  if(!molec)FATAL_BUG("No molecule passed to read_molecule_data!");

  /* check if there is a binary save file around which we can use */
  sprintf(instring,"%s.BIN",molec->filename);
  test = fopen(instring,"r");
  if( test ) bin_read_failed = read_bin_molecule(molec);
  else bin_read_failed = 1;

  if(!bin_read_failed) return;

  molec->num_atoms = 0;
  num_steps = 0;

  /******

    find and read out the number of atoms, checking for Walsh info along the way

  *******/
  skipcomments(infile,instring);
  upcase(instring);

  if( strstr(instring,"MOVIE-FILE") ) reading_movie = 1;
  else reading_movie = 0;

  while(!strstr(instring,"NUMBER OF ATOMS") && !strstr(instring,"WALSH INFO")){
    if(skipcomments(infile,instring)<0){
      error("End of File hit while looking for geometry data.\n");
      return;
    }
    upcase(instring);
  }
  if( strstr(instring,"NUMBER") ){
    skipcomments(infile,instring);
    sscanf(instring,"%d",&num_atoms);
  }else{
    skipcomments(infile,instring);
    skipcomments(infile,instring);
    sscanf(instring,"%d",&num_steps);
    /* if we hit walsh info, then we still need to find the number of atoms */
    while(!strstr(instring,"NUMBER OF ATOMS") ){
      skipcomments(infile,instring);
      upcase(instring);
    }
    skipcomments(infile,instring);
    sscanf(instring,"%d",&num_atoms);
  }

  if( reading_movie ){
    while(!strstr(instring,"NUMBER OF FRAMES")){
      if(skipcomments(infile,instring)<0){
        error("End of File hit while looking for frame data.\n");
        return;
      }
      upcase(instring);
    }
    skipcomments(infile,instring);
    sscanf(instring,"%d",&num_frames);
  }else{
    num_frames = 1;
  }


  /* get memory for the atoms */
  molec->atoms = (atom_type *)D_CALLOC(num_atoms*num_frames,sizeof(atom_type));
  if( !(molec->atoms) ){
    printf("Can't get memory for the molec atoms.");
    display("Ouch!");
    return;
  }


  which_step = 0;
  while( num_steps != 0  && which_step == 0){
    printf("There are %d walsh steps in the file.  Which would you like to see? ",
           num_steps);
    which_step = 1;
    readintparm("active step",&which_step);
    /*    scanf("%d\n",&which_step);*/
    if( which_step > num_steps || which_step <= 0){
      printf("%d is an invalid choice.  Try again.\n",which_step);
      which_step = 0;
    }
  }

  /* find the atomic positions */
  if( which_step == 0 ){
    which_step = 1;
/*
    skipcomments(infile,instring);
    upcase(instring);
*/
      while(!reading_movie && !strstr(instring,"ATOMIC POSITIONS") &&
            !strstr(instring,"WALSH STEP") &&
            !strstr(instring,"POSITIONS OF ATOMS FROM")){
      skipcomments(infile,instring);
      upcase(instring);
    }
  }else{

    for( i=0; i<which_step; i++ ){
      skipcomments(infile,instring);
      upcase(instring);
      /* find the proper walsh step */
      while(!strstr(instring,"WALSH_STEP")){
        skipcomments(infile,instring);
        upcase(instring);
      }
    }
    /* now find the positions of the atoms */
    while((!strstr(instring,"ATOMIC POSITIONS") ||
           !strstr(instring,"WALSH STEP")) &&
          !strstr(instring,"POSITIONS OF ATOMS FROM")){
      skipcomments(infile,instring);
      upcase(instring);
    }

  }

  if( reading_movie ){
    num_to_skip = 10;
    num_skipped = 0;
    printf("There are a total of %d frames.\n",num_frames);
    readintparm("number of frames to skip",&num_to_skip);
    num_frames = (int)floor((float)num_frames/(float)num_to_skip);
    printf("okay, I'll read %d frames\n",num_frames);
  }

  lines_read = 0;
  for(j=0;j<num_frames;j++){
    if( reading_movie ){
      num_skipped = 1;
      while( num_skipped != num_to_skip ){
        for(i=0;i<num_atoms;i++){
          lines_read++;
          if(skipcomments(infile,instring) < 0){
            printf("Problems with atom %d of frame %d after reading %d lines.\n",
                   i,j,lines_read);
            fatal("EOF hit reading movie!");
          }

        }
        num_skipped++;
      }
    }

    for(i=0;i<num_atoms;i++){
      /* make sure that we didn't hit end of file */
      if(skipcomments(infile,instring) < 0){
        error("Can't read in all the atoms!");
        display("Toooo bad.");
        return;
      }
      if( !reading_movie ){
        sscanf(instring,"%d %s %lf %lf %lf",&(molec->atoms[i].num),
                (molec->atoms[i].type),
               &(molec->atoms[i].loc.x),
               &(molec->atoms[i].loc.y),&(molec->atoms[i].loc.z));
        molec->atoms[i].num--;
      }else{
        sscanf(instring,"%s %lf %lf %lf",
                (molec->atoms[j*num_atoms+i].type),
               &(molec->atoms[j*num_atoms+i].loc.x),
               &(molec->atoms[j*num_atoms+i].loc.y),
               &(molec->atoms[j*num_atoms+i].loc.z));
        molec->atoms[j*num_atoms+i].num = j*num_atoms+i;
      }
    }
  }
  molec->num_atoms = num_atoms;
  molec->num_frames = num_frames;
  molec->current_frame = 0;


  /* fill in the atomic radii and colors */
  fill_atomic_info(num_atoms*num_frames,molec->atoms);

  /******

    go through all the atoms and find which ones should be connected
    by bonds in the output file.

  *******/
  molec->num_lines = (int *)D_CALLOC(molec->num_frames,sizeof(int));
  molec->bond_tol = DEF_BOND_TOL;
  molec->old_bond_tol = DEF_BOND_TOL;
#ifdef INCLUDE_BOND_VALENCE
  molec->bond_tol2 = DEF_BOND_TOL;
  molec->old_bond_tol2 = DEF_BOND_TOL;
#endif
  determine_connections(molec);

  /*****

    check to see if this is an extended system

  *****/
  extended_system = 0;
  eof_hit = skipcomments(infile,instring);
  upcase(instring);
  while( eof_hit >= 0 && !strstr(instring,"DIMENSION") ){
    eof_hit = skipcomments(infile,instring);
    upcase(instring);
  }
  if( eof_hit >= 0 ){
    /* !!! THIS NEEDS TO BE *INSIDE* THE WALSH LOOP */

    /* it is an extended system, so read out the dimensionality now */
    sscanf(instring,"%s %d",foostring,&(molec->num_dim));
    /* skip the next line */
    if(!reading_movie) skipcomments(infile,instring);

    /* read in the lattice vectors */
    molec->lattice_vect[0].x = 0.0;
    molec->lattice_vect[0].y = 0.0;
    molec->lattice_vect[0].z = 0.0;
    for(i=1;i<=molec->num_dim;i++){
      skipcomments(infile,instring);
      sscanf(instring,"%s %lf %lf %lf",foostring,&(molec->lattice_vect[i].x),
             &(molec->lattice_vect[i].y),&(molec->lattice_vect[i].z));
    }
    for(i=0;i<=molec->num_dim;i++){
      molec->orig_lattice[i].x = molec->lattice_vect[i].x;
      molec->orig_lattice[i].y = molec->lattice_vect[i].y;
      molec->orig_lattice[i].z = molec->lattice_vect[i].z;
    }

    molec->num_atoms_in_cell = num_atoms;
  }

  fprintf(stderr,"Done Reading Molecular Data.\n");
  display("GO!");

  /* that's that (I think) */
}

/****************************************************************************
 *
 *                   Procedure new_molecule
 *
 * Arguments: filename: pointer to type char
 *
 * Returns: none
 *
 * Action: does everything to get space for and read in a new molecule
 *
 ****************************************************************************/
void new_molecule(char *filename)
{
  char file_name[80];
  char *theinline;
  FILE *infile;


  /* set up a new object to hold the molecule */
  makenewobject();
  whichobj = head->obj;

  /* now build the molecule primitive */
  whichobj->prim = (prim_type *)D_CALLOC(1,sizeof(prim_type));
  if( !whichobj->prim )fatal("Can't get space for molecule primitive.");
  whichobj->prim->which = MOLECULE;

  whichobj->prim->molec = (molec_type *)D_CALLOC(1,sizeof(molec_type));
  if( !whichobj->prim->molec )
    fatal("Can't get space for molecule.");
#ifndef USING_THE_MAC
  if( !filename ){
    display("Look in the xterm...");
#ifndef USE_READLINE
    printf("Enter the file name containing the geometry data: ");
    scanf("%s\n",file_name);
#else
    theinline= readline("Enter the file name containing the geometry data: ");
    add_history(theinline);
    if( theinline ){
      sscanf(theinline,"%s",file_name);
      free(theinline);
    } else {
      error("Bad file name");
      file_name[0] = 0;
    }
#endif

  }else{
    strcpy(file_name,filename);
  }

  /* open the file */
  infile = fopen(file_name,"r");
  if(!infile){
    printf("Problems opening file: %s\n",file_name);
    display("oooooops!");
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
    return;
  }
#endif


  if( infile ){
    strcpy(whichobj->prim->molec->filename,file_name);
    read_molecule_data(infile,whichobj->prim->molec);
  }


  /* check to see if any atoms were actually read in.... */
  if(!whichobj->prim->molec->num_atoms){
    /* no... free the memory that we asked for */
    D_FREE(whichobj->prim->molec);
    D_FREE(whichobj->prim);
    D_FREE(whichobj);
    whichobj=0;
    head = head->next;
  }
  else{
    whichobj->scale.x=whichobj->scale.y=whichobj->scale.z=0.5;
    whichobj->trans.x=0;whichobj->trans.y=0;
    whichobj->trans.z=0;
    whichobj->cent.x = g_xmax/2;
    whichobj->cent.y = g_ymax/2;

#ifdef INTERACTIVE_USE
    /* create the button window */
    if(useButtons) {
      build_molec_button_window(&button_wins,whichobj->prim->molec);
      whichobj->prim->but_win = button_wins;
    }
#endif
    /* set up some default parameters */
    whichobj->prim->molec->draw_connectors = 1;
    whichobj->prim->molec->fancy_lines = 1;
    whichobj->prim->molec->outlines_on = 1;
    whichobj->prim->molec->shading_on = 1;
    whichobj->prim->molec->hydrogens_on = 1;
    whichobj->prim->molec->dummies_on = 1;
    whichobj->prim->molec->line_width = 1;
    whichobj->prim->molec->rad_mult = 1.0;
    whichobj->prim->molec->tubes_on = 1;
    whichobj->prim->molec->bond_rad = 0.05;
  }

}

