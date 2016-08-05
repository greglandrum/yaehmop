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

/*****************************************************************************
 *
 *               hierarchy.c
 *
 * This file contains the procedures for maintaining the object hierarchy tree
 *
 *****************************************************************************/
/***
  Recent Edit History:

  04.09.98 gL:
    added select_atoms_in_region
  09.09.98 gL:
    added invert_selected_atoms
  26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)
  06.04.99 gL:
     free_obj no longer frees the object itself, only its contents.
     This is the correct behavior.  Other functions changed to
     reflect this.
***/
#include "viewkel.h"

/****************************************************************************
 *
 *                   Procedure makenewobject
 *
 * Arguments: none
 * Returns: none
 *
 * Action:  allocates space for an object and puts it at the beginning of the
 *   head list
 *
 *****************************************************************************/
void makenewobject(void)
{
  head_type *new_head;
  object_type *new_object;

  /* allocate a new head node */
  new_head = (head_type *)D_CALLOC(1,sizeof(head_type));
  /* now allocate space for the new object */
  new_object = (object_type *)D_CALLOC(1,sizeof(object_type));
  if( !new_head || !new_object )fatal( "Memory allocation." );

  new_head->next = head;
  head = new_head;

  /* set the pointer to this object */
  new_head->obj=new_object;

  /* set up some of the information for the new object */
  new_object->scale.x=1.0;new_object->scale.y=1.0;new_object->scale.z=1.0;
  new_object->trans.x=0.0;new_object->trans.y=0.0;new_object->trans.z=camera->la.z;
}


/****************************************************************************
 *
 *                   Procedure addchild
 *
 * Arguments: parent, child: pointers to object_type
 * Returns: none
 *
 * Action:  adds the child to the parent by creating a new object node,
 *    copying the parents information into that node, and then pointing
 *    the parent node child pointers at the new node and at child.  Also
 *    sets the current pointers to child to equal zero.
 *
 *****************************************************************************/
void addchild( object_type *parent, object_type *child )
{
  int i;
  object_type *newobj;
  head_type *temp1, *temp2;

  if( parent == child ){
    display("No way bud.");
    return;
  }

  newobj=(object_type *)D_CALLOC(1,sizeof(object_type));
  if( !newobj ) fatal("Memory allocation");

  i=0;
  /* copy the parent's information into newobj */
  newobj->prim = parent->prim;
  parent->prim = 0;
  newobj->rot.x=parent->rot.x;newobj->rot.y=parent->rot.y;newobj->rot.z=parent->rot.z;
  newobj->scale.x=parent->scale.x;newobj->scale.y=parent->scale.y;newobj->scale.z=parent->scale.z;
  newobj->cent.x=parent->cent.x;newobj->cent.y=parent->cent.y;newobj->cent.z=parent->cent.z;
  while(parent->children[i]){
    newobj->children[i]=parent->children[i];
    parent->children[i]=0;
    i++;
  }
  /* there is no need to update the bounding box since the screen will be
     redrawn after this */

  /* now, the parents parms should be set */
  parent->scale.x=1;parent->scale.y=1;parent->scale.z=1;
  parent->rot.x=0;parent->rot.y=0;parent->rot.z=0;
  parent->cent.x=0;parent->cent.y=0;parent->cent.z=0;

  /* update the child's parameters */
  child->trans.x-=parent->trans.x;child->trans.y-=parent->trans.y;child->trans.z-=parent->trans.z;

  /* set the pointers to the new object and to child */
  parent->children[0]=newobj;
  parent->children[1]=child;

  /* now clear up the pointers to child */
  if( parenthead ){
    /* if the child comes off of a forest node, then kill that node */
    temp1=head;
    while( temp1 && temp1 != parenthead ){
      temp2=temp1;
      temp1=temp1->next;
    }
    if( !temp1 ){
      error("Something funny just happened.");
      return;
    }
    temp2->next=temp1->next;
    D_FREE(temp1);
  }
  else if( parentobj ){
    /* if the child is currently the child of another node, then fix that */
    i=0;
    /* find the child pointer in the array */
    while(parentobj->children[i]&&parentobj->children[i]!=child) i++;
    if(!parentobj->children[i]){
      error("Ooops, something broke.");
      return;
    }
    while(parentobj->children[i]){
      parentobj->children[i]=parentobj->children[i+1];
      i++;
    }
  }
  else error("What is going on here???");

  /* redraw the screen */
  redraw();
}


/****************************************************************************
 *
 *                   Procedure free_obj
 *
 * Arguments: obj: pointer to object_type
 * Returns: none
 *
 * Action:  frees up all the memory used by 'obj and its primitives.
 *
 *****************************************************************************/
void free_obj(object_type *obj)
{
  int i,j;

  if( obj->prim->molec ){
    molec_type *molec;
    molec = obj->prim->molec;

    /* the contents of the individual atoms need to be free'd too */
    for(i=0;i<molec->num_atoms;i++){
      if( molec->atoms[i].p_surf ){
        D_FREE(molec->atoms[i].p_surf);
      }
      if( molec->atoms[i].linesto ){
        D_FREE(molec->atoms[i].linesto);
      }
    }
    if( molec->atoms ) D_FREE(molec->atoms);
    if( molec->lines ) D_FREE(molec->lines);
    D_FREE(molec->num_lines);
    if( molec->polyhed_centers ) D_FREE(molec->polyhed_centers);
    if( molec->polyhed_rads ) D_FREE(molec->polyhed_rads);
    if( molec->polyhed_verts ) D_FREE(molec->polyhed_verts);
    if( molec->triangles ) D_FREE(molec->triangles);
    D_FREE(molec);
  }
  if( obj->prim->cont_plot ){
    contour_plot_type *cont_plot;
    iso_curve_type *curve1,*curve2;
    gnuplot_contour_type *cont1,*cont2;

    cont_plot = obj->prim->cont_plot;
    curve1 = cont_plot->data;
    while(curve1){
      if( curve1->points ) D_FREE(curve1->points);
      curve2 = curve1->next;
      D_FREE(curve1);
      curve1 = curve2;
    }
    curve1 = cont_plot->raw_data;
    while(curve1){
      if( curve1->points ) D_FREE(curve1->points);
      curve2 = curve1->next;
      D_FREE(curve1);
      curve1 = curve2;
    }
    for(i=0;i<cont_plot->num_curves;i++){
      cont1 = cont_plot->contours[i];
      while(cont1){
        D_FREE(cont1->coords);
        cont2 = cont1->next;
        D_FREE(cont1);
        cont1 = cont2;
      }
    }

    if( cont_plot->contours ) D_FREE(cont_plot->contours);

    if( cont_plot->levels_list) D_FREE(cont_plot->levels_list);
    D_FREE(cont_plot->styles);
    D_FREE(cont_plot->curves_to_display);
    D_FREE(cont_plot);
  }
  if( obj->prim->MO_surf ){
    molec_type *molec;
    MO_surface_type *MO_surf;
    MO_contours_type *cont1,*cont2;
    iso_curve_type *curve1,*curve2;

    MO_surf = obj->prim->MO_surf;
    molec = obj->prim->MO_surf->molec;

    if( MO_surf->characters) D_FREE(MO_surf->characters);
    if( MO_surf->kpoints) D_FREE(MO_surf->kpoints);
    if( MO_surf->MO_numbers) D_FREE(MO_surf->MO_numbers);
    for(i=0;i<MO_surf->num_centers;i++){
      D_FREE(MO_surf->MO_centers[i].AO_list);
    }
    D_FREE(MO_surf->MO_centers);
    for(i=0;i<MO_surf->num_centers_in_cell;i++){
      D_FREE(MO_surf->raw_MO_centers[i].AO_list);
    }
    if( MO_surf->raw_MO_centers != MO_surf->MO_centers )
      D_FREE(MO_surf->raw_MO_centers);
    for(i=0;i<MO_surf->num_unique;i++){
      for(j=0;j<MO_surf->unique_centers[i].num_AOs;j++){
        D_FREE(MO_surf->unique_centers[i].AO_list[j].rad_lookup_tbl->values);
        D_FREE(MO_surf->unique_centers[i].AO_list[j].rad_lookup_tbl);
        while(MO_surf->unique_centers[i].AO_list[j+1].rad_lookup_tbl ==
              MO_surf->unique_centers[i].AO_list[j].rad_lookup_tbl) j++;
      }
      D_FREE(MO_surf->unique_centers[i].AO_list);
      D_FREE(MO_surf->unique_centers[i].type);
    }
    D_FREE(MO_surf->unique_centers);

    if( MO_surf->MO_contours){
      cont1 = MO_surf->MO_contours->contours;
      while(cont1){
        D_FREE(cont1->coords);
        cont2 = cont1->next;
        D_FREE(cont1);
        cont1 = cont2;
      }
      curve1 = MO_surf->MO_contours->data;
      while(curve1){
        D_FREE(curve1->points);
        curve2 = curve1->next;
        D_FREE(curve1);
        curve1 = curve2;
      }
      if(MO_surf->MO_contours->levels_list)
        D_FREE(MO_surf->MO_contours->levels_list);
      D_FREE(MO_surf->MO_contours);
    }
    D_FREE(MO_surf);

    /* the contents of the individual atoms need to be free'd too */
    for(i=0;i<molec->num_atoms;i++){
      if( molec->atoms[i].p_surf ){
        D_FREE(molec->atoms[i].p_surf);
      }
      if( molec->atoms[i].linesto ){
        D_FREE(molec->atoms[i].linesto);
      }
    }
    if( molec->atoms ) D_FREE(molec->atoms);
    if( molec->lines ) D_FREE(molec->lines);
    D_FREE(molec->num_lines);
    if( molec->polyhed_centers ) D_FREE(molec->polyhed_centers);
    if( molec->polyhed_rads ) D_FREE(molec->polyhed_rads);
    if( molec->polyhed_verts ) D_FREE(molec->polyhed_verts);
    if( molec->triangles ) D_FREE(molec->triangles);
    D_FREE(molec);
  }
  if( obj->prim->graph ){
    D_FREE(obj->prim->graph->data);
    D_FREE(obj->prim->graph->raw_data);
    D_FREE(obj->prim->graph);
  }
  if( obj->prim->prop_graph ){
    D_FREE(obj->prim->prop_graph->the_data->data);
    D_FREE(obj->prim->prop_graph->the_data->raw_data);
    if(obj->prim->prop_graph->the_data->curve_names)
      D_FREE(obj->prim->prop_graph->the_data->curve_names);
    D_FREE(obj->prim->prop_graph->the_data->curves_to_display);
    D_FREE(obj->prim->prop_graph->the_data->styles);
    D_FREE(obj->prim->prop_graph->the_data->fills);
    D_FREE(obj->prim->prop_graph->the_data);

    if(obj->prim->prop_graph->the_integration->curve_names)
      D_FREE(obj->prim->prop_graph->the_integration->curve_names);
    D_FREE(obj->prim->prop_graph->the_integration->curves_to_display);
    D_FREE(obj->prim->prop_graph->the_integration->styles);
    D_FREE(obj->prim->prop_graph->the_integration->fills);
    D_FREE(obj->prim->prop_graph->the_integration);
    D_FREE(obj->prim->prop_graph->the_integration->data);
    D_FREE(obj->prim->prop_graph->the_integration->raw_data);

    D_FREE(obj->prim->prop_graph);
  }
  if( obj->prim->band_graph ){
    D_FREE(obj->prim->band_graph->special_points);
    D_FREE(obj->prim->band_graph->the_data->styles);
    D_FREE(obj->prim->band_graph->the_data->curves_to_display);
    D_FREE(obj->prim->band_graph->the_data->raw_data);
    D_FREE(obj->prim->band_graph->the_data->data);
    D_FREE(obj->prim->band_graph->the_data);
    D_FREE(obj->prim->band_graph);
  }
  if( obj->prim->walsh_graph ){
    D_FREE(obj->prim->walsh_graph->the_data->data);
    D_FREE(obj->prim->walsh_graph->the_data->raw_data);
    D_FREE(obj->prim->walsh_graph->the_data->styles);
    D_FREE(obj->prim->walsh_graph->the_data->curves_to_display);
    D_FREE(obj->prim->walsh_graph->the_data);
    D_FREE(obj->prim->walsh_graph->total_E->data);
    D_FREE(obj->prim->walsh_graph->total_E->raw_data);
    D_FREE(obj->prim->walsh_graph->total_E);
    D_FREE(obj->prim->walsh_graph);
  }
  if( obj->prim->p_surf ){
    D_FREE(obj->prim->p_surf->points);
    D_FREE(obj->prim->p_surf->colors);
    D_FREE(obj->prim->p_surf);
  }
  if( obj->prim->FMO_diagram ){
    FMO_level_type *lev1,*lev2;
    FMO_connect_type *conn1,*conn2;
    FMO_diagram_type *diagram;
    diagram = obj->prim->FMO_diagram;
    lev1 = diagram->levels;
    while(lev1){
      lev2 = lev1;
      lev1 = lev1->next;
      D_FREE(lev2);
    }
    conn1 = diagram->connects;
    while(conn1){
      conn2 = conn1;
      conn1 = conn1->next;
      D_FREE(conn2);
    }
    for(i=0;i<diagram->num_frags;i++){
      lev1 = diagram->frags[i].levels;
      while(lev1){
        lev2 = lev1;
        lev1 = lev1->next;
        D_FREE(lev2);
      }
      D_FREE(diagram->frags[i].raw_energies);
    }
    if(diagram->chg_mat) D_FREE(diagram->chg_mat);
    if( diagram->frags ) D_FREE(diagram->frags);
    if( diagram->raw_energies ) D_FREE(diagram->raw_energies);
    if( diagram->occups ) D_FREE(diagram->occups);
    D_FREE(diagram);
  }
  if( obj->prim->label ){
    if(obj->prim->label->atoms_to_label){
      D_FREE(obj->prim->label->atoms_to_label);
    }
    D_FREE(obj->prim->label);
  }

  /* free up the button window, if it exists */
  if( obj->prim->but_win && obj->prim->but_win->child )
    free_child_but_win(obj->prim->but_win->child);
  if( obj->prim->but_win ) free_but_win(obj->prim->but_win,&button_wins);

  D_FREE(obj->prim);
}






/****************************************************************************
 *
 *                   Procedure clearall
 *
 * Arguments: none
 * Returns: none
 *
 * Action:  clears the window by disposing of all the objects
 *
 *****************************************************************************/
void clearall(void)
{
  int i;
  object_type *obj1, *obj2;
  head_type *head1, *head2;

  /* make sure that there are some objects present */
  if( !head ) return;

  /* get rid of them */
  head1=head;
  while( head1 ){
    head2 = head1->next;
    obj1=head1->obj;
    /* go through the children */
    i=0;
    obj2=obj1->children[i];
    while( obj2 ){
      free_obj( obj2 );
      D_FREE(obj2);
      i++;
      obj2=obj1->children[i];
    }
    free_obj( obj1 );
    D_FREE(obj1);
    D_FREE(head1);
    head1 = head2;
  }
  head=0;
  whichobj=0;
  /* set the camera stuff back to their default values */
  camera->la.z=159.75;
  camera->la.x=0;
  camera->la.y=0;
  camera->vup.y=1;
  camera->vup.x=0;
  camera->vup.z=0;
  camera->hsize=10;
  camera->foclength=15.975;

  redrawgraph();
  display("Hierarchy tree purged");
}


/****************************************************************************
 *
 *                   Procedure clear_labels
 *
 * Arguments: none
 * Returns: none
 *
 * Action:  clears all the labels currently displayed
 *
 *****************************************************************************/
void clear_labels(void)
{
  int i;
  object_type *obj1, *obj2;
  head_type *head1, *head2;

  /* make sure that there are some objects present */
  if( !head ) return;

  /* get rid of them */
  head1=head;
  while( head1 ){
    if( head1->obj->prim->label){
      head2 = head1->next;
      obj1=head1->obj;
      /* go through the children */
      i=0;
      obj2=obj1->children[i];
      while( obj2 ){
        free_obj( obj2 );
        D_FREE(obj2);
        i++;
        obj2=obj1->children[i];
      }
      free_obj( obj1 );
      if ( head == head1 ) head = head2;
      D_FREE(obj1);
      D_FREE(head1);
      head1 = head2;
    } else{
      head1 = head1->next;
    }
  }
  whichobj=head->obj;

  redrawgraph();
  display("Labels purged");
}



/****************************************************************************
 *
 *                   Procedure select_object
 *
 * Arguments: xpos,ypos: integers
 *
 * Returns: none
 *
 * Action:   This finds the first object with a bounding box that
 *  includes the point 'xpos,'ypos.
 *    if no object is found then the current object (whichobj) is not changed.
 *
 *****************************************************************************/
void select_object(int xpos,int ypos)
{
  head_type *temphead;
  int done;

  /* step through all the objects */
  temphead = head;
  done = 0;
  while(temphead && !done ){
    if( xpos > (int)temphead->obj->bmin.x && xpos < (int)temphead->obj->bmax.x
       && ypos > (int)temphead->obj->bmin.y && ypos < (int)temphead->obj->bmax.y ){
      /* this is the one */
      whichobj = temphead->obj;
      display("Object changed.");
#if 0
fprintf(stderr,"pos: (%d,%d), bmin: (%f %f), bmax: (%f %f)\n",
        xpos,ypos,temphead->obj->bmin.x,temphead->obj->bmin.y,
        temphead->obj->bmax.x,temphead->obj->bmax.y);
#endif
      done = 1;
    }
    else temphead = temphead->next;
  }
  if( !done ) display("No object there.");
}



/****************************************************************************
 *
 *                   Function select_atom
 *
 * Arguments:        molec: pointer to molec_type
 *               xpos,ypos: integers
 *
 * Returns: int
 *
 * Action:   This finds the first atom that
 *    includes the point 'xpos,'ypos in its radius
 *    if no atom is found then nothing is done
 *
 *      Returns nonzero if an atom was selected
 *
 *****************************************************************************/
int select_atom(molec_type *molec,int xpos,int ypos)
{
  int i,j;
  int which;
  float dist;
  atom_type *atoms;

  if( !molec ) FATAL_BUG("Bogus molec passed to select_atom");

  if(molec->num_frames > 1)
    atoms = &(molec->atoms[(molec->current_frame%molec->num_frames)*molec->num_atoms]);
  else atoms = molec->atoms;

  /* step through the atoms */
  which = -1;
  for(i=0;i<molec->num_atoms;i++){
    /* don't bother with hidden atoms */
    if( !atoms[i].exclude ){
      dist = (xpos - atoms[i].screen_loc.x)*
        (xpos - atoms[i].screen_loc.x) +
          (ypos - atoms[i].screen_loc.y) *
            (ypos - atoms[i].screen_loc.y);
      if( dist <= (float)(atoms[i].screen_rad*
                          atoms[i].screen_rad)){
        /*******

          check the drawing order, if this atom was drawn
          after the currently selected one (i.e. if it's
          closer to the user in screen coords), then choose
          it instead.... this allows front atoms to have
          higher priority

        ********/
        if( which >= 0 ){
          if( atoms[which].draw_order >
             atoms[i].draw_order ){
            which = i;
          }
        } else{
          which = i;
        }
      }
    }
  }

  if( which >= 0 ){
    if( atoms[which].is_selected ){

      /*******

        removing selected atoms can be tricky if they
        are low numbered (we don't want to leave holes in
        the selected atoms array).

        so... loop over the other selected atoms and
        decrement their is_selected fields.
      *******/
      if( atoms[which].is_selected < num_selected){
        for(i=atoms[which].is_selected+1;i<=num_selected;i++){
          /* find this selected atom */
          for(j=0;j<molec->num_atoms;j++){
            if( atoms[j].is_selected == i ){
              atoms[j].is_selected--;
              break;
            }
          }
        }
      }
      atoms[which].is_selected = 0;
      num_selected--;
    }else{
      /*******

        keep track of the selection order, so that
        the handedness of stuff can be dealt with

        ********/
      num_selected++;
      atoms[which].is_selected = num_selected;
    }
    display("Atom selected.");
          return(1);
  }else{
    display("No atom found.");
    return(0);
  }
}

/****************************************************************************
 *
 *                   Function select_atoms_in_region
 *
 * Arguments:        molec: pointer to molec_type
 *               startx,starty: integers
 *               endx,endy: integers
 *
 * Returns: int
 *
 * Action:   This finds all of the atoms which
 *    fall inside the box defined by startx,starty,endx,endy
 *
 *    if no atoms are found then nothing is done
 *
 *      Returns nonzero if an atom was selected
 *
 *****************************************************************************/
int select_atoms_in_region(molec_type *molec,int startx,int starty,
                           int endx,int endy)
{
  int i;
  int which=0;
  atom_type *atoms;

  if( !molec ) FATAL_BUG("Bogus molec passed to select_atom");

  if(molec->num_frames > 1)
    atoms = &(molec->atoms[(molec->current_frame%molec->num_frames)*molec->num_atoms]);
  else atoms = molec->atoms;

  /* step through the atoms */
  which = -1;
  for(i=0;i<molec->num_atoms;i++){
    if( atoms[i].screen_loc.x >= startx &&
        atoms[i].screen_loc.x <= endx &&
        atoms[i].screen_loc.y >= starty &&
        atoms[i].screen_loc.y <= endy ){
      which = i;
      if( !atoms[which].is_selected ){
        num_selected++;
        atoms[which].is_selected = num_selected;
      }
    }
  }
  if(which ){
    display("Atom(s) selected.");
    return(1);
  }else{
    display("No atom found.");
    return(0);
  }
}

