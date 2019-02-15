/*******************************************************

Copyright (C) 2018 Greg Landrum
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

#include "bind.h"

// For booleans
#include "stdbool.h"


void main(int argc, char **argv){
  char err_string[240];
  int i, j;
  FILE *nullfile = fopen("nul","w");

  status_file = nullfile;
  output_file = nullfile;

  unit_cell = (cell_type *)calloc(1,sizeof(cell_type));
  details = (detail_type *)calloc(1,sizeof(detail_type));
  if(!unit_cell || !details) fatal("Can't allocate initial memory.");
  FILE *dest=nullfile;

  /*******
    initialize some variables
    ********/
  set_details_defaults(details);
  set_cell_defaults(unit_cell);

  safe_strcpy(details->title,"test job");

  // molecular calculation
  details->Execution_Mode = MOLECULAR;
  details->num_KPOINTS = 1;
  details->K_POINTS = (k_point_type *)calloc(1,sizeof(k_point_type));
  if( !details->K_POINTS ) fatal("Can't allocate the single k point.");
  details->K_POINTS[0].weight = 1.0;
  details->avg_props = 0;
  details->use_symmetry = 1;
  details->net_chg_PRT=1;
  details->ROP_mat_PRT=1;


  unit_cell->using_Zmat = 0;
  unit_cell->using_xtal_coords = 0;


  unit_cell->num_atoms = 3;
  unit_cell->atoms = (atom_type *)calloc(unit_cell->num_atoms,sizeof(atom_type));
  if(!unit_cell->atoms){
    sprintf(err_string,"Can't allocate memory for: %d atoms.",unit_cell->num_atoms);
    fatal("Can't allocate memory for the atoms.");
  }
  safe_strcpy(unit_cell->atoms[0].symb,"H");
  unit_cell->atoms[0].loc.x = 0.0;
  unit_cell->atoms[0].loc.y = 0.0;
  unit_cell->atoms[0].loc.z = 0.0;

  safe_strcpy(unit_cell->atoms[1].symb,"C");
  unit_cell->atoms[1].loc.x = 0.0;
  unit_cell->atoms[1].loc.y = 0.0;
  unit_cell->atoms[1].loc.z = 1.0;

  safe_strcpy(unit_cell->atoms[2].symb,"N");
  unit_cell->atoms[2].loc.x = 0.0;
  unit_cell->atoms[2].loc.y = 0.0;
  unit_cell->atoms[2].loc.z = 2.1;

  unit_cell->charge = 0.0;


  // shouldn't need to change below here
  fill_atomic_parms(unit_cell->atoms,unit_cell->num_atoms,NULL,NULL);
  unit_cell->num_raw_atoms = unit_cell->num_atoms;
  charge_to_num_electrons(unit_cell);
  build_orbital_lookup_table(unit_cell,&num_orbs,&orbital_lookup_table);


  /* install the sig_int handler */
  signal(SIGINT,handle_sigint);

  run_eht(dest);

  //pull properties
  for(i=0;i<unit_cell->num_atoms;i++){
   printf(">>>> Atom %d: %.2f\n",i+1,properties.net_chgs[i]);
  }

  for(i=0;i<unit_cell->num_atoms;i++){
  for(j=0;j<i;j++){
   printf(">>>> ROP %d-%d: %.2f\n",i+1,j+1,properties.ROP_mat[i*(i+1)/2 + j]);
  }
}

  free(unit_cell);
  free(details);
  exit(0);
}
