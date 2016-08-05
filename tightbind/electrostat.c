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
*   This is the stuff for dealing with the electrostatic potential terms
*
*  created:  greg landrum  January 1994
*
*****************************************************************************/
#include "bind.h"


/****************************************************************************
 *
 *    Function factorial
 *
 * Arguments: num: an integer
 *
 * Returns: an integer
 *
 * Action:  evaluates  num!  (obviously)
 *
 *
 ****************************************************************************/
int factorial(num)
  int num;
{
  int i;
  real accum;

  accum = 1;
  for(i=1;i<=num;i++){
    accum *= i;
  }

  return(accum);
}

/* this is the helper function for quicksorting the atomic energy level array */
int compare_energies(a,b)
  const void *a,*b;
{
  return( (int)( *(real *)a - *(real *)b) );
}

/****************************************************************************
 *
 *                   Procedure free_atomic_energy
 *
 * Arguments:  cell: pointer to cell type
 *           energy: pointer to real
 *       work_array: pointer to real
 *
 * Returns: none
 *
 * Action:  fills up the lowest num_electrons/2 free atomic levels to find
 *   the "free atomic energy" used in extended Hueckel binding calculations.
 *
 ****************************************************************************/
void free_atomic_energy(cell,energy,work_array)
  cell_type *cell;
  real *energy,*work_array;
{
  atom_type *atom;
  int i;
  int orbs_so_far;
  real temp,electrons_left;


  /*********
    to find the energy of the "free" atoms reasonably (like a 1 electron theory
    should) we have to fill the free levels (not by atoms).

    to do this build an array of the coulombic levels, sort it by energy, then
    fill it.
  *********/
  orbs_so_far = 0;
  for(i=0;i<cell->num_atoms;i++){
    atom = &(cell->atoms[i]);
    if( atom->ns != 0 ) work_array[orbs_so_far++] = atom->coul_s;
    if( atom->np != 0 ){
      /* three entries for p orbitals */
      work_array[orbs_so_far++] = atom->coul_p;
      work_array[orbs_so_far++] = atom->coul_p;
      work_array[orbs_so_far++] = atom->coul_p;
    }
    if( atom->nd != 0 ){
      work_array[orbs_so_far++] = atom->coul_d;
      work_array[orbs_so_far++] = atom->coul_d;
      work_array[orbs_so_far++] = atom->coul_d;
      work_array[orbs_so_far++] = atom->coul_d;
      work_array[orbs_so_far++] = atom->coul_d;
    }
  }

  /* now sort the array by energy */
  qsort((void *)work_array,orbs_so_far,sizeof(real),compare_energies);

  /* that's that, now loop until we use up the number of electrons */
  electrons_left = cell->num_electrons;
  temp = 0.0;
  orbs_so_far = 0;
  while( electrons_left > 0.0){
    if( electrons_left > 2.0 ){
      temp += 2.0 * work_array[orbs_so_far++];
      electrons_left -= 2.0;
    }
    else{
      temp += electrons_left * work_array[orbs_so_far++];
      electrons_left = 0.0;
    }
  }

  *energy = temp;
}

/****************************************************************************
 *
 *                   Procedure AO_occupations
 *
 * Arguments:  cell: pointer to cell type
 *         num_orbs: int
 *           OP_mat: pointer to real
 *     orbital_lookup_table: int
 *            accum: pointer to real
 *
 *
 *
 * Returns: none
 *
 * Action:  Determines the total number of electrons in each type of atomic
 *   orbital for each atom...  i.e. the number of electrons in the s orbital, the
 *   number in the p orbitals, etc.
 *
 *  the results are stored in 'accum which should be at least 'num_orbs long
 *
 ****************************************************************************/
void AO_occupations(cell,num_orbs,OP_mat,orbital_lookup_table,accum)
  cell_type *cell;
  int num_orbs;
  int *orbital_lookup_table;
  real *OP_mat;
  real *accum;
{
  atom_type *atom;
  static real *free_atom_occups=0;
  int orbs_so_far,orb_tab;
  real electrons_left;
  int i,j,k;


#if 0
/***************************************

  This code is for finding the real orbital populations....
  since there are problems with the way molecules dissociate in EHT, we
  have to use the occupations for the free atoms until the code to do the
  electrostatics has been written "for real"

***************************************/


  orbs_so_far = 0;
  for(i=0;i<cell->num_atoms;i++){
    orb_tab = orbital_lookup_table[i];
    atom = &(cell->atoms[i]);

    /*******
      note that the diagonal element(s) for each orbital are hit twice since we add up
       .5 of each contribution.
    *******/
    if( atom->ns != 0 ){
      /* sweep down */
      for(j=orb_tab;j<num_orbs;j++){
        accum[orbs_so_far] += .5*OP_mat[j*num_orbs + orb_tab];
      }
      /* sweep across */
      for(j=0;j<=orb_tab;j++){
        accum[orbs_so_far] += .5*OP_mat[orb_tab*num_orbs + j];
      }
fprintf(stderr,"Atom %d s orbital occupation: %lg\n",i,accum[orbs_so_far]);
      orbs_so_far++;
    }

    if( atom->np != 0 ){
      for(j=orb_tab + BEGIN_P;j <= orb_tab + END_P;j++){
        for( k=j;k<num_orbs;k++){
          accum[orbs_so_far] += .5*OP_mat[k*num_orbs+j];
        }
        for( k=0;k<=j;k++){
          accum[orbs_so_far] += .5*OP_mat[j*num_orbs+k];
        }
      }

fprintf(stderr,"Atom %d p orbital occupation: %lg\n",i,accum[orbs_so_far]);
      orbs_so_far++;
    }
    if( atom->nd != 0 ){
      for(j=orb_tab + BEGIN_D;j <= orb_tab + END_D;j++){
        for(k=j;k<num_orbs;k++){
          accum[orbs_so_far] += .5*OP_mat[k*num_orbs+j];
        }
        for(k=0;k<=j;k++){
          accum[orbs_so_far] += .5*OP_mat[k*num_orbs+j];
        }
      }
      orbs_so_far++;
    }
  }
#endif

  /********
    figure out the free atom orbital occupations, we only have to do
    this once.
  ********/

  /* get space to store the free occupations, if we need it */
  if( !free_atom_occups ){
    free_atom_occups = (real *)calloc(num_orbs,sizeof(real));
    if( !free_atom_occups)
      fatal("Can't allocate space for free atoms occupations in AO_occupations");

    /******
      fill the array of free atom occupations

       this isn't set up for d electrons yet.
    *******/
    orbs_so_far = 0;
    for( i=0;i<cell->num_atoms;i++){
      atom = &(cell->atoms[i]);
      electrons_left = atom->num_valence;

      if( atom->ns ){
        if( electrons_left >= 2.0 ){
          free_atom_occups[orbs_so_far++] = 2.0;
          electrons_left -= 2.0;
        }
        else{
          free_atom_occups[orbs_so_far++] = electrons_left;
          electrons_left = 0.0;
        }
      }
      if( atom->np ){
        if( electrons_left >= 6.0 ){
          free_atom_occups[orbs_so_far++] = 6.0;
          electrons_left -= 6.0;
        }
        else{
          free_atom_occups[orbs_so_far++] = electrons_left;
          electrons_left = 0.0;
        }
      }
    }
  }

  /*****
    copy the free atom_occupation numbers into the accum array
    which was passed in
  ******/
  bcopy((char *)free_atom_occups,(char *)accum,num_orbs*sizeof(real));
}

/****************************************************************************
 *
 *                   Procedure eval_electrostat
 *
 * Arguments:  cell: pointer to cell type
 *         num_orbs: int
 *         eigenset: eigenset_type
 *      occupations: pointer to real
 *           OP_mat: pointer to real
 *     orbital_lookup_table: int
 * electrostat_term,total_E: pointers to real
 *         eHMO_term: pointer to real
 *            accum: pointer to real
 *
 *
 *
 * Returns: none
 *
 * Action:   Evaluates the electrostatic repulsion and total energy for the
 *    molecule.
 *
 *   'accum is used to build the orbital occupation array.  it should be at least
 *     num_orbs long
 *
 * The equation for the electrostatic repulsion term comes from Gion's first
 *   EHMO-ASED paper:    J. Phys. Chem. _93_ 5366 (1989)  eqns 17 and the appendix
 *
 *  I've tried to keep the notation in the code similar.
 *
 ****************************************************************************/
void eval_electrostatics(cell,num_orbs,eigenset,occupations,OP_mat,orbital_lookup_table,
                         electrostat_term,eHMO_term,total_E,accum,net_chgs)
  cell_type *cell;
  int num_orbs;
  int *orbital_lookup_table;
  real *occupations,*OP_mat;
  eigenset_type eigenset;
  real *electrostat_term,*eHMO_term,*total_E;
  real *accum,*net_chgs;
{
  static real *atomic_energy=0;
  atom_type *atomA,*atomB;
  int i,orbs_so_far,orb_tab;
  int n,l,p;
  int numA,numB;
  real rhoA_term,rhoB_term;
  real ZA,ZB;
  real R,zeta;
  real p_sum;
  real electro_accum;
  int electrons_left;

  /**************

    make sure that we have memory to accumulate the "free" atomic energies

  ***************/
  if( !atomic_energy ){
    atomic_energy = (real *)calloc(cell->num_atoms,sizeof(real));
    if( !atomic_energy )
      fatal("Cannot allocate memory for atomic energy array in eval_electrostatics.");
  }

  /*******
    zero out the accumulator array
  ********/
  bzero((char *)accum,num_orbs*sizeof(real));

  /*******
    build the orbital occupation array
  ********/
  AO_occupations(cell,num_orbs,OP_mat,orbital_lookup_table,accum);


#if 0
  /******
    NOTE:  Due to the fact that I want to get some quick results, this is only
    tested on diatomics at the moment.
    ******/
  numA = 0;
  numB = 1;
#endif

  /******
    find the energies of the "Free" atoms

    once again this is not dealing with d orbitals
  *******/
  for(i=0;i<cell->num_atoms;i++){
    atomA = &(cell->atoms[i]);
    electrons_left = atomA->num_valence;
    atomic_energy[i] = 0.0;
    if( atomA->ns != 0 ){
      if(electrons_left > 1){
        atomic_energy[i] += atomA->coul_s * 2.0;
        electrons_left -= 2;
      }
      else{
        atomic_energy[i] += atomA->coul_s;
        electrons_left = 0;
      }
    }
    if( atomA->np != 0 ){
      atomic_energy[i] += atomA->coul_p * (real)electrons_left;
    }
  }

  /*****

    find the extended Hueckel energy

  ******/
  *eHMO_term = 0.0;
  for(i=0;i<num_orbs;i++){
    *eHMO_term += EIGENVAL(eigenset,i)*occupations[i];
  }


  /**********

    loop over all the atoms

  **********/
  *electrostat_term = 0.0;
  for(numA=0;numA<cell->num_atoms;numA++){
    for(numB=numA+1; numB<cell->num_atoms; numB++){
      atomA = &(cell->atoms[numA]);
      atomB = &(cell->atoms[numB]);

      /************
        This isn't set up particularly efficiently...
        if this ever becomes a bottleneck, the code is easily improved.
      *************/


      /*******
        get the distance from the symmetric distance matrix and convert it to bohrs...
        (this reference assumes that numA < numB)
      ********/
      R = cell->distance_mat[ numB*(numB+1)/2 + numA ]/BOHR;

      /******

        evaluate the electrostatic term due to center A

      *******/
      rhoA_term = 0.0;
      orbs_so_far = orbital_lookup_table[numA];
      if( atomA->ns != 0){
        zeta = atomA->exp_s;
        n = atomA->ns;
        p_sum = 0;
        for(p=1; p <= 2*n; p++){
          p_sum += pow((2.0*R*zeta),(2*n-p)) *   (real)p / (real)factorial(2*n-p);
        }
        p_sum *= exp(-2.0*R*zeta)/(2.0*n*R);

        rhoA_term += accum[orbs_so_far] * (1/R - p_sum);
        orbs_so_far++;
      }
      if( atomA->np != 0){
        zeta = atomA->exp_p;
        n = atomA->np;
        p_sum = 0;
        for(p=1; p <= 2*n; p++){
          p_sum += pow((2.0*R*zeta),(2*n-p)) *   (real)p / (real)factorial(2*n-p);
        }
        p_sum *= exp(-2.0*R*zeta)/(2.0*n*R);

        rhoA_term += accum[orbs_so_far] * (1/R - p_sum);
        orbs_so_far++;
      }
      if( atomA->nd != 0){
        /********
          d orbitals don't work yet, since I haven't gotten around to figuring out
          how to deal with double zeta wave functions yet.
          ********/
        error("You can't optimize structures with d orbitals yet... sorry.");
      }

      /******
        evaluate the term due to center B
      *******/
      rhoB_term = 0.0;
      orbs_so_far = orbital_lookup_table[numB];
      if( atomB->ns != 0){
        zeta = atomB->exp_s;
        n = atomB->ns;
        p_sum = 0;
        for(p=1; p <= 2*n; p++){
          p_sum += pow((2.0*R*zeta),(2*n-p)) *   (real)p / (real)factorial(2*n-p);
        }
        p_sum *= exp(-2.0*R*zeta)/(2.0*n*R);

        rhoB_term += accum[orbs_so_far] * (1/R - p_sum);
        orbs_so_far++;
      }
      if( atomB->np != 0){
        zeta = atomB->exp_p;
        n = atomB->np;
        p_sum = 0;
        for(p=1; p <= 2*n; p++){
          p_sum += pow((2.0*R*zeta),(2*n-p)) *   (real)p / (real)factorial(2*n-p);
        }
        p_sum *= exp(-2.0*R*zeta)/(2.0*n*R);

        rhoB_term += accum[orbs_so_far] * (1/R - p_sum);
        orbs_so_far++;
      }
      if( atomB->nd != 0){
        /********
          d orbitals don't work yet, since I haven't gotten around to figuring out
          how to deal with double zeta wave functions yet.
          ********/
        error("You can't optimize structures with d orbitals yet... sorry.");
      }

      /*********
        figure out the nuclear charges...

        This uses Slater's rules i.e. each electron in the next shell in shields for .8,
        anything farther in shields for 1. We're not taking the valence orbitals
        into account here (they come into play in the rho terms calculated above).
      *********/
      ZA = atomA->num_valence;
      ZB = atomB->num_valence;

      /* things are wierd for d orbitals (and they aren't implemented yet anyway) */
#if 0
      if( atomA->nd == 0 ){
        /********
          we can use the quantum number of the s and p orbitals to figure out the number
          of electrons in the next shell in.
          *********/
#define MULTI .15
        switch(atomA->ns){
        case 2:
          ZA += NEXT_SHELL_IN * 2;
          break;
        case 3:
          ZA += NEXT_SHELL_IN * 6;
          break;
        case 4:
          ZA += NEXT_SHELL_IN * 16;
          break;
        }
      }


      if( atomB->nd == 0 ){
        /********
          we can use the quantum number of the s and p orbitals to figure out the number
          of electrons in the next shell in.
          *********/
        switch(atomB->ns){
        case 2:
          ZB += NEXT_SHELL_IN * 2;
          break;
        case 3:
          ZB += NEXT_SHELL_IN * 6;
          break;
        case 4:
          ZB += NEXT_SHELL_IN * 16;
          break;
        }
      }
#endif

      /*******
        now do the total electrostatic energy (this distance needs to be in bohr).
      *******/

      electro_accum = ZA*ZB/R;
      fprintf(stderr, "Za:  %lf Zb: %lf  R: %lf  ZaZb/R: %lf Za*rhoB: %lf Zb*rhoA %lf\n",
              ZA,ZB,R,electro_accum,ZA*rhoB_term,ZB*rhoA_term);

      electro_accum -= .5*(ZA*rhoB_term + ZB*rhoA_term);

      electro_accum *= 27.2;

      *electrostat_term += electro_accum;
      fprintf(stderr,"total electrostat thus far: %lf\n",*electrostat_term);
    }

    /*******
      subtract off the energy of atom i
    ********/
#if 0
    free_atomic_energy(cell,&atomic_energy,accum);
#endif

    *eHMO_term -= atomic_energy[numA];
  }
  /*****
    set the total energy and then return
  *****/
  *total_E = *electrostat_term + *eHMO_term;

  return;
}
