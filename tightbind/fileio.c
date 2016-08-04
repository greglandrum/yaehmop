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
 *     this file contains stuff for dealing with file input and output
 *      i.e. reading and writing data
 *
 *  created:  greg landrum  August 1993
 *
 *****************************************************************************/

/***
  Edit History:

  March '98: WG
    - write f orbital parameters to .out [write_atom_parms()]
    - read in f orbital parameters, normalize f zeta coeff's [fill_atomic_parms()]
    - COHP keyword definitions [read_inputfile()]
  01.04.99: gL
    - changed the reading of projected DOS contributions so that it is now
      possible to put a '\' at the end of long lines and continue them
      on the next line.  This is a quick "solution" which should be
      helped along by some more sophisticated parsing of the projections.
  07.07.99: gL
    - fixed up the stupid problem with charges not being filled in when
      the charge keyword appeared below the parameters specification of
      an input file.
    - added some anal-retentive changes to string dimensioning as insurance.
      This is almost 100% guaranteed to be unnecessary, but the changes were
      trivial.
***/

#include "bind.h"

#ifndef SYM_OPS_DEFINED
#include "symmetry.h"
#endif

#ifndef EHT_PARM_FILE
#define EHT_PARM_FILE "eht_parms.dat"
#endif

/* hopefully this will be way more than enough */
#define MAX_CUSTOM_ATOMS 40
atom_type custom_atoms[MAX_CUSTOM_ATOMS];



/****************************************************************************
 *
 *                   Procedure read_geom_frag
 *
 * Arguments: infile: pointer to type FILE
 *         geom_frag: pointer to geom_frag_type
 *
 * Returns: none
 *
 * Action:  This reads in the coordinates of a geometry fragment
 *           from 'infile.  It does *not* take care of evaluating
 *           z matrix positions or anything of the sort.
 *
 *****************************************************************************/
void read_geom_frag(FILE *infile,geom_frag_type *geom_frag)
{
  int i;
  int foo_int,which;
  char instring[400];
  int *values;
  int num_values;
  geom_frag_type *new_frag;

  /* read out the number(s) of the atom(s) this fragment will replace */
  skipcomments(infile,instring,FATAL);
  parse_integer_string(instring,&values,&num_values);
  geom_frag->which = values[0]-1;

  /* now get the number of atoms */
  skipcomments(infile,instring,FATAL);
  sscanf(instring,"%d",&(geom_frag->num_atoms));

  /* allocate space for the atoms */
  geom_frag->atoms = (atom_type *)calloc(geom_frag->num_atoms,sizeof(atom_type));
  if(!geom_frag->atoms){
    fatal("Can't allocate memory for geom_frag atoms.");
  }

  /* read'em out! */
  for(i=0;i<geom_frag->num_atoms;i++){
    skipcomments(infile,instring,FATAL);
    /*  figure out which atom this is */
    sscanf(instring,"%d",&which);
    which--;

    /* a modicum of error checking */
    if( which < 0 ){
      fatal("Atom number specified which is less than zero.");
    }

    if( which >= geom_frag->num_atoms ){
      fatal("Atom number specified which is larger than num_atoms.");
    }

    if( !geom_frag->using_Z_mat ){
      sscanf(instring,"%d %s  %lf %lf %lf",&foo_int,
             geom_frag->atoms[which].symb,
             &(geom_frag->atoms[which].loc.x),&(geom_frag->atoms[which].loc.y),
             &(geom_frag->atoms[which].loc.z));
    }
    else{
      sscanf(instring,"%d %s %d %lf %d %lf %d %lf",&foo_int,
             geom_frag->atoms[which].symb,
             &(geom_frag->atoms[which].Zmat_loc.ref1),
             &(geom_frag->atoms[which].Zmat_loc.bond_length),
             &(geom_frag->atoms[which].Zmat_loc.ref2),
             &(geom_frag->atoms[which].Zmat_loc.alpha),
             &(geom_frag->atoms[which].Zmat_loc.ref3),
             &(geom_frag->atoms[which].Zmat_loc.beta));
      geom_frag->atoms[which].Zmat_loc.ref1--;
      geom_frag->atoms[which].Zmat_loc.ref2--;
      geom_frag->atoms[which].Zmat_loc.ref3--;
    }
    geom_frag->atoms[which].symb[2] = 0;

    /*********
      'which_atom is used to store where the atom was positioned in the
      geometry specification.  This is necessary in order to be able
      to correctly construct the Z matrix (if it's being used)
      *********/
    geom_frag->atoms[which].which_atom = i;
  }

  /* now copy the specification into geom_frags for other replacements */
  for(i=1;i<num_values;i++){
    new_frag = (geom_frag_type *)calloc(1,sizeof(geom_frag_type));
    if( !new_frag ) fatal("Can't get memory for new_frag in read_geom_frag");
    bcopy(geom_frag,new_frag,sizeof(geom_frag_type));

    new_frag->which = values[i]-1;
    new_frag->atoms = (atom_type *)calloc(geom_frag->num_atoms,
                                          sizeof(atom_type));
    if(!new_frag->atoms)fatal("can't get new_frag->atoms in read_geom_frag");
    bcopy(geom_frag->atoms,new_frag->atoms,geom_frag->num_atoms*sizeof(atom_type));

    geom_frag->next = new_frag;
    geom_frag = new_frag;
  }
}







/****************************************************************************
 *
 *                   Procedure write_atom_coords
 *
 * Arguments:  atoms: pointer to atom_type
 *         num_atoms: int
 *        using_Zmat: char
 *        using_xtal_coords: char
 *
 * Returns: none
 *
 * Action:
 *   writes all the atomic positions to the output file.
 *
 *****************************************************************************/
void write_atom_coords(atom_type *atoms,int num_atoms,
                       char using_Zmat,char using_xtal_coords)
{
  int i,j;
  char found_this_one;

  fprintf(output_file,"; ********* Atoms within the unit cell:  *********\n");
  fprintf(output_file,"# NUMBER OF ATOMS: \n\t%d\n",num_atoms);
  if( using_Zmat ){
    fprintf(output_file,"# Z-Matrix orientation of atoms:\n");
  }
  else if(using_xtal_coords){
    fprintf(output_file,"# Crystallographic coordinates of atoms:\n");
  }
  else{
    fprintf(output_file,"# ATOMIC POSITIONS\n");
  }

  for(i=0;i<num_atoms;i++){
    /* write the parameters */
    if( !using_Zmat ){
      fprintf(output_file,"% 4d % 4s % -8.6f % -8.6f % -8.6f\n",
              i+1,atoms[i].symb,atoms[i].loc.x,atoms[i].loc.y,atoms[i].loc.z);
    }
    else{
      fprintf(output_file,"% 4d % 4s % 4d % -8.6f % 4d % -8.6f % 4d % 8.6f\n",
              i+1,atoms[i].symb,
              atoms[i].Zmat_loc.ref1+1,atoms[i].Zmat_loc.bond_length,
              atoms[i].Zmat_loc.ref2+1,atoms[i].Zmat_loc.alpha,
              atoms[i].Zmat_loc.ref3+1,atoms[i].Zmat_loc.beta);
    }

  }
}

/****************************************************************************
 *
 *                   Procedure write_atom_parms
 *
 * Arguments: details: pointer to detail_type
 *              atoms: pointer to atom_type
 *          num_atoms: int
 *         print_them: char
 *
 * Returns: none
 *
 * Action:
 *   writes all the atomic parms to the output file.
 *
 *****************************************************************************/
void write_atom_parms(detail_type *details,atom_type *atoms,int num_atoms,
                      char print_them)
{
  static int first_call=1,max_num_unique=0;
  atom_type *atom;
  int i,j;
  char found_this_one;

  if( first_call ){
    unique_atoms = (atom_type *)calloc(num_atoms,sizeof(atom_type));
    if(!unique_atoms)fatal("Can't allocate unique atom list");
    num_unique_atoms = 0;
    max_num_unique = num_atoms;
  }

  /********

    on the first call we want to find the unique atoms.

    if we're doing either zeta variation or charge iteration,
    this means any atom which is being varied.  That way we
    print out the values of atoms which are unique, but don't print
    out twenty billion things that are the same.

    **********/
  for( i=0; i<num_atoms; i++ ){
    found_this_one = 0;
    if( !details->vary_zeta &&
       !(details->do_chg_it || atoms[i].chg_it_vary)  &&
       !(details->do_muller_it || atoms[i].chg_it_vary) ){
      for(j=0;j<num_unique_atoms&&!found_this_one;j++){
        if(!strncmp(atoms[i].symb,unique_atoms[j].symb,2)) found_this_one=1;
      }
    }
    /* if we didn't find it... */
    if( !found_this_one ){
      /* just copy the data into the unique_atoms array */
      bcopy(&(atoms[i]),&(unique_atoms[num_unique_atoms]),
            sizeof(atom_type));
      num_unique_atoms++;
      if( num_unique_atoms == max_num_unique ){
        max_num_unique += num_atoms;
        unique_atoms = (atom_type *)my_realloc((int *)unique_atoms,
                                            max_num_unique*sizeof(atom_type));
        if( !unique_atoms ) fatal("Can't realloc unique_atoms.");
      }
    }
    first_call = 0;
  }

  if( print_them ){
    /* now write out Hueckel parameters for the unique atoms */
    fprintf(output_file,"\n\n# ******** Extended Hueckel Parameters ********\n");
    fprintf(output_file,";  FORMAT  quantum number orbital: Hii, <c1>, exponent1,\
 <c2>, <exponent2>\n\n");

    for(i=0;i<num_unique_atoms;i++){
      if( details->vary_zeta ||
         ((details->do_chg_it || details->do_muller_it) &&
          atoms[unique_atoms[i].which_atom].chg_it_vary)){

        fprintf(output_file,
                "ATOM: %s(%d)   Atomic number: %d  # Valence Electrons: %d\n",
                unique_atoms[i].symb,unique_atoms[i].which_atom+1,
                unique_atoms[i].at_number,
                unique_atoms[i].num_valence);
        atom = &atoms[unique_atoms[i].which_atom];
      }else{
        fprintf(output_file,
                "ATOM: %s   Atomic number: %d  # Valence Electrons: %d\n",
                unique_atoms[i].symb,unique_atoms[i].at_number,
                unique_atoms[i].num_valence);
        atom = &unique_atoms[i];
      }
      if( atom->ns != 0){
        fprintf(output_file,"\t%dS: %10.4lf %10.4lf\n",atom->ns,
                atom->coul_s,atom->exp_s);
      }
      if( atom->np != 0){
        fprintf(output_file,"\t%dP: %10.4lf %10.4lf\n",atom->np,
                atom->coul_p,atom->exp_p);
      }
      if( atom->nd != 0){
        fprintf(output_file,"\t%dD: %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf\n",
                atom->nd,
                atom->coul_d,atom->coeff_d1,
                atom->exp_d,atom->coeff_d2,
                atom->exp_d2);
      }
      /* parameters for an atom with f orbitals */

      if( atom->nf != 0){
        fprintf(output_file,"\t%dF: %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf\n",
                atom->nf,
                atom->coul_f,atom->coeff_f1,
                atom->exp_f,atom->coeff_f2,
                atom->exp_f2);
      }
    }
  }
}


/****************************************************************************
 *
 *                   Procedure fill_chg_it_parms
 *
 * Arguments:  atoms: pointer to atom_type
 *         num_atoms: int
 *         num_lines: int
 *            infile: pointer to file type
 *
 * Returns: none
 *
 * Action:
 *   sets the atoms up with charge iteration parameters
 *   atomic wavefunction parameters not supplied by the input file are read
 *   out of a parameter file with the name indicated in the #define for
 *   EHT_PARM_FILE given above (or in the makefile)
 *
 *****************************************************************************/
void fill_chg_it_parms(atoms,num_atoms,num_lines,infile)
  atom_type *atoms;
  int num_atoms,num_lines;
  FILE *infile;
{
  int i,j,num_read;
  char symb[10],instring[240];
  real s_A,s_B,s_C,p_A,p_B,p_C,d_A,d_B,d_C;

  for(i=0;i<num_lines;i++){
    skipcomments(infile,instring,FATAL);
    upcase(instring);
    num_read = sscanf(instring,"%s %lf %lf %lf %lf %lf %lf %lf %lf %lf",symb,
                      &s_A,&s_B,&s_C,&p_A,&p_B,&p_C,&d_A,&d_B,&d_C);
    /* check to see how much info we actually read in */
    if( num_read < 7 ){
      p_A = p_B = p_C = 0.0;
    }
    if( num_read < 10 ){
      d_A = d_B = d_C = 0.0;
    }

    /* now loop over all the atoms and pop these parameters in */
    for( j=0;j<num_atoms;j++){
      upcase(atoms[j].symb);

      if( !strncmp(atoms[j].symb,symb,2) ){
        atoms[j].s_A = s_A; atoms[j].s_B = s_B, atoms[j].s_C = s_C;
        atoms[j].p_A = p_A; atoms[j].p_B = p_B, atoms[j].p_C = p_C;
        atoms[j].d_A = d_A; atoms[j].d_B = d_B, atoms[j].d_C = d_C;
      }
    }
  }
}

/****************************************************************************
 *
 *                   Procedure fill_atomic_parms
 *
 * Arguments:  atoms: pointer to atom_type
 *         num_atoms: int
 *            infile: pointer to file type
 *
 * Returns: none
 *
 * Action:
 *   sets all the atoms up with extended hueckel parameters
 *   atomic wavefunction parameters not supplied by the input file are read
 *   out of a parameter file with the name indicated in the #define for
 *   EHT_PARM_FILE given above (or in the makefile)
 *
 *****************************************************************************/
void fill_atomic_parms(atoms,num_atoms,infile)
  atom_type *atoms;
  int num_atoms;
  FILE *infile;
{
  char err_string[240],instring[80];
  char *parm_file_name;
  point_type saveloc;
  Z_mat_type saveZloc;
  FILE *parmfile;
  int i,j;
  int num_read;
  int save_which;
  char found;
  real temp;

  int num_custom;

  char ang[2],tstring[10];
  int atnum,nzeta,nquant,nval;
  real Hii,exp1,exp2,c1,c2;



  /* open the parameter file */
#ifndef USING_THE_MAC
  parm_file_name = (char *)getenv("BIND_PARM_FILE");
#else
  parm_file_name = 0;
#endif
  if( !parm_file_name ){
    parm_file_name = (char *)calloc(240,sizeof(char));
    if(!parm_file_name)fatal("can't get memory for parm_file_name");
    safe_strcpy(parm_file_name,EHT_PARM_FILE);
  }
  parmfile = fopen(parm_file_name,"r");

  bzero(custom_atoms,MAX_CUSTOM_ATOMS*sizeof(atom_type));

  /* make sure that it opened, but don't exit if not... */
#ifndef USING_THE_MAC
  if(!parmfile){
    safe_strcpy(err_string,"Can't open parameter file: ");
    strcat(err_string,parm_file_name);
    strcat(err_string," continuing anyway <watch out!>.");
    error(err_string);
  }
#else
#include "Mac_Fopen.h"
  if( !parmfile ){
          error("Can't open parm file, please specify a name");
          parmfile = choose_mac_file(parm_file_name,MAC_FOPEN_OPEN_NOCD);
          if( !parmfile ) fatal("Still can't do it.");
  }
#endif
  /* loop over the atoms and get the parameters */
  num_custom = 0;
  for(i=0;i<num_atoms;i++){

#if 0
    /******
      if the first character in the symbol is a space, then move
      the second character to the first position, this makes searching
      for parameters much easier
      ******/
    if( atoms[i].symb[0] == ' ' || atoms[i].symb[0] == '\t'){
      atoms[i].symb[0] = atoms[i].symb[1];
      atoms[i].symb[1] = 0;
    }
    /* remove tab characters and spaces */
    if(atoms[i].symb[1] == '\t' || atoms[i].symb[1] == ' '){
      atoms[i].symb[1] = 0;
    }
    /* remove carriage returns */
    if(atoms[i].symb[1] == '\n'){
      atoms[i].symb[1] = 0;
    }

#endif
    /*******
      if it's a dummy atom, then we don't have to really do anything
      ********/
    if( atoms[i].symb[0] == '&' ){
      atoms[i].at_number = -1;
    }


    /*******
      if it's a special atom appearing for the first time, then get
      its parameters from the input file
      ********/
    else if(atoms[i].symb[0] == '*'){
      skipcomments(infile,instring,FATAL);
      num_read = sscanf(instring,"%s %d %d %d %lf %lf %d %lf %lf %d %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf",
                        custom_atoms[num_custom].symb,
                        &(custom_atoms[num_custom].at_number),
                        &(custom_atoms[num_custom].num_valence),
                        &(custom_atoms[num_custom].ns),&(custom_atoms[num_custom].exp_s),
                        &(custom_atoms[num_custom].coul_s),
                        &(custom_atoms[num_custom].np),&(custom_atoms[num_custom].exp_p),
                        &(custom_atoms[num_custom].coul_p),
                        &(custom_atoms[num_custom].nd),&(custom_atoms[num_custom].exp_d),
                        &(custom_atoms[num_custom].coul_d),
                        &(custom_atoms[num_custom].coeff_d1),
                        &(custom_atoms[num_custom].exp_d2),
                        &(custom_atoms[num_custom].coeff_d2),
                        &(custom_atoms[num_custom].nf),&(custom_atoms[num_custom].exp_f),
                        &(custom_atoms[num_custom].coul_f),
                        &(custom_atoms[num_custom].coeff_f1),
                        &(custom_atoms[num_custom].exp_f2),
                        &(custom_atoms[num_custom].coeff_f2));

#if 0
      /* patch up the custom name */
      if(custom_atoms[num_custom].symb[0]==' '||custom_atoms[num_custom].symb[0] =='\t'){
        custom_atoms[num_custom].symb[0] = custom_atoms[num_custom].symb[1];
        custom_atoms[num_custom].symb[1] = 0;
      }
      if(custom_atoms[num_custom].symb[1] == '\t'){
        custom_atoms[num_custom].symb[1] = 0;
      }

#endif
      upcase(custom_atoms[num_custom].symb);
      /***********
        we're not guaranteed to have gotten information on all the orbitals,
        (the user probably won't have anything for f orbitals when giving
        alternate parameters for hydrogen...).  To make sure that the
        program doesn't have bogus values for unspecified parameters, check
        to see how many we actually got out of the file.
        ***********/
      switch(num_read){
      case 6:
        custom_atoms[num_custom].np=0;
        custom_atoms[num_custom].nd=0;
        custom_atoms[num_custom].nf=0;
        break;
      case 9:
        custom_atoms[num_custom].nd=0;
        custom_atoms[num_custom].nf=0;
        break;
      case 15:
        custom_atoms[num_custom].nf=0;
        break;
      }

      /*******
        Normalize the d-coefficients.
        ******/
      if(custom_atoms[num_custom].coeff_d2 != 0){
        temp = 4.0*(custom_atoms[num_custom].exp_d*custom_atoms[num_custom].exp_d2/
                  pow(custom_atoms[num_custom].exp_d+custom_atoms[num_custom].exp_d2,2.0));
        temp = pow(temp,((real)custom_atoms[num_custom].nd+.5));

        temp = sqrt(custom_atoms[num_custom].coeff_d1*custom_atoms[num_custom].coeff_d1+
                    custom_atoms[num_custom].coeff_d2*custom_atoms[num_custom].coeff_d2+
                    2.0*temp*custom_atoms[num_custom].coeff_d1*
                    custom_atoms[num_custom].coeff_d2);
        temp = 1.0/temp;

        custom_atoms[num_custom].coeff_d1 *= temp;
        custom_atoms[num_custom].coeff_d2 *= temp;
      }

      /*******
        Normalize the f-coefficients.
        ******/
      if(custom_atoms[num_custom].coeff_f2 != 0){
        temp = 4.0*(custom_atoms[num_custom].exp_f*custom_atoms[num_custom].exp_f2/
                  pow(custom_atoms[num_custom].exp_f+custom_atoms[num_custom].exp_f2,2.0));
        temp = pow(temp,((real)custom_atoms[num_custom].nf+.5));

        temp = sqrt(custom_atoms[num_custom].coeff_f1*custom_atoms[num_custom].coeff_f1+
                    custom_atoms[num_custom].coeff_f2*custom_atoms[num_custom].coeff_f2+
                    2.0*temp*custom_atoms[num_custom].coeff_f1*
                    custom_atoms[num_custom].coeff_f2);
        temp = 1.0/temp;

        custom_atoms[num_custom].coeff_f1 *= temp;
        custom_atoms[num_custom].coeff_f2 *= temp;
      }


      /************

        make sure that the atomic list has these parameters
        since I copy the whole block at once, I need to save the
        atomic position and then re-insert it into the atom structure

        **************/
      bcopy((char *)&(atoms[i].loc),(char *)&(saveloc),sizeof(point_type));
      bcopy((char *)&(atoms[i].Zmat_loc),(char *)&saveZloc,sizeof(Z_mat_type));
      save_which = atoms[i].which_atom;
      bzero((char *)&(atoms[i]),sizeof(atom_type));
      bcopy((char *)&(custom_atoms[num_custom]),(char *)&(atoms[i]),
            sizeof(atom_type));

      atoms[i].which_atom = save_which;
      bcopy((char *)&saveloc,(char *)&(atoms[i].loc),sizeof(point_type));
      bcopy((char *)&saveZloc,(char *)&(atoms[i].Zmat_loc),sizeof(Z_mat_type));

      num_custom++;
    }
    else{
      found = 0;

      /*****
        check to see if this is a custom atom that has already been hit
        ******/
          upcase(atoms[i].symb);
      for(j=0;j<num_custom;j++){

                upcase(custom_atoms[j].symb);
        if(!strncmp(atoms[i].symb,custom_atoms[j].symb,2)){
          /* BINGO! copy the data and save the location. */
          bcopy((char *)&(atoms[i].loc),(char *)&(saveloc),sizeof(point_type));
          bcopy((char *)&(atoms[i].Zmat_loc),(char *)&saveZloc,sizeof(Z_mat_type));
          save_which = atoms[i].which_atom;

          bcopy((char *)&(custom_atoms[j]),(char *)&(atoms[i]),
                sizeof(atom_type));
          atoms[i].which_atom = save_which;
          bcopy((char *)&saveloc,(char *)&(atoms[i].loc),sizeof(point_type));
          bcopy((char *)&saveZloc,(char *)&(atoms[i].Zmat_loc),sizeof(Z_mat_type));

          found = 1;
        }
      }

      /******
        look for the parameters in the param file
        *******/
      rewind(parmfile);

      while(!found && skipcomments(parmfile,instring,IGNORE)>=0 ){
        /* compare two characters and see if this is the right atom */
        sscanf(instring,"%s",tstring);
        upcase(tstring);
        while(!strncmp(atoms[i].symb,tstring,2)){
          /* it is... read out the data */
          sscanf(instring,"%s %d %d %d %d %s %lf %lf %lf %lf %lf",
                 tstring,&atnum,&nval,&nzeta,&nquant,ang,&Hii,&exp1,&exp2,&c1,&c2);

          atoms[i].at_number = atnum;
          atoms[i].num_valence = nval;
          /* figure out which type of orbital this is */
          switch(ang[0]){
          case 's':
          case 'S':
            if( Hii != 0.0 )
              atoms[i].ns = nquant;
            else atoms[i].ns = 0;

            atoms[i].exp_s = exp1;
            atoms[i].coul_s = Hii;
            break;
          case 'p':
          case 'P':
            if( Hii != 0.0 )
              atoms[i].np = nquant;
            else atoms[i].np = 0;
            atoms[i].exp_p = exp1;
            atoms[i].coul_p = Hii;
            break;
          case 'd':
          case 'D':
            if( Hii != 0.0 )
              atoms[i].nd = nquant;
            else atoms[i].nd = 0;
            atoms[i].exp_d = exp1;
            atoms[i].coul_d = Hii;
            atoms[i].coeff_d1 = c1;
            if( nzeta == 2 ){
              atoms[i].exp_d2 = exp2;
              atoms[i].coeff_d2 = c2;

              /*******
                this is some kind of wierd coefficient adjustment business that
                they do in the original source... I'm not sure why...
                ******/
              temp = 4.0*(exp1*exp2/pow(exp1+exp2,2.0));
              temp = pow(temp,(real)nquant+.5);

              temp = sqrt(c1*c1+c2*c2+2*temp*c1*c2);
              temp = 1.0/temp;

              atoms[i].coeff_d1 *= temp;
              atoms[i].coeff_d2 *= temp;
            }
            break;
          case 'f':
          case 'F':
            if( Hii != 0.0 )
              atoms[i].nf = nquant;
            else atoms[i].nf = 0;
            atoms[i].exp_f = exp1;
            atoms[i].coul_f = Hii;
            atoms[i].coeff_f1 = c1;
            if( nzeta == 2 ){
              atoms[i].exp_f2 = exp2;
              atoms[i].coeff_f2 = c2;

              /*******
                this is some kind of wierd coefficient adjustment business that
                they do in the original source... I'm not sure why...
                ******/
              temp = 4.0*(exp1*exp2/pow(exp1+exp2,2.0));
              temp = pow(temp,(real)nquant+.5);

              temp = sqrt(c1*c1+c2*c2+2*temp*c1*c2);
              temp = 1.0/temp;

              atoms[i].coeff_f1 *= temp;
              atoms[i].coeff_f2 *= temp;
            }
            break;
          }
          /********
            now read in the next line to make sure that we have all the oribtals
            for this atom.
            ********/
          skipcomments(parmfile,instring,IGNORE);
          sscanf(instring,"%s",tstring);
          found = 1;
        }
      }
      if(!found){
        fprintf(stderr,"Can't find parameters for atom: %s\n",atoms[i].symb);
        fatal("Parameter acquisition failure :-(");
      }
    }
  }
  if( parmfile ) fclose(parmfile);
}



/****************************************************************************
 *
 *                   Procedure parse_printing_options
 *
 * Arguments: infile: pointer to type FILE
 *           details: pointer to detail_type
 *              cell: pointer to cell_type
 *
 * Returns: none
 *
 * Action:  This parses all the printing options that were given...
 *
 *****************************************************************************/
void parse_printing_options(infile,details,cell)
  FILE *infile;
  detail_type *details;
  cell_type *cell;
{
  char instring[240];
  char type_string[40];
  BOOLEAN *which_option;
  int EOF_hit;
  int which;
  printing_info_type *print_info;

  /**********

    NOTE: when adding keywords to this, remember that strstr searches for
    substrings.  Therefore it is important that longer strings come before
    shorter strings with a common substring.
    e.g. REDUCED OVERLAP POP  must come before OVERLAP POP in the if..else
    construct or else the matching will be done incorrectly
    (strstr("REDUCED OVERLAP POP", "OVERLAP POP") returns a value, whereas:
    strstr("OVERLAP POP","REDUCED OVERLAP POP") returns zero.

    **********/

  /* read until we hit either the EOF or the keyword END_PRINT */
  EOF_hit = 0;
  while( EOF_hit > -1 ){
    EOF_hit = skipcomments(infile,instring,IGNORE);
    upcase(instring);
    if( EOF_hit > -1 ){

      /*----------------------------------------------------------------------*/
      if(strstr(instring,"END_PRINT") || strstr(instring,"END PRINT")){
        EOF_hit = -1;
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"DIST")){
        details->distance_mat_PRT = 1;
        which_option = &(details->distance_mat_PRT);
        which = PRT_DIST;
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"AVERAGE REDUCED OVERLAP POP")){
        details->avg_ROP_mat_PRT = 1;
        which_option = &(details->avg_ROP_mat_PRT);
        which = PRT_AVG_ROP;
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"AVERAGE OVERLAP POP")){
        details->avg_OP_mat_PRT = 1;
        which_option = &(details->avg_OP_mat_PRT);
        which = PRT_AVG_OP;
      }
      else if(strstr(instring,"MOD REDUCED OVERL")){
        details->mod_ROP_mat_PRT = 1;
        which_option = &(details->mod_ROP_mat_PRT);
        which = PRT_ROP;
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"MOD OVERLAP POP")){
        details->mod_OP_mat_PRT = 1;
        which_option = &(details->mod_OP_mat_PRT);
        which = PRT_OP;
      }

      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"REDUCED OVERLAP POP")){
        details->ROP_mat_PRT = 1;
        which_option = &(details->ROP_mat_PRT);
        which = PRT_ROP;
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"OVERLAP POP")){
        details->OP_mat_PRT = 1;
        which_option = &(details->OP_mat_PRT);
        which = PRT_OP;
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"OVERLAP")){
        details->overlap_mat_PRT = 1;
        which_option = &(details->overlap_mat_PRT);
        which = PRT_OVERLAP;
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"REDUCED CHARGE MAT")){
        details->Rchg_mat_PRT = 1;
        which_option = &(details->Rchg_mat_PRT);
        which = PRT_RCHG_MAT;
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"CHARGE MAT")){
        details->chg_mat_PRT = 1;
        which_option = &(details->chg_mat_PRT);
        which = PRT_CHG_MAT;
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"WAVE FUNC") || strstr(instring,"WAVEFUNC")){
        details->wave_fn_PRT = 1;
        which_option = &(details->wave_fn_PRT);
        which = PRT_WAVE_FUNC;
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"MOD NET CHARGE")){
        details->mod_net_chg_PRT = 1;
        which_option = &(details->mod_net_chg_PRT);
        which = PRT_NET_CHG;
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"NET CHARGE")){
        details->net_chg_PRT = 1;
        which_option = &(details->net_chg_PRT);
        which = PRT_NET_CHG;
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"ELECTROSTAT")){
        details->electrostat_PRT = 1;
        which_option = &(details->electrostat_PRT);
        which = PRT_ELECTROSTAT;
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"FERMI")){
        details->fermi_PRT = 1;
        which_option = &(details->fermi_PRT);
        which = PRT_FERMI;
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"LEVELS")){
        details->levels_PRT = 1;
        which_option = &(details->levels_PRT);
        which = PRT_LEVELS;
      }

      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"HAMIL")){
        details->hamil_PRT = 1;
        which_option = &(details->hamil_PRT);
        which = PRT_HAMIL;
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"ORBITAL") && strstr(instring,"ENERG")){
        details->energies_PRT = 1;
        which_option = &(details->energies_PRT);
        which = PRT_ORB_ENERGY;
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"ORBITAL") && strstr(instring,"COEFF")){
        details->energies_PRT = 1;
        which_option = &(details->energies_PRT);
        which = PRT_ORB_COEFF;
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"ENERG")){
        details->energies_PRT = 1;
        which_option = &(details->energies_PRT);
        which = PRT_ENERGIES;
      }
      else if( strstr(instring,"ORBITAL MAP") ){
        details->orbital_mapping_PRT = 1;
      }
      else{
        fprintf(stderr,"Invalid Printing Option: %s in input file.\n",
                instring);
        fprintf(stderr,"Did you forget to include end_print?\n");
        error("Bad print option.");
      }

      /********

        Okay, now we've parsed the keyword. Check to see if there are
        any other options

        ********/
      if(strstr(instring,"WALSH")){
        /* they want to follow this variable along the Walsh diagram... */

        /* get space for the printing option */
        print_info = (printing_info_type *)calloc(1,sizeof(printing_info_type));
        if( !print_info ) fatal("Can't get memory to store printing options.");

        /* put this structure at the head of the printing options linked list */
        print_info->next = details->step_print_options;
        details->step_print_options = print_info;

        /* we know the kind of the printing option already (it was set above) */
        print_info->which_to_print = which;

        /* now read out the parameters */
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%s %d %d",type_string,&(print_info->contrib1),
               &(print_info->contrib2));

        upcase(type_string);

        /* account for the C indexing thing */
        print_info->contrib1--;
        print_info->contrib2--;


        /*******
          some printing options only take one argument, to make sure these don't alert
           the error checking, adjust those now.
          ********/
        if( which == PRT_NET_CHG ){
          print_info->contrib2 = 0;
        }
        else if(which == PRT_ELECTROSTAT ){
          print_info->contrib1 = 0;
          print_info->contrib2 = 0;
        }

        /* figure out the type */
        if(strstr(type_string,"ORB")){
          print_info->type = P_DOS_ORB;
          /******
            error checking
            *******/
          if( print_info->contrib1 < 0 || print_info->contrib2 < 0 ){
            fatal("Projected orbital number less than 1 in walsh printing region.");
          }
          /*****
            unfortunately, it's not yet possible to check to see if they have projected
            out an invalid orbital, since we don't yet know how many orbitals there are.
            this will have to be checked later.
            *****/
        }
        else if(strstr(type_string,"ATOM")){
          print_info->type = P_DOS_ATOM;

          /******
            error checking
            *******/
          if( print_info->contrib1 < 0 || print_info->contrib2 < 0 ){
            fatal("Projected atom number less than 1 in walsh printing region.");
          }
          if( print_info->contrib1 > cell->num_atoms-1 ||
             print_info->contrib2 > cell->num_atoms-1 ){
            fatal("Projected atom number greater than num_atoms in walsh printing region.");
          }
        }
        else fatal("Invalid type specified for a printing option.");
      } else if( strstr(instring,"TRANSP") ){
        /* they want the transpose of the quantity they are printing */
        *which_option = *which_option | PRT_TRANSPOSE_FLAG;
      }
    }
  }
}


/****************************************************************************
 *
 *                   Procedure parse_equiv_atoms
 *
 * Arguments: infile: pointer to type FILE
 *           details: pointer to detail_type
 *              cell: pointer to cell_type
 *
 * Returns: none
 *
 * Action:  This parses the equivalent atom list given by the user
 *
 *****************************************************************************/
void parse_equiv_atoms(infile,details,cell)
  FILE *infile;
  detail_type *details;
  cell_type *cell;
{
  char instring[400];
  int num_equiv;
  int num_read,*values_read=0;
  int i,j;
  equiv_atom_type *equiv_list;

  /* read out the number of lines giving lists of equivalent atoms */
  skipcomments(infile,instring,FATAL);
  sscanf(instring,"%d",&num_equiv);

  for(i=0;i<num_equiv;i++){
    /* get space for the new equiv atoms list entry */
    equiv_list = (equiv_atom_type *)calloc(1,sizeof(equiv_atom_type));
    if( !equiv_list ) fatal("Can't get memory for equiv_list entry");

    /* insert it in the list */
    equiv_list->next = cell->equiv_atoms;
    cell->equiv_atoms = equiv_list;

    skipcomments(infile,instring,FATAL);
    parse_integer_string(instring,&values_read,&num_read);
    for(j=0;j<num_read;j++) values_read[j]--;

    /* get space to store the equivalent atoms */
    equiv_list->num_equiv = num_read;
    equiv_list->equiv_atoms = (int *)calloc(num_read,sizeof(int));
    if( !equiv_list->equiv_atoms ) fatal("Can't get memory for equiv_atoms array");

    /* now add the entries to the array */
    bcopy(values_read,equiv_list->equiv_atoms,num_read*sizeof(int));

    /* that's that. */
  }

#ifdef DEBUG
  /* write out the equivalent atoms */
  fprintf(output_file,"; Equivalent atoms for Muller iteration\n");
  i=1;
  equiv_list = cell->equiv_atoms;
  while(equiv_list){
    fprintf(output_file,"%d ",i++);
    for(j=0;j<equiv_list->num_equiv;j++){
      fprintf(output_file,"%d,",equiv_list->equiv_atoms[j]);
    }
    fprintf(output_file,"\n");
    equiv_list = equiv_list->next;
  }
#endif
}
/****************************************************************************
 *
 *                   Procedure parse_muller_parms
 *
 * Arguments: infile: pointer to type FILE
 *           details: pointer to detail_type
 *              cell: pointer to cell_type
 *
 * Returns: none
 *
 * Action:  This parses the muller iteration parameters that were given.
 *
 *  this is currently set up to read stuff from the file format that
 *    Edgar provided, so some of the ordering is mixed up due to
 *    the change I made in the notation used.
 *
 *****************************************************************************/
void parse_muller_parms(infile,details,cell)
  FILE *infile;
  detail_type *details;
  cell_type *cell;
{
  char instring[400],foostring[80];
  int i,j;
  int num_parms;
  int *numbers_read=0,num_read;
  real E_vals[7],Z_vals[7];
  atom_type *atom,*atom2;

  /* read out the number of params that are going to be given */
  skipcomments(infile,instring,FATAL);
  sscanf(instring,"%d",&num_parms);

  for(i=0;i<num_parms;i++){
    /* read out the line with the atom numbers */
    skipcomments(infile,instring,FATAL);
    /* parse the atom numbers */
    parse_integer_string(instring,&numbers_read,&num_read);

    for(j=0;j<num_read;j++) numbers_read[j]--;

    /* now get the actual parameters */
    if(num_read != 0 ){
      atom = &cell->atoms[numbers_read[0]];
      if(atom->nd != 0 ){
        skipcomments(infile,instring,FATAL);
        if( sscanf(instring,"%s %lf %lf %lf %lf %lf %lf %lf %lf",
               foostring,&atom->init_d_occup,
               &E_vals[0],&E_vals[1],
               &E_vals[2],&E_vals[3],
               &E_vals[4],&E_vals[5],
               &E_vals[6]) != 9 ) fatal("Bad line in Muller parms section");
        skipcomments(infile,instring,FATAL);
        if( sscanf(instring,"%lf %lf %lf %lf %lf %lf %lf",
               &Z_vals[0],&Z_vals[1],
               &Z_vals[2],&Z_vals[3],
               &Z_vals[4],&Z_vals[5],
               &Z_vals[6]) != 7) fatal("Bad line in Muller parms section");

        /* now copy in the parameters in the order I use */
        atom->muller_d_E[0] = E_vals[6];
        atom->muller_d_E[1] = E_vals[3];
        atom->muller_d_E[2] = E_vals[4];
        atom->muller_d_E[3] = E_vals[5];
        atom->muller_d_E[4] = E_vals[2];
        atom->muller_d_E[5] = E_vals[1];
        atom->muller_d_E[6] = E_vals[0];
        atom->muller_d_Z[0] = Z_vals[6];
        atom->muller_d_Z[1] = Z_vals[3];
        atom->muller_d_Z[2] = Z_vals[4];
        atom->muller_d_Z[3] = Z_vals[5];

        /********

          repeat the process for the s and p parameters

        *********/
        skipcomments(infile,instring,FATAL);
        if( sscanf(instring,"%s %lf %lf %lf %lf %lf %lf %lf %lf",
               foostring,&atom->init_s_occup,
               &E_vals[0],&E_vals[1],
               &E_vals[2],&E_vals[3],
               &E_vals[4],&E_vals[5],
               &E_vals[6]) != 9) fatal("Bad line in Muller parms section");
        skipcomments(infile,instring,FATAL);
        if( sscanf(instring,"%lf %lf %lf %lf %lf %lf %lf",
               &Z_vals[0],&Z_vals[1],
               &Z_vals[2],&Z_vals[3],
               &Z_vals[4],&Z_vals[5],
               &Z_vals[6]) != 7 ) fatal("Bad line in Muller parms section");
        atom->muller_s_E[0] = E_vals[6];
        atom->muller_s_E[1] = E_vals[3];
        atom->muller_s_E[2] = E_vals[4];
        atom->muller_s_E[3] = E_vals[5];
        atom->muller_s_E[4] = E_vals[2];
        atom->muller_s_E[5] = E_vals[1];
        atom->muller_s_E[6] = E_vals[0];
        atom->muller_s_Z[0] = Z_vals[6];
        atom->muller_s_Z[1] = Z_vals[3];
        atom->muller_s_Z[2] = Z_vals[4];
        atom->muller_s_Z[3] = Z_vals[5];

        skipcomments(infile,instring,FATAL);
        if(sscanf(instring,"%s %lf %lf %lf %lf %lf %lf %lf %lf",
               foostring,&atom->init_p_occup,
               &E_vals[0],&E_vals[1],
               &E_vals[2],&E_vals[3],
               &E_vals[4],&E_vals[5],
               &E_vals[6])!= 9) fatal("Bad line in Muller parms section");
        skipcomments(infile,instring,FATAL);
        if( sscanf(instring,"%lf %lf %lf %lf %lf %lf %lf",
               &Z_vals[0],&Z_vals[1],
               &Z_vals[2],&Z_vals[3],
               &Z_vals[4],&Z_vals[5],
               &Z_vals[6]) != 7) fatal("Bad line in Muller parms section");
        atom->muller_p_E[0] = E_vals[6];
        atom->muller_p_E[1] = E_vals[3];
        atom->muller_p_E[2] = E_vals[4];
        atom->muller_p_E[3] = E_vals[5];
        atom->muller_p_E[4] = E_vals[2];
        atom->muller_p_E[5] = E_vals[1];
        atom->muller_p_E[6] = E_vals[0];
        atom->muller_p_Z[0] = Z_vals[6];
        atom->muller_p_Z[1] = Z_vals[3];
        atom->muller_p_Z[2] = Z_vals[4];
        atom->muller_p_Z[3] = Z_vals[5];


      }
      else if(atom->np != 0 ){
        skipcomments(infile,instring,FATAL);
         if( sscanf(instring,"%s %lf %lf %lf %lf %lf",
               foostring,&atom->init_s_occup,
               &E_vals[0],&E_vals[1],
               &E_vals[2],&E_vals[3]) != 6) fatal("Bad line in Muller parms section");
        skipcomments(infile,instring,FATAL);
        if( sscanf(instring,"%lf %lf %lf %lf",
               &Z_vals[0],&Z_vals[1],
               &Z_vals[2],&Z_vals[3]) != 4) fatal("Bad line in Muller parms section");
        /* copy them in */
        atom->muller_s_E[0] = E_vals[3];
        atom->muller_s_E[1] = E_vals[1];
        atom->muller_s_E[2] = E_vals[2];
        atom->muller_s_E[3] = E_vals[0];
        atom->muller_s_Z[0] = Z_vals[3];
        atom->muller_s_Z[1] = Z_vals[1];
        atom->muller_s_Z[2] = Z_vals[2];
        skipcomments(infile,instring,FATAL);
         if( sscanf(instring,"%s %lf %lf %lf %lf %lf",
               foostring,&atom->init_p_occup,
               &E_vals[0],&E_vals[1],
               &E_vals[2],&E_vals[3]) != 6)fatal("Bad line in Muller parms section");
        skipcomments(infile,instring,FATAL);
        if( sscanf(instring,"%lf %lf %lf %lf",
               &Z_vals[0],&Z_vals[1],
               &Z_vals[2],&Z_vals[3]) != 4 )fatal("Bad line in Muller parms section");
        /* copy them in */
        atom->muller_p_E[0] = E_vals[3];
        atom->muller_p_E[1] = E_vals[1];
        atom->muller_p_E[2] = E_vals[2];
        atom->muller_p_E[3] = E_vals[0];
        atom->muller_p_Z[0] = Z_vals[3];
        atom->muller_p_Z[1] = Z_vals[1];
        atom->muller_p_Z[2] = Z_vals[2];

      }
      else if(atom->ns != 0 ){
        skipcomments(infile,instring,FATAL);
        if( sscanf(instring,"%s %lf %lf %lf %lf",
               foostring,&atom->init_s_occup,
               &E_vals[0],&E_vals[1],&E_vals[2])!= 5)
          fatal("Bad line in Muller parms section");
        skipcomments(infile,instring,FATAL);
        if( sscanf(instring,"%lf %lf %lf",
                   &Z_vals[0],&Z_vals[1],&Z_vals[2]) != 3)
          fatal("Bad line in Muller parms section");
        atom->muller_s_E[0] = E_vals[2];
        atom->muller_s_E[1] = E_vals[1];
        atom->muller_s_E[2] = E_vals[0];
        atom->muller_s_Z[0] = Z_vals[2];
        atom->muller_s_Z[1] = Z_vals[1];
        atom->muller_s_Z[2] = Z_vals[0];
      }

      /*******

        figure out the initial parameters based on the initial occupations

      *******/
      calc_muller_init_parms(atom);


      /* now copy the parameters into the other atoms of the same type */
      atom->chg_it_vary = 1;
      for(j=1;j<num_read;j++){
        atom2 = &cell->atoms[numbers_read[j]];
        bcopy(atom->muller_s_E,atom2->muller_s_E,7*sizeof(real));
        bcopy(atom->muller_s_Z,atom2->muller_s_Z,4*sizeof(real));
        bcopy(atom->muller_p_E,atom2->muller_p_E,7*sizeof(real));
        bcopy(atom->muller_p_Z,atom2->muller_p_Z,4*sizeof(real));
        bcopy(atom->muller_d_E,atom2->muller_d_E,7*sizeof(real));
        bcopy(atom->muller_d_Z,atom2->muller_d_Z,4*sizeof(real));
        atom2->chg_it_vary = 1;
        atom2->coul_s = atom->coul_s;atom2->exp_s = atom->exp_s;
        atom2->coul_p = atom->coul_p;atom2->exp_p = atom->exp_p;
        atom2->coul_d = atom->coul_d;atom2->exp_d = atom->exp_d;
        atom2->coeff_d1 = atom->coeff_d1;atom2->coeff_d2 = atom->coeff_d2;
        atom2->exp_d2 = atom->exp_d2;
      }
    }else{
      fatal("Atom numbers not specified in muller parms section.");
    }
  }
  fprintf(output_file,"; Muller Iteration being carried out.\n");
  fprintf(output_file,"; Here are the initial parameters:\n");
  write_atom_parms(details,cell->atoms,cell->num_atoms,1);
}
/****************************************************************************
 *
 *                   Procedure parse_charge_iteration
 *
 * Arguments: infile: pointer to type FILE
 *           details: pointer to detail_type
 *              cell: pointer to cell_type
 *
 * Returns: none
 *
 * Action:  This parses all the charge iteration options that were given...
 *
 *****************************************************************************/
void parse_charge_iteration(infile,details,cell)
  FILE *infile;
  detail_type *details;
  cell_type *cell;
{
  char instring[240];
  char tempstring[240];
  BOOLEAN *which_option;
  int i,j;
  int EOF_hit;
  int which;
  int *numbers_read;
  int num_read,num_to_vary;
  int num_parm_lines;
  chg_it_parm_type *parms;


  /*********

      NOTE: before adding keywords to this, make sure to read the
        warning in parse_printing_options.

  ***********/

  parms = &(details->chg_it_parms);
  /* read until we hit either the EOF or the keyword END CHARGE */
  EOF_hit = 0;
  while( EOF_hit > -1 ){
    EOF_hit = skipcomments(infile,instring,IGNORE);
    upcase(instring);
    if( EOF_hit > -1 ){
      /*----------------------------------------------------------------------*/
      if(strstr(instring,"END_CHARGE") || strstr(instring,"END CHARGE")){
        EOF_hit = -1;
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring, "VARY")){
        skipcomments(infile,instring,FATAL);
        parse_integer_string(instring,&numbers_read,&num_read);
        /* loop through and tag each atom that's going to be iterated */
        for(j=0;j<num_read;j++) cell->atoms[numbers_read[j]-1].chg_it_vary = 1;
        free(numbers_read);
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring, "ADJUST")){
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%lf",&parms->adjust);
      }

      /*----------------------------------------------------------------------*/
      else if(strstr(instring, "MAX IT")){
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%d",&parms->max_it);
      }

      /*----------------------------------------------------------------------*/
      else if(strstr(instring, "LAMBDA")){
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%lf",&parms->lambda);
      }

      /*----------------------------------------------------------------------*/
      else if(strstr(instring, "DAMP1")){
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%lf",&parms->damp1);
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring, "DAMP2")){
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%lf",&parms->damp2);
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring, "DAMP3")){
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%lf",&parms->damp3);
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring, "LAMPRI")){
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%lf",&parms->lampri);
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring, "TOLER")){
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%lf",&parms->tolerance);
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring, "VARIABLE")){
        parms->variable_step = 1;
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring, "PARAM")){
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%d",&num_parm_lines);
        fill_chg_it_parms(cell->atoms,cell->num_atoms,num_parm_lines,infile);
      }
      else{
        safe_strcpy(tempstring,"Bad charge iteration section line: ");
        strcat(tempstring,instring);
        error(tempstring);
      }
    }
  }
}





/****************************************************************************
 *
 *                   Procedure read_inputfile
 *
 * Arguments:  cell: pointer to cell_type
 *          details: pointer to detail_type
 *             name: pointer to type char
 *         num_orbs: pointer to int
 *     orbital_lookup_table: pointer to pointer to int
 *
 * Returns: none
 *
 * Action: reads all the data out of the file indicated by 'name...
 *   the information is read into the two structures 'cell and 'details
 *
 *****************************************************************************/
void read_inputfile(cell,details,name,num_orbs,orbital_lookup_table,the_file)
  cell_type *cell;
  detail_type *details;
  char *name;
  int *num_orbs;
  int **orbital_lookup_table;
  FILE *the_file;
{
  walsh_details_type *walsh;
  band_info_type *band_info;
  geom_frag_type *geom_frag;
  FMO_frag_type *FMO_frag;
  char err_string[MAX_STR_LEN];
  char instring[MAX_STR_LEN],tempstring[MAX_STR_LEN];
  char string1[MAX_STR_LEN],string2[MAX_STR_LEN],string3[MAX_STR_LEN];
  char numstring[MAX_STR_LEN];
  FILE *infile;
  k_point_type *points;
  int num_k_points,max_k_points;
  real *real_ptr;
  int *int_ptr,*iarray;
  int reading_atom_nums;
  int idle;
  real ecut,eerr;
  char got_geom;
  char found,got_params,got_muller_params,placed,done;
  char Zmat;
  int i,j,k;
  int itab,temp,found_end;
  int which,foo_int,tot_num_so_far;
  int num_COOPS;
  COOP_type *COOP_ptr,*new_COOP;
  real weight,tot_K_weight;
  int max_p_DOS;
  p_DOS_type *p_DOS;

  real begin_walsh,end_walsh;


  /* open the input file */
  if( !the_file ){
          infile = fopen(name,"r");
  } else{
          infile = the_file;
  }

  /* make sure that it opened */
  if(!infile){
    safe_strcpy(err_string,"Can't open input file: ");
    strcat(err_string,name);
    fatal(err_string);
  }

  /*******
    initialize some variables
    ********/
  details->walsh_details.num_steps = 1;
  details->walsh_details.num_vars = 0;
  details->use_symmetry = 0;
  details->find_princ_axes = 0;
  details->vary_zeta = 0;
  details->avg_props = 0;
  details->just_geom = 0;
  details->dump_overlap = 0;
  details->dump_hamil = 0;
  details->sparsify_value = 0.0;
  got_params = 0;
  got_muller_params = 0;
  got_geom = 0;
  details->Execution_Mode = FAT;
  details->the_const = THE_CONST;
  details->weighted_Hij = 1;
  details->eval_electrostat = 0;
  details->close_nn_contact = NN_DISTANCE;
  details->symm_tol = SYMM_TOL;
  details->muller_mix = MULLER_MIX_DEF;
  details->muller_E_tol = MULLER_E_TOL_DEF;
  details->muller_Z_tol = MULLER_Z_TOL_DEF;
  details->num_moments = 4;
  details->line_width = 80;
  details->k_offset = K_OFFSET;
  cell->equiv_atoms = 0;

  cell->charge = -1000.0;

  /***************

    As of version 2.0, new3 style file i/o is no longer supported, so
    we don't need the not new3! thing anymore... if it's
    there, we'll ignore it.

  ***************/
#ifdef SUPPORT_NEW3_FILEIO
  /* start reading that puppy in! */
  skipcomments(infile,instring,FATAL);


  /*********
    check to see what kind of file this is...
    if the first non-comment line doesn't start out with the string:
    NOT NEW3!
    then we'll assume that it is a new3 input file and use a separate
    file I/O routine.
    **********/
  safe_strcpy(numstring,instring);
  upcase(numstring);
  if( !strstr(numstring,"NOT N") ){
    /******
      this is a new3 input file, therefore this line is the title and
      we should call the new3 input routine.
      ******/
    safe_strcpy(details->title,instring);
    read_NEW3file(cell,details,infile);

  }
  else{
#endif
    /* read the title */
    skipcomments(infile,instring,FATAL);

#ifndef SUPPORT_NEW3_FILEIO
    /* check to see if this line is "Not new3!" */
    safe_strcpy(numstring,instring);
    upcase(numstring);
    if( strstr(numstring,"NOT N") ) skipcomments(infile,instring,FATAL);
#endif

    safe_strcpy(details->title,instring);
    fprintf(output_file,"#JOB_TITLE: %s\n", details->title);

    /**********

      now loop through the keywords until we hit the end of the file.


      NOTE: before adding keywords to this, make sure to read the
        warning in parse_printing_options.

    ***********/
    while(skipcomments(infile,instring,IGNORE)>-1){
      /* convert the string to upper case */
      upcase(instring);

      /*----------------------------------------------------------------------*/
      if(strstr(instring,"MOLECUL")){
        fprintf(status_file,"Doing a molecular calculation.\n");
        details->Execution_Mode = MOLECULAR;
        details->num_KPOINTS = 1;
        details->K_POINTS = (k_point_type *)calloc(1,sizeof(k_point_type));
        if( !details->K_POINTS ) fatal("Can't allocate the single k point.");
        details->K_POINTS[0].weight = 1.0;
      }

      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"THIN")){
        fprintf(status_file,"Doing a THIN mode extended calculation.\n");
        details->Execution_Mode = THIN;
      }

      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"FAT")){
        fprintf(status_file,"Doing a FAT mode extended calculation.\n");
        details->Execution_Mode = FAT;
      }
#ifdef INCLUDE_NETCDF_SUPPORT
      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"NETCDF")){
        fprintf(status_file,"I'll be writing a netCDF file.");
        details->do_netCDF = 1;
      }
#endif
      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"ELECTROSTAT")){
        fprintf(status_file,"Electrostatic contributions to the energy will be\
 evaluated.\n");
        details->eval_electrostat = 1;
      }

      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"NONWEIGHTED")){
        fprintf(status_file,"The nonweighted Hij form will be used.\n");
        details->weighted_Hij = 0;
      }

      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"THE CONST")){
        fprintf(status_file,"Reading in a new value for the constant K.\n");
        if( sscanf(instring,"%s %s %lf",string1,string2,
                   &(details->the_const)) != 3 ){
          skipcomments(infile,instring,FATAL);
          sscanf(instring,"%lf",&(details->the_const));
        }
        fprintf(status_file,"K = %lf\n",details->the_const);
      }

      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"NEAREST NEIGHBOR")){
        if( sscanf(instring,"%s %s %lf",string1,string2,
                   &(details->close_nn_contact)) != 3 ){
          skipcomments(infile,instring,FATAL);
          sscanf(instring,"%lf",&(details->close_nn_contact));
        }
      }

      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"ZETA") ){
        fprintf(status_file,"Self consistant variation of orbital exponents will be \
done.\n");
        details->vary_zeta = 1;
      }

      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"JUST GE") ){
        fprintf(status_file,"Only the geometries will be generated.\n");
        details->just_geom = 1;
      }

      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"AVERAGE PROP") ){
        fprintf(status_file,"Average Propeties will be calculated.\n");
        details->avg_props = 1;
      }

      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"GEOM FRAG") ){
        if( got_geom )
          fatal("You must specify the Geometry keyword after the Geom Frag keyword.");

        if( details->do_muller_it )
          fatal("Muller iteration and Geom Frags are incompatible (for now).");
        if( details->do_chg_it )
          fatal("Charge iteration and Geom Frags are incompatible (for now).");

        /* get memory for the geom_frag and slap it in the list */
        geom_frag = (geom_frag_type *)calloc(1,sizeof(geom_frag_type));
        if( !geom_frag ) fatal("Can't get memory for geom_frag\n");
        geom_frag->next = cell->geom_frags;
        cell->geom_frags = geom_frag;

        if( strstr(instring,"Z MATRIX") ){
          geom_frag->using_Z_mat = 1;
        }

        read_geom_frag(infile,geom_frag);
      }

      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"GEOM")){
        /* check to see if we are doing a Z matrix input file */
        if( strstr(instring,"Z MATRIX") || strstr(instring,"Z-MATRIX")){
          Zmat=1;
          cell->using_Zmat = 1;
        }
        else{
          Zmat=0;
          cell->using_Zmat = 0;

          /* check to see if we're using crystallographic coords */
          if( strstr(instring,"CRYSTALLO")){
            cell->using_xtal_coords = 1;
          }else{
            cell->using_xtal_coords = 0;
          }
        }
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%d",&(cell->num_atoms));

        /* allocate space for the atoms */
        cell->atoms = (atom_type *)calloc(cell->num_atoms,sizeof(atom_type));
        if(!cell->atoms){
          sprintf(err_string,"Can't allocate memory for: %d atoms.",cell->num_atoms);
          fatal("Can't allocate memory for the atoms.");
        }

        skipcomments(infile,instring,FATAL);
        /**************

          check to see if we are going to be reading numbers at the beginning
          of each line

        **************/
        if( !cell->using_Zmat ){
          foo_int = sscanf(instring,"%d %s %lf %lf %lf",&which,
                            cell->atoms[0].symb,&cell->atoms[0].loc.x,
                            &cell->atoms[0].loc.y,&cell->atoms[0].loc.z);
          if( foo_int == 5 ) reading_atom_nums = 1;
          else reading_atom_nums = 0;
        }else{
          reading_atom_nums = 1;
        }

        /*************

          read in the individual atomic locations

        **************/
        for(i=0;i<cell->num_atoms;i++){

          /*  figure out which atom this is */
          if( reading_atom_nums ){
            sscanf(instring,"%d",&which);
            which--;
          } else{
            which = i;
          }

          /* a modicum of error checking */
          if( which < 0 ){
            fatal("Atom number specified which is less than zero.");
          }

          if( which >= cell->num_atoms ){
            fatal("Atom number specified which is larger than num_atoms.");
          }

          if( !Zmat ){
            if( reading_atom_nums ){
              sscanf(instring,"%d %s  %lf %lf %lf",&foo_int,
                     cell->atoms[which].symb,
                     &(cell->atoms[which].loc.x),&(cell->atoms[which].loc.y),
                     &(cell->atoms[which].loc.z));
            } else{
              sscanf(instring,"%s  %lf %lf %lf",
                     cell->atoms[which].symb,
                     &(cell->atoms[which].loc.x),&(cell->atoms[which].loc.y),
                     &(cell->atoms[which].loc.z));
            }
          }
          else{
            sscanf(instring,"%d %s %d %lf %d %lf %d %lf",&foo_int,
                   cell->atoms[which].symb,
                   &(cell->atoms[which].Zmat_loc.ref1),&(cell->atoms[which].Zmat_loc.bond_length),
                   &(cell->atoms[which].Zmat_loc.ref2),&(cell->atoms[which].Zmat_loc.alpha),
                   &(cell->atoms[which].Zmat_loc.ref3),&(cell->atoms[which].Zmat_loc.beta));
            cell->atoms[which].Zmat_loc.ref1--;
            cell->atoms[which].Zmat_loc.ref2--;
            cell->atoms[which].Zmat_loc.ref3--;

          }
          cell->atoms[which].symb[2] = 0;

          /*********
            'which_atom is used to store where the atom was positioned in the
            geometry specification.  This is necessary in order to be able
            to correctly construct the Z matrix (if it's being used)
            *********/
          cell->atoms[which].which_atom = i;

          if(i+1 < cell->num_atoms) skipcomments(infile,instring,FATAL);
        }
        fprintf(status_file,"Read: %d atoms\n",cell->num_atoms);

        cell->num_raw_atoms = cell->num_atoms;

        got_geom = 1;
      } /* end of keyword GEOMETRY */

      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"ELECTR")){
        if( cell->charge != -1000.0 ){
          fprintf(stderr,
                  "Both Electrons and Charge keywords specified.  Charge ignored.\n");
          cell->charge = -1000.0;
        }
        if( sscanf(instring,"%s %lf",string1,&(cell->num_electrons)) < 2 ){
          skipcomments(infile,instring,FATAL);
          sscanf(instring,"%lf",&(cell->num_electrons));
        }
      }

      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"CHARGE ITER") ){
        if( cell->geom_frags )
          fatal("Charge iteration and Geom Frags are incompatible (for now).");

        details->do_chg_it = 1;
        details->chg_it_parms.lambda = .1;
        details->chg_it_parms.damp1 = .1;
        details->chg_it_parms.damp2 = .25;
        details->chg_it_parms.damp3 = 0;
        details->chg_it_parms.lampri = .25;
        details->chg_it_parms.max_it = 100;
        details->avg_props = 1;
        parse_charge_iteration(infile,details,cell);
      }


      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"CHARGE")){
        if( cell->num_electrons != 0.0 ){
          fprintf(stderr,
                  "Both Electrons and Charge keywords specified.  Electrons ignored.\n");
          cell->num_electrons = 0.0;
        }
        if( sscanf(instring,"%s %lf",string1,&(cell->charge)) != 2 ){
          skipcomments(infile,instring,FATAL);
          sscanf(instring,"%lf",&(cell->charge));
        }
        if( got_params ){
          charge_to_num_electrons(cell);
        }
      }

      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"SYMM TOL")){
        if( sscanf(instring,"%s %s %lf",string1,string2,&(details->symm_tol)) != 3 ){
          skipcomments(infile,instring,FATAL);
          sscanf(instring,"%lf",&(details->symm_tol));
        }
        fprintf(status_file,"Symmetry Tolerance reset to %lf\n",
                details->symm_tol);
      }


      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"PRINC")){
        fprintf(status_file,"The principle axes of the molecule will be found.\n");
        details->find_princ_axes = 1;
      }
      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"WALSH")){
        walsh = &(details->walsh_details);
        /*************

          deal with a walsh diagram

          **************/
        /* read in the number of variables and steps */
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%d %d", &(walsh->num_vars),&(walsh->num_steps));

        /* get space to store the variables' values */
        walsh->values = (real *)calloc(walsh->num_vars*walsh->num_steps,
                                        sizeof(real));
        if(!walsh->values) fatal("Can't get memory for Walsh value array.");

        /* now read in the individual values */
        for(i=0;i<walsh->num_vars;i++){
          itab = i*walsh->num_steps;
          skipcomments(infile,instring,FATAL);

          /* check to see if we are doing AUTO-WALSH(tm) */
          if( instring[0] == '!' || instring[1] == '!' ){
            sscanf(instring,"!%lf,%lf",&begin_walsh,&end_walsh);
            auto_walsh(&(walsh->values[itab]),walsh->num_steps,begin_walsh,
                       end_walsh);
          }
          else{
            /* we're doing it the old way... */

            /*  use strtok to get the first comma delimited number */
            safe_strcpy(numstring,(char *)strtok(instring,(const char *)","));
            for(j=0;j<walsh->num_steps;j++){
              /*  use strtok to get the next comma delimited number */
              sscanf(numstring,"%lf",&(walsh->values[itab+j]));
              safe_strcpy(numstring,(char *)strtok(0,","));
            }
          }
        }
        /* put the walsh information in the output file */
        fprintf(output_file,"\n# Walsh Information:\n");
        fprintf(output_file,"%d Variables will be changed with\n",walsh->num_vars);
        fprintf(output_file,"%d steps each.\n\n",walsh->num_steps);

        /* symmetry will always be used for Walsh diagrams, so turn it on now. */
        fprintf(status_file,"Symmetry will be used in the calculations.\n");
        details->use_symmetry = 1;
      }/* end of keyword WALSH */

      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"PARAM") ){
        /* make sure we have the geomtery already */
        if( !got_geom ) fatal("You must specify the Geometry before the Parameters.");

        /*************

          fill in the atomic parameters

          *************/
        fill_atomic_parms(cell->atoms,cell->num_atoms,infile);

        /* do the geom frags, if need be */
        geom_frag = cell->geom_frags;
        while(geom_frag){
          fill_atomic_parms(geom_frag->atoms,geom_frag->num_atoms,infile);
          geom_frag = geom_frag->next;
        }

        if( cell->charge != -1000.0 ){
          charge_to_num_electrons(cell);
        }


        got_params = 1;
        fprintf(status_file,"Completed parameter acquisition.\n");
      } /* end of keyword PARAMETERS */

      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"LATTICE") ){
        if( details->Execution_Mode == MOLECULAR ){
          fatal("Lattice parameters should not be included in molecular calculations.");
        }
        else{
          if(!got_geom) fatal("the Lattice keyword must come after Geometry");

          fprintf(status_file,"Processing lattice parameters.\n");
          skipcomments(infile,instring,FATAL);
          sscanf(instring,"%d",&(cell->dim));

          skipcomments(infile,instring,FATAL);
          sscanf(instring,"%d %d %d",
                 &(cell->overlaps[0]),&(cell->overlaps[1]),&(cell->overlaps[2]));

          /* zero out un-needed overlaps (just in case) */
          if( cell->dim < 3) {
            cell->overlaps[2] = 0;
          }
          if(cell->dim < 2){
            cell->overlaps[1] = 0;
          }
          /* a mote of error checking */
          if(cell->dim < 1){
            fatal("Silly dimensionality given in the lattice parameter section.");
          }

          skipcomments(infile,instring,FATAL);
          sscanf(instring,"%d %d",
                 &(cell->tvects[0].begin),&(cell->tvects[0].end));

          if(cell->tvects[0].end <= cell->num_atoms - cell->dim )
            fatal("Please make the ends of the lattice vectors the HIGHEST numbered atoms.");

          if( cell->dim > 1){
            skipcomments(infile,instring,FATAL);
            sscanf(instring,"%d %d",
                   &(cell->tvects[1].begin),&(cell->tvects[1].end));
            if(cell->tvects[1].end <= cell->num_atoms - cell->dim )
              fatal("Please make the ends of the lattice vectors the HIGHEST numbered atoms.");
          }
          if(cell->dim > 2){
            skipcomments(infile,instring,FATAL);
            sscanf(instring,"%d %d",
                   &(cell->tvects[2].begin),&(cell->tvects[2].end));

            if(cell->tvects[2].end <= cell->num_atoms - cell->dim )
              fatal("Please make the ends of the lattice vectors the HIGHEST numbered atoms.");
          }

          /*********
            now check to see if the order of any of the end points
            needs to be changed
            *********/
          for(j=0;j<3 && j<cell->dim;j++){
            /* first decrement the tabs since we index arrays from 0 */
            cell->tvects[j].end--;
            cell->tvects[j].begin--;

            if(cell->tvects[j].end < cell->tvects[j].begin){
              temp = cell->tvects[j].end;
              cell->tvects[j].end = cell->tvects[j].begin;
              cell->tvects[j].begin = temp;
            }
            else if(cell->tvects[j].end == cell->tvects[j].begin &&
                    cell->tvects[j].end ){
              /* both ends of the vector are the same */
              fatal("Begin and end of a translation vector are the same.");
            }
          }
          cell->num_atoms -= cell->dim;
          cell->num_raw_atoms = cell->num_atoms;

        }
      } /* end of keyword LATTICE */

      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"K POINTS AUTO") ||
              strstr(instring,"K-POINTS AUTO") ||
              strstr(instring,"KPOINTS AUTO")){
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%d %d %d",&(details->points_per_axis[0]),
               &(details->points_per_axis[1]),&(details->points_per_axis[2]));
        details->use_automatic_kpoints = 1;
                if(details->points_per_axis[0] == 0 && details->points_per_axis[1] == 0
                         && details->points_per_axis[2] == 0 ){
                         fatal("You forgot to specify the mesh for the automatic k points.");
                         }

        fprintf(status_file,"An automatic k-point mesh (%dx%dx%d) will be used.\n",
                &(details->points_per_axis[0]),&(details->points_per_axis[1]),
                &(details->points_per_axis[2]));
      }
      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"HIGH SYMM")){
        details->use_high_symm_p=1;
      }
      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"K POINTS") ||
              strstr(instring,"K-POINTS") ||
              strstr(instring,"KPOINTS") ){
        /* a little bit of error checking */
        if( details->Execution_Mode == MOLECULAR ){
          fprintf(status_file,"Warning: a K point set was specified for a molecular \
calculation.\n");
        }

        /* read out the number of k points. */
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%d",&(details->num_KPOINTS));

        /* get space for the k points */
        points = (k_point_type *)calloc(details->num_KPOINTS,sizeof(k_point_type));
        if(!points)fatal("Can't allocate memory for k point set.");

        fprintf(status_file,"Reading %d K points.\n",details->num_KPOINTS);

        /* now read the points */
        for(i=0;i<details->num_KPOINTS;i++){
          skipcomments(infile,instring,FATAL);
          sscanf(instring,"%lf %lf %lf %lf",
                 &(points[i].loc.x),&(points[i].loc.y),
                 &(points[i].loc.z),
                 &(points[i].weight));
        }
        /* set the pointer in the details structure */
        details->K_POINTS = points;
      } /* end of keyword K POINTS */

      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"K OFFSET") ){
        if(sscanf(instring,"%s %s %lf",string1,string2,&(details->k_offset)) != 3){
          skipcomments(infile,instring,FATAL);
          sscanf(instring,"%lf",&(details->k_offset));
        }
      }
      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"ZERO OVER") ){
        /* read out the number of overlaps to zero */
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%d",&(details->num_overlaps_off));

        /* get space to store the zeroed overlaps */
        details->overlaps_off = (overlap_cancel_type *)
          calloc(details->num_overlaps_off,sizeof(overlap_cancel_type));
        if(!details->overlaps_off)
          fatal("Can't get space for zeroed overlaps");

        /* read out the overlaps to zero */
        for(i=0;i<details->num_overlaps_off;i++){
          skipcomments(infile,instring,FATAL);
          upcase(instring);
          sscanf(instring,"%s %d %d %s",numstring,
                 &(details->overlaps_off[i].which1),
                 &(details->overlaps_off[i].which2),
                 err_string);
          /* set the type */
          if( strstr(numstring,"ATOM") ){
            details->overlaps_off[i].type = P_DOS_ATOM;
          } else if( strstr(numstring,"ORB") ){
            details->overlaps_off[i].type = P_DOS_ORB;
          } else error("Invalid type of overlap to zero.");

          /* is this an intercell overlap being turned off? */
          if(strstr(err_string,"INTER")){
            details->overlaps_off[i].inter_cell = 1;
          } else details->overlaps_off[i].inter_cell = 0;

          /* now decrement the contributions */
          details->overlaps_off[i].which1--;
          details->overlaps_off[i].which2--;
        }
      } /* end of keyword "Zero Overlap" */
      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"SYM")){
        fprintf(status_file,"Symmetry will be used in the calculations.\n");
        details->use_symmetry = 1;
      }
      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"NO TOTAL DOS") ){
        details->no_total_DOS_PRT = 1;
      }
      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"PROJECT") ){
        /* read out the number of projections */
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%d",&(details->num_proj_DOS));

        /* get space for the projections */
        details->proj_DOS = (p_DOS_type *)calloc(details->num_proj_DOS,
                                                 sizeof(p_DOS_type));
        if( !details->proj_DOS )
          fatal("Can't get memory to hold projected DOS info.");

        /* now read the individual projections */
        for(i=0;i<details->num_proj_DOS;i++){
          skipcomments(infile,instring,FATAL);

          p_DOS = &(details->proj_DOS[i]);

          /* first get memory for the contributions to this projection */
          max_p_DOS = 5;
          p_DOS->weights = (real *)calloc(max_p_DOS,sizeof(real));
          if( !p_DOS->weights )fatal("Can't get space to store p_DOS->weights");
          p_DOS->contributions = (int *)calloc(max_p_DOS,sizeof(int));
          if( !p_DOS->contributions )
            fatal("Can't allocate memory to hold projected DOS contributions.");

          /******
            use strtok to read out comma separated projections
          *******/
          safe_strcpy(tempstring,(char *)strtok(instring,",\n"));
          sscanf(instring,"%s %d %lf",numstring,&p_DOS->contributions[0],
                 &p_DOS->weights[0]);

          /*******
            subtract one from the specified number, since arrays are indexed from
            zero in C
          *******/
          p_DOS->contributions[0]--;
          p_DOS->weight_sum = p_DOS->weights[0];

          p_DOS->num_contributions = 1;

          /* check to see the type of contribution */
          upcase(numstring);
          if( strstr(numstring,"ORB") ){
            p_DOS->type = P_DOS_ORB;
          }
          else if(strstr(numstring,"ATOM")){
            p_DOS->type = P_DOS_ATOM;
          }
          else if(strstr(numstring,"FMO")){
            p_DOS->type = P_DOS_FMO;
          }
          else fatal("Invalid projected DOS type in input file.");

          /***
            loop until all of the contributions have been read out
          ****/
          done = 0;
          while(!done){
            /* get the next contribution */
            if( tempstring[strlen(tempstring)-1] == '\\' ){
              skipcomments(infile,instring,FATAL);
              safe_strcpy(tempstring,(char *)strtok(instring,",\n"));
            } else{
              safe_strcpy(tempstring,(char *)strtok(0,",\n"));
            }
            if( !tempstring || tempstring[0] == 0){
              done = 1;
            }
            else{
              sscanf(tempstring,"%d %lf",
                     &p_DOS->contributions[p_DOS->num_contributions],
                     &p_DOS->weights[p_DOS->num_contributions]);

              p_DOS->weight_sum += p_DOS->weights[p_DOS->num_contributions];

              /*******
                subtract one from the specified number, since arrays are indexed from
                zero in C
              *******/
              p_DOS->contributions[p_DOS->num_contributions]--;
              p_DOS->num_contributions++;
              /* check to see if more memory is required */
              if( p_DOS->num_contributions == max_p_DOS ){
                max_p_DOS += 5;

                real_ptr = p_DOS->weights;
                p_DOS->weights = (real *)calloc(max_p_DOS,sizeof(real));
                if( !p_DOS->weights )fatal("Can't get space to enlarge p_DOS->weights");

                /* copy over the old data */
                bcopy((char *)real_ptr,(char *)p_DOS->weights,p_DOS->num_contributions*
                      sizeof(real));

                int_ptr = p_DOS->contributions;
                p_DOS->contributions = (int *)calloc(max_p_DOS,sizeof(int));
                if( !p_DOS->contributions )
                  fatal("Can't allocate memory to enlarge projected DOS contributions.");

                /* copy over the old data */
                bcopy((char *)int_ptr,(char *)p_DOS->contributions,
                      p_DOS->num_contributions*sizeof(int));
              }
            }
          }
        }
      } /* end of keyword PROJECT */

      else if( strstr(instring,"COOP") ){
        /* read out the total number of COOPs in the file */
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%d",&num_COOPS);

        /* loop until we've read in all of the COOP's in the file */
        for(i=0;i<num_COOPS;i++){
          /* get space for the new COOP */
          new_COOP = (COOP_type *)calloc(1,sizeof(COOP_type));
          if( !(new_COOP) ) fatal("Can't get space for a COOP.");

          /* read in the information */
          skipcomments(infile,instring,FATAL);
          sscanf(instring,"%s %d %d %d %lf %lf %lf",
                 numstring,&(new_COOP->which),&(new_COOP->contrib1),
                 &(new_COOP->contrib2),&(new_COOP->cell.x),
                 &(new_COOP->cell.y),&(new_COOP->cell.z));

          /* deal with the array indexing thing */
          new_COOP->contrib1--;
          new_COOP->contrib2--;

          /* figure out what type of COOP it is */
          upcase(numstring);
          if(strstr(numstring,"H-ATOM")){
            new_COOP->type = P_DOS_ATOM;new_COOP->energy_weight = TRUE;
          }
          else if(strstr(numstring,"ATOM")){
            new_COOP->type = P_DOS_ATOM;new_COOP->energy_weight = FALSE;
          }
          else if(strstr(numstring,"H-ORB")){
            new_COOP->type = P_DOS_ORB;new_COOP->energy_weight = TRUE;
          }
          else if(strstr(numstring,"ORB")){
            new_COOP->type = P_DOS_ORB;new_COOP->energy_weight= FALSE;
          }
          else if(strstr(numstring,"H-FMO")){
            new_COOP->type = P_DOS_FMO;new_COOP->energy_weight= TRUE;
          }
          else if(strstr(numstring,"FMO")){
            new_COOP->type = P_DOS_FMO;new_COOP->energy_weight= FALSE;
          }
          else fatal("Invalid COOP type in input file.");

          /*******
            now put this COOP into the proper place in the linked list.
            Start out by finding the right type column, then put this element
            at the head of the column.
          ********/
          COOP_ptr = details->the_COOPS;
          placed = 0;
          while(COOP_ptr && !placed){
            if( COOP_ptr->which == new_COOP->which ){
              new_COOP->next_to_avg = COOP_ptr->next_to_avg;
              COOP_ptr->next_to_avg = new_COOP;
              placed = 1;
            }
            COOP_ptr = COOP_ptr->next_type;
          }
          if( !placed ){
            /* put the COOP at the end of the list */
            COOP_ptr = details->the_COOPS;

            if( !COOP_ptr ){
              /* this is the first element */
              details->the_COOPS = new_COOP;
            }
            else{

              /* find the end of the list */
              while(COOP_ptr->next_type) COOP_ptr = COOP_ptr->next_type;

              COOP_ptr->next_type = new_COOP;
            }
          }
        }
      } /* end of keyword COOP */

      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"MO PRINT") ){
        /* deal with printing MO's */

        /* first read out the number of MOs to print */
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%d",&(details->num_MOs_to_print));

        /* get memory */
        details->MOs_to_print = (int *)calloc(details->num_MOs_to_print,
                                              sizeof(int));
        if(!details->MOs_to_print)fatal("Can't get memory for MOs_to_print.");

        /* now read out the MOs that will be printed */
        for(i=0;i<details->num_MOs_to_print;i++){
          skipcomments(infile,instring,FATAL);
          sscanf(instring,"%d",&(details->MOs_to_print[i]));
          details->MOs_to_print[i]--;
        }
      } /* end of keyword MO PRINT */
      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"PRINT") ){
        /*****
          dealing with parsing the printing options is hairy enough for them
          to get their own function.
          *****/
        parse_printing_options(infile,details,cell);

      } /* end of keyword PRINT */

      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"BAND") ){
        /* first get space for the band info storage */
        details->band_info =
          (band_info_type *)calloc(1,sizeof(band_info_type));
              if( !details->band_info )fatal("Can't get space for the band info.");
        band_info = details->band_info;

        /* read in the number of k points per symmetry line*/
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%d",&(band_info->points_per_line));

        /* now get the number of special points */
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%d",&(band_info->num_special_points));

        /* get space for the special points */
        band_info->special_points = (special_point_type *)
          calloc(band_info->num_special_points,sizeof(special_point_type));
        if(!band_info->special_points)
          fatal("Can't get space for special point storage.");

        /* get space for the K points */
        band_info->lines = (k_point_type *)
          calloc((band_info->num_special_points-1)*
                  band_info->points_per_line+1,
                 sizeof(special_point_type));
        if(!band_info->lines)
          fatal("Can't get space for storage of symmetry lines.");

        /* read in the special points */
        for(i=0;i<band_info->num_special_points;i++){
          skipcomments(infile,instring,FATAL);
          sscanf(instring,"%s %lf %lf %lf",&band_info->special_points[i].name,
                 &band_info->special_points[i].loc.x,
                 &band_info->special_points[i].loc.y,
                 &band_info->special_points[i].loc.z);
        }

        /* okay, got'em all, now generate the symmetry lines */
        gen_symm_lines(band_info);

        /* that's all we need to do. */

      } /* end of keyword BAND */
      else if( strstr(instring,"DIAGWO") ){
        details->diag_wo_overlap=1;
      }
      else if( strstr(instring,"NUM MOMENTS") ){
        if( sscanf(instring,"%s %s %d",string1,string2,&(details->num_moments)) != 3 ){
          skipcomments(infile,instring,FATAL);
          sscanf(instring,"%d",&(details->num_moments));
        }
      }
      else if( strstr(instring,"MOMENTS") ){
        details->do_moments=1;
      }
      else if( strstr(instring,"DUMP HAM") ){
        details->dump_hamil=1;
      }
      else if( strstr(instring,"DUMP OV") ){
        details->dump_overlap=1;
      }
      else if( strstr(instring,"DUMP SPARSE") ){
        details->dump_sparse_mats=1;
      }
      else if( strstr(instring,"DUMP DIST") ){
        details->dump_dist_mat=1;
      }
      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"FMO") ){
        if( details->num_FCO_frags ){
          error("Both FMO and FCO analysis specified.  Ignoring the FCO stuff.");
          if( details->FMO_frags ) free(details->FMO_frags);
          if( details->FMO_props ) free(details->FMO_props);
          details->num_FCO_frags = 0;
        }

        /* read out the number of fragments that they want to do */
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%d",&(details->num_FMO_frags));

        /* allocate memory to store the fragments */
        details->FMO_frags = (FMO_frag_type *)calloc(details->num_FMO_frags,
                                                     sizeof(FMO_frag_type));
        if( !details->FMO_frags ) fatal("Can't get memory for details->FMO_frags.");

        /* read out the number of electrons per fragment using strtok */
        skipcomments(infile,instring,FATAL);
        safe_strcpy(tempstring,(char *)strtok(instring,",\n"));
        sscanf(tempstring,"%lf",&(details->FMO_frags[0].num_electrons));
        for( i=1;i<details->num_FMO_frags;i++){
          safe_strcpy(tempstring,(char *)strtok(0,",\n"));

          /* error checking */
          if( !tempstring || tempstring[0] == 0 ){
            fatal("Number of electrons unspecified for some fragments.");
          }

          sscanf(tempstring,"%lf",&(details->FMO_frags[i].num_electrons));
        }

        /* now read out the fragments */
        tot_num_so_far = 0;
        for(i=0;i<details->num_FMO_frags;i++){
          FMO_frag = &(details->FMO_frags[i]);
          skipcomments(infile,instring,FATAL);

          /* use the special integer parsing routine */
          parse_integer_string(instring,&(FMO_frag->atoms_in_frag),
                               &(FMO_frag->num_atoms));

          for( j=0; j<FMO_frag->num_atoms; j++){
            FMO_frag->atoms_in_frag[j]--;

            /* some error checking */
            if( FMO_frag->atoms_in_frag[j] >= cell->num_atoms ){
              fatal("Atom number larger than num_atoms in FMO specification.");
            }
          }
          tot_num_so_far += FMO_frag->num_atoms;
        }
        if( tot_num_so_far != cell->num_atoms ){
          fprintf(stderr,"Number of atoms in FMO specification: %d not equal to num_atoms: %d\n",
                  tot_num_so_far,cell->num_atoms);
        }
      } /*** end of keyword FMO ***/

      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"FCO") ){
        /* a bit of error checking */
        if( details->num_FMO_frags ){
          error("Both FMO and FCO analysis specified.  Ignoring the FMO stuff.");
          if( details->FMO_frags ) free(details->FMO_frags);
          if( details->FMO_props ) free(details->FMO_props);
          details->num_FMO_frags = 0;
        }

        /* read out the number of fragments that they want to do */
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%d",&(details->num_FCO_frags));

        /* allocate memory to store the fragments */
        details->FMO_frags = (FMO_frag_type *)calloc(details->num_FCO_frags,
                                                     sizeof(FMO_frag_type));
        if( !details->FMO_frags ) fatal("Can't get memory for details->FMO_frags.");

        /* read out the number of electrons per fragment using strtok */
        skipcomments(infile,instring,FATAL);
        safe_strcpy(tempstring,(char *)strtok(instring,",\n"));
        sscanf(tempstring,"%lf",&(details->FMO_frags[0].num_electrons));
        for( i=1;i<details->num_FCO_frags;i++){
          safe_strcpy(tempstring,(char *)strtok(0,",\n"));

          /* error checking */
          if( !tempstring || tempstring[0] == 0 ){
            fatal("Number of electrons unspecified for some fragments.");
          }

          sscanf(tempstring,"%lf",&(details->FMO_frags[i].num_electrons));
        }

        /* now read out the fragments */
        tot_num_so_far = 0;
        for(i=0;i<details->num_FCO_frags;i++){
          FMO_frag = &(details->FMO_frags[i]);
          skipcomments(infile,instring,FATAL);

          /* use the special integer parsing routine */
          parse_integer_string(instring,&(FMO_frag->atoms_in_frag),
                               &(FMO_frag->num_atoms));

          for( j=0; j<FMO_frag->num_atoms; j++){
            FMO_frag->atoms_in_frag[j]--;

            /* some error checking */
            if( FMO_frag->atoms_in_frag[j] >= cell->num_atoms ){
              fatal("Atom number larger than num_atoms in FCO specification.");
            }
          }
          tot_num_so_far += FMO_frag->num_atoms;
        }
        if( tot_num_so_far != cell->num_atoms ){
          fprintf(stderr,"Number of atoms in FCO specification: %d not equal to num_atoms: %d\n",
                  tot_num_so_far,cell->num_atoms);
        }
      } /*** end of keyword FMO ***/


      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"CRYSTAL SPEC") ){
        fprintf(status_file,"Processing crystal specifications.\n");
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%lf %lf %lf",&(cell->xtal_defn.axis_lengths[0]),
               &(cell->xtal_defn.axis_lengths[1]),
               &(cell->xtal_defn.axis_lengths[2]));
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%lf %lf %lf",&(cell->xtal_defn.angles[0]),
               &(cell->xtal_defn.angles[1]),
               &(cell->xtal_defn.angles[2]));

      } /*** end of keyword CRYSTAL SPEC ***/


      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"ORBITAL OCCUP") ){
        if( details->Execution_Mode != MOLECULAR ){
          fatal("Orbital Occupations cannot be specified for extended systems");
        }
        else{
          skipcomments(infile,instring,FATAL);
          sscanf(instring,"%d",&details->num_orbital_occups);
          /* get the memory */
          details->orbital_occups =
            (orbital_occup_type *)calloc(details->num_orbital_occups,sizeof(orbital_occup_type));
          if( !details->orbital_occups ) fatal("Can't get space for orbital occupations");

          for(i=0;i<details->num_orbital_occups;i++){
            skipcomments(infile,instring,FATAL);
            sscanf(instring,"%d %lf",&(details->orbital_occups[i].orb),
                   &(details->orbital_occups[i].occup));
            details->orbital_occups[i].orb--;
          }
        }
      } /*** end of keyword ORBITAL OCCUP ***/


      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"SPARSIFY") ){
        if( sscanf(instring,"%s %lf",string1,&(details->sparsify_value)) != 2 ){
          skipcomments(infile,instring,FATAL);
          sscanf(instring,"%lf",&details->sparsify_value);
        }
      }

      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"RHO") ){
        if( sscanf(instring,"%s %lf",string1,&(details->rho)) != 2 ){
          skipcomments(infile,instring,FATAL);
          sscanf(instring,"%lf",&details->rho);
        }
      }

      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"JUST AVERA") ){
        details->just_avgE = 1;
      }

      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"JUST MATR") ){
        details->just_matrices = 1;
      }

      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"ALTERNATE O") ){
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%d",&details->num_occup_AVG);
        skipcomments(infile,instring,FATAL);
        sscanf(instring,"%lf",&details->occup_AVG_step);
      }

      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"MULLER IT") ){
        if( cell->geom_frags )
          fatal("Muller iteration and Geom Frags are incompatible (for now).");

        details->do_muller_it = 1;
        details->avg_props = 1;
      }
      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"MULLER E TOL") ){
        if( sscanf(instring,"%s %s %s %lf",string1,string2,string3,
                   &(details->muller_E_tol)) != 4 ){
          skipcomments(infile,instring,FATAL);
          sscanf(instring,"%lf",&details->muller_E_tol);
        }
      }
      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"MULLER Z TOL") ){
        if( sscanf(instring,"%s %s %s %lf",string1,string2,string3,
                   &(details->muller_Z_tol)) != 4 ){
          skipcomments(infile,instring,FATAL);
          sscanf(instring,"%lf",&details->muller_Z_tol);
        }
      }
      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"MULLER MIX") ){
        if( sscanf(instring,"%s %s %lf",string1,string2,
                   &(details->muller_mix)) != 3 ){
          skipcomments(infile,instring,FATAL);
          sscanf(instring,"%lf",&details->muller_mix);
        }
      }
      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"MULLER PARMS") ){
        /********
          we need the atomic parameters before we can deal with the
          muller iteration parameters, get them now if we need to.
        ********/
        if( !got_params ){
          fill_atomic_parms(cell->atoms,cell->num_atoms,infile);
          got_params = 1;
        }
        parse_muller_parms(infile,details,cell);
        got_muller_params = 1;
      }
      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"EQUIV ATOMS") ){
        parse_equiv_atoms(infile,details,cell);
      }
      /*----------------------------------------------------------------------*/
      else if( strstr(instring,"LINE WIDTH") ){
        if( sscanf(instring,"%s %s %d",string1,string2,
                   &(details->line_width)) != 3 ){
          skipcomments(infile,instring,FATAL);
          sscanf(instring,"%d",&details->line_width);
        }
      }

      /*----------------------------------------------------------------------*/
      /* hmmm, we shouldn't have gotten here. spew some error messages */
      else{
        safe_strcpy(tempstring,"Bad input line: ");
        strcat(tempstring,instring);
        error(tempstring);
      }
    }  /* end of keyword loop */
    /********
      make sure that we have parameters for the atoms by this point.
      If there was no PARAM keyword in the file, then we won't have,
      so read them in now
    ********/
    if(!got_params){
      fill_atomic_parms(cell->atoms,cell->num_atoms,infile);

      /* do the geom_frags if we need to */
      geom_frag = cell->geom_frags;
      while(geom_frag){
        fill_atomic_parms(geom_frag->atoms,geom_frag->num_atoms,infile);
        geom_frag = geom_frag->next;
      }

      /******

        if just the charge was specified, then convert that
        to a number of electrons

      ******/
      if( cell->charge != -1000.0 ){
        charge_to_num_electrons(cell);
      }

      fprintf(status_file,"Completed parameter acquisition.\n");
    }

    /*******

      make sure that we have the muller parameters if we're doing
      muller iteration....

    *******/
    if( details->do_muller_it && !got_muller_params ){
      fatal("You forgot to specify the Muller parms.");
    }

    /* now write out the atomic parameters and positions */
    write_atom_coords(cell->atoms,cell->num_raw_atoms,cell->using_Zmat,
                     cell->using_xtal_coords);

    geom_frag = cell->geom_frags;
    while(geom_frag){
      fprintf(output_file,"\n\n;------ Geom Frag Data\n");
      fprintf(output_file,"#Fragment Centered at atom %d\n",geom_frag->which);
      write_atom_coords(geom_frag->atoms,geom_frag->num_atoms,geom_frag->using_Z_mat,
                        0);
      write_atom_parms(details,geom_frag->atoms,geom_frag->num_atoms,0);
      geom_frag = geom_frag->next;
    }

    write_atom_parms(details,cell->atoms,cell->num_raw_atoms,1);

    /* build the orbital lookup table */
    build_orbital_lookup_table(cell,num_orbs,orbital_lookup_table);

#ifdef SUPPORT_NEW3_FILEIO
  }
#endif
}

/**************************************************************************************
   NOTE:  AUTO-WALSH is a registered trademark of Sneaky Weasel Software and should
    not be used without permission or even mentioned without the appropriate reverence
    and respect.

   The preceding message was a big joke.
***************************************************************************************/


