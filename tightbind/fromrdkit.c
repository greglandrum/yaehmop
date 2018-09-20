/****************************************************************************
 *
 *                   Procedure fromrdkit
 *
 * Arguments:  cell: pointer to cell_type
 *          details: pointer to detail_type
 *             name: pointer to type char
 *         num_orbs: pointer to int
 *     orbital_lookup_table: pointer to pointer to int
 *
 * Returns: none
 *
 * Action: convert rdkit data into
 *   the two structures 'cell and 'details
 *
 *****************************************************************************/
void fromrdkit(cell,details,name,num_orbs,orbital_lookup_table,Mol)
  cell_type *cell;
  detail_type *details;
  char *name;
  int *num_orbs;
  int **orbital_lookup_table;
  RDKit::Mol *MOL;
{
  walsh_details_type *walsh;
  band_info_type *band_info;
  geom_frag_type *geom_frag;
  FMO_frag_type *FMO_frag;
  char err_string[MAX_STR_LEN];
  char instring[MAX_STR_LEN],tempstring[MAX_STR_LEN];
  char string1[MAX_STR_LEN],string2[MAX_STR_LEN],string3[MAX_STR_LEN];
  char numstring[MAX_STR_LEN];
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


% list of all possible parameters
  details->Execution_Mode = MOLECULAR;
  details->num_KPOINTS = 1;
  details->K_POINTS = (k_point_type *)calloc(1,sizeof(k_point_type));
  if( !details->K_POINTS ) fatal("Can't allocate the single k point.");
	details->K_POINTS[0].weight = 1.0;

    //details->Execution_Mode = THIN|FAT;
    //details->do_netCDF = 1;
    details->eval_electrostat = 1;
    //details->weighted_Hij = 0;
    //details->close_nn_contact;
    //details->the_const;
    //details->vary_zeta = 1; // ZETA
    //details->just_geom = 1;
    //details->avg_props = 1;

    

    // passing Zmat to cell
    cell->using_xtal_coords = 0;
    cell->atoms = (atom_type *)calloc(cell->num_atoms,sizeof(atom_type));
    if(!cell->atoms){
	  sprintf(err_string,"Can't allocate memory for: %d atoms.",cell->num_atoms);
	  fatal("Can't allocate memory for the atoms.");
	}
    Zmat=0;
	cell->using_Zmat = 0;
    // cell->atoms[0].symb,&cell->atoms[0].loc.x,&cell->atoms[0].loc.y,&cell->atoms[0].loc.z)
	reading_atom_nums = 1;




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
      } 
    /* end of keyword GEOMETRY */
      

    
      /*----------------------------------------------------------------------*/
      else if(strstr(instring,"ELECTR")){
	if( cell->charge != -1000.0 ){
	  fprintf(stderr,
		  "Both Electrons and Charge keywords specified.  Charge ignored.\n");
	  cell->charge = -1000.0;
	}
	
    /// must specify only one of the those two parameters
    cell->num_electrons = 19;
    cell->charge = 1;
    charge_to_num_electrons(cell);
    details->symm_tol;

    details->find_princ_axes = 1;

    // remove walsh calls

    //details->use_symmetry = 1;
    //details->use_high_symm_p=1;

    //details->num_MOs_to_print
    details->MOs_to_print = (int *)calloc(details->num_MOs_to_print,
					      sizeof(int));
	if(!details->MOs_to_print)fatal("Can't get memory for MOs_to_print.");

	/* now read out the MOs that will be printed */
	for(i=0;i<details->num_MOs_to_print;i++){
	  skipcomments(infile,instring,FATAL);
	  sscanf(instring,"%d",&(details->MOs_to_print[i]));
	  details->MOs_to_print[i]--;
	}


    // check this method to get back data to rdkit
	parse_printing_options(infile,details,cell);
	

    //details->sparsify_value;
    //details->rho;
    //details->just_avgE = 1;
    //details->just_matrices = 1;


    /********
      make sure that we have parameters for the atoms by this point.
      If there was no PARAM keyword in the file, then we won't have,
      so read them in now
    ********/
    if(!got_params){
      fill_atomic_parms(cell->atoms,cell->num_atoms,infile);

      /******

	if just the charge was specified, then convert that
	to a number of electrons

      ******/
      if( cell->charge != -1000.0 ){
	charge_to_num_electrons(cell);
      }
    }      

    /* now write out the atomic parameters and positions */
    write_atom_coords(cell->atoms,cell->num_raw_atoms,cell->using_Zmat,
		     cell->using_xtal_coords);
    
    write_atom_parms(details,cell->atoms,cell->num_raw_atoms,1);

    /* build the orbital lookup table */
    build_orbital_lookup_table(cell,num_orbs,orbital_lookup_table);

      
      
    /*======================= list of other parameters required ???? ==================
      
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

=========================*/

