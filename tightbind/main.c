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
*     this is the main procedure for the program bind
*
*  read the file README.bind to find out more
*
*  created:  greg landrum  August 1993
*
*****************************************************************************/
/***
  Edit History:

  March '98: WG
    - print FMO-FMO COOP's to .out

  07.05.98 gL:
    added netCDF support.
   14.02.2004 gL:
     added a "10th anniversary" message :-)

***/

#include "bind.h"
FILE *COHP_file;
const char greetings[]="Welcome to the 10th Anniversary edition of YAeHMOP!\n";

#ifdef USING_THE_MAC
extern FILE *choose_mac_file(char *,char);
#include "Mac_Fopen.h"
#endif

#ifndef USING_THE_MAC
void main(argc, argv)
  int argc;
  char **argv;
#else
void main()
#endif
{
  FILE *temp_file;
  char file_name[80],err_string[240];
char test_string[80];
  int zeta_converged,Hii_converged;
  real new_num_electrons;
  COOP_type *COOP_ptr;
  int i;
  FILE *the_file=0;
  int walsh_step;
#ifdef USING_THE_MAC
  int argc;
  char argv[4][80];

        /* set up some stuff for Sioux */
        //SIOUXSettings.standalone = FALSE;
        SIOUXSettings.asktosaveonclose = FALSE;
        SIOUXSettings.autocloseonquit = FALSE;
        printf("Starting bind.\n");

  the_file = choose_mac_file(argv[1],MAC_FOPEN_OPEN_CD);
  if( !the_file ) {
          fatal("User cancelled intial file open");
  } else{
          argc = 2;
          strcpy(file_name,argv[1]);
  }

  /* get the command line arguments */
//  argc = ccommand(&argv);

#endif
  /* make sure the program was called with the right arguments */
  if( argc < 2){
    fprintf(stderr,"Usage: bind <inputfile>\n");
    exit(666);
  }

  /* install the sig_int handler */
  signal(SIGINT,handle_sigint);

  fprintf(stderr,greetings);


  /************

    get space for the atomic data and experimental (numerical experiments!! YES!)
    parameters...   (and make sure we got it)

  **************/
  unit_cell = (cell_type *)calloc(1,sizeof(cell_type));
  details = (detail_type *)calloc(1,sizeof(detail_type));
  if(!unit_cell || !details) fatal("Can't allocate initial memory.");

  /* check to see if we can open the input file */
  temp_file = fopen(argv[1],"r");
  if( !temp_file ){
    sprintf(err_string,"Can't open input file: %s",argv[1]);
    fatal(err_string);
  }
  fclose(temp_file);

  /* open the file that will be used to dump progress reports */
  strcpy(file_name,argv[1]);
  strcat(file_name,".status");
  status_file = fopen(file_name,"w+");
  if(!status_file)fatal("Can't open status file!");

  /* open the file that will be used for results */
  strcpy(file_name,argv[1]);
  strcat(file_name,".out");
  output_file = fopen(file_name,"w+");
  if(!output_file)fatal("Can't open results file!");
  fprintf(output_file,"#BIND_OUTPUT version: %s\n\n",VERSION_STRING);



  /********

    read in the data

  *********/
  read_inputfile(unit_cell,details,argv[1],&num_orbs,&orbital_lookup_table,the_file);

  /* copy the file name into the details structure */
  strcpy(details->filename,argv[1]);

  /* open the file that will be used for walsh output (if we need one) */
  if(details->walsh_details.num_vars != 0){
    strcpy(file_name,argv[1]);
    strcat(file_name,".walsh");
    walsh_file = fopen(file_name,"w+");
    if(!walsh_file)fatal("Can't open Walsh results file!");
  }

  /**********

    In order to be able to have the number of symmetry elements change
      with the distortions in a Walsh diagram, determine what elements
      are present throughout the distortion.

  **********/
  if(details->walsh_details.num_vars != 0){
    walsh_update(unit_cell,details,0,0);
    find_walsh_sym_ops(unit_cell,details);
  }else{
    if( unit_cell->using_Zmat )
      eval_Zmat_locs(unit_cell->atoms,unit_cell->num_atoms,unit_cell->dim,1);
    else if( unit_cell->using_xtal_coords ) eval_xtal_coord_locs(unit_cell,1);
    if( unit_cell->geom_frags && !details->walsh_details.num_vars )
      process_geom_frags(unit_cell);

    if(details->use_symmetry) find_sym_ops(details,unit_cell);
  }


  /*****
    if we are using automagic k-points, at this point
     we may need to do a walsh_update to get a set of atomic
     positions, then we can calculate that atomic positions
     to determine the lattice constants so that a valid
     reciprocal lattice can be built.
  ******/
  if(details->use_automatic_kpoints){
    if(details->walsh_details.num_vars != 0 ) walsh_update(unit_cell,details,0,0);

    automagic_k_points(details,unit_cell);
  }

  /********

    allocate space for the various arrays that are going to be needed

  *********/
  allocate_matrices(unit_cell,details,&Hamil_R,&Overlap_R,
                    &Hamil_K,&Overlap_K,&cmplx_hamil,&cmplx_overlap,
                    &eigenset,&work1,&work2,
                    &work3,&cmplx_work,
                    &properties,&avg_prop_info,num_orbs,
                    &tot_overlaps,orbital_lookup_table,&orbital_ordering);

  /*********

    Put any error checking that needs to be done into this
     function

  *********/
  check_for_errors(unit_cell,details,num_orbs);


  /******

    dump some useful information into the output file.

  ******/
  fprintf(output_file,"\n; Number of orbitals\n");
  fprintf(output_file,"#Num_Orbitals: %d\n",num_orbs);

  if( details->orbital_mapping_PRT ){
    fprintf(output_file,"\n; Orbital Mapping\n");
    for(i=0;i<num_orbs;i++){

      map_orb_num_to_name(test_string,i,orbital_lookup_table,num_orbs,
                          unit_cell->atoms,unit_cell->num_atoms);
      fprintf(output_file,"%d \t %s\n",i+1,test_string);
    }
  }


  /**********

    now do the calculation (loop over walsh diagram points)

    NOTE: this loops always gets executed at least once, since we
      assume that walsh_details.num_steps has been set to 1.

  ***********/
  for( walsh_step=0; walsh_step<details->walsh_details.num_steps; walsh_step++){

    /* open the file that will be used for band output (if we need one) */
    if(details->band_info){
      if( details->walsh_details.num_steps > 1 )
        sprintf(file_name,"%s.step%d.band",argv[1],walsh_step+1);
      else
        sprintf(file_name,"%s.band",argv[1]);
      if( band_file ) fclose(band_file);
      band_file = fopen(file_name,"w+");
      if(!band_file)fatal("Can't open band results file!");
    }

    /* open the file that will be used for FMO output (if we need one) */
    if(details->num_FMO_frags){
      if( details->walsh_details.num_steps > 1 )
        sprintf(file_name,"%s.step%d.FMO",argv[1],walsh_step+1);
      else
        sprintf(file_name,"%s.FMO",argv[1]);
      if(FMO_file) fclose(FMO_file);
      FMO_file = fopen(file_name,"w+");
      if(!FMO_file)fatal("Can't open FMO results file!");
      /******

        put the header into the FMO file

      *******/
      init_FMO_file(details,num_orbs,unit_cell->num_electrons);
    }


    /*******

      if we are doing a walsh diagram, set up the positions now

      *******/
    if( details->walsh_details.num_vars != 0 ){
      walsh_update(unit_cell,details,walsh_step,1);

      /* reset the charges in the update_zetas procedure */
      update_zetas(unit_cell,properties.net_chgs,(real)unit_cell->num_atoms*ZETA_TOL,
                   &zeta_converged,RESET);
    }

    if( details->Execution_Mode != MOLECULAR ){
      display_lattice_parms(unit_cell);
    }

    /**********

      generate the distance_matrix

      ***********/
    build_distance_matrix(unit_cell,details);


    /* check to see if any calculations are necessary */
    if(!(details->just_geom)){

      /**********

        loop until the zeta values converge (if we are doing either a self
        consistent calculation or charge iteration,  otherwise this
        just gets executed once)

        ***********/
      zeta_converged = 0;
      Hii_converged = 0;
      while( !zeta_converged || !Hii_converged ){
        /*************

          if we evaluate all of the overlaps once, then do it now...

          **************/
        if( (details->Execution_Mode == FAT && details->store_R_overlaps ) ||
           details->Execution_Mode == MOLECULAR ){
          /* build the R space overlap matrix */
          R_space_overlap_matrix(unit_cell,details,Overlap_R,num_orbs,
                                 tot_overlaps,orbital_lookup_table,0);

          /***********

            if we're doing FMO, do all the work for it now...

            for standard FMO (projecting molecular orbitals) we only
            need to do this once.

            ***********/
          if( details->num_FMO_frags || (details->num_FCO_frags &&
                                         details->Execution_Mode == MOLECULAR) ){

            /******

              build the hamiltonian... this is molecular type hamiltonian, so
              it's built differently here when we are doing an extended system.

              *******/
            full_R_space_Hamiltonian(unit_cell,details,Overlap_R,Hamil_R,num_orbs,
                                     orbital_lookup_table,1);

            /* first build the matrices */
            build_FMO_overlap(details,num_orbs,unit_cell->num_atoms,Overlap_R,
                              orbital_lookup_table);
            build_FMO_hamil(details,num_orbs,unit_cell->num_atoms,Hamil_R,
                            orbital_lookup_table);

            /* now diagonalize them */
            diagonalize_FMO(details,work1,work2,work3,cmplx_hamil,cmplx_overlap,cmplx_work);

            /* generate the transform matrices */
            gen_FMO_tform_matrices(details);
          }
          /* build the real hamiltonian */
          full_R_space_Hamiltonian(unit_cell,details,Overlap_R,Hamil_R,
                                   num_orbs,orbital_lookup_table,0);


        }
        else if( details->Execution_Mode == FAT &&
                !details->store_R_overlaps ){
          fprintf(stderr,"Storing S(k) instead of S(R)\n");
          fprintf(status_file,
                  "Storing the %d S(k)'s instead of the %d S(R)'s to save memory.\n",
                  details->num_KPOINTS,tot_overlaps);

          build_all_K_overlaps(unit_cell,details,Overlap_R,Overlap_K,
                               num_orbs,tot_overlaps,orbital_lookup_table);

          /* build the real hamiltonian */
          full_R_space_Hamiltonian(unit_cell,details,Overlap_R,Hamil_R,
                                   num_orbs,orbital_lookup_table,0);

        }
        else if( details->Execution_Mode == THIN ){
          /* just find the diagonal elements now... */
          R_space_Hamiltonian(unit_cell,details,Overlap_R,Hamil_R,num_orbs,
                              orbital_lookup_table);
        }

        /***********

          time to actually do the real work, i.e. do all the K points
          (or just diagonalize the matrices for the molecular case).


          work2 comes back holding the occupation numbers.
          work3 has the reduced overlap matrix.

          ************/
        loop_over_k_points(unit_cell,details,Overlap_R,Hamil_R,Overlap_K,
                           Hamil_K,cmplx_hamil,cmplx_overlap,
                           eigenset,work1,work2,work3,cmplx_work,
                           &properties,
                           avg_prop_info,num_orbs,orbital_lookup_table);


        if( !details->just_matrices ){
          /******

          evaluate the electrostatic term for molecular optimizations

          work2 is used to pass in the occupation numbers, and
          work3 is filled (in the
          function) with the orbital occupation numbers.
          *******/
          if( details->Execution_Mode == MOLECULAR && details->eval_electrostat ){
            eval_electrostatics(unit_cell,num_orbs,eigenset,work2,
                                properties.OP_mat,
                                orbital_lookup_table,&electrostatic_term,
                                &eHMO_term,
                                &total_energy,work3,properties.net_chgs);

            /* display the results */
            fprintf(output_file,"\n; Energy Partitioning:\n");
            fprintf(output_file,"\t        extended Hueckel Energy: %lg eV\n",
                    eHMO_term);
            fprintf(output_file,"\t Electrostatic Repulsion Energy: %lg eV\n",
                    electrostatic_term);
            fprintf(output_file,"\t                   Total Energy: %lg eV\n",
                    total_energy);

            fprintf(stderr,"%lg %lg %lg %lg\n", unit_cell->distance_mat[1],
                    eHMO_term,electrostatic_term,total_energy);
          }

          /*********
          do the average properties calculations
          *********/
          if( details->avg_props ){
            sort_avg_prop_info(details,num_orbs,avg_prop_info,orbital_ordering);

            find_crystal_occupations(details,unit_cell->num_electrons,num_orbs,
                                     orbital_ordering,&(properties.Fermi_E));

            /******
              now determine net charges, and orbital occupations

              the AO occupations come back in work2 in case anything else needs to
              be done with them.
              *******/
            calc_avg_occups(details,unit_cell,num_orbs,orbital_ordering,
                            avg_prop_info,&properties,work2);



            /******

              update the zeta values for the self-consistent procedure....

              for the moment this only works for molecular calculations.

              *******/
            if( details->Execution_Mode == MOLECULAR && details->vary_zeta ){
              update_zetas(unit_cell,properties.net_chgs,
                           (real)unit_cell->num_atoms*ZETA_TOL,&zeta_converged,NORMAL);
            }
            else{
              zeta_converged = 1;
            }



            /********

              now that we have the average occupations, we can go on and
              do charge iteration if it is required.

              *********/
            if( details->do_chg_it ){
              update_chg_it_parms(details,unit_cell,work2,&Hii_converged,num_orbs,
                                  orbital_lookup_table);

            }else if( details->do_muller_it ){

              print_avg_occups(details,unit_cell,num_orbs,orbital_ordering,
                               avg_prop_info,properties,work2);

              update_muller_it_parms(details,unit_cell,work2,
                                     &Hii_converged,num_orbs,
                                     orbital_lookup_table);

            } else{
              Hii_converged = 1;
            }


            if( Hii_converged && zeta_converged){
              /* write the parms used to obtain that data. */
              write_atom_parms(details,unit_cell->atoms,unit_cell->num_atoms,1);
              /******
                find the average overlap population and reduced OP matrices
                as well as the the average OP's that the user requested.
                *******/
              if( details->avg_OP_mat_PRT || details->avg_ROP_mat_PRT ){
                calc_avg_OP(details,unit_cell,num_orbs,orbital_ordering,
                            avg_prop_info,Overlap_R,properties);
              }

#ifdef INCLUDE_NETCDF_SUPPORT
              if( details->do_netCDF ){
                netCDF_init_file(details,unit_cell,num_orbs,
                                 walsh_step);
                netCDF_write_Es(details,num_orbs,avg_prop_info);
                netCDF_write_MOs(details,num_orbs,avg_prop_info);
              }
#endif
              /* Density of States */
              if( !details->no_total_DOS_PRT || !details->just_avgE ){
                gen_total_DOS(details,unit_cell,num_orbs,avg_prop_info,orbital_ordering);

                /* Projected Density of States */
                if( details->num_proj_DOS ){
                  gen_projected_DOS(details,unit_cell,num_orbs,avg_prop_info,
                                    orbital_ordering,orbital_lookup_table);
                }
                fprintf(output_file,"# END OF DOS\n\n");
              }
              /* check to see if we need to do a COOP */
              if( details->the_COOPS ){
                gen_COOP(details,unit_cell,num_orbs,avg_prop_info,Overlap_R,
                         orbital_ordering,orbital_lookup_table);
              }

              /*************

                print out the stuff that will appear at the bottom of the file.

                *************/

              /* Fermi level */
              fprintf(output_file,"\n;  The Fermi Level was determined for %d K points based on\n",
                      details->num_KPOINTS);
              fprintf(output_file,";     an ordering of %d crystal orbitals occupied by %lf electrons\n",
                      num_orbs*details->num_KPOINTS,unit_cell->num_electrons);
              fprintf(output_file,";      in the unit cell (%lf electrons total)\n",
                      unit_cell->num_electrons*(real)details->num_KPOINTS);
              fprintf(output_file,"#Fermi_Energy:  %lf\n",properties.Fermi_E);

              /* print the moments if we generated them */
              if( details->do_moments && details->moments ){
                fprintf(output_file,"; Moments Analysis\n");
                fprintf(output_file,";  Moments are normalized by the 0th moment and\n");
                fprintf(output_file,";   referenced to the 1st\n");
                fprintf(output_file,"Moment \t Value\n");
                for(i=0;i<=details->num_moments;i++){
                  fprintf(output_file,"%d \t %8.4lg\n",i,details->moments[i]);
                }
              }

              print_avg_occups(details,unit_cell,num_orbs,orbital_ordering,
                               avg_prop_info,properties,work2);

              /* deal with multiple occupations, if there are any */
              if( details->num_occup_AVG ){
                new_num_electrons = unit_cell->num_electrons;
                fprintf(output_file,
                        "\n;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
                fprintf(output_file,"\n# ALTERNATE OCCUPATION ANALYSIS\n");
                fprintf(output_file,"%d Alternate Occupations were done\n",
                        details->num_occup_AVG);
                for(i=0;i<=details->num_occup_AVG;i++){
                  if( i ) new_num_electrons += details->occup_AVG_step;
                  find_crystal_occupations(details,new_num_electrons,
                                           num_orbs,orbital_ordering,
                                           &(properties.Fermi_E));
                  calc_avg_occups(details,unit_cell,num_orbs,orbital_ordering,
                                  avg_prop_info,&properties,work2);

                  fprintf(output_file,"\n#NUM_ELECTRONS_PER_CELL: %lf\n",
                          new_num_electrons);
                  fprintf(output_file,"#Fermi_Energy:  %lf\n",
                          properties.Fermi_E);
                  print_avg_occups(details,unit_cell,num_orbs,orbital_ordering,
                                   avg_prop_info,properties,work2);
                  if( details->the_COOPS ){
                    gen_avg_COOPs(details,unit_cell,num_orbs,avg_prop_info,Overlap_R,
                                  orbital_ordering,orbital_lookup_table);
                  }
                  if( details->num_FMO_frags ){
                    calc_avg_FMO_occups(details,num_orbs,orbital_ordering,
                                        avg_prop_info,work2);
                  }

                  /********

                    now print out the average values of the COOP's that
                    the user asked for

                  **********/
                  COOP_ptr = details->the_COOPS;

                  if( COOP_ptr ){
                    fprintf(output_file,"\n; Average Values of COOP's\n");
                  }
                  while(COOP_ptr){
                    if( COOP_ptr->type == P_DOS_ATOM ){
                      fprintf(output_file,"(%d) Between atoms %s%d and %s%d: %lf\n",
                              COOP_ptr->which,
                              unit_cell->atoms[COOP_ptr->contrib1].symb,
                              COOP_ptr->contrib1+1,
                              unit_cell->atoms[COOP_ptr->contrib2].symb,
                              COOP_ptr->contrib2+1,
                              COOP_ptr->avg_value);
                    }
                    else if(COOP_ptr->type == P_DOS_ORB){
                      fprintf(output_file,"(%d) Between orbitals %d and %d: %lf\n",
                              COOP_ptr->which,
                              COOP_ptr->contrib1+1,
                              COOP_ptr->contrib2+1,
                              COOP_ptr->avg_value);
                    }
                    else if(COOP_ptr->type == P_DOS_FMO){
                      fprintf(output_file,"(%d) Between fmo's %d and %d: %lf\n",
                              COOP_ptr->which,
                              COOP_ptr->contrib1+1,
                              COOP_ptr->contrib2+1,
                              COOP_ptr->avg_value);
                    }
                    COOP_ptr = COOP_ptr->next_type;
                  }

                }

                /* to be safe, redo the original occupation stuff */
                find_crystal_occupations(details,unit_cell->num_electrons,
                                         num_orbs,orbital_ordering,
                                         &(properties.Fermi_E));
                calc_avg_occups(details,unit_cell,num_orbs,orbital_ordering,
                                avg_prop_info,&properties,work2);
              }


              if( !details->just_avgE && !details->num_occup_AVG ){
                if( details->num_FMO_frags ){
                  calc_avg_FMO_occups(details,num_orbs,orbital_ordering,
                                      avg_prop_info,work2);
                }

                /********

                  now print out the average values of the COOP's that
                  the user asked for

                  **********/
                COOP_ptr = details->the_COOPS;

                if( COOP_ptr ){
                  fprintf(output_file,"\n; Average Values of COOP's\n");
                }
                while(COOP_ptr){
                  if( COOP_ptr->type == P_DOS_ATOM ){
                    fprintf(output_file,"(%d) Between atoms %s%d and %s%d: %lf\n",
                            COOP_ptr->which,
                            unit_cell->atoms[COOP_ptr->contrib1].symb,
                            COOP_ptr->contrib1+1,
                            unit_cell->atoms[COOP_ptr->contrib2].symb,
                            COOP_ptr->contrib2+1,
                            COOP_ptr->avg_value);
                  }
                  else if(COOP_ptr->type == P_DOS_ORB){
                    fprintf(output_file,"(%d) Between orbitals %d and %d: %lf\n",
                            COOP_ptr->which,
                            COOP_ptr->contrib1+1,
                            COOP_ptr->contrib2+1,
                            COOP_ptr->avg_value);
                  }
                  else if(COOP_ptr->type == P_DOS_FMO){
                    fprintf(output_file,"(%d) Between fmo's %d and %d: %lf\n",
                            COOP_ptr->which,
                            COOP_ptr->contrib1+1,
                            COOP_ptr->contrib2+1,
                            COOP_ptr->avg_value);
                  }
                  COOP_ptr = COOP_ptr->next_type;
                }
              } /* end of if(!details->just_avgE) */
            } /* end of if(Hii_converged && zeta_converged) */
          } /* end of if(details->avg_props) */
          else{
            /* we need to set convergence stuff here */
            Hii_converged = 1;
            zeta_converged = 1;
          }
        } /* end of if(!details->just_matrices */
        else{
          Hii_converged = 1;
          zeta_converged = 1;
        }
      }/* end of convergence loop */

      /*********
        check to see if a band structure is being done.
        if so deal with it.
        **********/
      if( details->band_info ){
        construct_band_structure(unit_cell,details,Overlap_R,
                                 Hamil_R,Overlap_K,
                                 Hamil_K,cmplx_hamil,cmplx_overlap,
                                 eigenset,work1,work2,work3,cmplx_work,
                                 num_orbs,orbital_lookup_table);

      }
    }

    /******
      print out any important walsh results
      ******/
    if( details->walsh_details.num_vars != 0 ){
      if( details->Execution_Mode != MOLECULAR ){
        walsh_output( details,unit_cell,num_orbs,eigenset,Overlap_K,Hamil_K,
                     properties,orbital_lookup_table,walsh_step);
      }else{
        walsh_output( details,unit_cell,num_orbs,eigenset,Overlap_R,Hamil_R,
                     properties,orbital_lookup_table,walsh_step);
      }
    }
  } /* end of walsh loop */

  fprintf(status_file,"Done!\n");
  fprintf(stdout,"Done!\n");

  /* close the files and exit */
#ifdef INCLUDE_NETCDF_SUPPORT
  if( details->do_netCDF ){
    netCDF_close_file(details);
  }
#endif

  fclose(status_file);
  fclose(output_file);
  exit(0);
}
