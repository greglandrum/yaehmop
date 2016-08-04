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
*     this has the routines needed for netCDF support in bind
*
*  created:  greg landrum  May 7 1998
*
*****************************************************************************/
#include "bind.h"

#ifdef INCLUDE_NETCDF_SUPPORT

/* the max number of parameters */
#define PARAM_LEN 20

#define TRY_CDF(__fcall__)\
{int errval; errval=__fcall__;netCDF_handle_error(errval);}

void netCDF_handle_error(int status)
{
  char errstring[240];

  if (status != NC_NOERR) {
    sprintf(errstring,"netCDF error: %s",nc_strerror(status));
    fatal(errstring);
  }
}


/****************************************************************************
 *
 *                   Procedure netCDF_init_file
 *
 * Arguments: details: pointer to detail_type
 *              cell: pointer to cell type
 *          num_orbs: int
 *          walsh_step: int
 *
 * Returns: none
 *
 * Action:  Writes header information into a netCDF file.
 *
 ****************************************************************************/
void netCDF_init_file(detail_type *details,cell_type *cell,int num_orbs,
                      int walsh_step)
{
  int i,j,itab,jtab;
  char outname[240];
  int fileid;
  int dim_array[10],count_array[10];
  int symb_var_ID,pos_var_ID,param_var_ID;
  int kpoints_ID,kpoint_weights_ID,num_electrons_ID;
  int lattice_ID;
  real temp_vect[PARAM_LEN];

  if( !details->do_netCDF )
    FATAL_BUG("init_netCDF called without CDF writing on.");

  /* get space for the netCDF_info structure */
  details->netCDF_info = (netCDF_info_type *)my_calloc(1,sizeof(netCDF_info_type));

  if( details->walsh_details.num_steps > 1){
    sprintf(outname,"%s.step%d.nCDF",walsh_step,details->filename);
  }
  else{
    sprintf(outname,"%s.nCDF",details->filename);
  }

  /* create (or overwrite) the CDF file */
  TRY_CDF(nc_create(outname,NC_CLOBBER,&fileid));

  /* set the pointer in details */
  details->netCDF_info->file_ID = fileid;


  /* write the definitions */
  TRY_CDF(nc_def_dim(fileid,"num_orbs",
                     (size_t)num_orbs,
                     &details->netCDF_info->num_orbs_ID));
  TRY_CDF(nc_def_dim(fileid,"num_kpoints",
                     (size_t)details->num_KPOINTS,
                     &details->netCDF_info->num_KPOINTS_ID));
  TRY_CDF(nc_def_dim(fileid,"tot_num_orbs",
                     (size_t)details->num_KPOINTS*num_orbs,
                     &details->netCDF_info->tot_num_orbs_ID));
  TRY_CDF(nc_def_dim(fileid,"num_atoms",
                     (size_t)cell->num_atoms,
                     &details->netCDF_info->num_atoms_ID));
  TRY_CDF(nc_def_dim(fileid,"num_dim",
                     (size_t)cell->dim,
                     &details->netCDF_info->num_dim_ID));
  TRY_CDF(nc_def_dim(fileid,"name_length",
                     (size_t)ATOM_SYMB_LEN,
                     &details->netCDF_info->name_len_ID));
  TRY_CDF(nc_def_dim(fileid,"param_length",
                     (size_t)PARAM_LEN,
                     &details->netCDF_info->param_len_ID));

  /* it's kind of dumb that we have to do this */
  TRY_CDF(nc_def_dim(fileid,"space_dim",
                     (size_t)3,
                     &details->netCDF_info->space_dim_ID));


  /* start with the atoms */
  dim_array[0] = details->netCDF_info->num_atoms_ID;
  dim_array[1] = details->netCDF_info->name_len_ID;
  TRY_CDF(nc_def_var(details->netCDF_info->file_ID,
                     "Atom_Symbols",NC_CHAR,2,(int *)dim_array,
                     &symb_var_ID));
  dim_array[0] = details->netCDF_info->num_atoms_ID;
  dim_array[1] = details->netCDF_info->param_len_ID;
  TRY_CDF(nc_def_var(details->netCDF_info->file_ID,
                     "Atomic_Parameters",NC_DOUBLE,
                     2,(int *)dim_array,&param_var_ID));
  dim_array[0] = details->netCDF_info->num_atoms_ID;
  dim_array[1] = details->netCDF_info->space_dim_ID;
  TRY_CDF(nc_def_var(details->netCDF_info->file_ID,
                     "Atomic_Positions",NC_DOUBLE,
                     2,(int *)dim_array,&pos_var_ID));
  dim_array[0] = details->netCDF_info->num_dim_ID;
  dim_array[1] = details->netCDF_info->space_dim_ID;
  TRY_CDF(nc_def_var(details->netCDF_info->file_ID,
                     "Lattice_Vects",NC_DOUBLE,
                     2,(int *)dim_array,&lattice_ID));

  dim_array[0] = details->netCDF_info->num_KPOINTS_ID;
  dim_array[1] = details->netCDF_info->space_dim_ID;
  TRY_CDF(nc_def_var(details->netCDF_info->file_ID,
                     "Kpoint_locs",NC_DOUBLE,
                     2,(int *)dim_array,&kpoints_ID));
  dim_array[0] = details->netCDF_info->num_KPOINTS_ID;
  TRY_CDF(nc_def_var(details->netCDF_info->file_ID,
                     "Kpoint_weights",NC_DOUBLE,
                     1,(int *)dim_array,&kpoint_weights_ID));

  /* number of electrons */
  TRY_CDF(nc_def_var(details->netCDF_info->file_ID,
                     "Num_Electrons",NC_DOUBLE,
                     0,(int *)dim_array,&num_electrons_ID));


  /* close define mode */
  TRY_CDF(nc_enddef(fileid));

  /****
    write header data
  ****/

  /* loop over atoms */
  for(i=0;i<cell->num_atoms;i++){

    /* write the symbol */
    dim_array[0] = i;
    dim_array[1] = 0;
    count_array[0] = 1;
    count_array[1] = ATOM_SYMB_LEN;
    TRY_CDF(nc_put_vara_text(details->netCDF_info->file_ID,symb_var_ID,
                             (size_t *)dim_array,(size_t *)count_array,
                             cell->atoms[i].symb));

    /* then the position */
    temp_vect[0] = cell->atoms[i].loc.x;
    temp_vect[1] = cell->atoms[i].loc.y;
    temp_vect[2] = cell->atoms[i].loc.z;
    dim_array[0] = i;
    dim_array[1] = 0;
    count_array[0] = 1;
    count_array[1] = 3;
    TRY_CDF(nc_put_vara_double(details->netCDF_info->file_ID,pos_var_ID,
                               (size_t *)dim_array,(size_t *)count_array,
                               (double *)temp_vect));

    /* build an array of the parameters */
    temp_vect[0] = (double)cell->atoms[i].at_number;
    temp_vect[1] = (double)cell->atoms[i].num_valence;
    temp_vect[2] = (double)cell->atoms[i].ns;
    temp_vect[3] = (double)cell->atoms[i].coul_s;
    temp_vect[4] = (double)cell->atoms[i].exp_s;
    temp_vect[5] = (double)cell->atoms[i].np;
    temp_vect[6] = (double)cell->atoms[i].coul_p;
    temp_vect[7] = (double)cell->atoms[i].exp_p;
    temp_vect[8] = (double)cell->atoms[i].nd;
    temp_vect[9] = (double)cell->atoms[i].coul_d;
    temp_vect[10] = (double)cell->atoms[i].exp_d;
    temp_vect[11] = (double)cell->atoms[i].coeff_d1;
    temp_vect[10] = (double)cell->atoms[i].exp_d2;
    temp_vect[13] = (double)cell->atoms[i].coeff_d2;
    temp_vect[14] = (double)cell->atoms[i].nf;
    temp_vect[15] = (double)cell->atoms[i].coul_f;
    temp_vect[16] = (double)cell->atoms[i].exp_f;
    temp_vect[17] = (double)cell->atoms[i].coeff_f1;
    temp_vect[18] = (double)cell->atoms[i].exp_f2;
    temp_vect[19] = (double)cell->atoms[i].coeff_f2;

    /* then write them */
    dim_array[0] = i;
    dim_array[1] = 0;
    count_array[0] = 1;
    count_array[1] = PARAM_LEN;
    TRY_CDF(nc_put_vara_double(details->netCDF_info->file_ID,param_var_ID,
                               (size_t *)dim_array,(size_t *)count_array,
                               (double *)temp_vect));
  }

  /* slap in the translation vectors */
  for(i=0;i<cell->dim;i++){
    itab = cell->tvects[i].begin;
    jtab = cell->tvects[i].end;
    temp_vect[0] = cell->atoms[jtab].loc.x-cell->atoms[itab].loc.x;
    temp_vect[1] = cell->atoms[jtab].loc.y-cell->atoms[itab].loc.y;
    temp_vect[2] = cell->atoms[jtab].loc.z-cell->atoms[itab].loc.z;

    dim_array[0] = i;
    dim_array[1] = 0;
    count_array[0] = 1;
    count_array[1] = 3;
    TRY_CDF(nc_put_vara_double(details->netCDF_info->file_ID,lattice_ID,
                               (size_t *)dim_array,(size_t *)count_array,
                               (double *)temp_vect));
  }

  /* and the k points */
  for(i=0;i<details->num_KPOINTS;i++){
    temp_vect[0] = details->K_POINTS[i].loc.x;
    temp_vect[1] = details->K_POINTS[i].loc.y;
    temp_vect[2] = details->K_POINTS[i].loc.z;
    dim_array[0] = i;
    dim_array[1] = 0;
    count_array[0] = 1;
    count_array[1] = 3;
    TRY_CDF(nc_put_vara_double(details->netCDF_info->file_ID,kpoints_ID,
                               (size_t *)dim_array,(size_t *)count_array,
                               (double *)temp_vect));
    count_array[0] = i;
    if( cell->dim > 0 ) temp_vect[0] = details->K_POINTS[i].weight;
    else temp_vect[0] = 1;
    TRY_CDF(nc_put_var1_double(details->netCDF_info->file_ID,
                               kpoint_weights_ID,
                               (size_t *)count_array,(double *)temp_vect));
  }

  TRY_CDF(nc_put_var1_double(details->netCDF_info->file_ID,
                             num_electrons_ID,
                             (size_t *)0,(double *)&cell->num_electrons));

  /* it's safe to write the other data now */
}


/****************************************************************************
 *
 *                   Procedure netCDF_write_MOs
 *
 * Arguments: details: pointer to detail_type
 *           num_orbs: int
 *          prop_info: pointer to avg_prop_info_type
 *
 * Returns: none
 *
 * Action:  Writes orbital energies to netCDF file
 *
 ****************************************************************************/
void netCDF_write_MOs(detail_type *details,int num_orbs,
                      avg_prop_info_type *prop_info)
{
  int orbR_ID;
  int orbI_ID;
  int dim_array[10],count_array[10];
  size_t cell_dim;

  /* read out the dimension */
  TRY_CDF(nc_inq_dimlen(details->netCDF_info->file_ID,
                        details->netCDF_info->num_dim_ID,
                        (size_t *)&cell_dim));
  /****
    "allocate" the variables
  ****/

  /* switch to definition mode */
  TRY_CDF(nc_redef(details->netCDF_info->file_ID));

  /* coefficient arrays */
  dim_array[0] = details->netCDF_info->num_KPOINTS_ID;
  dim_array[1] = details->netCDF_info->num_orbs_ID;
  dim_array[2] = details->netCDF_info->num_orbs_ID;
  TRY_CDF(nc_def_var(details->netCDF_info->file_ID,
                     "Crystal_OrbitalsR",NC_FLOAT,
                     3,(int *)dim_array,&orbR_ID));
  /* only write imaginary coefficients if they are required */
  if( cell_dim > 0 ){
    TRY_CDF(nc_def_var(details->netCDF_info->file_ID,
                       "Crystal_OrbitalsI",NC_FLOAT,
                       3,(int *)dim_array,&orbI_ID));
  }
  /* close define mode */
  TRY_CDF(nc_enddef(details->netCDF_info->file_ID));

  /* now write the variables */
  dim_array[0] = 0;
  dim_array[1] = 0;
  dim_array[2] = 0;
  count_array[0] = details->num_KPOINTS;
  count_array[1] = num_orbs;
  count_array[2] = num_orbs;
  TRY_CDF(nc_put_vara_float(details->netCDF_info->file_ID,orbR_ID,
                             (size_t *)dim_array,(size_t *)count_array,
                             (float *)prop_info->orbs));
  if( cell_dim > 0 ){
    TRY_CDF(nc_put_vara_float(details->netCDF_info->file_ID,orbI_ID,
                              (size_t *)dim_array,(size_t *)count_array,
                              (float *)prop_info->orbsI));
  }

}

/****************************************************************************
 *
 *                   Procedure netCDF_write_Es
 *
 * Arguments: details: pointer to detail_type
 *           num_orbs: int
 *          prop_info: pointer to avg_prop_info_type
 *
 * Returns: none
 *
 * Action:  Writes header information into a netCDF file.
 *
 ****************************************************************************/
void netCDF_write_Es(detail_type *details,int num_orbs,
                     avg_prop_info_type *prop_info)
{
  int E_var_ID;
  int dim_array[10],count_array[10];

  TRY_CDF(nc_redef(details->netCDF_info->file_ID));
  dim_array[0] = details->netCDF_info->num_KPOINTS_ID;
  dim_array[1] = details->netCDF_info->num_orbs_ID;
  TRY_CDF(nc_def_var(details->netCDF_info->file_ID,
                     "Orbital_Energies",NC_FLOAT,
                     2,(int *)dim_array,&E_var_ID));
  /* close define mode */
  TRY_CDF(nc_enddef(details->netCDF_info->file_ID));

  dim_array[0] = 0;
  dim_array[1] = 0;
  count_array[0] = details->num_KPOINTS;
  count_array[1] = num_orbs;
  TRY_CDF(nc_put_vara_float(details->netCDF_info->file_ID,E_var_ID,
                             (size_t *)dim_array,(size_t *)count_array,
                             (float *)prop_info->energies));

}



/****************************************************************************
 *
 *                   Procedure netCDF_close_file
 *
 * Arguments: details: pointer to detail_type
 *
 * Returns: none
 *
 * Action:  closes the net CDF file in 'details
 *
 ****************************************************************************/
void netCDF_close_file(detail_type *details)
{
  if( !details->do_netCDF )
    FATAL_BUG("close_netCDF called without CDF writing on.");

  if( details->netCDF_info->file_ID ){
   TRY_CDF(nc_close(details->netCDF_info->file_ID));
  }
}

#endif
