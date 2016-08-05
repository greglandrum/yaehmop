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
 *   These are the things needed for generating COOP data
 *
 *  created:  greg landrum  April 1994
 *
 ******************************************************************************/

 /***
   Edit History:

   March '98: WG
     - COHP analysis added [eval_COOP()]
     - FMO COOP's added [eval_COOP()]

   January '99: WG
     - function 'intercell_COOP_check()' added
 ***/
 #include "bind.h"


 /****************************************************************************
  *
  *                   Function eval_COOP
  *
  * Arguments:   COOP: pointer to COOP_type
  *            details: pointer to detail_type
  *              cell: pointer to cell_type
  *          num_orbs: int
  *         prop_info: pointer to avg_prop_info_type
  *        R_overlaps: pointer to hermetian_matrix_type
  *  orbital_ordering: pointer to K_orb_ptr_type
  * orbital_lookup_table: pointer to int
  *
  *
  * Returns: real
  *
  * Action:   This evaluates and returns the actual COOP.
  *
  ****************************************************************************/
 real eval_COOP(COOP,details,cell,num_orbs,prop_info,R_overlaps,orbital_ordering,
                 orbital_lookup_table)
   COOP_type *COOP;
   detail_type *details;
   cell_type *cell;
   int num_orbs;
   avg_prop_info_type *prop_info;
   hermetian_matrix_type R_overlaps;
   K_orb_ptr_type *orbital_ordering;
   int *orbital_lookup_table;
 {
   static point_type *cell_dim=0;
   static real *overlap_store=0;
   static int last_overlap=-1;
   int i,ii,j,itab,jtab;
   real accum,accumI,temp,temp2;
   real answer, answer_contrib,Hii_1, Hii_2;
   int begin1,begin2,end1,end2;
   int found1,found2,offset1,offset2;
   int orbital_sum,increment,set1,set2;
   int *FMO_map1,*FMO_map2,frag_1,frag_2,fmo_index1,fmo_index2;
   int which_overlap;
   hermetian_matrix_type overlap;
   real phaseR,phaseI;
   float *MO_ptr,*MO_ptrI,*FMO_ptr,*FMO_ptrI;
   real *FMO_AOptr1,*FMO_AOptr2;
   real ao_term;
   point_type k,R;
   real kdotR;
   char doing_unit_cell;

   overlap.dim = num_orbs;

   /* check to see if we need to fill in the dimensions of the unit cell */
   if( !cell_dim ){
     cell_dim = (point_type *)calloc(3,sizeof(point_type));
     if( !cell_dim ) fatal("Can't get space for cell_dim in eval_COOP.");

     /* fill in the dimensions */
     for(i=0;i<3;i++){
       itab = cell->tvects[i].begin;
       jtab = cell->tvects[i].end;
       cell_dim[i].x = cell->atoms[jtab].loc.x-cell->atoms[itab].loc.x;
       cell_dim[i].y = cell->atoms[jtab].loc.y-cell->atoms[itab].loc.y;
       cell_dim[i].z = cell->atoms[jtab].loc.z-cell->atoms[itab].loc.z;
     }
   }

   /* check to see if we need space to store overlap matrices */
   if( !overlap_store ){
     overlap_store = (real *)calloc(num_orbs*num_orbs,sizeof(real));
     if( !overlap_store ) fatal("Can't get space for overlap_store.");
   }

   /* get pointers to the orbital information */
   MO_ptr = &(prop_info->orbs[orbital_ordering->MO*num_orbs]);
   MO_ptrI = &(prop_info->orbsI[orbital_ordering->MO*num_orbs]);

   /* figure out which overlap matrix and phase factor we should be using */
   which_overlap = overlap_tab_from_vect(&(COOP->cell),cell);
   overlap.mat = &(R_overlaps.mat[which_overlap*num_orbs*num_orbs]);

   /* check to see if we are doing a COOP w/in the unit cell */
   if( COOP->cell.x == 0.0 && COOP->cell.y == 0.0 && COOP->cell.z == 0.0 ){
     doing_unit_cell = 1;
   } else{
     doing_unit_cell = 0;
   }

   k.x = details->K_POINTS[orbital_ordering->Kpoint].loc.x;
   k.y = details->K_POINTS[orbital_ordering->Kpoint].loc.y;
   k.z = details->K_POINTS[orbital_ordering->Kpoint].loc.z;
   R.x = COOP->cell.x;
   R.y = COOP->cell.y;
   R.z = COOP->cell.z;

   /* evaluate the dot product */
   kdotR = TWOPI*(k.x*R.x + k.y*R.y + k.z*R.z);

   /* determine the phase */
   phaseR = cos(kdotR);
   phaseI = sin(kdotR);

   /*********

     In order to get the correct overlap matrix for COOP's
     outside of the unit cell, we need to undo the
     original symmetrization that was done when the overlap
     matrices were built.  This is taken care of as part of evaluating
     the COOP.  This is why the overlap terms appear in such a strange
     manner.

   *********/
   accum = 0;
   switch(COOP->type){
   case P_DOS_ORB:
     /* this is the sum of the coefficients */
     accum = ((real)MO_ptr[COOP->contrib1] * MO_ptr[COOP->contrib2] +
              (real)MO_ptrI[COOP->contrib1] * MO_ptrI[COOP->contrib2]) /
                ((real)MULTIPLIER*(real)MULTIPLIER);
     accumI = ((real)MO_ptr[COOP->contrib1] * MO_ptrI[COOP->contrib2] -
               (real)MO_ptrI[COOP->contrib1] * MO_ptr[COOP->contrib2]) /
                 ((real)MULTIPLIER*(real)MULTIPLIER);
     /*****

       now multiply on the overlap matrix and phase factor

       This is hairy because of the symmetrization thing

     ******/
     if(COOP->contrib1 != COOP->contrib2){
       if( details->Execution_Mode != MOLECULAR ){
         if( !doing_unit_cell ){
           if( COOP->contrib2 > COOP->contrib1 ){
             answer = 2.0*(accum*phaseR - accumI*phaseI)*
               (overlap.mat[COOP->contrib2*num_orbs+COOP->contrib1]+
                overlap.mat[COOP->contrib1*num_orbs+COOP->contrib2]);
           } else{
             answer = 2.0*(accum*phaseR - accumI*phaseI)*
               (overlap.mat[COOP->contrib2*num_orbs+COOP->contrib1]-
                overlap.mat[COOP->contrib1*num_orbs+COOP->contrib2]);
           }
         } else{
           if( COOP->contrib1 > COOP->contrib2 ){
             answer = 4.0*overlap.mat[COOP->contrib2*num_orbs+COOP->contrib1]*
               (accum*phaseR - accumI*phaseI);
           } else{
             answer = 4.0*overlap.mat[COOP->contrib1*num_orbs+COOP->contrib2]*
               (accum*phaseR - accumI*phaseI);
           }
         }
       } else{
         answer = 4.0*overlap.mat[COOP->contrib1*num_orbs+COOP->contrib2]*
           (accum*phaseR - accumI*phaseI);
       }
     }
     else{
       if( details->Execution_Mode != MOLECULAR ){
         if( !doing_unit_cell ){
           answer = 2.0*overlap.mat[COOP->contrib2*num_orbs+COOP->contrib1]*
             (accum*phaseR + accumI*phaseI);
         } else{
           answer = 2.0*overlap.mat[COOP->contrib2*num_orbs+COOP->contrib1]*
             (accum*phaseR + accumI*phaseI);
         }
       } else{
         answer = 2.0*overlap.mat[COOP->contrib1*num_orbs+COOP->contrib2]*
           (accum*phaseR + accumI*phaseI);
       }
     }

     /**************************************
         energy weighting stuff for COOP's
      **************************************/

     if(COOP->energy_weight){
       found1=0;found2=0;offset1=0;

       for(ii=0; ii < cell->num_atoms; ii++)
         {
           find_atoms_orbs(num_orbs,cell->num_atoms,ii,orbital_lookup_table,&begin1,&end1);

           if((COOP->contrib1>=begin1)&&(COOP->contrib1<end1)){

             found1=1;
             offset1=COOP->contrib1 - begin1;

             if(offset1<1){
               Hii_1=cell->atoms[ii].coul_s;
             }
             if(offset1>0 && offset1<4){
               Hii_1=cell->atoms[ii].coul_p;
             }
             if(offset1>3 && offset1<9){
               Hii_1=cell->atoms[ii].coul_d;
             }
             if(offset1>8){
               Hii_1=cell->atoms[ii].coul_f;
             }
           }

           if((COOP->contrib2>=begin1)&&(COOP->contrib2<end1)){

             found2=1;
             offset1=COOP->contrib2 - begin1;

             if(offset1<1){
               Hii_2=cell->atoms[ii].coul_s;
             }
             if(offset1>0 && offset1<4){
               Hii_2=cell->atoms[ii].coul_p;
             }
             if(offset1>3 && offset1<9){
               Hii_2=cell->atoms[ii].coul_d;
             }
             if(offset1>8){
               Hii_2=cell->atoms[ii].coul_f;
             }
           }
         }
       if(!(found1 && found2)){
         FATAL_BUG("Can't find Hii's");
       }

       if(doing_unit_cell && (COOP->contrib1 == COOP->contrib2)){
         answer *= Hii_1;
       }
       else{
         if(details->weighted_Hij){

           /* stuff for weighted constant */
           temp = Hii_1 + Hii_2;
           temp2 = (Hii_1 - Hii_2)/temp;
           temp2 *= temp2;
           temp2 = 0.5*(details->the_const + temp2 + temp2*temp2*(1-details->the_const));
           answer *= temp2*temp;
         }
         else{
           answer *=(0.5*details->the_const*(Hii_1+Hii_2));
         }
       }
     }
     break;

   case P_DOS_ATOM:
     answer=0;answer_contrib=0;

     /* find the orbitals */
     find_atoms_orbs(num_orbs,cell->num_atoms,COOP->contrib1,orbital_lookup_table,
                     &begin1,&end1);
     if( begin1 >= 0 ){
       find_atoms_orbs(num_orbs,cell->num_atoms,COOP->contrib2,orbital_lookup_table,
                       &begin2,&end2);

       if( begin2 >= 0 ){

         offset1=0; /* set counter for orbitals of 1st atom */
         for(i=begin1; i<end1; i++){

           if(COOP->energy_weight){
             if(offset1<1){
               Hii_1=cell->atoms[COOP->contrib1].coul_s;
             }
             if(offset1>0 && offset1<4){
               Hii_1=cell->atoms[COOP->contrib1].coul_p;
             }
             if(offset1>3 && offset1<9){
               Hii_1=cell->atoms[COOP->contrib1].coul_d;
             }
             if(offset1>8){
               Hii_1=cell->atoms[COOP->contrib1].coul_f;
             }
           }

           offset2=0; /* reset counter for orbitals of 2nd atom */
           for(j=begin2; j<end2; j++){

             if(COOP->energy_weight){
               if(offset2<1){
                 Hii_2=cell->atoms[COOP->contrib2].coul_s;
               }
               if(offset2>0 && offset2<4){
                 Hii_2=cell->atoms[COOP->contrib2].coul_p;
               }
               if(offset2>3 && offset2<9){
                 Hii_2=cell->atoms[COOP->contrib2].coul_d;
               }
               if(offset2>8){
                 Hii_2=cell->atoms[COOP->contrib2].coul_f;
               }
             }
             accum = ((real)MO_ptr[i] * MO_ptr[j] +
                      (real)MO_ptrI[i] * MO_ptrI[j]) /
                        ((real)MULTIPLIER*(real)MULTIPLIER);
             accumI = ((real)MO_ptr[i] * MO_ptrI[j] -
                       (real)MO_ptrI[i] * MO_ptr[j]) /
                         ((real)MULTIPLIER*(real)MULTIPLIER);

             if( i != j ){
               if( details->Execution_Mode != MOLECULAR ){
                 if( !doing_unit_cell ){

                   if( j > i ){
                     answer_contrib = 2.0*(overlap.mat[j*num_orbs+i] +
                                    overlap.mat[i*num_orbs+j])*
                                      (accum*phaseR -  accumI*phaseI);
                   } else{
                     answer_contrib = 2.0*(overlap.mat[j*num_orbs+i] -
                                    overlap.mat[i*num_orbs+j])*
                                      (accum*phaseR -  accumI*phaseI);
                   }

                 } else{
                   if( i > j ){
                     answer_contrib = 4.0*overlap.mat[j*num_orbs+i] *
                       (accum*phaseR -  accumI*phaseI);
                   } else{
                     answer_contrib = 4.0*overlap.mat[i*num_orbs+j] *
                       (accum*phaseR -  accumI*phaseI);
                   }
                 }
               }else{
                 if( i > j ){
                   answer_contrib = 4.0*overlap.mat[j*num_orbs+i] *
                     (accum*phaseR -  accumI*phaseI);
                 } else{
                   answer_contrib = 4.0*overlap.mat[i*num_orbs+j] *
                     (accum*phaseR -  accumI*phaseI);
                 }
               }
             }
             else{
               if( details->Execution_Mode != MOLECULAR ){
                 answer_contrib = 2.0*overlap.mat[j*num_orbs+i] *
                   (accum*phaseR -  accumI*phaseI);
               }else{
                 answer_contrib = 2.0*overlap.mat[i*num_orbs+j] *
                   (accum*phaseR -  accumI*phaseI);
               }
             }

             /* add energy weighting factor (if appropriate) */
             if(COOP->energy_weight){

               if(doing_unit_cell && (COOP->contrib1 == COOP->contrib2)){
                 answer_contrib *= Hii_1;
               }
               else{
                 if(details->weighted_Hij){

                   temp = Hii_1 + Hii_2;
                   temp2 = (Hii_1 - Hii_2)/temp;
                   temp2 *= temp2;
                   temp2 = 0.5*(details->the_const + temp2 + temp2*temp2*(1-details->the_const));
                   answer_contrib *= temp2*temp;
                 }
                 else{
                   answer_contrib *= 0.5*details->the_const*(Hii_1+Hii_2);
                 }
               }
             }
             answer += answer_contrib; answer_contrib=0;
             offset2++;
           }
           offset1++;
         }
       }
     }
     break;

   case P_DOS_FMO:
     answer=0;

     /* set pointers to FMO coefficients */
     FMO_ptr = &(prop_info->FMO_orbs[orbital_ordering->MO*num_orbs]);
     FMO_ptrI = &(prop_info->FMO_orbsI[orbital_ordering->MO*num_orbs]);

     /* eval sum of coefficients in FMO expansion */
     accum = ((real)FMO_ptr[COOP->contrib1]*FMO_ptr[COOP->contrib2]
              + (real)FMO_ptrI[COOP->contrib1]*FMO_ptrI[COOP->contrib2]) /
                ((real)MULTIPLIER*(real)MULTIPLIER);

     accumI = ((real)FMO_ptr[COOP->contrib1]*FMO_ptrI[COOP->contrib2]
               - (real)FMO_ptrI[COOP->contrib1]*FMO_ptr[COOP->contrib2]) /
                 ((real)MULTIPLIER*(real)MULTIPLIER);

     /* identify fragments for FMO to AO projection of overlap */

     orbital_sum=0; increment=0; set1=0; set2=0;

     for(i=0;i<details->num_FMO_frags;i++){
       increment = details->FMO_frags[i].num_orbs;
       if((COOP->contrib1 < (orbital_sum + increment))&&(!set1)){
         frag_1=i; set1=1;
         fmo_index1=(COOP->contrib1 - orbital_sum);
       }
       if((COOP->contrib2 < (orbital_sum + increment))&&(!set2)){
         frag_2=i; set2=1;
         fmo_index2=(COOP->contrib2 - orbital_sum);
       }
       orbital_sum += increment;
     }

     /* set pointers to FMO's in AO basis */
     FMO_AOptr1 = details->FMO_frags[frag_1].eigenset.vectR;
     FMO_AOptr2 = details->FMO_frags[frag_2].eigenset.vectR;

     /* build FMO to AO map for fragments */
     FMO_map1 = calloc(details->FMO_frags[frag_1].num_orbs,sizeof(int));

     /* loop over atoms in fragment */
     increment=0;
     for(i=0;i<details->FMO_frags[frag_1].num_atoms;i++){

       find_atoms_orbs(num_orbs,cell->num_atoms,details->FMO_frags[frag_1].atoms_in_frag[i],
                       orbital_lookup_table,&begin1,&end1);

       while(begin1<end1)
         {
           FMO_map1[increment]=begin1;
           begin1++; increment++;
         }

       if(increment > details->FMO_frags[frag_1].num_orbs){
         FATAL_BUG("FMO to AO map exceeds number of orbitals in fragment");
       }
     }

     if(frag_1 == frag_2){
       FMO_map2 = FMO_map1;
     }else{
       /* set up FMO to AO map for frag2 */
       FMO_map2 = calloc(details->FMO_frags[frag_2].num_orbs,sizeof(int));
       increment=0;

       /* loop over atoms in fragment 2 */
       for(i=0;i<details->FMO_frags[frag_2].num_atoms;i++){
         find_atoms_orbs(num_orbs,cell->num_atoms,details->FMO_frags[frag_2].atoms_in_frag[i],
                         orbital_lookup_table,&begin2,&end2);
         while(begin2<end2)
           {
             FMO_map2[increment]=begin2;
             begin2++; increment++;
           }

         if(increment > details->FMO_frags[frag_2].num_orbs){
           FATAL_BUG("FMO to AO map exceeds number of fragment orbitals");
         }
       }
     }

     /* loop over AO's in FMO's and eval COOP */
     for(i=0;i<details->FMO_frags[frag_1].num_orbs;i++){
       for(j=0;j<details->FMO_frags[frag_2].num_orbs;j++){

         /* identify AO coefficients within each FMO */
         ao_term = FMO_AOptr1[fmo_index1*details->FMO_frags[frag_1].eigenset.dim + i]*
           FMO_AOptr2[fmo_index2*details->FMO_frags[frag_2].eigenset.dim + j];

         if(ao_term){
           if(FMO_map1[i] != FMO_map2[j]){
             if(details->Execution_Mode != MOLECULAR){
               if(!doing_unit_cell){
                 if(FMO_map2[j] > FMO_map1[i]){
                   answer_contrib = 2.0*(accum*phaseR - accumI*phaseI)*ao_term*
                     (overlap.mat[FMO_map2[j]*num_orbs + FMO_map1[i]]+
                      overlap.mat[FMO_map1[i]*num_orbs + FMO_map2[j]]);
                 }
                 else{
                   answer_contrib = 2.0*(accum*phaseR - accumI*phaseI)*ao_term*
                     (overlap.mat[FMO_map2[j]*num_orbs + FMO_map1[i]]-
                      overlap.mat[FMO_map1[i]*num_orbs + FMO_map2[j]]);
                 }
               }
               else{
                 if(FMO_map1[i] > FMO_map2[j]){
                   answer_contrib = 2.0*(accum*phaseR - accumI*phaseI)*ao_term*
                     overlap.mat[FMO_map2[j]*num_orbs + FMO_map1[i]];
                 }
                 else{
                   answer_contrib = 2.0*(accum*phaseR - accumI*phaseI)*ao_term*
                     overlap.mat[FMO_map1[i]*num_orbs + FMO_map2[j]];
                 }
                 if(details->num_FMO_frags){
                   answer_contrib *= 2.0;
                 }
               }
             }
             else{
               answer_contrib = 4.0*(accum*phaseR - accumI*phaseI)*ao_term;

               if((FMO_map1[i] < FMO_map2[j])||(FMO_map1 == FMO_map2)){
                 answer_contrib *= overlap.mat[FMO_map1[i]*num_orbs + FMO_map2[j]];
               }
               else{
                 answer_contrib *= overlap.mat[FMO_map2[j]*num_orbs + FMO_map1[i]];
               }
             }
           }
           else{
             if(details->Execution_Mode != MOLECULAR){
               if(!doing_unit_cell){
                 answer_contrib = 2.0*(accum*phaseR + accumI*phaseI)*ao_term*
                   overlap.mat[FMO_map2[j]*num_orbs + FMO_map1[i]];
               }
               else{
                 answer_contrib = 2.0*(accum*phaseR + accumI*phaseI)*ao_term*
                   overlap.mat[FMO_map2[j]*num_orbs + FMO_map1[i]];
               }
             }
             else{
               answer_contrib = 2.0*(accum*phaseR + accumI*phaseI)*ao_term*
                 overlap.mat[FMO_map1[i]*num_orbs + FMO_map2[j]];
             }
           }

           if(COOP->energy_weight){
             found1=0; found2=0; ii=0; offset1=0;

             for(ii=0; ii < cell->num_atoms ; ii++)
               {
                 find_atoms_orbs(num_orbs,cell->num_atoms,ii,orbital_lookup_table,&begin1,&end1);

                 if((FMO_map1[i] >= begin1) && (FMO_map1[i] < end1)){
                   found1=1;
                   offset1=FMO_map1[i] - begin1;

                   if(offset1<1){
                     Hii_1=cell->atoms[ii].coul_s;
                   }
                   if(offset1>0 && offset1<4){
                     Hii_1=cell->atoms[ii].coul_p;
                   }
                   if(offset1>3 && offset1<9){
                     Hii_1=cell->atoms[ii].coul_d;
                   }
                   if(offset1>8){
                     Hii_1=cell->atoms[ii].coul_f;
                   }
                 }
                 if((FMO_map2[j] >= begin1) && (FMO_map2[j] < end1)){
                   found2=1;
                   offset1=FMO_map2[j] - begin1;

                   if(offset1<1){
                     Hii_2=cell->atoms[ii].coul_s;
                   }
                   if(offset1>0 && offset1<4){
                     Hii_2=cell->atoms[ii].coul_p;
                   }
                   if(offset1>3 && offset1<9){
                     Hii_2=cell->atoms[ii].coul_d;
                   }
                   if(offset1>8){
                     Hii_2=cell->atoms[ii].coul_f;
                   }
                 }
               }
             if(!(found1 && found2)){
               FATAL_BUG("Can't find Hii's in P_DOS_FMO block");
             }

             if(doing_unit_cell && (FMO_map1[i] == FMO_map2[j])){
               answer_contrib *= Hii_1;
             }
             else{

               if(details->weighted_Hij){

                 /* weighted constant */
                 temp = Hii_1 + Hii_2;
                 temp2 = (Hii_1 - Hii_2)/temp;
                 temp2 *= temp2;
                 temp2 = 0.5*(details->the_const + temp2 + temp2*temp2*(1 - details->the_const));
                 answer_contrib *= temp2*temp;
               }
               else{
                 answer_contrib *= (0.5*details->the_const*(Hii_1 + Hii_2));
               }
             }
           }
           answer += answer_contrib; /* update COOP value */
         }
       }
     }
     /* free FMO_map memory blocks */
     free(FMO_map1);
     if(FMO_map2 != FMO_map1) free(FMO_map2);
   }
   /* we're done, return the answer that we have calculated */
   return(answer);
 }


 /****************************************************************************
  *
  *                   Procedure gen_COOP
  *
  * Arguments: details: pointer to detail_type
  *              cell: pointer to cell type
  *          num_orbs: int
  *     avg_prop_info: pointer to avg_prop_info_type
  *        R_overlaps: pointer to hermetian_matrix_type
  *  orbital_ordering: pointer to K_orb_ptr_type
  *  orbital_lookup_table: pointer to int.
  *
  *
  * Returns: none
  *
  * Action:  Generates a COOP and dumps it to the output file.
  *
  *
  *  The definition used for the COOP (P) is:
  *    P_uv = 2*C_u*C_v*S_uv
  *  The terms are summed for Atom-Atom COOPs.
  *
  *  This expression is evaluated for every Crystal orbital, weighted by the
  *    weighting factor for the K point of the orbitals.
  *
  *  When evaluating COOPs between cells, there is an additional phase factor
  *    present in the expression for P_uv.
  *
  ****************************************************************************/
 void gen_COOP(details,cell,num_orbs,avg_prop_info,R_overlaps,orbital_ordering,
               orbital_lookup_table)
   detail_type *details;
   cell_type *cell;
   int num_orbs;
   avg_prop_info_type *avg_prop_info;
   hermetian_matrix_type R_overlaps;
   K_orb_ptr_type *orbital_ordering;
   int *orbital_lookup_table;
 {
   int i,j;
   int tot_num_orbs,num_COOPS;
   int num_to_avg;
   real COOP_accum,temp;
   real this_E,diff;
   real tot_num_K;
   real num_occup_bands;
   COOP_type *COOP_ptr1,*COOP_ptr2;

   tot_num_orbs = num_orbs * details->num_KPOINTS;

   fprintf(output_file,"# COOP (Crystal Orbital Overlap Population) results\n");

   /* count up the number of COOP curves */
   COOP_ptr1 = details->the_COOPS;
   num_COOPS = 0;
   while(COOP_ptr1){
     num_COOPS++;
     COOP_ptr1 = COOP_ptr1->next_type;
   }
   fprintf(output_file,"%d curves will be generated.\n",num_COOPS);

   /*******

     first determine some values that are needed for averaging purposes

   *******/

   /* the total k weighting that is used */
   tot_num_K = 0.0;
   num_occup_bands = 0.0;
   for( i=0; i<details->num_KPOINTS; i++){
     tot_num_K += details->K_POINTS[i].weight*details->K_POINTS[i].num_filled_bands;
     num_occup_bands += details->K_POINTS[i].num_filled_bands;
   }

   /******

     go through the list of COOPS and loop through the crystal orbitals for each one

   *******/
   COOP_ptr1 = details->the_COOPS;
   while(COOP_ptr1){

     /* initialize the average value */
     COOP_ptr1->avg_value = 0.0;

     /* print out the contributions to this COOP */
     fprintf(output_file,"; Contributions to this COOP are: \n");
     COOP_ptr2 = COOP_ptr1;
     num_to_avg = 0;
     while(COOP_ptr2){
       num_to_avg++;
       fprintf(output_file,"; COOP between");
       if( COOP_ptr2->type == P_DOS_ORB ){
         fprintf(output_file," orbitals ");
       }
       else if( COOP_ptr2->type == P_DOS_ATOM ){
         fprintf(output_file," atoms ");
       }
       else if( COOP_ptr2->type == P_DOS_FMO ){
         fprintf(output_file," FMO's ");
       }
       fprintf(output_file,"%d and %d in cell: ",COOP_ptr2->contrib1+1,
               COOP_ptr2->contrib2+1);
       fprintf(output_file,"( %lg %lg %lg ).\n",COOP_ptr2->cell.x,COOP_ptr2->cell.y,
               COOP_ptr2->cell.z);
       COOP_ptr2 = COOP_ptr2->next_to_avg;
     }

     fprintf(output_file,"#BEGIN CURVE\n");

     /* check for inversion of intercell vector */
     intercell_COOP_check(COOP_ptr1);

     /*******
       Add up all the contributions to the COOP that lie within DOS_DEGEN_TOL
       of each other. (This is just like how the DOS is generated).
     ********/
     i=0;
     while(i<tot_num_orbs){
       COOP_ptr2 = COOP_ptr1;
       COOP_accum = 0.0;
       while(COOP_ptr2){
         temp = details->K_POINTS[orbital_ordering[i].Kpoint].weight *
           eval_COOP(COOP_ptr2,details,cell,num_orbs,
                     &(avg_prop_info[orbital_ordering[i].Kpoint]),
                     R_overlaps,
                     &(orbital_ordering[i]),orbital_lookup_table);
         COOP_accum += temp;

         /*****
           accumulate the average value.  Since the calculation done in
           eval_COOP assumes that every orbital has 2 electrons in it, we
           need to divide temp by 2 and multiply by the number of
           electrons actually in the orbital.
         ******/
         COOP_ptr1->avg_value += orbital_ordering[i].occup * temp / 2.0;

         COOP_ptr2 = COOP_ptr2->next_to_avg;
       }
       this_E = (real)*(orbital_ordering[i].energy);
       j=i+1;

       if( j < tot_num_orbs ){
         diff = (real)*(orbital_ordering[i].energy) - (real)*(orbital_ordering[j].energy);
         while(fabs(diff) < DOS_DEGEN_TOL && j<tot_num_orbs){
           COOP_ptr2 = COOP_ptr1;
           while(COOP_ptr2){
             temp = details->K_POINTS[orbital_ordering[j].Kpoint].weight *
               eval_COOP(COOP_ptr2,details,cell,num_orbs,
                         &(avg_prop_info[orbital_ordering[j].Kpoint]),
                         R_overlaps,
                         &(orbital_ordering[j]),orbital_lookup_table);
             COOP_accum += temp;

             COOP_ptr1->avg_value += orbital_ordering[j].occup * temp / 2.0;

             COOP_ptr2 = COOP_ptr2->next_to_avg;
           }
           j++;

           if( j < tot_num_orbs ){
             diff = (real)*(orbital_ordering[i].energy) - (real)*(orbital_ordering[j].energy);
           }
         }
       }
       i = j;
       /* write out the result */
       fprintf(output_file,"%lg %lg\n",COOP_accum*num_occup_bands /
               (num_to_avg*details->num_KPOINTS*tot_num_K),
               this_E);
     }

     /* now correct the accumulated value to make it the AVERAGE value */
     COOP_ptr1->avg_value = COOP_ptr1->avg_value * num_occup_bands /
       (num_to_avg * details->num_KPOINTS * tot_num_K);

     fprintf(output_file,"#END CURVE\n");
     COOP_ptr1 = COOP_ptr1->next_type;
   }
 }




 /****************************************************************************
  *
  *                   Procedure gen_avg_COOPs
  *
  * Arguments: details: pointer to detail_type
  *              cell: pointer to cell type
  *          num_orbs: int
  *     avg_prop_info: pointer to avg_prop_info_type
  *        R_overlaps: pointer to hermetian_matrix_type
  *  orbital_ordering: pointer to K_orb_ptr_type
  *  orbital_lookup_table: pointer to int.
  *
  *
  * Returns: none
  *
  * Action:  Generates just the average COOPs and dumps them to the output file.
  *   This is for the alternate occupation analysis.
  *   Everything here is done the same way as in gen_COOP
  *
  ****************************************************************************/
 void gen_avg_COOPs(details,cell,num_orbs,avg_prop_info,R_overlaps,orbital_ordering,
                    orbital_lookup_table)
   detail_type *details;
   cell_type *cell;
   int num_orbs;
   avg_prop_info_type *avg_prop_info;
   hermetian_matrix_type R_overlaps;
   K_orb_ptr_type *orbital_ordering;
   int *orbital_lookup_table;
 {
   int i,j;
   int tot_num_orbs,num_COOPS;
   int num_to_avg;
   real COOP_accum,temp;
   real this_E,diff;
   real tot_num_K;
   real num_occup_bands;
   COOP_type *COOP_ptr1,*COOP_ptr2;

   tot_num_orbs = num_orbs * details->num_KPOINTS;

   /*******

     first determine some values that are needed for averaging purposes

   *******/

   /* the total k weighting that is used */
   tot_num_K = 0.0;
   num_occup_bands = 0.0;
   for( i=0; i<details->num_KPOINTS; i++){
     tot_num_K += details->K_POINTS[i].weight*details->K_POINTS[i].num_filled_bands;
     num_occup_bands += details->K_POINTS[i].num_filled_bands;
   }

   /******

     go through the list of COOPS and loop through the crystal orbitals for each one

   *******/
   COOP_ptr1 = details->the_COOPS;
   while(COOP_ptr1){


     COOP_ptr2 = COOP_ptr1;
     num_to_avg = 0;
     while(COOP_ptr2){
       num_to_avg++;
       COOP_ptr2 = COOP_ptr2->next_to_avg;
     }

     /* initialize the average value */
     COOP_ptr1->avg_value = 0.0;

     /* check for inversion of intercell vector */
     intercell_COOP_check(COOP_ptr1);

     /*******
       Add up all the contributions to the COOP that lie within DOS_DEGEN_TOL
       of each other. (This is just like how the DOS is generated).
     ********/
     i=0;
     while(i<tot_num_orbs){
       COOP_ptr2 = COOP_ptr1;
       COOP_accum = 0.0;
       while(COOP_ptr2){
         temp = details->K_POINTS[orbital_ordering[i].Kpoint].weight *
           eval_COOP(COOP_ptr2,details,cell,num_orbs,
                     &(avg_prop_info[orbital_ordering[i].Kpoint]),
                     R_overlaps,
                     &(orbital_ordering[i]),orbital_lookup_table);
         COOP_accum += temp;

         /*****
           accumulate the average value.  Since the calculation done in
           eval_COOP assumes that every orbital has 2 electrons in it, we
           need to divide temp by 2 and multiply by the number of
           electrons actually in the orbital.
         ******/
         COOP_ptr1->avg_value += orbital_ordering[i].occup * temp / 2.0;

         COOP_ptr2 = COOP_ptr2->next_to_avg;
       }
       this_E = (real)*(orbital_ordering[i].energy);
       j=i+1;

       if( j < tot_num_orbs ){
         diff = (real)*(orbital_ordering[i].energy) - (real)*(orbital_ordering[j].energy);
         while(fabs(diff) < DOS_DEGEN_TOL && j<tot_num_orbs){
           COOP_ptr2 = COOP_ptr1;
           while(COOP_ptr2){
             temp = details->K_POINTS[orbital_ordering[j].Kpoint].weight *
               eval_COOP(COOP_ptr2,details,cell,num_orbs,
                         &(avg_prop_info[orbital_ordering[j].Kpoint]),
                         R_overlaps,
                         &(orbital_ordering[j]),orbital_lookup_table);
             COOP_accum += temp;

             COOP_ptr1->avg_value += orbital_ordering[j].occup * temp / 2.0;

             COOP_ptr2 = COOP_ptr2->next_to_avg;
           }
           j++;

           if( j < tot_num_orbs ){
             diff = (real)*(orbital_ordering[i].energy) - (real)*(orbital_ordering[j].energy);
           }
         }
       }
       i = j;
     }

     /* now correct the accumulated value to make it the AVERAGE value */
     COOP_ptr1->avg_value = COOP_ptr1->avg_value * num_occup_bands /
       (num_to_avg * details->num_KPOINTS * tot_num_K);
     COOP_ptr1 = COOP_ptr1->next_type;
   }
 }

/***********************************************
*                                              *
* Function : intercell_COOP_check              *
*                                              *
* Arguments: COOP  ptr. to COOP_type           *
*                                              *
* Action   : sets intercell COOP contributions *
*            and determines intercell vector   *
*                                              *
***********************************************/

void intercell_COOP_check(COOP)

COOP_type *COOP;
{
  int temp;

  /* do we need to invert the cell vector ? */
  if((COOP->cell.z == 0.0)&&(COOP->cell.y <= 0.0)){
    if(!((COOP->cell.y == 0.0)&&(COOP->cell.x >= 0.0))){
      /* reverse COOP contrib's */
      temp = COOP->contrib1;
      COOP->contrib1 = COOP->contrib2;
      COOP->contrib2 = temp;
      /* invert xy plane cell vector components */
      COOP->cell.x *= -1.0;
      COOP->cell.y *= -1.0;
    }
  }
  else if(COOP->cell.z < 0){
    /* reverse COOP contrib's */
    temp = COOP->contrib1;
    COOP->contrib1 = COOP->contrib2;
    COOP->contrib2 = temp;
    /* invert cell vector */
    COOP->cell.x *= -1.0;
    COOP->cell.y *= -1.0;
    COOP->cell.z *= -1.0;
  }
}

