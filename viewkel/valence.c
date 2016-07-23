/*******************************************************
*      Copyright (C) 1997, 1999 Greg Landrum
*
*  This file is part of yaehmop.
*
*   This is free software.
* 
*  Permission is granted to modify, or otherwise fold, spindle, and mutilate this
*    code provided all copyright notices are left intact.
*
*  This code may be distributed to your heart's content, in whatever form,
*    provided no fee is charged for the distribution, all copyright notices are
*    left intact, and the source is distributed (without fee) along with any
*    binaries to anyone who requests it.
*
*  There are, of course, no warranties at all on this program.
*
********************************************************************/


/******************************

  valence.c

  This file contains the implementation of the
  bond length->bond valence scheme based upon

  This file implements the bond valence method using
   the parameters and equations given in:
   M. O'Keeffe and N. E. Brese, JACS 113, 3226-3229 (1991)

  Created by gL November 1997

*******************************/
#include "viewkel.h"
#include "valence.h"

#ifdef INCLUDE_BOND_VALENCE
float calc_R0(float c1,float r1,float c2,float r2)
{
  float result;
  
  result = sqrt(c1)-sqrt(c2);
  result *= result;
  result *= r1*r2;
  result /= c1*r1 + c2*r2;
  result += r1+r2;
  
  return(result);
}

/****************************************************************************
 *
 *                   Procedure bond_length_to_bond_valence
 *
 * Arguments:  atom1,atom2: pointers to atom_type
 *                 length: float
 *        R0_val,valence: pointers to float
 *            
 * Returns: none
 *
 * Action: Calculates the bond valence between atoms 'atom1 and 'atom2
 *   using the bond length 'length.
 *
 *  If *R0_val is nonzero, it will be used to calculate the valence, which
 *   is returned in 'valence, otherwise the ci and ri parameters from the 
 *   atoms are used.  
 *
 *  On Return: 'RO_val contains the R0 value (calculated or passed in) 
 *   and 'valence contains the calculated value of the valence.
 *  
 ****************************************************************************/
void bond_length_to_bond_valence(atom_type *atom1,atom_type *atom2,
				 float length,float *R0_val,
				 float *valence)
{
  
  /* do we need to calculate R0? */
  if(*R0_val == 0.0){
   
    if( atom1->ci == 0.0 || atom2->ci == 0.0 ) *R0_val = 0.0;
    else{
      *R0_val = calc_R0(atom1->ci,atom1->ri,atom2->ci,atom2->ri);
    }
  } 


  /******

    if R0 is zero at this point, there is a problem.
    return a valence of 0

  *******/
  if(*R0_val == 0.0 ){
    *valence = 0.0;
    return;
  }
  *valence = (float)exp((double)(*R0_val-length)/(double)BVAL);
    
}
#endif

