/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2009-06-05 16:26:09 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials-x.y.z/src/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * MEX implementation for CVODES Matlab interface.
 * -----------------------------------------------------------------
 * -----------------------------------------------------------------
 * Modified for cvodesTB++ by: Alexandre Donzé @ imag.fr
 * -----------------------------------------------------------------
 */


#include <string.h>
#include <stdlib.h>
#include "cvm.h"
#include "nvm.h"

/*
 * ---------------------------------------------------------------------------------
 * Definitions for global variables (declared in cvm.h)
 * ---------------------------------------------------------------------------------
 */

booleantype cvm_quad;      /* Forward quadratures? */
booleantype cvm_quadB;     /* Backward quadratures? */
booleantype cvm_asa;       /* Adjoint sensitivity? */
booleantype cvm_fsa;       /* Forward sensitivity? */
booleantype cvm_mon;       /* Forward monitoring? */ 
booleantype cvm_monB;      /* Backward monitoring? */ 

cvm_CVODESdata cvm_Cdata = NULL;  /* CVODES data */
cvm_MATLABdata cvm_Mdata = NULL;  /* MATLAB data */
 
/*
 * ---------------------------------------------------------------------------------
 * Main entry point
 * ---------------------------------------------------------------------------------
 */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[] )
{
  int mode;
  /* 
     Modes:
     
     1 - initialize CVODES solver
     2 - initialize forward sensitivity calculations
     3 - initialize adjoint sensitivity calculations
     4 - initialize backward solver

    11 - reinitialize CVODES solver
    12 - reinitialize forward sensitivity calculations
    13 - toggle FSA off
    14 - reinitialize backward solver

    20 - solve problem
    21 - solve backward problem

    30 - get integrator stats
    31 - get backward integrator stats
    32 - extract data from cvode_mem
    33 - set one optional input at a time (NYI)

    40 - finalize
  */

  mode = (int)mxGetScalar(prhs[0]);

  mexUnlock();

  switch(mode) {
  case 1:
    if (cvm_Cdata != NULL) {
      /* a previous pb was initialized, we must clear  memory */
      CVM_Free(nlhs, plhs, nrhs-1, &prhs[1]);
      CVM_final();
    }
    CVM_init();
    CVM_Initialization(0, nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 2:
    CVM_SensInitialization(0, nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 3:
    CVM_AdjInitialization(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 4:
    CVM_InitializationB(0, nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 11:
    CVM_Initialization(1, nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 12:
    CVM_SensInitialization(1, nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 13:
    CVM_SensToggleOff(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 14:
    CVM_InitializationB(1, nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 20:
    CVM_Solve(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 21:
    CVM_SolveB(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 30:
    CVM_Stats(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 31:
    CVM_StatsB(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 32:
    CVM_Get(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 33:
    CVM_Set(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 40:
    CVM_Free(nlhs, plhs, nrhs-1, &prhs[1]);
    CVM_final();
    return;
  }
  
  /* do not call CVM_makePersistent after free */
  if (mode != 40) {
    CVM_makePersistent();
    mexLock();
  }
   return;
}
