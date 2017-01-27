/*
 * -----------------------------------------------------------------
 * $Revision: 1.11 $
 * $Date: 2010-10-26 10:35:36 $
 * -----------------------------------------------------------------
 * Programmer: donze @ eecs berkeley edu
 * -----------------------------------------------------------------
 * Copyright ? uh.  
 * -----------------------------------------------------------------
 * MEX implementation for CVODES Matlab interface.
 * -----------------------------------------------------------------
 */



#include <string.h>
#include <stdlib.h>
#include "cvm.h"
#include "nvm.h"

#include "dynamics.h"
#include "param_set.h"

#ifndef SKIP_CMAC
#include "reach.h"
#include "td.h"
#else 
#include "traj.h"
#endif

/*
 * ---------------------------------------------------------------------------------
 * Definitions for global variables (declared in cvm.h)
 * ---------------------------------------------------------------------------------
 */

extern booleantype cvm_fsa;       /* Forward sensitivity? */
extern cvm_CVODESdata cvm_Cdata;  /* CVODES data */

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

    50 - get f values
    51 - get jacobian values
    61 - compute simple trajectories 
    93 - compute trajectories with sensitivities
    90 - compute only expansion, doesn't store trajectories
    92 - compute trajectories with expansion

    63 - compute delta-spaced trajectories
    91 - compute delta-spaced trajectories with expansion

    60 - compute trajectories with costs
    62 - temporal differences routine

    95 - refine parameter sets
    96 - reachability using expansion

    40 - finalize
  */

  mode = (int)mxGetScalar(prhs[0]);

  int Ns;
  mexUnlock();

  switch(mode) {
  case 1:
    if (cvm_Cdata) {
      /* a previous pb was initialized, we must clear  memory */
      CVM_Free(nlhs, plhs, nrhs-1, &prhs[1]);
      CVM_final();
    }
    CVM_init();
    f_data = new Fdata(); 
    CVM_Initialization(0, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 2:
    if (cvm_Cdata == NULL)
      mexErrMsgTxt("Init Solver first.");
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
    
  case 50:
    if (cvm_Cdata == NULL)
      mexErrMsgTxt("Init solver first.");
    CVM_GetF(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 51:
    if (cvm_Cdata == NULL)
      mexErrMsgTxt("Init solver first.");
    CVM_GetJf(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 61:    
    if (cvm_Cdata == NULL)
      mexErrMsgTxt("Init solver first.");
    CVM_ComputeTraj(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 93:
    if (cvm_Cdata == NULL)
      mexErrMsgTxt("Init solver first.");
    if (!cvm_fsa)
      mexErrMsgTxt("Init sensitivity first.");

    ((Fdata*) f_data)->xsJump = N_VCloneVectorArray_Serial(Ns,y);           

    CVM_ComputeTrajSensi(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 90:
    if (cvm_Cdata == NULL)
      mexErrMsgTxt("Init solver first.");
    if (!cvm_fsa)
      mexErrMsgTxt("Init sensitivity first.");

    /* stores sensitivity jump */
 
   ((Fdata*) f_data)->xp = (realtype *) mxMalloc(DIMX*sizeof(realtype)); 
   ((Fdata*) f_data)->delta = (realtype *) mxMalloc(DIMX*sizeof(realtype));
  
    for (int i=0; i<DIMX ; i++) {      
      ((Fdata*) f_data)->xp[i] = 0.;
      ((Fdata*) f_data)->delta[i] = 0.;
    }
  
    ((Fdata*) f_data)->xsJump = N_VCloneVectorArray_Serial(Ns,y);        

    CVM_ComputeExpansion(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 92:
    if (cvm_Cdata == NULL)
      mexErrMsgTxt("Init solver first.");
    if (!cvm_fsa)
      mexErrMsgTxt("Init sensitivity first.");

    ((Fdata*) f_data)->xsJump = N_VCloneVectorArray_Serial(Ns,y);        

    CVM_ComputeExpTraj(nlhs, plhs, nrhs-1, &prhs[1]);
    break;


#ifndef SKIP_CMAC

  case 63:    
    if (cvm_Cdata == NULL)
      mexErrMsgTxt("Init solver first.");

    ((Fdata*) f_data)->xp = (realtype *) mxMalloc(DIMX*sizeof(realtype)); 
    ((Fdata*) f_data)->delta = (realtype *) mxMalloc(DIMX*sizeof(realtype));
    
     for (int i=0; i<DIMX ; i++) {
    
      ((Fdata*) f_data)->xp[i] = 0.;
      ((Fdata*) f_data)->delta[i] = 0.;
      
     }
    CVM_CmacComputeTraj(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

    
  case 91:
    if (cvm_Cdata == NULL)
      mexErrMsgTxt("Init solver first.");
    if (!cvm_fsa)
      mexErrMsgTxt("Init sensitivity first.");
    
    /* stores sensitivity jump */

    ((Fdata*) f_data)->xsJump = N_VCloneVectorArray_Serial(Ns,y);        
    ((Fdata*) f_data)->xp = (realtype *) mxMalloc(DIMX*sizeof(realtype)); 
    ((Fdata*) f_data)->delta = (realtype *) mxMalloc(DIMX*sizeof(realtype));
  
    for (int i=0; i<DIMX ; i++) {      
      ((Fdata*) f_data)->xp[i] = 0.;
      ((Fdata*) f_data)->delta[i] = 0.;
    }
    
    CVM_CmacComputeExpTraj(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 60:
    if (cvm_Cdata == NULL)
      mexErrMsgTxt("Init solver first.");
    CVM_ComputeCosts(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
    
  case 62:    
    if (cvm_Cdata == NULL)
      mexErrMsgTxt("Init solver first.");
    CVM_TD(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

    
    
  case 95:
    CVM_Refine(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 96:
    if (cvm_Cdata == NULL)
      mexErrMsgTxt("Init solver first.");
    if (!cvm_fsa)
      mexErrMsgTxt("Init sensitivity first.");
    
    ((Fdata*) f_data)->xsJump = N_VCloneVectorArray_Serial(Ns,y);        

    CVM_Reach(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

#endif

  /* ---------------------------------- */

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

void FreeFdata(void * &_f_data, mxArray * &_mxData){

  Fdata* data = (Fdata*) _f_data;
  data->p = NULL;
  data->u = NULL;
  
#ifndef SKIP_CMAC
  data->C = NULL;
#endif
  /* free memory related to sensitivity */

  if (cvm_fsa) {
    N_VDestroyVectorArray(data->xsJump,Ns);    
  }

  delete data;

}

#ifndef _GETU
void GetU(void *_f_data, void *&_data_dest){}
#endif

/* root function for delta-spaced trajectories */

int g_delta(realtype _t, N_Vector _y, realtype *_gout, void *_g_data) {

  Fdata* _data = (Fdata*) _g_data;
  
  if (_data->dimg>0) {    
    g(_t,_y,_gout,_g_data);
  }
  for(int i=0 ; i< DIMX; i++) {
    _gout[i] = abs(NV_Ith_S(_y,i)-_data->xp[i])-_data->delta[i];
    
    if ((_gout[i]>0)&&(_gout[i]< (-_data->delta[i]/100.)))
      _gout[i+_data->dimg]=0;
  }
  
  return 0;

}



