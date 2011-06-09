#ifndef TRAJ_H
#define TRAJ_H

#include "dynamics.h"

/*
 * ---------------------------------------------------------------------------------
 * Redability replacements for cvodes global variables
 * ---------------------------------------------------------------------------------
 */

#define cvode_mem   (cvm_Cdata->cvode_mem)
#define bp_data     (cvm_Cdata->bp_data) 
#define y           (cvm_Cdata->y) 
#define yQ          (cvm_Cdata->yQ) 
#define yS          (cvm_Cdata->yS) 
#define N           (cvm_Cdata->N) 
#define Nq          (cvm_Cdata->Nq) 
#define Ng          (cvm_Cdata->Ng) 
#define Ns          (cvm_Cdata->Ns) 
#define Nd          (cvm_Cdata->Nd) 
#define Nc          (cvm_Cdata->Nc) 
#define ls          (cvm_Cdata->ls) 
#define pm          (cvm_Cdata->pm) 
#define ism         (cvm_Cdata->ism) 
#define cvadj_mem   (cvm_Cdata->cvadj_mem) 
#define interp      (cvm_Cdata->interp) 
#define yB          (cvm_Cdata->yB) 
#define yQB         (cvm_Cdata->yQB) 
#define NB          (cvm_Cdata->NB) 
#define NqB         (cvm_Cdata->NqB) 
#define lsB         (cvm_Cdata->lsB) 
#define pmB         (cvm_Cdata->pmB) 
#define errmsg      (cvm_Cdata->errmsg)
#define f_data      (cvm_Cdata->f_data)
#define uround      (cvm_Cdata->uround)

#define itol (cvm_Cdata->itol)
#define reltol (cvm_Cdata->reltol)
#define Sabstol (cvm_Cdata->Sabstol)
#define NV_abstol (cvm_Cdata->NV_abstol)
#define yS0  (cvm_Cdata->yS0) 

int CVM_ComputeExpansion(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
int CVM_CmacComputeExpTraj(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
int CVM_ComputeExpTraj(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

int CVM_ComputeTrajSensi(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
int CVM_ComputeTraj(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

int CVM_GetF(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
int CVM_GetJf(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

void ComputeExpansion(Array2D pts, Array2D XS0, Array2D epsi, Array1D tspan, Array2D &Xf, Array2D &XSf, Array2D &ExpaMax);

void ComputeExpTraj(Array2D pts, Array2D XS0, Array2D epsi, Array1D tspan, Array<trajectory*,1> &trajArray, Array2D &Xf, Array2D &XSf, Array2D &ExpaMax);


#endif
