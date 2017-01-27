#ifndef BREACH_H
#define BREACH_H

#define SKIP_CMAC 

#include <cvodes/cvodes.h>           /* main integrator header file */
#include <cvodes/cvodes_dense.h>     /* use CVDENSE linear solver */
#include <cvodes/cvodes_band.h>      /* use CVBAND linear solver */
#include <cvodes/cvodes_diag.h>      /* use CVDIAG linear solver */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fct. and macros */
#define Ith(x,i) NV_Ith_S(x,i)
#include <sundials/sundials_types.h> /* definition of realtype */
#include <sundials/sundials_math.h>  /* contains the macros ABS, SQR, and EXP */
#include <math.h>
#include "mex.h"
#include "cvm.h"
#include "nvm.h"
#include "mextools.h"
#include "random/uniform.h"
#include "random/normal.h"
#include "param_set.h"

#ifndef SKIP_CMAC
  #include "CmacVF.h"
#endif

using namespace blitz;
using namespace ranlib;

typedef enum {TRAJ_RAS, TRAJ_TEST_NOT_PASSED, TRAJ_TEST_PASSED, TRAJ_REACHED_SSTATE, TRAJ_OSCILLATES} traj_status; 

class trajectory {

 public:
  
  Array1D *current_x;
  Array1D *p0;
  Array1D *time;
  Array2D *X;
  Array2D *U;
  Array2D *Expa;
  Array2D *XS; 
  int length; 
  int status; // O by default

  /* Used for RL */ 

  Array1D *BellRes;
  Array1D *iV;
  Array1D *Vtilde;
  double perf;
 
  /* Constructor */
  
  trajectory() {
    length = 0;   
    time = new Array1D(0);
    X = new Array2D(0,0);
    U = new Array2D(0,0);
    Expa = new Array2D(0,0);
    XS = new Array2D(0,0);
    BellRes = new Array1D(0);
    iV = new Array1D(0);
    Vtilde = new Array1D(0);
    perf =0;
    status = 0; 
  };

  void ComputeTraj(Array1D& tspan);
  void ComputeTraj(Array1D& tspan, int (trajectory::*test_function)());
  void ComputeTraj(Array1D& tspan, vector<int>& idx_u, Array1D& tin, Array2D& uval);
  void ComputeTrajSensi(Array1D& tspan); 
  void ComputeExpa(const Array1D& epsi, Array1D &ExpaMax);

  void affiche();
  void write_mxTraj(mxArray* & mxTraj);
  int test1();

  ~trajectory() { 
    delete time;
    delete X;
    delete U;
    delete Expa;
    delete XS;
    delete BellRes;
    delete iV;
    delete Vtilde;
  };
};

ostream& operator<<(ostream& os, const trajectory& traj);

class FdataCommon {
 public:

  /* Used for f */
  
  realtype* u;
  realtype* p;
  int dimu;
  int dimx;
  int dimp;  // nb of parameters
  int dimg;  // nb of root functions
  
  /* Used for g */

  realtype *xp;
  realtype *delta;

  /* Used for sensitivity related computation  */

  int ns; // number of sensitivities to compute
  N_Vector *xsJump;


  /* used for RL */
#ifndef SKIP_CMAC
  CmacVF *C;
#endif
  
};

mxArray* Traj2mxStruct(Array<trajectory*,1> trajArray);
int g_delta(realtype t, N_Vector y, realtype *gout, void *g_data);
void FreeFdata(void * &f_data, mxArray * &mxData);
void GetU(void * f_data, void *&res);

int giDQ(int ig,realtype t, N_Vector x, realtype* gx, realtype *dgi,void *g_data);
int ComputeSensiJump(int ig, realtype t, N_Vector x, N_Vector *xS, void* f_data);

int Jac(long int N, DenseMat J, realtype t,
               N_Vector x, N_Vector fx, void *jac_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);



#endif
