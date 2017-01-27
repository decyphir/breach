#ifndef TD_H
#define TD_H


#include "reach.h"


#define S_G      RCONST(1.0)
#define MAXCOST  RCONST(2.0)

struct td_options {
  double lambda;
  double eta;
  double gamma;
  Array1D delta;
  
  td_options() {
    lambda=-1.;
    eta= .1;
    gamma= .9;
    delta = Array1D(1);
    delta = 0;
  }
};

struct td_stats {
  int     Av_nb_new;
  double  Av_Perf;
  double  Av_BellRes;
 
  td_stats() {
    Av_nb_new=0;
    Av_Perf=0.;
    Av_BellRes=0.;
  }
};

int CVM_ComputeCosts(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
int CVM_TD(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

void ComputeCosts(Array1D tspan, Array2D pts, Array2D &Xf, Array1D &Costs,Array<trajectory*,1>  &trajArray);
void TD(Array1D tspan, Array2D pts);

#endif

