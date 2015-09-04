#include "mex.h"
#include "math.h"
#define EPS (1e-12)

void mexFunction(int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ) 
{
    
  double *f_in = (double *)mxGetPr(prhs[0]);
  int m = mxGetM(prhs[0]);
  int n=  mxGetN(prhs[0]);

  double dt = mxGetScalar(prhs[1]);
  double *f_out;
  int i,j;
  double cum_sum = EPS;
  int is_pos;
  
  plhs[0] = mxCreateDoubleMatrix(m,n, mxREAL);
  f_out =  mxGetPr(plhs[0]);
    
  if (f_in[0] < 0)
    is_pos = 0;
  else 
    is_pos =1;
    
  for(i = 1; i<n; i++){    
    switch (is_pos) {
    case 1:
      if ( f_in[i] < 0) {
	cum_sum = -EPS;
	is_pos=0;
      }
      else
	cum_sum++;
      break;
    case 0:
      if ( f_in[i] < 0)
	cum_sum--;
      else {
	cum_sum =EPS;
	is_pos = 1;
      }     
      break;
    default:
      break;
    }
    f_out[i] = cum_sum*dt;
  } 
}
