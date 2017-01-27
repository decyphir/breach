#include "mex.h"
#include "math.h"


void mexFunction(int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ) 
{
    
  double *f_in = (double *)mxGetPr(prhs[0]);
  int m = mxGetM(prhs[0]);
  int n=  mxGetN(prhs[0]);
  double *f_out;
  double * min_current = (double *) mxMalloc( m*sizeof(double));
  int i,j;
  
  plhs[0] = mxCreateDoubleMatrix(m,n, mxREAL);
  f_out =  mxGetPr(plhs[0]);
  

  for(i=0; i<m;i++)
      min_current[i] = mxGetInf();
  
  for(i = 0; i<m; i++){
    for(j = 0; j<n; j++) {
      if (f_in[j*m+i]<min_current[i])
	min_current[i] = f_in[j*m+i];
      f_out[j*m+i]= min_current[i];    
    }
  } 

  mxFree(min_current);
}


