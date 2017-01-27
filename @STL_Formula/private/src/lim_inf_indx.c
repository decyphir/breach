#include "mex.h"
#include "math.h"


void mexFunction(int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ) 
{
    
  double *f_in = (double *)mxGetPr(prhs[0]);
  int m = mxGetM(prhs[0]);
  int n=  mxGetN(prhs[0]);
  double *f_out, *indx_out;

  double *min_current = (double*) mxMalloc(m*sizeof(double));
  double *indx_min_current= (double*) mxMalloc(m*sizeof(double)); 

  int i,j;
  
  plhs[0] = mxCreateDoubleMatrix(m,n, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(m,n, mxREAL);

  f_out =  mxGetPr(plhs[0]);
  indx_out = mxGetPr(plhs[1]);

  for(i=0; i<m;i++)
      min_current[i] = mxGetInf();
  
  for(i = 0; i<m; i++){
    for(j = n-1; j>=0; j--) {
      if (f_in[j*m+i]<min_current[i]){
	min_current[i] = f_in[j*m+i];
	indx_min_current[i] = (double) j;
      }
      f_out[j*m+i] = min_current[i];    
      indx_out[j*m+i] = indx_min_current[i]+1;
    }
  } 
  mxFree(min_current);
  mxFree(indx_min_current);

}


