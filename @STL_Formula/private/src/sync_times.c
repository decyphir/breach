#include "mex.h"
#include "math.h"

int sync_times( double * t1, const int n1, double *t2,  const int n2, double * t_out);

void mexFunction(int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ) 
{
    
  double *t1 = (double *)mxGetPr(prhs[0]);
  int n1 = mxGetN(prhs[0]);
  double *t2 = (double *)mxGetPr(prhs[1]);
  int n2=  mxGetN(prhs[1]);
  double *t_out = (double*) mxMalloc((n1+n2)*sizeof(double)); 
  int nout;
  nout = sync_times(t1, n1, t2, n2, t_out);
  int i; 
  double * mx_tout;

  plhs[0] = mxCreateDoubleMatrix(1,nout, mxREAL);
  mx_tout = mxGetPr(plhs[0]);

  for(i=0;i<nout; i++ ) {
    mx_tout[i] = t_out[i];
  }
  
}

int sync_times( double * t1, const int n1, double *t2,  const int n2, double * t_out) {
  
  int i,i1, i2, nout;

  i =0;  i1 = 0; i2 = 0;  nout =0;

  while(1) {
    
    if (t1[i1]<t2[i2]) {
      t_out[i] =t1[i1];
      i1++;
      i++;nout++;
    }
    else 
      if (t1[i1]>t2[i2]) {
	t_out[i] =t2[i2];
	i2++;
	i++;nout++;
      }
      else
	{
	  t_out[i] =t2[i2];
	  i2++;
	  i1++;
	  i++;nout++;
	}
    
    if (i1+1 == n1) {
      for (; i2<n2;i2++) {
	t_out[i] = t2[i2];
	i++; nout++;
      }
      return nout;
    }

    if (i2+1 == n2) {
      for (; i1<n1;i1++) {
	t_out[i] = t1[i1];
	i++; nout++;
      }
      return nout;
    }          
  }
} 
