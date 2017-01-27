#include "mex.h"
#include "math.h"

#define MIN(a,b) (a<=b?a:b)
#define MAX(a,b) (a>=b?a:b)

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[] ) 
{
    
    double *phi1 = (double *)mxGetPr(prhs[0]);
    double *phi2 = (double *)mxGetPr(prhs[1]);
    double *minphi1_win;
    int i_inf = (int) mxGetScalar(prhs[2]); 
    int i_sup = (int) mxGetScalar(prhs[3]);  /* -1 means +\infty*/
  
    int n1 =  mxGetN(prhs[0]);
    int n = mxGetN(prhs[1]);
    double * f_out;
    int i,win,j;
    double a_max,a;
    double *a_buffer;
    int a_current_ind, jc, size_buf;
    plhs[0] = mxCreateDoubleMatrix(1,n1, mxREAL);
    f_out =  mxGetPr(plhs[0]);
  
 
    if (i_sup==-1) {

        if (i_inf==0) {
            if (n1 != n)
                mexErrMsgTxt("inputs should have same dimensions");

            a_max = MIN(phi1[n-1],phi2[n-1]);
            i = 1;
            f_out[n-i] = a_max;
      
            for (i=2 ; i<=n; i++ ){
                a = MIN(a_max, phi1[n-i]);
                a_max = MAX(MIN(phi1[n-i],phi2[n-i]), a);
                f_out[n-i] = a_max;    
            } 
        }
        else /* [ti, \infty] */
            {
                /* 	mexPrintf("\n [ti infty] \n\n "); */
                minphi1_win = (double *) mxGetPr(prhs[4]);
                a_max = phi1[n1-1];
                i=1; 
                for (i=2;i<i_inf;i++) {
                    a_max = MIN(a_max, phi1[n1-i]);
                    f_out[n1-i] = a_max;
                    /* 	  mexPrintf("%g ", a_max); */
                    /* 	  if (i%15==0) */
                    /* 	    mexPrintf("\n"); */
                }
                /* 	mexPrintf("\n begin real stuff...\n\n "); */
                for (i=1; i<=n; i++ ){
                    a = MIN(a_max, phi1[n-i]);
                    a_max = MAX(MIN(minphi1_win[n-i], phi2[n-i]), a);
                    f_out[n-i] = a_max;    
                    /* 	  mexPrintf("%g ", a_max); */
                }	
                /* 	mexPrintf("\n"); */
            }
    }
    else /* [ti, tf] : this might be improvable .... */
        {
            /*   mexPrintf("\n [ti tf] \n\n "); */
            minphi1_win = (double *) mxGetPr(prhs[4]);
            win = i_sup-i_inf+1;
            a_buffer = (double *) mxCalloc(win, (win)* sizeof(double));
            a_max = phi1[n1-1];
            i=1; 	

            /* deal with tail smaller than ti */ 
            for (i=2;i<i_inf;i++) {
                a_max = MIN(a_max, phi1[n1-i]);
                f_out[n1-i] = a_max;
                /* 	mexPrintf("%g ", a_max); */
                /* 	if (i%15==0) */
                /* 	  mexPrintf("\n"); */
            }
      
            /* fill buffer */
            a_current_ind =0;
            size_buf =0;
            /*       mexPrintf("\n filling buffer...\n\n "); */
            for (i=1; i<=n; i++ ){ 
	
                a_max = MIN(minphi1_win[n-i], phi2[n-i]);
                a_buffer[a_current_ind]= a_max;
                size_buf= MIN(size_buf+1, win);
	
                for (j=1; j<size_buf; j++) { /* damn, here I lose */ 
                    jc= (j+a_current_ind)%(size_buf);
                    a_buffer[jc] = MIN(a_buffer[jc], phi1[n-i]);
                    if (a_buffer[jc]>=a_max)
                        a_max = a_buffer[jc];
                }
	
                a_current_ind =  (a_current_ind+1)%(win); 	  	  	
                f_out[n-i] = a_max;    
                /* 	mexPrintf("%g ", a_max); */
                /* 	if (i%15==0) */
                /* 	  mexPrintf("\n"); */
            }
            /*        mexPrintf("\n");      */
        }    
}
