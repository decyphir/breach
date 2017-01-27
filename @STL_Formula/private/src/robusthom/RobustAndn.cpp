#include "robustness.h"
#include "mex_routines.h"

/*
 * RobustAndn(Tn,Vn)
 *
 * Tn is a cell array of time steps
 * Vn is a cell array of signal values
 *
 * Tn(i) and Vn(i) must have the same size and this size must be != 1
 */
void mexFunction(int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[])
{
	/* read matlab inputs */

#ifdef DEBUG__
	cout << "Entering RobustAndn" << endl;
#endif

	int num_y = mxGetNumberOfElements(prhs[0]); // total number of signals
	// std::cout << "num_y:" <<mxGetNumberOfElements(prhs[0])  << std::endl;

	/* check for empty signal */
	int i_init; // index of first non empty signal
	for(i_init=0 ; i_init<num_y && mxGetN(mxGetCell(prhs[0],i_init))==0 ; ++i_init)
		;

	if(i_init==num_y) // all signals are empty
	{
		plhs[0] = mxCreateDoubleMatrix(1, 0, mxREAL);
		plhs[1] = mxCreateDoubleMatrix(1, 0, mxREAL);
		return;
	}

	if(mxGetN(mxGetCell(prhs[0],i_init))!=mxGetN(mxGetCell(prhs[1],i_init)))
	{
		mexWarnMsgTxt("RobustAndn: lengths of time steps and signal are different.");
	}

	/* initialize the first signal */
	Signal *y, *ytmp, *yandn;
	yandn = new Signal(mxGetPr(mxGetCell(prhs[0],i_init)), mxGetPr(mxGetCell(prhs[1],i_init)), mxGetN(mxGetCell(prhs[0],i_init)));
	// std::cout << "\n\n y:" << *y << std::endl;

	/* compute and robustness */
	for(int i=i_init+1;i<num_y;++i)
	{
		if(mxGetN(mxGetCell(prhs[0],i))!=0)
		{
			if(mxGetN(mxGetCell(prhs[0],i))!=mxGetN(mxGetCell(prhs[1],i)))
			{
				mexWarnMsgTxt("RobustAndn: lengths of time steps and signal are different.");
			}

			y = new Signal(mxGetPr(mxGetCell(prhs[0],i)), mxGetPr(mxGetCell(prhs[1],i)), mxGetN(mxGetCell(prhs[0],i)));
			// std::cout << "\n\n y:" << *y << std::endl;
			ytmp = yandn;
			yandn = computeAnd(ytmp, y);
			// std::cout << "\n\n yandn in loop:" << *yandn << std::endl;
			delete y;
            delete ytmp;
		}
	}
	// std::cout << "yandn after computation:" << *yandn << std::endl;

	int N = yandn->size();
	plhs[0] = mxCreateDoubleMatrix(1, N+1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1, N+1, mxREAL);

	writeSignal(*yandn, plhs[0], plhs[1]);
	delete yandn;
}
