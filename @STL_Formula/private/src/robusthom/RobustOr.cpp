#include "robustness.h"
#include "mex_routines.h"

/*
 * RobustOr(T1,V1,T2,V2)
 *
 * T1 is the time steps for first signal
 * V1 is values of the first signal
 * T2 is the time steps for second signal
 * V2 is values of the second signal
 *
 * size of T1 and V1 must be the same and must be != 1
 * size of T2 and V2 must be the same and must be != 1
 */
void mexFunction(int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[])
{
	/* read matlab inputs */

#ifdef DEBUG__
	cout << "Entering RobustEv" << endl;
#endif

	// manage empty signal
	if(mxGetN(prhs[0])==0) // if the first signal is empty
	{
		plhs[0] = mxDuplicateArray(prhs[2]); // may be empty
		plhs[1] = mxDuplicateArray(prhs[3]);
		return;
	}
	if(mxGetN(prhs[2])==0) // if the second signal is empty
	{
		plhs[0] = mxDuplicateArray(prhs[0]);
		plhs[1] = mxDuplicateArray(prhs[1]);
		return;
	}

	// check data coherence
	if(mxGetN(prhs[0])!=mxGetN(prhs[1]) || mxGetN(prhs[2])!=mxGetN(prhs[3]))
	{
		mexWarnMsgTxt("RobustOr: lengths of time steps and signal are different.");
	}

	// signals
	Signal *y1, *y2;

	// std::cout << "N:" <<mxGetN(prhs[0])  << std::endl;

	y1 = new Signal(mxGetPr(prhs[0]), mxGetPr(prhs[1]), mxGetN(prhs[0]));
	y2 = new Signal(mxGetPr(prhs[2]), mxGetPr(prhs[3]), mxGetN(prhs[2]));

	// std::cout << "\n\n y2:" << *y2 << std::endl;

	/* compute or robustness */
	// std::cout << "\n\n Computing or :" << std::endl;

	Signal *yor = computeOr(y1, y2);
	// std::cout << "yor:" << *yor << std::endl;

	int N = yor->size();
	plhs[0] = mxCreateDoubleMatrix(1, N+1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1, N+1, mxREAL);

	writeSignal(*yor, plhs[0], plhs[1]);
	delete y1;
    delete y2;
    delete yor;
}
