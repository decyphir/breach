#include "robustness.h"
#include "mex_routines.h"
#define INF std::numeric_limits<double>::infinity()

/*
 * RobustUntil(...)
 *
 * TOWRITE
 */
void mexFunction(int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[])
{
	/* read matlab inputs */

#ifdef DEBUG__
    mexPrintf("Entering RobustUntil\n");
#endif

	// manage empty input
	if(mxGetN(prhs[0])==0 || mxGetN(prhs[2])==0) // if the first signal is empty
	{
		//TODO manage meaningfully empty inputs
		mexWarnMsgTxt("RobustUntil: got an empty signal! Will answer empty signal.");
		plhs[0] = mxCreateDoubleMatrix(1, 0, mxREAL);
		plhs[1] = mxCreateDoubleMatrix(1, 0, mxREAL);
		return;
	}

	// check data coherence
	if(mxGetN(prhs[0])!=mxGetN(prhs[1]))
	{
		mexWarnMsgTxt("RobustUntil: lengths of time steps and signal are different.");
	}

	// signals
	Signal *y1, *y2;

	// std::cout << "N:" <<mxGetN(prhs[0])  << std::endl;
	y1 = new Signal(mxGetPr(prhs[0]), mxGetPr(prhs[1]), mxGetN(prhs[0]));
	y2 = new Signal(mxGetPr(prhs[2]), mxGetPr(prhs[3]), mxGetN(prhs[2]));

	// interval
	double *I = mxGetPr(prhs[4]);
	double a = I[0];
	double b = I[1];

	/* compute until robustness */
	//std::cout << "Computing until_[" << a << ":" << b << "] y" << std::endl;

	Signal *yunt;

	if(b==INF)
	{
		if (a==0)
			yunt = computeUntil(y1,y2);
		else
		{
			Signal * yalw1 = computeBoundedGlobally(y1, a);
			Signal * yuntmp = computeUntil(y1, y2);
			yuntmp->shift(-a);
			yunt = computeAnd(yalw1, yuntmp);
			delete yalw1;
            delete yuntmp;
		}
	}
	else
		yunt = computeTimedUntil(y1, y2, a, b );

	int N = yunt->size();
	plhs[0] = mxCreateDoubleMatrix(1, N+1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1, N+1, mxREAL);

	writeSignal(*yunt, plhs[0], plhs[1]);
	delete y1;
    delete y2;
    delete yunt;
}
