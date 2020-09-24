#include "robustness.h"
#include "mex_routines.h"
#define INF std::numeric_limits<double>::infinity()

using namespace std;

/*
 * RobustAvEv(T, V, I)
 *
 * T is time steps
 * V is values at time steps t
 * I is the eventually interval
 *
 * size of T and v must be the same and their size must be != 1
 */
void mexFunction(int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[])
{

#ifdef DEBUG__
	std::cout << "Entering RobustAvEv" << endl;
#endif

	/* read matlab inputs */
	if(mxGetN(prhs[0])==0) // manage empty signal
	{
		plhs[0] = mxCreateDoubleMatrix(1, 0, mxREAL);
		plhs[1] = mxCreateDoubleMatrix(1, 0, mxREAL);
		return;
	}

	// check data coherence
	if(mxGetN(prhs[0])!=mxGetN(prhs[1]))
	{
		mexWarnMsgTxt("RobustAvEvRight: lengths of time steps and signal are different.");
	}

	// signal
	Signal *y;
	y = new Signal(mxGetPr(prhs[0]), mxGetPr(prhs[1]), mxGetN(prhs[0]));
	//std::cout << "y created" << endl;
	// interval
	double *I = mxGetPr(prhs[2]);
	double a = I[0];
	double b = I[1];

	/* compute ev robustness */
    //  std::cout << "Computing ev_[" << a << ":" << b << "] y" << std::endl;

	if (a>0)
		y->shift(-a);

	Signal *yev;
	//double area;

	if (b-a < 0)
	{
		mexWarnMsgTxt("RobustAvEvRight: in interval [a, b], b-a < 0 ! Treated as 0.");
		yev = y;
	}
	else
	{
		if (b-a == 0)
			yev = y;
		else
		{

			if (b-a >= (y->endTime- y->beginTime))
			{
				yev = computeEventually(y);
			}
			else
			{
				yev = computeRightAvEventually(y, b-a);
        //     	area = computeAverageEventuallyArea(y, b-a);
			}
		}
	}
	// std::cout << "yev:" << *yev << " size:" << yev->size() << std::endl;

	int N = yev->size();
	plhs[0] = mxCreateDoubleMatrix(1, N+1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1, N+1, mxREAL);
    //plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);

    //double *mxArea = mxGetPr(plhs[2]);
    //mxArea[0] = area;
	writeSignal(*yev, plhs[0], plhs[1]);
	delete y;
	delete yev;
}
