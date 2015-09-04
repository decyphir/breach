#include "robustness.h"
#include <cstdlib>

double random_double()
{
      double r = (double)rand()/(double)RAND_MAX;
      return r;
}

void randomValues(double * V, double delta, int n) {
	for(int i=0; i < n; i++) {
		V[i]=random_double()*delta;
	}
}
void randomTimes(double * T, double delta, int n) {
	double prev=0;	
	for(int i=0; i < n; i++) {
		T[i]=random_double()*delta + .02 * delta + prev;
		prev=T[i];
	}
}


int main() {
	Signal *x, *y, *z;
	double *xV, *xT, *yV, *yT;

	int n=10; //size of test
	std::cout << "début calcul, avec n=" << n << std::endl;

	//parameters for random generator
	const double delta_x=10;
	const double delta_y=20;
	const double delta_t=1;


	//INIT
	xV=new double[n];
	yV=new double[n];
	xT=new double[n];
	yT=new double[n];

	randomValues(xV, delta_x, n);
	randomValues(yV, delta_y, n);
	randomTimes(xT, delta_t, n);
	randomTimes(yT, delta_t, n);

	x=new Signal(xT, xV, n);
	y=new Signal(yT, yV, n);

	delete[] xV;
    delete[] xT;
    delete[] yV;
    delete[] yT;

	std::cout << "signal x: " << *x;
	std::cout << "signal y: " << *y;


	z=new Signal(*x);
	z->shift(-100);
	
	std::cout << "signal z: " << *z;
	std::cout << "signal z défini entre " << z->beginTime << " et " << z->endTime << std::endl;
}
