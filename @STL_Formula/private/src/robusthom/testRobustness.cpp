#include "robustness.h"
#include "mex_routines.h"

/*
 * To compile, go into breach/@STL_Formula/private/src/robusthom and run
 * the following command line (replace .o by .obj if using Windows):
 *  mex testRobustness.cpp robustness.o signal.o mex_routines.o
 *  
 * Then launch:
 *  testRobustness
 *
 */

//int main() {
void mexFunction(int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[])
{
	Signal *x, *nx, *y, *ny, *z, *nz;

// TEST AND (robustAnd still does not manage singular value)
    double t[3] = {1,2,3};
    double v[3];
    v[0] = 39; v[1] = 1; v[2] = 100;
	x=new Signal(t,v,3);  // Signal(*T,*V,n)
    v[0] = 42; v[1] = 20; v[2] = 50;
	y=new Signal(t,v,3);
	z=computeAnd(x, y);
	std::cout << "calcul x AND y (test 1)" << std::endl;
	std::cout << *z << std::endl;
	
	delete x;
    delete y;

	x=new Signal(1,39,100); // Signal(T,V,n)
    x->push_back(Sample(99,39,0)); // Sample(time,value,derivative)
    x->endTime = 99;
	y=new Signal(1,37.5,100);
	y->push_back(Sample(12,37.5,0.5));
	y->push_back(Sample(14,38.5,0.5));
	y->push_back(Sample(17,40,-0.5));
    y->endTime = 17;

	std::cout << "calcul x AND y (test 2)" << std::endl;
	z=computeAnd(y, x);
	std::cout << *z << std::endl;
	delete z;

// TEST OR
	std::cout << "calcul x OR y par De Morgan" << std::endl;

	nx=computeNot(x);
	ny=computeNot(y);
	nz=computeAnd(nx, ny);
	z=computeNot(nz);
	std::cout << *z;
	delete z;

	std::cout << "calcul x OR y direct" << std::endl;
	z=computeOr(y, x);
	std::cout << *z << std::endl;
	delete z;


// TEST EVENTUALLY
	x=new Signal(1,37.5,100);
	x->push_back(Sample(10,50,-5));
	x->push_back(Sample(14,38.5,0.5));
	x->push_back(Sample(18,36,-0.5));
	x->push_back(Sample(80,20,0));
    x->endTime = 80;

	std::cout << "calcul E x" << std::endl;
	z=computeEventually(x);
	std::cout << *z << std::endl;

	delete z;

	std::cout << "calcul G x" << std::endl;
	nx=computeNot(x);
	nz=computeEventually(nx);
	z=computeNot(nz);
	std::cout << *z << std::endl;


// TEST1 UNTIL
	x=new Signal(1,37.5,100);
	x->push_back(Sample(10,50,-5));
	x->push_back(Sample(14,38.5,0.5));
	x->push_back(Sample(18,36,-0.5));
	x->push_back(Sample(80,20,0));
    x->endTime = 80;

	y=new Signal(10,39.5,30);
	y->push_back(Sample(12,40,0.5));
	y->push_back(Sample(26,30,0));
    y->endTime = 26;

	z=computeUntil(y,x);
	std::cout << *z << std::endl;


// TEST2 UNTIL
	x=new Signal();
	x->beginTime=0;
	x->endTime=12;
	x->push_back(Sample(0,3,3));
	x->push_back(Sample(0.5,4.5,-1));
	x->push_back(Sample(1,4,1));
	x->push_back(Sample(1.5,4.5,-1));
	x->push_back(Sample(2.5,4,1));
	x->push_back(Sample(4.5,6,0));
	x->push_back(Sample(5.5,6,-1./3.));
	x->push_back(Sample(8.5,5,0));
	x->push_back(Sample(10,7,-1));

	y=new Signal();
	y->beginTime=0;
	y->endTime=12;

	y->push_back(Sample(0,4,0));
	y->push_back(Sample(2,3,1.6));
	y->push_back(Sample(4.5,7,-2));
	y->push_back(Sample(5.5,5,2));
	y->push_back(Sample(6.5,6,1));
	y->push_back(Sample(7,6.5,-1));
	y->push_back(Sample(7.5,6,0));
	y->push_back(Sample(8,4,2));
	y->push_back(Sample(9,6,0));

	z=computeUntil(x,y);
	std::cout << *z << std::endl;	


// TEST BOUNDED EVENTUALLY AND BOUNDED GLOBALLY
	x=new Signal();
	x->beginTime=0;
	x->endTime=12;
	x->push_back(Sample(0,1,6));
	x->push_back(Sample(0.5,2,6));
	x->push_back(Sample(1,5,-0.25));
	x->push_back(Sample(3,4,0));
	x->push_back(Sample(5,4.5,2));
	x->push_back(Sample(5.5,6.5,-1));
	x->push_back(Sample(7,5,2));
	x->push_back(Sample(7.5,5.5,0.25));
	x->push_back(Sample(9.5,6,-1));
	x->push_back(Sample(11,3,1));
	x->push_back(Sample(11.5,2,1));

	y=computeBoundedEventually(x,1);
	std::cout << "Eventually x \n" << *y;

	z=computeBoundedGlobally(x,1);
	std::cout << "Globally x \n" << *z;	

	delete z;

	nx=computeNot(x);
	nz=computeBoundedEventually(nx,1);
	z=computeNot(nz);
	std::cout << "Globally x (par De Morgan) \n" << *z;


//TEST TIMED UNTIL
	x=new Signal();
	x->beginTime=0;
	x->endTime=12;
	x->push_back(Sample(0,3,3));
	x->push_back(Sample(0.5,4.5,-1));
	x->push_back(Sample(1,4,1));
	x->push_back(Sample(1.5,4.5,-1));
	x->push_back(Sample(2.5,4,1));
	x->push_back(Sample(4.5,6,0));
	x->push_back(Sample(5.5,6,-1./3.));
	x->push_back(Sample(8.5,5,0));
	x->push_back(Sample(10,7,-1));

	y=new Signal();
	y->beginTime=0;
	y->endTime=12;

	y->push_back(Sample(0,4,0));
	y->push_back(Sample(2,3,1.6));
	y->push_back(Sample(4.5,7,-2));
	y->push_back(Sample(5.5,5,2));
	y->push_back(Sample(6.5,6,1));
	y->push_back(Sample(7,6.5,-1));
	y->push_back(Sample(7.5,6,0));
	y->push_back(Sample(8,4,2));
	y->push_back(Sample(9,6,0));

	z=computeTimedUntil(x,y,1,2);
	std::cout << "x U_[1:2] y \n" << *z;	

}
