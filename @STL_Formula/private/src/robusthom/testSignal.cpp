#include "signal.h"

int main() {
	Signal *y, *z;
	Signal::iterator i;

	y=new Signal(37.5,1.,100);
	y->push_back(Sample(12,37.5,0.5));
	y->push_back(Sample(14,38.5,0.5));

	z=new Signal(39.5,10,30);

	i=y->begin();
	while(i->time != 14 && i != y->end()) i++;
	
	std::cout << "intersection de y avec z au temps: " << i->timeIntersect(z->front()) << std::endl;

	//std::cout << y->begin()->valueAt(12) << std::endl;
	
	std::cout << *y;
	y->simplify();
	std::cout << *y;

	//std::cout << *z;	
}
