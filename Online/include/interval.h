#ifndef __INTERVAL_H
#define __INTERVAL_H

#include <iostream>
#include <string>
#include "tools.h"

using namespace std;

namespace CPSGrader {

class interval {
public:
	double begin;
	double end;
	string begin_str;
	string end_str;

	interval();
	interval(const interval& that);
	interval operator=(const interval& that);
	void copy(const interval& that);

	// tries to convert string into doubles, if fails, keep the string assuming this is a parameter or expression
	interval(string b, string e);
	interval(double b, double e): begin(b), end(e) {};

	string to_string() const;

	void print() const {
		cout << to_string();
	};

	void print(std::ostream &os) const {
	    os << to_string();
	};
}
;

std::ostream& operator<<(std::ostream& os, const interval& I);

}

#endif
