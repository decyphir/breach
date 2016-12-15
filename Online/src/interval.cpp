#include "stdafx.h"
#include <interval.h>
#include "tools.h"
#include "transducer.h"
using namespace std;

namespace CPSGrader {

interval::interval() {
	begin = -1;
	end = -1;
}

// tries to convert string into doubles, if fails, keep the string assuming this is a parameter or expression
interval::interval(string b, string e) {
	if (!s_to_d(b, begin))
		begin_str = b;
	if (!s_to_d(e, end))
		end_str = e;
}

interval::interval(const interval& that) {
	copy(that);
}
void interval::copy(const interval& that) {
	begin_str= that.begin_str;
	end_str = that.end_str;
	begin = that.begin;
	end = that.end;
}

interval interval::operator=(const interval& that) {
	if (this != &that)
		copy(that);
	return *this;
}

string interval::to_string() const {
	ostringstream o;
	o << "[";
	if (begin_str.empty())
		o << begin;
	else
		o << begin_str;
	o << ",";
	if (end_str.empty())
		o << end;
	else
		o << end_str;
	o << "]";
	return o.str();
}

std::ostream& operator<<(std::ostream& os, const interval& I) {
	I.print(os);
	return os;
}
;

}
