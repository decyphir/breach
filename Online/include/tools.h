/*
 * tools.h
 *
 *  Created on: Mar 31, 2014
 *      Author: alex
 */

#ifndef TOOLS_H_
#define TOOLS_H_
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;
typedef vector<vector<double> > trace_data;

/**  string to double conversion */
inline bool s_to_d(std::string const& s,
                             double &x )
{
  std::istringstream i(s);
  x =0;
  char c;
  if (!(i >> x) || ( i.get(c)))
    return false;
  else
	return true;
}

/** double to string conversion */
inline string d_to_s(
	double const &x)
{
	std::ostringstream o;
	o << x;
	return o.str();
}

/** Returns a m x n random trace */
trace_data rand_trace_data( int m, int n);

const std::string current_date_time();

void print(const trace_data&);
bool read_trace(const string &trace_file_name, trace_data &data);

#endif /* TOOLS_H_ */
