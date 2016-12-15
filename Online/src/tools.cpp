/*
 * tools.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: alex
 */
#include "stdafx.h"

#include <cstdlib>
#include <vector>
#include <string>
#include <ctime>
#include <stdio.h>
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using namespace std;

double random_double()
{
      double r = (double)rand()/(double)RAND_MAX;
      return r;
}

trace_data rand_trace_data( int m, int n) {
	trace_data rand_trace;

	double dt= 0.1;
	for(int i=0; i<m; i++) {
		vector<double> v;
		v.push_back(i*dt);
		for(int j=1; j<n; j++) {
			v.push_back(random_double());
		}
		rand_trace.push_back(v);
	}
	return rand_trace;

}

const std::string current_date_time() {
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%Y-%m-%d-%H-%M-%S", &tstruct);

	return buf;
}

void print(const trace_data & data) {

	vector<double>::const_iterator iter_value;
	vector<vector<double>>::const_iterator iter_sample;

	for(iter_sample=data.begin();iter_sample != data.end(); iter_sample++) {
		for(iter_value = iter_sample->begin(); iter_value != iter_sample->end(); iter_value++ ) {
			cout << *iter_value << " ";
		}
		cout << endl;
	}
}


bool read_trace(const string &trace_file_name, vector< vector<double> > &data)
{

	std::ifstream file;
	file.open(trace_file_name.c_str());

	if (!file.is_open()) {
		cout << "trace file " << trace_file_name << " not found" << endl;
		return false;
	}

	string line;
	while (!std::getline(file, line, '\n').eof()) {
		istringstream reader(line);
		vector<double> lineData;
		while (!reader.eof()) {
			double val;
			reader >> val;
			lineData.push_back(val);
		}
		data.push_back(lineData);
	}
	file.close();
	return true;
}
