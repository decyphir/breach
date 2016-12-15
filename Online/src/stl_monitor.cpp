#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <limits>
#include "tools.h"
#include "transducer.h"
#include "signal.h"
#include "stl_driver.h"
#include <GetOpt.h>

using namespace std;
using namespace stl_transducer;

void read_trace(const string &trace_file_name, trace_data &data) {

	std::ifstream file;
	file.open(trace_file_name.c_str());

	if (!file.is_open()) {
		cout << "trace file " << trace_file_name << " not found" << endl;
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
}

int main(int argc, char** argv) {

	// default trace file (use -t filename to override)
	string sol_path =
			"/Users/alex/workspace/MOOC149/InterestingSolutions/Correct/Garvit/";
	string env = "env1";
	string trace_file_name = sol_path + env + "/trace.txt";

	/* Command line options:
	 *
	 * 	 -p sets parser debugger on
	 *   -s sets scanning debugger on
	 *   -t filename uses filename as trace file
	 */

	bool debug_parsing(false), debug_scanning(false);
	bool random_trace(false); int len_trace(0);

	for (int ai = 1; ai < argc; ++ai) {
		if (argv[ai] == std::string("-p")) {
			debug_parsing = true;
		} else if (argv[ai] == std::string("-s")) {
			debug_scanning = true;
		} else if (argv[ai] == std::string("-t")) {
			ai++;
			trace_file_name = string(argv[ai]);
		} else if (argv[ai] == std::string("-r")) {
			random_trace = true;
			ai++;
			len_trace = atoi(argv[ai]);
		}
	}

	// instantiates the Driver
	trace_data data;

	if (!random_trace)
		read_trace(trace_file_name, data);
	else
		data = rand_trace_data(len_trace, 23);

	Driver stl_driver = Driver(data);
	stl_driver.trace_parsing = debug_parsing;
	stl_driver.trace_scanning = debug_scanning;

	// read a file with expressions (last command line argument is stl file)

	std::fstream infile(argv[argc - 1]);
	if (!infile.good()) {
		cout << "Could not open file: " << argv[argc - 1] << endl;
		return 0;
	}

	 transducer * phi;
	 double rob;
	 Signal z;

	 bool result = stl_driver.parse_stream(infile, argv[argc - 1]);
	 //cout << "file parsed, result=" << result <<  endl;

	if (result) {
		string phi_name;
		double rob_ref;

		// loop over all formulas defined in the file
		for (auto elem: stl_driver.formula_map ) {
			phi_name = elem.first;
			phi = elem.second->clone();
			rob= phi->compute_robustness();
			z =  phi->z;
			rob_ref = stl_driver.param_map["rob_"+phi_name];
			cout << endl << phi_name << endl << "rob=" << rob << "  rob_ref=" << rob_ref;

			if (fabs(rob-rob_ref)> numeric_limits<double>::epsilon())
				cout << " -> PROBLEM" << endl;
			else
				cout << endl;
		}
	}
	cout << "done." << endl;
	return 0;

}

