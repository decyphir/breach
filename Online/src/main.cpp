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
#include "stl_driver.h"

using namespace std;
using namespace CPSGrader;

/** This is CPSGrader default executable. It reads a test plan and output a report.
* For each test in the test plan, a data file is read, and a set of stl tests are
* executed. */
int main(int argc, char** argv) {

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
		} else if (argv[ai] == std::string("-r")) {
			random_trace = true;
			ai++;
			len_trace = atoi(argv[ai]);
		}
	}

	// instantiates the STLDriver
	STLDriver stl_driver = STLDriver();
	stl_driver.trace_parsing = debug_parsing;
	stl_driver.trace_scanning = debug_scanning;

	// read a file with expressions (last command line argument is stl file)
	std::fstream infile(argv[argc - 1]);
	if (!infile.good()) {
		cout << "Could not open file: " << argv[argc - 1] << endl;
		return 0;
	}

	bool result = stl_driver.parse_stream(infile, argv[argc - 1]);

	if (result) {
		stl_driver.print();
		while (!stl_driver.trace_test_queue_empty()) {
			string trace_file_name = stl_driver.get_next_trace_test_env();
			bool success = stl_driver.read_trace_file(trace_file_name);
			if (success)
				stl_driver.run_next_trace_test();
			else {
				stl_driver.report = stl_driver.report + "Skipped test " + stl_driver.get_next_trace_test()->id +"\n";
				stl_driver.pop_next_trace_test();
			}
		}
		cout << "Report:" << endl;
		cout << stl_driver.report << endl;
		cout << "Test log:" << endl;
		cout << stl_driver.test_log << endl;
	}
	else
		cout << "Problem parsing test_plan file." << endl;

	cout << "done." << endl;

	return 0;

}

