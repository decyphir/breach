// $Id: stl_driver.h 176 2016-04-05 03:50:15Z donze $ 	
/** \file stl_driver.h Declaration of the stl_transducer::Driver class. */

#ifndef STL_DRIVER_H
#define STL_DRIVER_H

#include <string>
#include <vector>
#include <map>
#include <deque>
#include "transducer.h"
#include "stl_parser.h"
#include "tools.h"

typedef CPSGrader::Parser::token token;
typedef CPSGrader::Parser::token_type token_type;

namespace CPSGrader {

/** A class encapsulating an STL formula plus feedback messages to be written in
 * the report based on its satisfaction or violation status.  */
class stl_test {
public:
	string test_id;
	map<string, double> param_map;
	transducer* formula;
	string report_positive;
	string report_negative;
	bool is_error;

	stl_test(){
		is_error = false;
		formula = nullptr;
	};

	stl_test(const string & _test_id,const map<string, double> &_param_map, transducer* _formula, const string &_report_pos,
			const string &_report_neg, bool _is_err):
				test_id(_test_id),
				param_map(_param_map),
				formula(_formula),
				report_positive(_report_pos),
				report_negative(_report_neg),
				is_error(_is_err)
	{
	};

	~stl_test() {
		delete formula;
	}

};

struct trace_test {
	string id;
	string env;
	double sim_time;
	bool visu;
	deque<stl_test> tests;
};

/** The Driver class brings together all components. It creates an instance of
 * the Parser and Scanner classes and connects them. Then the input stream is
 * fed into the scanner object and the parser gets it's token
 * sequence. Furthermore the driver object is available in the grammar rules as
 * a parameter. Therefore the driver class contains a reference to the
 * structure into which the parsed data is saved. */
class STLDriver {
public:

	/// DATA

	/** enable debug output in the flex scanner */
	bool trace_scanning;

	/** enable debug output in the bison parser */
	bool trace_parsing;

	/** stream name (file or input stream) used for error messages. */
	std::string streamname;

	/** Reserved words for STL parser */
	static map<string, token_type> reserved;

	/** parameters in formulas */
	map<string, double> param_map;

    /** signals in formulas */
	map<string, int> signal_map;

	/** formulas defined by the driver */
	map<string, transducer*> formula_map;
    
	/** tests (sets of formulas)  */
	map<string, stl_test> stl_test_map;

	/** data array - time is first column */
	trace_data data;
	deque<trace_test> trace_test_queue;

	string report;
	string test_log;
	int nb_test_pos;
	int nb_test_total;
	bool error_flag;

	/// CONSTRUCTORS

	/** construct a new parser driver context */
	STLDriver();

	/** construct a parser driver with data */
	STLDriver(trace_data _trace);

	~STLDriver() {};

	/** Create a new trace_tests structure and add it to the queue */
	void add_trace_test(const string &test_id, const string &trace_cfg, double sim_time, bool visu);

	/** adds an stl_test from the stl_test_map to the current trace_test with the local params value
	 * does nothing if the stl_test does not exist */
	void add_stl_test(const string &test_id);
	void add_stl_test(const string &test_id,  const map<string,double> &local_param_map);

	/** adds an stl_test to the current trace_tests - does nothing if the queue is empty */
	void add_stl_test(const string &test_id,  const map<string,double> &local_param_map, transducer* formula, const string &report_pos,
			const string &report_neg, bool is_err);

	/// PARSER

	/** Invoke the scanner and parser for a stream.
	 * @param in	input stream
	 * @param sname	stream name for error messages
	 * @return		true if successfully parsed
	 */
	bool parse_stream(std::istream& in, const std::string& sname =
			"stream input");

	/** Invoke the scanner and parser on an input string.
	 * @param input	input string
	 * @param sname	stream name for error messages
	 * @return		true if successfully parsed
	 */
	bool parse_string(const std::string& input, const std::string& sname =
			"string stream");

	/** Invoke the scanner and parser on a file. Use parse_stream with a
	 * std::ifstream if detection of file reading errors is required.
	 * @param filename	input file name
	 * @return		true if successfully parsed
	 */
	bool parse_file(const std::string& filename);

	// To demonstrate pure handling of parse errors, instead of
	// simply dumping them on the standard error output, we will pass
	// them to the driver using the following two member functions.

	/** Error handling with associated line number. This can be modified to
	 * output the error e.g. to a dialog box. */
	void error(const class location& l, const std::string& m);

	/** General error handling. This can be modified to output the error
	 * e.g. to a dialog box. */
	void error(const std::string& m);

	/** Pointer to the current lexer instance, this is used to connect the
	 * parser to the scanner. It is used in the yylex macro. */
	class Scanner* lexer;

	/// API

	/** runs the tests of the front of trace_tests_queue and pop them
	 * assumes that data has been filled with the appropriate trace data */
	bool run_tests();

	/** monitor a single formula requires data is not empty */
	double test_formula(const string &);

	/** gets next formula to test */
	inline transducer * get_next_formula() const {

		if (!trace_test_queue.front().tests.empty())
			return trace_test_queue.front().tests.front().formula;
		else
			return nullptr;
	}

	inline const stl_test * get_next_stl_test() const {

		if (!trace_test_queue.front().tests.empty())
			return &(trace_test_queue.front().tests.front());
		else
			return nullptr;
	}

	inline const trace_test * get_next_trace_test() const {
		if (!trace_test_queue.empty())
			return &(trace_test_queue.front());
		else
			return nullptr;
	}

	inline const bool trace_test_queue_empty() {
		return trace_test_queue.empty();
	}

	inline const string get_next_trace_test_env() {
		if (!trace_test_queue_empty())
			return trace_test_queue.front().env;
		else
			return string();
	}

	inline void pop_next_trace_test() {
		if (!trace_test_queue_empty())
			trace_test_queue.pop_front();
	}

	// TODO rid of run_tests altogether
	inline void run_next_trace_test() {
		run_tests();
	}

	/// UTILITY FUNCTIONS

	/** clear assigned formulas and trace_test_queue */
	void clear();

	/** Prints driver to stream output */
	void print(ostream &out) const;

	/** Prints driver to standard output */
	inline void print() const { print(cout); };


	/** Read a trace file */
	inline bool read_trace_file(string trace_file_name) {
		data.clear();
		return read_trace(trace_file_name, data);
	}

	/** dump all assigned formulas satisfaction function to a file */
	void dump();
	void dump_trace_file(const string&);
	bool dump_test_log_file(const string&);
    void print_trace(ostream &os);
    void print_trace();
}
;

} // namespace CPSGrader

#endif // STL_DRIVER_H
