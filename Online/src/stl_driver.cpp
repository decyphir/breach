// $Id: stl_driver.cpp 378 2014-04-27 17:47:40Z donze $
/** \file driver.cc Implementation of the example::Driver class. */

#include "stdafx.h"
#include <fstream>
#include <sstream>

#include "stl_driver.h"
#include "transducer.h"
#include "stl_scanner.h"

namespace CPSGrader {

map<string, token_type> STLDriver::reserved = map<string,token_type>(); // filled in stl_scanner.lpp

    //map<string, short> transducer::signal_map = map<string,short>();

    //map<string, double> transducer::param_map = map<string, double>();

double Signal::BigM = 100.;
//vector<vector<double>> *transducer::trace_data_ptr = nullptr;

STLDriver::STLDriver()
: trace_scanning(false),
  trace_parsing(false),
  report(""),
  lexer(nullptr),
  nb_test_total(0),
  nb_test_pos(0),
  error_flag(false)
{
}

STLDriver::STLDriver(trace_data _trace)
: trace_scanning(false),
  trace_parsing(false),
  report(""),
  lexer(nullptr),
  nb_test_total(0),
  nb_test_pos(0),
  error_flag(false)
{
	data = _trace;
};

bool STLDriver::parse_stream(std::istream& in, const std::string& sname)
{
	streamname = sname;
	Scanner scanner(&in);
	scanner.set_debug(trace_scanning);
	this->lexer = &scanner;

	Parser parser(*this);
	parser.set_debug_level(trace_parsing);

	return (parser.parse() == 0);

}

bool STLDriver::parse_file(const std::string &filename)
{
	std::ifstream in(filename.c_str());
	if (!in.good()) return false;
	return parse_stream(in, filename);
}

bool STLDriver::parse_string(const std::string &input, const std::string& sname)
{
	std::istringstream iss(input);
	return parse_stream(iss, sname);
}

/** Error handling with associated line number. This can be modified to
 * output the error e.g. to a dialog box. */
void STLDriver::error(const class location& l,
		const std::string& m)
{
	std::cerr << l << ": " << m << std::endl;
}

/** General error handling. This can be modified to output the error
 * e.g. to a dialog box. */
void STLDriver::error(const std::string& m)
{
	std::cerr << m << std::endl;
}

/** clear assigned formulas and trace_test_queue */
void STLDriver::clear() {

	for (auto formula = formula_map.begin(); formula != formula_map.end(); formula++) {
		if (formula->second != 0) {
			delete formula->second;
			formula->second = 0;
		}
	}

	for (auto stest = stl_test_map.begin(); stest != stl_test_map.end(); stest++) {
		if (stest->second.formula != 0) {
			delete stest->second.formula;
			stest->second.formula = 0;
		}
	}


	formula_map.clear();
	param_map.clear();
	stl_test_map.clear();

	while (!trace_test_queue.empty()) {
		transducer * formula = get_next_formula();
		while (formula!=nullptr) {
			delete formula;
			trace_test_queue.front().tests.pop_front();
			formula = get_next_formula();
		}
		trace_test_queue.pop_front();
	}
	report.clear();
	test_log.clear();
	error_flag = false;

	nb_test_total =0;
	nb_test_pos = 0;

}



void STLDriver::print(ostream &out) const {

	out << "\nAssigned formulas:" << endl;
	out << "-------------------" << endl;

		for (auto formula = formula_map.begin(); formula != formula_map.end(); formula++){
			out << formula->first << ":" << endl;
			out << *(formula->second) << endl;
		}

	out << "\nDefault Parameters:" << endl;
	out << "---------------------" << endl;

	for (auto param = param_map.begin(); param != param_map.end(); param++){
		out << param->first << "=";
		out << param->second << endl;
	}

	out << "\nTrace tests:" <<  endl;
	out << "---------------" << endl;

	string indent = "    ";
	for (auto it= trace_test_queue.begin(); it!= trace_test_queue.end(); it++) {

		out << (*it).id << ": ";
		out << (*it).env <<", ";
		out << "SimTime: " << (*it).sim_time << " Visu:" << (*it).visu << endl;
		for (auto its=  (*it).tests.begin(); its!= (*it).tests.end(); its++) {
			// transducer->param_map =  param_map; FIXME probably needs a set_param for transducers
			out << indent << its->test_id << endl;
			for (auto elem = its->param_map.begin(); elem != its->param_map.end(); elem++)
				out << indent << indent << elem->first << "=" <<  elem->second << endl;

			out << indent << indent <<  *(*its).formula << endl;
		}
		out << endl;
	}
	out << endl << endl;
}

/** compute robustness for all formulas defined in the driver and write results in files */
void STLDriver::dump() {

	//transducer::param_map = param_map; FIXME

	for (auto formula = formula_map.begin(); formula != formula_map.end(); formula++){
		formula->second->trace_data_ptr = &data;
		formula->second->init_horizon();
		formula->second->compute_robustness();
		(formula->second->z).dumpFile(formula->first + ".out");
	}
}

/** run all stl tests for the next trace_test and write results in report */
bool STLDriver::run_tests(){

	error_flag = false;
	string indent("    ");
	if (!trace_test_queue.empty()) {

		const stl_test * curr_test= get_next_stl_test();
		transducer * formula = get_next_formula();

		report.append("stl_test " + trace_test_queue.front().id + "\n");

		// Initialize global parameters
		//transducer->param_map =  param_map; FIXME 

		while (formula != nullptr) {

			// initialize local parameters FIXME 
			//for (auto elem = curr_test->param_map.begin(); elem != curr_test->param_map.end(); elem++){
			//	transducer->param_map[elem->first] = elem->second;
			//}

			formula->trace_data_ptr = &data;
			formula->init_horizon();
			double rob = formula->compute_robustness();

			test_log.append(curr_test->test_id + ": " + d_to_s(rob) + "\n");
			report.append(indent + curr_test->test_id + ": " + d_to_s(rob));

			if (rob > 0.){
				nb_test_pos++;
				report += "  --> PB\n";
				if (!curr_test->report_positive.empty())
					report += indent + curr_test->report_positive + "\n";
				if (curr_test->is_error)
					error_flag = true;
			}
			else {
				report += "  --> OK\n";
				if (!curr_test->report_negative.empty())
					report += indent + curr_test->report_negative + "\n";
			}
			
			//delete formula;
			trace_test_queue.front().tests.pop_front();
			formula = get_next_formula();
			curr_test = get_next_stl_test();
		}
	}
	report += "\n";
	trace_test_queue.pop_front();
	return error_flag;
}


void STLDriver::add_trace_test(const string &test_id, const string &trace_cfg, double sim_time, bool visu) {
	trace_test T;
	T.id = test_id;
	T.env = trace_cfg;
	T.sim_time = sim_time;
	T.visu = visu;
	trace_test_queue.push_back(T);
}

/** Note: test inserted in the stl_test_map have no local parameters */
void STLDriver::add_stl_test(const string & test_id) {

	nb_test_total++;
	if (!trace_test_queue.empty()) {
		map<string, stl_test>::iterator elem;
		if ((elem = stl_test_map.find(test_id)) != stl_test_map.end()) {
			trace_test_queue.back().tests.push_back(elem->second);
			trace_test_queue.back().tests.back().formula = (elem->second).formula->clone();
			trace_test_queue.back().tests.back().param_map.clear();

		}
	}
}

void STLDriver::add_stl_test(const string & test_id,const map<string, double> &local_param_map) {

	nb_test_total++;
	if (!trace_test_queue.empty()) {
		map<string, stl_test>::iterator elem;
		if ((elem = stl_test_map.find(test_id)) != stl_test_map.end()) {

//			cout << "Test " << elem->first;
//			if (!local_param_map.empty()) {
//				cout << " with local param:" << endl;
//				for (auto elemp = local_param_map.begin(); elemp != local_param_map.end(); elemp++)
//					cout << elemp->first << "=" <<  elemp->second << endl;
//			}
//			else
//				cout << endl;

			trace_test_queue.back().tests.push_back(elem->second);
			trace_test_queue.back().tests.back().formula = (elem->second).formula->clone();
			trace_test_queue.back().tests.back().param_map = local_param_map;

		}
	}
}

void STLDriver::add_stl_test(const string & test_id,const map<string, double> &local_param_map, transducer* formula, const string &report_pos,
		const string &report_neg, bool is_err) {

	nb_test_total++;
	if (!trace_test_queue.empty()) {
//		cout << "Inserting test " << test_id;
//		if (!local_param_map.empty()) {
//			cout << " with local param:" << endl;
//			for (auto elemp = local_param_map.begin(); elemp != local_param_map.end(); elemp++)
//				cout << elemp->first << "=" <<  elemp->second << endl;
//		}
//		else
//			cout << endl;

		stl_test_map[test_id] = stl_test(test_id, local_param_map, formula->clone(), report_pos,report_neg, is_err);
		stl_test_map[test_id].param_map = local_param_map;
		stl_test_map[test_id].formula = formula->clone();
		trace_test_queue.back().tests.push_back(stl_test_map[test_id]);
		trace_test_queue.back().tests.back().formula = formula->clone();

	}
}

double STLDriver::test_formula(const string & phi_in) {
	//transducer->param_map = param_map;

	if (data.empty()){
		cout << "Empty data" << endl;
		return 0.;
	}
	string funky_name = "f_u_n_k_y_p_h__i_n_a_m_e"; // seriously?
	string str_to_parse = funky_name + ":=" + phi_in;
	if (parse_string(str_to_parse)){
		transducer * phi = formula_map[funky_name]->clone();
		formula_map.erase(funky_name);
		phi->trace_data_ptr = &data;
		phi->init_horizon();
		double rob = phi->compute_robustness();
		delete phi;
		return rob;
	}
	else {
		cout << "Couldn't parse formula: " << phi_in << endl;
		return 0.;
	}
}


void STLDriver::print_trace(ostream &os){
		for (auto ii = data.begin(); ii != data.end(); ii++){
			for (auto jj = (*ii).begin(); jj != (*ii).end(); jj++) {
				os << *jj << " ";
			}
			os << endl;
        }
}

void STLDriver::print_trace(){
    print_trace(cout);
}

void STLDriver::dump_trace_file(const string & filename){
	fstream ofs;
	ofs.open(filename.c_str(), std::ofstream::out);
	if (ofs.is_open()) {
		for (auto ii = data.begin(); ii != data.end(); ii++){
			for (auto jj = (*ii).begin(); jj != (*ii).end(); jj++) {
				ofs << *jj << " ";
			}
			ofs << endl;
		}
		ofs.close();
	}
	else {
		cout << "Couldn't open file " << filename.c_str() << " for writing signal" << endl; // TODO implement exception
	}
}

bool STLDriver::dump_test_log_file(const string & filename){

	fstream ofs;
	ofs.open(filename.c_str(), std::ofstream::out);
	if (ofs.is_open()) {
		ofs << test_log << endl;
		ofs << "Number of positive STL tests: " << nb_test_pos << "/" << nb_test_total;
		ofs.close();
		return true;
	}
	else {
		cout << "Couldn't open file " << filename.c_str() << " for writing signal" << endl; // TODO implement exception
		return false;
	}

}


} // namespace CPSGrader
