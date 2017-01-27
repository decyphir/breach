#ifndef SIGNAL_H
#define SIGNAL_H

#include <deque>
#include <iostream>
#include <fstream>
#include <limits>
#include "mex.h"
//#include <math.h>
//#include <algorithm>
// TODO maybe make TOP and BOTTOM static signal attribute ...
#define BOTTOM (-Signal::BigM) //-std::numeric_limits<double>::infinity()
#define TOP (Signal::BigM) //std::numeric_limits<double>::infinity()

using namespace std;

/*
double fmin(double a,double b)
{
  return (a<b)?a:b;
}

double fmax(double a,double b)
{
  return (b<a)?a:b;
}
*/
namespace CPSGrader {

class Point {
public:
	double time;
	double value;

	Point(): time(0.), value(0.) { }
	Point(double t, double v) : time(t), value(v) { }

	friend std::ostream & operator<<(std::ostream &, const Point &);
};

class Sample : public Point {
public:
	double derivative;

	Sample() {derivative =0;}
	Sample(double t, double v, double d) : Point(t,v), derivative(d) { }	

	double valueAt(const double &) const;
	double timeIntersect(const Sample &) const;

	friend std::ostream & operator<<(std::ostream &, const Sample &);
};

inline
double Sample::valueAt(const double & t) const {
	return value + derivative * (t - time);
}

inline
double Sample::timeIntersect(const Sample & point) const {
	return (value - point.value + (point.derivative * point.time) - (derivative * time)) / (point.derivative - derivative);
}

//piecewise-constant, right-continuous signals
class Sequence : public std::deque<Point> { 
public:
	Sequence() { }

	friend std::ostream & operator<<(std::ostream &, const Sequence &);
};

//piecewise-linear, right-continuous signals
class Signal : public std::deque<Sample> {
public:

	static double BigM;

	double beginTime;
	double endTime;
	Signal(): beginTime(0.), endTime(0.) { };
	Signal(double, double, int); 
	Signal(double *, double *, int); //create continuous signal from array of sampling points (time, value) with linear interpolation
	void appendSample(double, double);
	void appendSample(double, double, double);
    void appendSignal(Signal);
    void simplify(); //remove sampling points where (y,dy) is continuous.
	void resize(double, double, double); //restricts/extends the signal to [s,t) with default value v where not defined
	void shift(double); //shifts the signal of delta_t time units
    void addLastSample(); // add a sample at endTime. 
	void removeInf();


	// write signal to file
	void dumpFile(const string filename) const {
		fstream ofs;
		ofs.open(filename.c_str(), std::ofstream::out);
		if (ofs.is_open()) {
			for (auto i = begin(); i != end(); i++) {
				ofs << *i << std::endl;
			}
			ofs.close();
		}
		else {
			cout << "Couldn't open file " << filename.c_str() << " for writing signal" << endl; // TODO implement exception
		}
	}

    /** Robustness computation functions - all these methods first clear content */

    /// Signal becomes neg of signal argument
    void compute_not(const Signal&); 

    /// Signal becomes conjunction (min) of the two arguments
    void compute_and(const Signal&, const Signal&);

    /// Signal becomes disjunction (max) of the two arguments
    void compute_or(const Signal&, const Signal&);

    /// Signal becomes implication (max of not x and y) of the two arguments
    void compute_implies(const Signal&, const Signal&);

    /// Signal becomes robust eventually of argument
    void compute_eventually(const Signal&); 

    /// ev_[0, a) 
    void compute_bounded_eventually(const Signal&, double);
    
    /// ev_[a,b]
    void compute_timed_eventually(const Signal&, double, double);

    /// Signal becomes robust globally of argument
    void compute_globally(const Signal&); 

    /// alw_[0, a) 
    void compute_bounded_globally(const Signal&, double);
    
    /// alw_[a,b]
    void compute_timed_globally(const Signal&, double, double);
    
    /// Signal becomes robust untimed until 
    void compute_until(const Signal &, const Signal &);
    
    /// until_[a,b]
    void compute_timed_until(const Signal&, const Signal&, double, double );

    /** TODO Auxiliary routines (should be privatized probably) */
    
    /// computes max lemire window of size a of argument signal 
    void compute_plateau_max(const Signal &, double);
    void compute_plateau_min(const Signal &, double);


	friend std::ostream & operator<<(std::ostream &, const Signal &);
};
}
#endif
