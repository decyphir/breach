#ifndef SIGNAL_H
#define SIGNAL_H

#include <deque>
#include <iostream>
#include <limits>
#include "mex.h"
#define BOTTOM -std::numeric_limits<double>::infinity()
#define TOP std::numeric_limits<double>::infinity()

//#define DEBUG__

using namespace std;

inline 
double fmin(double a,double b)
{
  return (a<b)?a:b;
}

inline 
double fmax(double a,double b)
{
  return (b<a)?a:b;
}


class Point {
public:
	double time;
	double value;

	Point() { }
	Point(double t, double v) : time(t), value(v) { }

	friend std::ostream & operator<<(std::ostream &, const Point &);
};


class Sample : public Point {
public:
	double derivative;

	Sample() { }
	Sample(double t, double v, double d) : Point(t,v), derivative(d) { }	
    Sample constant() const;
   
	double valueAt(const double &) const;
	double timeIntersect(const Sample &) const;
	double area(const double &) const;

	friend std::ostream & operator<<(std::ostream &, const Sample &);
};

inline
Sample Sample::constant() const {
  return Sample(time, value, 0);
}

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
	double beginTime;
	double endTime;

	Signal() { }	
	Signal(Sequence);        
    Signal(double, double, int); 
	Signal(double *, double *, int); //create continuous signal from array of sampling points (time, value) with linear interpolation
    
	void simplify(); //remove sampling points where (y,dy) is continuous.
	void resize(double, double, double); //restricts/extends the signal to [s,t) with default value v where not defined
	void shift(double); //shifts the signal of delta_t time units
    Signal* reverse();
    int push_front(Sample);
    void print() const;
	friend std::ostream & operator<<(std::ostream &, const Signal &);

};

inline
double Sample::area(const double & t) const {
  if (t > time) {
    return (value + valueAt(t)) * (t - time) / 2;
  }else{
    return 0;
  }
}

#endif
