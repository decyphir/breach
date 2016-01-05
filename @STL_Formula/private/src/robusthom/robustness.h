#include <cmath>
#include "signal.h"


/*
 * MAIN PROCEDURES
 */

// these functions allocate a new Signal which address is returned
Signal * computeNot(Signal *);
Signal * computeAnd(Signal *, Signal *);
Signal * computeOr(Signal *, Signal *);
Signal * computeEventually(Signal *);
Signal * computeUntil(Signal *, Signal *);
Signal * computeTimedUntil(Signal *, Signal *, double, double);


/*
Modified for AvMTL
*/
Signal * constZero(double beginT, double endT);
Signal * computeLeftAvEventually(Signal *, double);
Signal * computeRightAvEventually(Signal *, double);
Signal * computeLeftAvGrobally(Signal *, double);
Signal * computeRightAvGrobally(Signal *, double);

/* 
 * SUBROUTINES
 */

//adds to front of the signal passed as pointer the robustness of a conjuction of a linear segment, represented as (Sample, double) with an arbiratry signal 
void computeSegmentAnd(Signal *, const Sample &, double, Signal::reverse_iterator &); 

//specialization of the former for constant segment, represented as (double, double, double)
void computeConstantSegmentAnd(Signal *, double, double, double, Signal::reverse_iterator &);

//adds to front of the signal passed as pointer the robustness of a conjuction for arbitrary signals passed as iterators between times passed as doubles
void computePartialAnd(Signal *, Signal::reverse_iterator &, Signal::reverse_iterator &, double, double); 


//disjunction: identical to conjunction, operator "<" switched with operator ">"
void computeSegmentOr(Signal *, const Sample &, Signal::reverse_iterator &, double, double);
void computeConstantSegmentOr(Signal *, double, Signal::reverse_iterator &, double, double);
void computePartialOr(Signal *, Signal::reverse_iterator &, Signal::reverse_iterator &, double, double);


//adds to front of the signal passed as pointer the robustness of an eventually for an aribrary signal passed as iterator between times passed as doubles
void computePartialEventually(Signal *, Signal::reverse_iterator &, double, double);


//idem computeSegmentAnd, with an extra parameter for the last known value of the Until
void computeSegmentUntil(Signal *, const Sample &, double, Signal::reverse_iterator &, double);


// these functions work on the whole length of arbitrary signals and allocate a new Signal which address is returned.
// the second argument is a temorization constraint "a", such that we compute the robustness of E_[0;a] or G_[0;a] respectively
Signal * plateauMax(Signal *, double);
Signal * computeBoundedEventually(Signal *, double);

Signal * plateauMin(Signal *, double);
Signal * computeBoundedGlobally(Signal *, double);

// For AvMTL
//void addPoint(Signal* M, Sample newSample, double, double &);
Signal * plateauMaxAux(Signal *, double);
Signal * plateauMinAux(Signal *, double);

