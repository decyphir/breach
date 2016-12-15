#include "stdafx.h"
#include <transducer.h>
#include <string>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include "tools.h"
#include <exception>
#include "signal_expr.h"
#include <cmath>

using namespace std;

namespace CPSGrader {

/* Unary operators on signals */

double abs_transducer::compute_robustness() {

	// update child robustness
	child->compute_robustness();

	auto iter = child->z.begin();
	
	for (; iter != child->z.end();
			iter++) {
		double t= (*iter).time;
		double v = (*iter).value;
		z.appendSample(t, fabs(v));
	}
	return z.front().value;
}

/* Binary signal operators - all operators assume that left and right operand have the
 * same number of samples.
 */
double plus_transducer::compute_robustness() {

	// update children robustness
	childL->compute_robustness();
	childR->compute_robustness();

	auto iterL = childL->z.begin();
	auto iterR = childR->z.begin();

	for (; iterL != childL->z.end();
			iterL++, iterR++) {
		double t= (*iterL).time;
		double vL = (*iterL).value;
		double vR = (*iterR).value;
		z.appendSample(t, vL+vR);
	}
	return z.front().value;
}

double minus_transducer::compute_robustness() {

	childL->compute_robustness();
	childR->compute_robustness();

	auto iterL = childL->z.begin();
	auto iterR = childR->z.begin();

	for (; iterL != childL->z.end();
			iterL++, iterR++) {
		double t= (*iterL).time;
		double vL = (*iterL).value;
		double vR = (*iterR).value;
		z.appendSample(t, vL-vR);
	}
	return z.front().value;
}

double mult_transducer::compute_robustness() {

	childL->compute_robustness();
	childR->compute_robustness();

	auto iterL = childL->z.begin();
	auto iterR = childR->z.begin();

	for (; iterL != childL->z.end();
			iterL++, iterR++) {
		double t= (*iterL).time;
		double vL = (*iterL).value;
		double vR = (*iterR).value;
		z.appendSample(t, vL*vR);
	}
	return z.front().value;
}

}
