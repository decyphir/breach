#include "stdafx.h"
#include "robustness.h"
#include <list>
#include <algorithm>

// needed sometimes when compiling under windows 
#undef min
#undef max

#define fmin(a,b) ((a<b)?a:b) 
#define fmax(a,b) ((b<a)?a:b) 

//#define DEBUG__

using namespace std;

namespace CPSGrader { 
void Signal::compute_not(const Signal &y) {
	clear();
	beginTime = y.beginTime;
	endTime = y.endTime;
	for (auto i = y.begin(); i != y.end(); i++) {
		push_back(Sample(i->time, -(i->value), -(i->derivative)));
	}
}

void Signal::compute_and(const Signal &x, const Signal &y) {

	clear();
	auto i = x.rbegin();
	auto j = y.rbegin();

	beginTime = fmax(x.beginTime, y.beginTime);
	endTime = fmin(x.endTime, y.endTime);

	while ((i->time >= endTime) && i != x.rend()-1)
		i++;
	while ((j->time >= endTime) && j != y.rend()-1)
		j++;

	computePartialAnd(this, i, j, this->beginTime, this->endTime);
	simplify();

}

void Signal::compute_or(const Signal &x, const Signal &y) {

	clear();
	Signal::const_reverse_iterator i = x.rbegin();
	Signal::const_reverse_iterator j = y.rbegin();

	beginTime = fmax(x.beginTime, y.beginTime);
	endTime = fmin(x.endTime, y.endTime);

	while ((i->time >= endTime) && i != x.rend()-1)
		i++;
	while ((j->time >= endTime) && j != y.rend()-1)
		j++;

	computePartialOr(this, i, j, beginTime, endTime);
	simplify();

}

void Signal::compute_implies(const Signal &x, const Signal &y) {

	clear();
	Signal not_x;
	not_x.compute_not(x);
	compute_or(not_x, y);

}

void Signal::compute_eventually(const Signal &x) {

	beginTime = x.beginTime;
	endTime = x.endTime;
	Signal::const_reverse_iterator i = x.rbegin();
	computePartialEventually(this, i, beginTime, endTime);

}

void Signal::compute_globally(const Signal  &x) {
#ifdef DEBUG__
	cout << ">  Signal::compute_globally:                   IN." << endl;
#endif

	Signal  z1, z2;
	z1.compute_not(x);
	z2.compute_eventually(z1);
	compute_not(z2);

#ifdef DEBUG__
	cout << "<  Signal::compute_globally:                   OUT." << endl;
#endif
}

void Signal::compute_until(const Signal &x, const Signal &y) {
#ifdef DEBUG__
	cout << ">>>Signal::compute_until                      IN" << endl;
#endif

	if (x.empty()||y.empty()) {
		this->clear();
		return;
	}

	double s, t;
	double z_max = BOTTOM;

	Signal::const_reverse_iterator i = x.rbegin();
	Signal::const_reverse_iterator j = y.rbegin();

	beginTime = fmax(x.beginTime, y.beginTime);
	endTime = fmin(x.endTime, y.endTime);

	while (i->time >= endTime)
		i++;
	while (j->time >= endTime)
		j++;

	s = beginTime;
	t = endTime;

	while (i->time > s) {

		computeSegmentUntil(this, *i, t, j, z_max);
		z_max = front().value;

		if (j->time == i->time)
			j++;
		t = i->time;
		i++;
	}

	if (i->time == s)
		computeSegmentUntil(this, *i, t, j, z_max);
	else
		computeSegmentUntil(this, Sample(s, i->valueAt(s), i->derivative), t, j,
				z_max);

	simplify();

#ifdef DEBUG__
	cout << "<<<Signal::compute_until                      OUT" << endl;
#endif

}

void Signal::compute_bounded_eventually(const Signal& x, double a) {
	//TODO fix compute_eventually, compute_bounded_eventually, etc...
  #ifdef DEBUG__
	cout << ">  Signal::compute_bounded_eventually         IN." << endl;
  #endif

	if (x.empty()) {
		this->clear();
		return;
	}

	Signal z1, z2, z3;

	z1.compute_plateau_max(x, a);
	z2.resize(x.beginTime + a, x.endTime + a, BOTTOM);
	z2.shift(-a);
	z3.compute_or(z2, z1);

	compute_or(x, z3);
	simplify();

  #ifdef DEBUG__
	cout << "<  Signal::compute_bounded_eventually         OUT." << endl;
  #endif


}

void Signal::compute_timed_eventually(const Signal& x, double a, double b) {
	//TODO fix compute_eventually, compute_bounded_eventually, etc...
  #ifdef DEBUG__
	cout << ">  Signal::compute_timed_eventually           IN." << endl;
  #endif

	if (x.empty()) {
		this->clear();
		return;
	}

	Signal y = x;
	Signal *rob= nullptr;

	if (a>0)
		y.shift(-a);

	if (b-a <= 0)
	{
		*this = y;
		return;
	}
	else
	{
		if (b-a >= (y.endTime- y.beginTime))
		{
			rob = computeEventually(&y);
		}
		else
		{
			rob = computeBoundedEventually(&y, b-a);
		}
	}
	*this = *rob;
	delete rob;

  #ifdef DEBUG__
	cout << ">  Signal::compute_timed_eventually           OUT." << endl;
  #endif
}

void Signal::compute_bounded_globally(const Signal& x, double a) {
#ifdef DEBUG__
	cout << "Entering compute_bounded_eventually" << endl;
#endif

	if (x.empty()) {
		this->clear();
		return;
	}

	Signal z1, z2, z3;

	z1.compute_plateau_min(x, a);
	z2.resize(x.beginTime + a, x.endTime + a, TOP);
	z2.shift(-a);
	z3.compute_and(z2, z1);

	compute_and(x, z3);
	simplify();

}

void Signal::compute_timed_globally(const Signal& x, double a, double b) {
#ifdef DEBUG__
	cout << "Signal::compute_timed_globally:               IN." << endl;
	cout << "IN: " << x << endl;
#endif

	if (x.empty()) {
		this->clear();
		return;
	}

	Signal y=x;

	Signal *rob = nullptr;

	if (a>0)
		y.shift(-a);

	if (b-a <= 0)
	{
		*this = y;
		return;
	}
	else
	{
		if (b-a >= (y.endTime- y.beginTime))
		{
			rob = computeGlobally(&y);
		}
		else
		{
			rob=computeBoundedGlobally(&y, b-a);
		}
	}

	*this = *rob;
	delete rob;

#ifdef DEBUG__
	cout << "OUT:" << *this << endl; 
	cout << "compute_timed_globally:                       OUT." << endl;
#endif
}


void Signal::compute_timed_until(const Signal& x , const Signal& y, double a, double b) {

	Signal y1 = x;
	Signal y2 = y;
	Signal *yunt;

	if (x.empty()||y.empty()) {
		this->clear();
		return;
	}


	if(b > std::min(y1.endTime, y2.endTime))
	{
		if (a==0)
			yunt = computeUntil(&y1,&y2);
		else
			if (a > y1.endTime) {
				yunt = computeGlobally(&y1);
			}
			else
			{
				Signal * yalw1 = computeBoundedGlobally(&y1, a);
				Signal * yuntmp = computeUntil(&y1,&y2);
				yuntmp->shift(-a);
				yunt = computeAnd(yalw1, yuntmp);
				delete yalw1;
				delete yuntmp;
			}
	}
	else
		yunt = computeTimedUntil(&y1,&y2, a, b );

	*this = *yunt;
	delete yunt;

}

// TODO fix and use this version of plateau max...
void Signal::compute_plateau_max(const Signal &x, double a) {
   #ifdef DEBUG__
	cout << "Entering compute_plateau_max" << endl;
   #endif
	bool new_candidate = true, end_candidate = true;
	double t, v;
	Signal::const_reverse_iterator j;
	Sequence M; //sorted in ascending times from front to back
	Sequence y; //maximum of x(t) and x(t-) at discontinuity points of x
	Sequence::iterator i;

	clear();

	beginTime = x.beginTime;
	endTime = x.endTime;

	//PRECOMPUTATION: list the local maximums of x
	t = x.endTime;
	v = BOTTOM;
	for (j = x.rbegin(); j != x.rend(); j++) {
		y.push_front(Point(t, fmax(v, j->valueAt(t))));
		//    cout << "y.time:" << y.front().time << " y->value:" << y.front().value << endl;
		t = j->time;
		v = j->value;
	}
	y.push_front(Point(x.front().time, x.front().value));

	//INIT: read values in [0, a)
	i = y.begin();

	while (i->time < x.beginTime + a) {
		while (!M.empty() && i->value >= M.back().value) {
			M.pop_back();
		}
		//    cout << "i->time:" << i->time << " i->value:" << i->value << endl;

		M.push_back(*i);
		i++;
	}

	if (i->time == x.beginTime + a)
		new_candidate = true;
	else
		new_candidate = false;

	end_candidate = false;
	t = x.beginTime;

	// push first element in M
	push_back(Sample(t, M.front().value, 0));

	//STEP
	bool cont = true;
	while (cont) {
		//UPDATE OF CANDIDATE LIST
		//candidate crosses t: remove it from list
		if (end_candidate) {
			//      cout << "end_candidate: M.pop front" << endl;
			M.pop_front();
		}

		//sample crosses t+a: add it to candidate list
		if (new_candidate) {
			//cout << "new_candidate: doing stuff" << endl;
			while (!M.empty() && i->value >= M.back().value) {
				M.pop_back();
			}
			//detect if new maximum is found
			if (!M.empty()) {
				//if M non empty then t + a does not generate new maximum
				new_candidate = false;
			}
			//add candidate
			M.push_back(*i);

			//increment iterator
			i++;
		}

		//OUTPUT OF NEW MAXIMUM

		//next best candidate
		if (!M.empty()&&(t> back().time)) {
			push_back(Sample(t, M.front().value, 0));
		}


		//NEXT EVENT DETECTION
		if (!M.empty()) {
			if (i != y.end()) {
				if (i->time - a == M.front().time) {
					t = M.front().time;
					new_candidate = true;
					end_candidate = true;
				} else if (i->time - a < M.front().time) {
					t = i->time - a;
					new_candidate = true;
					end_candidate = false;
				} else { //M.back().time < i->time - a
					t = M.front().time;
					new_candidate = false;
					end_candidate = true;

				}
			} else {
				t = M.front().time;
				new_candidate = false;
				end_candidate = true;
			}
		} else {
			if (i != y.end()) {
				t = i->time - a;
				new_candidate = true;
				end_candidate = false;
			} else {
				new_candidate = false;
				end_candidate = false;
			}
		}
		cont = (end_candidate || new_candidate);
	}

}


//copy of plateauMax, operator < switched with operator >, TOP replaces BOTTOM
void Signal::compute_plateau_min(const Signal &x, double a) {
	bool new_candidate, end_candidate;
	double t, v;
	Signal::const_reverse_iterator j;
	Sequence M; //sorted in ascending times from front to back
	Sequence y; //minimum of x(t) and x(t-) at discontinuity points of x
	Sequence::iterator i;

	clear();
	beginTime = x.beginTime;
	endTime = x.endTime;

	//PRECOMPUTATION: list the local minimums of x
	t = x.endTime;
	v = TOP;
	for (j = x.rbegin(); j != x.rend(); j++) {
		y.push_front(Point(t, fmin(v, j->valueAt(t))));
		t = j->time;
		v = j->value;
	}
	y.push_front(Point(x.front().time, x.front().value));

	//INIT: read values in [0, a)
	i = y.begin();
	while (i->time < x.beginTime + a) {
		while (!M.empty() && i->value <= M.back().value) {
			M.pop_back();
		}
		M.push_back(*i);
		i++;
	}
	if (i->time == x.beginTime + a)
		new_candidate = true;
	else
		new_candidate = false;

	end_candidate = false;
	t = x.beginTime;

	// push first element in M
	push_back(Sample(t, M.front().value, 0));

	//STEP
	bool cont = true;
	while (cont) {
		//UPDATE OF CANDIDATE LIST
		//candidate crosses t: remove it from list
		if (end_candidate) {
			M.pop_front();
		}

		//sample crosses t+a: add it to candidate list
		if (new_candidate) {
			while (!M.empty() && i->value <= M.back().value) {
				M.pop_back();
			}
			//detect if new minimum is found
			if (!M.empty()) {
				//if M non empty then t + a does not generate new minimum
				new_candidate = false;
			}
			//add candidate
			M.push_back(*i);

			//increment iterator
			i++;
		}

		//OUTPUT OF NEW MINIMUM
		//no candidate
		//			if (M.empty()) {
		//				push_back(Sample(t, TOP, 0));
		//			}
		//next best candidate
		if (!M.empty()&&(t> back().time)) {
			push_back(Sample(t, M.front().value, 0));
		}

		//NEXT EVENT DETECTION
		if (!M.empty()) {
			if (i != y.end()) {
				if (i->time - a == M.front().time) {
					t = M.front().time;
					new_candidate = true;
					end_candidate = true;
				} else if (i->time - a < M.front().time) {
					t = i->time - a;
					new_candidate = true;
					end_candidate = false;
				} else { //M.back().time < i->time - a
					t = M.front().time;
					new_candidate = false;
					end_candidate = true;

				}
			} else {
				t = M.front().time;
				new_candidate = false;
				end_candidate = true;
			}
		} else {
			if (i != y.end()) {
				t = i->time - a;
				new_candidate = true;
				end_candidate = false;
			} else {
				new_candidate = false;
				end_candidate = false;
			}
		}
		cont = (new_candidate || end_candidate);
	}

}


/*---------------------------------------------------------------------------*
 *          MAIN FUNCTIONS (TO BE DEPRECATED SOME DAY)                       *
 *---------------------------------------------------------------------------*/

Signal * computeNot(Signal * y) {
    Signal * z = new Signal();
    Signal::iterator i;
    
	z->beginTime = y->beginTime;
	z->endTime = y->endTime;
    
	for (i = y->begin(); i != y->end(); i++) {
		z->push_back(Sample(i->time, -(i->value), -(i->derivative)));
	}
	return z;
}

Signal * computeAnd(Signal * x, Signal * y) {

	Signal::const_reverse_iterator i = x->rbegin();
	Signal::const_reverse_iterator j = y->rbegin();

	Signal * z = new Signal();

	z->beginTime = fmax(x->beginTime, y->beginTime);
	z->endTime = fmin(x->endTime, y->endTime);

	while (i->time >= z->endTime)
		i++;
	while (j->time >= z->endTime)
		j++;

	computePartialAnd(z, i, j, z->beginTime, z->endTime);
	z->simplify();

	return z;
}


Signal * computeOr(Signal * x, Signal * y) {

#ifdef DEBUG__
	cout << ">  computeOr:                                 IN." << endl;
#endif

	Signal::const_reverse_iterator i = x->rbegin();
	Signal::const_reverse_iterator j = y->rbegin();

	Signal * z = new Signal();

	z->beginTime = fmax(x->beginTime, y->beginTime);
	z->endTime = fmin(x->endTime, y->endTime);

	while (i->time >= z->endTime)
		i++;
	while (j->time >= z->endTime)
		j++;

	computePartialOr(z, i, j, z->beginTime, z->endTime);
	z->simplify();

#ifdef DEBUG__
	cout << "<  computeOr:                                 OUT." << endl;
#endif
	return z;

}

Signal * computeImplies(Signal * x, Signal * y) {

	Signal *not_x = computeNot(x);
	Signal *z = computeOr(not_x, y);

	delete not_x;
	return z;
}

Signal * computeEventually(Signal * x) {

#ifdef DEBUG__
	cout << ">  computeEventually:                         IN." << endl;
#endif

	Signal * z = new Signal();
	z->beginTime = x->beginTime;
	z->endTime = x->endTime;

	Signal::const_reverse_iterator i = x->rbegin();
	computePartialEventually(z, i, z->beginTime, z->endTime);

#ifdef DEBUG__
	cout << "IN: x:" << *x << endl; 
	cout << "OUT: z: " << *z << endl;
	cout << "<  computeEventually:                         OUT." << endl;
#endif

return z;
}


Signal * computeGlobally(Signal * x) {
#ifdef DEBUG__
	cout << ">  computeGlobally:                            IN." << endl;
#endif

	Signal * z1 = computeNot(x);
	Signal * z2 = computeEventually(z1);
	Signal *z = computeNot(z2);

	delete z1;
	delete z2;

#ifdef DEBUG__
	cout << "<  computeGlobally:                            OUT." << endl;
#endif
	return z;
}


Signal * computeUntil(Signal * x, Signal * y) {
#ifdef DEBUG__
	cout << ">  computeUntil:                              IN." << endl;
#endif

	double s, t;
	double z_max = BOTTOM;
	Signal * z = new Signal();

	Signal::const_reverse_iterator i = x->rbegin();
	Signal::const_reverse_iterator j = y->rbegin();

	z->beginTime = fmax(x->beginTime, y->beginTime);
	z->endTime = fmin(x->endTime, y->endTime);

	while (i->time >= z->endTime)
		i++;
	while (j->time >= z->endTime)
		j++;

	s = z->beginTime;
	t = z->endTime;

	while (i->time > s) {

		computeSegmentUntil(z, *i, t, j, z_max);
		z_max = z->front().value;

		if (j->time == i->time)
			j++;
		t = i->time;
		i++;
	}

	if (i->time == s)
		computeSegmentUntil(z, *i, t, j, z_max);
	else
		computeSegmentUntil(z, Sample(s, i->valueAt(s), i->derivative), t, j,
				z_max);

	z->simplify();

#ifdef DEBUG__
	cout << "<  computeUntil:                              OUT." << endl;
#endif
	return z;
}

Signal * computeBoundedEventually(Signal * x, double a) {
#ifdef DEBUG__
	cout << ">  computeBoundedEventually:                  IN." << endl;
	cout << "IN: " << *x << endl;
#endif

	Signal *z, *z1, *z2, *z3;

	z1 = plateauMax(x, a);
	z2 = new Signal(*x);
	z2->resize(x->beginTime + a, x->endTime + a, BOTTOM);
	z2->shift(-a);

	z3 = computeOr(z2, z1);

	delete z1;
	delete z2;

	z = computeOr(x, z3);

	delete z3;
	z->simplify();

#ifdef DEBUG__
	cout << "OUT: " << *z << endl;
	cout << "<  computeBoundedEventually:                  OUT." << endl;
#endif
	return z;

}

Signal * computeBoundedGlobally(Signal * x, double a) {

	Signal *z, *z1, *z2, *z3;

	//cout << "computeBoundedGlobally with a=" << a << endl;
	//cout << "x="  << *x << endl;
	z1 = plateauMin(x, a);

	z2 = new Signal(*x);
	//cout << *z2 << endl;
	z2->resize(x->beginTime + a, x->endTime + a, TOP);
	z2->shift(-a);

	//cout << *z1 << endl;
	//cout << *z2 << endl;

	z3 = computeAnd(z2, z1);

	delete z1;
	delete z2;


	z = computeAnd(x, z3);
	delete z3;
	z->simplify();

	return z;
}


Signal * computeTimedUntil(Signal * x, Signal * y, double a, double b) {
	Signal *z, *z2, *z3, *z4;

	z = new Signal();

	z2 = computeBoundedEventually(y, b - a);

	z3 = computeUntil(x, y);
	z4 = computeAnd(z2, z3);

	delete z2;
	delete z3;

	if (a > 0) {
		Signal * z1 = computeBoundedGlobally(x, a);
		z4->shift(-a);
		z = computeAnd(z1, z4);
		delete z1;
	} else
		z = computeAnd(x, z4);

	delete z4;
	return z;
}

/*---------------------------------------------------------------------------*
 *                         CONJUNCTION SUBROUTINES                           *
 *---------------------------------------------------------------------------*/

//PRECONDITIONS: j->time < t, i.time < t.
//POSTCONDITIONS: j->time <= i.time.
void computeSegmentAnd(Signal * z, const Sample & i, double t,
		Signal::const_reverse_iterator & j) {
	bool continued = false;
	double s = j->time;

	// for every sample *j in (i.time, t)
	while (s > i.time) {
		if (i.valueAt(t) < j->valueAt(t)) {
			if (i.valueAt(s) > j->value) {
				t = i.timeIntersect(*j);
				z->push_front(Sample(t, i.valueAt(t), i.derivative));
				z->push_front(Sample(s, j->value, j->derivative));
				continued = false;
			} else
				continued = true;
		} else if (i.valueAt(t) == j->valueAt(t)) {
			if (i.valueAt(s) > j->value) {
				if (continued) {
					z->push_front(Sample(t, i.valueAt(t), i.derivative));
					continued = false;
				}
				z->push_front(Sample(s, j->value, j->derivative));
			} else
				continued = true;
		} else {
			if (i.valueAt(s) < j->value) {
				if (continued) {
					z->push_front(Sample(t, i.valueAt(t), i.derivative));
				}
				t = i.timeIntersect(*j);
				z->push_front(Sample(t, j->valueAt(t), j->derivative));
				continued = true;
			} else {
				if (continued) {
					z->push_front(Sample(t, i.valueAt(t), i.derivative));
					continued = false;
				}
				z->push_front(Sample(s, j->value, j->derivative));
			}
		}

		//increment reverse iterator j
		t = s;
		j++;
		s = j->time;
	}

	//here we may have j->time < i.time
	// "i" values of z are no longer "continued"
	s = i.time;
	if (i.valueAt(t) < j->valueAt(t)) {
		if (i.value > j->valueAt(s)) {
			t = i.timeIntersect(*j);
			z->push_front(Sample(t, i.valueAt(t), i.derivative));
			z->push_front(Sample(s, j->valueAt(s), j->derivative));
		} else {
			z->push_front(i);
		}
	} else if (i.valueAt(t) == j->valueAt(t)) {
		if (i.value > j->valueAt(s)) {
			if (continued) {
				z->push_front(Sample(t, i.valueAt(t), i.derivative));
			}
			z->push_front(Sample(s, j->valueAt(s), j->derivative));
		} else {
			z->push_front(i);
		}
	} else {
		if (i.value < j->valueAt(s)) {
			if (continued) {
				z->push_front(Sample(t, i.valueAt(t), i.derivative));
			}
			t = i.timeIntersect(*j);
			z->push_front(Sample(t, j->valueAt(t), j->derivative));
			z->push_front(i);
		} else {
			if (continued) {
				z->push_front(Sample(t, i.valueAt(t), i.derivative));
			}
			z->push_front(Sample(s, j->valueAt(s), j->derivative));
		}
	}

}


//TODO  this routine is unconvential but safe memory-wise (no new/delete) so fix not urgent
void computePartialAnd(Signal * z, Signal::const_reverse_iterator & i,
		Signal::const_reverse_iterator & j, double s, double t) {

	while (i->time > s) {
		computeSegmentAnd(z, *i, t, j);
		if (j->time == i->time)
			j++;
		t = i->time;
		i++;
	}

	if (i->time == s)
		computeSegmentAnd(z, *i, t, j);
	else
		computeSegmentAnd(z, Sample(s, i->valueAt(s), i->derivative), t, j);

}

/*---------------------------------------------------------------------------*
 *                         DISJUNCTION SUBROUTINES                           *
 *---------------------------------------------------------------------------*/

//TODO  this routine is unconvential but safe memory-wise (no new/delete) so fix not urgent
// copy of computeSegmentAnd, operator "<" switched with operator ">"
void computeSegmentOr(Signal * z, const Sample & i, double t,
		Signal::const_reverse_iterator & j) {
	bool continued = false;
	double s = j->time;

	// for every sample *j in (i.time, t)
	while (s > i.time) {
		if (i.valueAt(t) > j->valueAt(t)) {
			if (i.valueAt(s) < j->value) {
				t = i.timeIntersect(*j);
				z->push_front(Sample(t, i.valueAt(t), i.derivative));
				z->push_front(Sample(s, j->value, j->derivative));
				continued = false;
			} else
				continued = true;
		} else if (i.valueAt(t) == j->valueAt(t)) {
			if (i.valueAt(s) < j->value) {
				if (continued) {
					z->push_front(Sample(t, i.valueAt(t), i.derivative));
					continued = false;
				}
				z->push_front(Sample(s, j->value, j->derivative));
			} else
				continued = true;
		} else {
			if (i.valueAt(s) > j->value) {
				if (continued) {
					z->push_front(Sample(t, i.valueAt(t), i.derivative));
				}
				t = i.timeIntersect(*j);
				z->push_front(Sample(t, j->valueAt(t), j->derivative));
				continued = true;
			} else {
				if (continued) {
					z->push_front(Sample(t, i.valueAt(t), i.derivative));
					continued = false;
				}
				z->push_front(Sample(s, j->value, j->derivative));
			}
		}

		//increment reverse iterator j
		t = s;
		j++;
		s = j->time;
	}

	//here we may have j->time < i.time
	// "i" values of z are no longer "continued"
	s = i.time;
	if (i.valueAt(t) > j->valueAt(t)) {
		if (i.value < j->valueAt(s)) {
			t = i.timeIntersect(*j);
			z->push_front(Sample(t, i.valueAt(t), i.derivative));
			z->push_front(Sample(s, j->valueAt(s), j->derivative));
		} else {
			z->push_front(i);
		}
	} else if (i.valueAt(t) == j->valueAt(t)) {
		if (i.value < j->valueAt(s)) {
			if (continued) {
				z->push_front(Sample(t, i.valueAt(t), i.derivative));
			}
			z->push_front(Sample(s, j->valueAt(s), j->derivative));
		} else {
			z->push_front(i);
		}
	} else {
		if (i.value > j->valueAt(s)) {
			if (continued) {
				z->push_front(Sample(t, i.valueAt(t), i.derivative));
			}
			t = i.timeIntersect(*j);
			z->push_front(Sample(t, j->valueAt(t), j->derivative));
			z->push_front(i);
		} else {
			if (continued) {
				z->push_front(Sample(t, i.valueAt(t), i.derivative));
			}
			z->push_front(Sample(s, j->valueAt(s), j->derivative));
		}
	}

}

//TODO  this routine is unconvential but safe memory-wise (no new/delete) so fix not urgent
void computePartialOr(Signal * z, Signal::const_reverse_iterator & i,
		Signal::const_reverse_iterator & j, double s, double t) {

	while (i->time > s) {
		computeSegmentOr(z, *i, t, j);
		if (j->time == i->time)
			j++;
		t = i->time;
		i++;
	}

	if (i->time == s)
		computeSegmentOr(z, *i, t, j);
	else
		computeSegmentOr(z, Sample(s, i->valueAt(s), i->derivative), t, j);

}

/*---------------------------------------------------------------------------*
 *                              UNTIL SUBROUTINES                            *
 *---------------------------------------------------------------------------*/

//TODO  this routine is unconvential but safe memory-wise (no new/delete) so fix not urgent
void computePartialEventually(Signal* z, Signal::const_reverse_iterator & i, double s, double t) {
	bool continued = false;
	double z_max = BOTTOM;
	while (i->time > s) {
		if (i->derivative >= 0) {
			if (z_max < i->valueAt(t)) {
				if (continued) {
					z->push_front(Sample(t, z_max, 0));
				}
				z_max = i->valueAt(t);
			}
			continued = true;
			//z->push_front(Sample(i->time, z_max, 0));
		} else if (i->valueAt(t) >= z_max) {
			if (continued) {
				z->push_front(Sample(t, z_max, 0));
				continued = false;
			}
			z_max = i->value;
			z->push_front(*i);
		} else if (z_max >= i->value) {
			continued = true;
			//z->push_front(Sample(i->time, z_max, 0));
		} else {
			z->push_front(
					Sample(i->time + (z_max - i->value) / i->derivative, z_max,
							0)); //time at which y reaches value next_z
			z->push_front(*i);
			z_max = i->value;
			continued = false;
		}

		t = i->time;
		i++;
	}

	//leftmost sample *i may not be on s
	//"z_max" values of z are not longer "continued".
	if (i->derivative >= 0) {
		if (z_max < i->valueAt(t)) {
			if (continued) {
				z->push_front(Sample(t, z_max, 0));
			}
			z_max = i->valueAt(t);
		}
		z->push_front(Sample(s, z_max, 0));
	} else if (i->valueAt(t) >= z_max) {
		if (continued) {
			z->push_front(Sample(t, z_max, 0));
		}
		z->push_front(Sample(s, i->valueAt(s), i->derivative));
	} else if (z_max >= i->valueAt(s)) {
		z->push_front(Sample(s, z_max, 0));
	} else {
		z->push_front(Sample(s + (z_max - i->value) / i->derivative, z_max, 0)); //time at which y reaches value next_z
		z->push_front(Sample(s, i->valueAt(s), i->derivative));
	}

}

//TODO  this routine is unconvential but safe memory-wise (no new/delete) so fix not urgent
void computeSegmentUntil(Signal * z, const Sample & i, double t,
		Signal::const_reverse_iterator & j, double z_max) {
	Signal *z1, *z2, *z3;
	Signal::const_reverse_iterator k, l;
	double s = i.time;

	if (i.derivative <= 0) {
		z1 = new Signal();
		computeSegmentAnd(z1, i, t, j);

		z2 = new Signal();
		k = z1->rbegin();
		computePartialEventually(z2, k, s, t);
		delete z1;

		l = z2->rbegin();
		computeSegmentOr(z, Sample(s, fmin(z_max, i.valueAt(t)), 0), t, l);
		delete z2;
	} else {
		z1 = new Signal();
		computePartialEventually(z1, j, s, t);

		z2 = new Signal();
		k = z1->rbegin();
		computeSegmentAnd(z2, i, t, k);
		delete z1;

		z1 = new Signal();
		z3 = new Signal();
		z3->push_front(Sample(s, z_max, 0));
		k = z3->rbegin();
		computeSegmentAnd(z1, i, t, k);
		delete z3;

		k = z1->rbegin();
		l = z2->rbegin();
		computePartialOr(z, k, l, s, t);
		delete z1;
		delete z2;
	}
}

/*---------------------------------------------------------------------------*
 *                           TIMED UNTIL SUBROUTINES                         *
 *---------------------------------------------------------------------------*/

//computation of the max induced by discontinuity points of x in t+(0,a] 
//a is assumed to be smaller than length of x.
//based on Lemire algorithm for "streaming min-max filter"
Signal * plateauMax(Signal * x, double a) {
    #ifdef DEBUG___
	cout << ">>>plateauMax:                                IN.====================================================================================" << endl;
    cout << "IN: a=" << a << " x=" << *x << endl;
 	#endif
	bool new_candidate = true, end_candidate = true;
	double t, v;
	Signal::const_reverse_iterator j;
	Sequence M; //sorted in ascending times from front to back
	Sequence y; //maximum of x(t) and x(t-) at discontinuity points of x
	Sequence::iterator i;

	Signal * z = new Signal();
	z->beginTime = x->beginTime;
	z->endTime = x->endTime;

	//PRECOMPUTATION: list the local maximums of x
	t = x->endTime;
	v = BOTTOM;
	for (j = x->rbegin(); j != x->rend(); j++) {
		y.push_front(Point(t, fmax(v, j->valueAt(t))));
		t = j->time;
		v = j->value;
	}
	y.push_front(Point(x->front().time, x->front().value));

	#ifdef DEBUG___
	cout << "y: " << y << endl;
	#endif

	//INIT: read values in [0, a)
	i = y.begin();

	while (i->time < x->beginTime + a) {
		while (!M.empty() && i->value >= M.back().value) {
			M.pop_back();
		}
		//    cout << "i->time:" << i->time << " i->value:" << i->value << endl;

		M.push_back(*i);
		i++;
	}

	#ifdef DEBUG___
	cout << "M:" << M << endl;
	#endif
	
	if (i->time == x->beginTime + a)
		new_candidate = true;
	else
		new_candidate = false;

	end_candidate = false;
	t = x->beginTime;

	// push first element in M
	//z->push_back(Sample(t, M.front().value, 0));

	//cout << "z:" << *z << endl;
	//STEP
	bool cont = true;
	while (cont) {
		//UPDATE OF CANDIDATE LIST
		#ifdef DEBUG___
		cout << "++++++++++++++++++++++++++++++++++++++\nUPDATE OF CANDIDATE LIST" << endl;
		#endif
		//candidate crosses t: remove it from list
		if (end_candidate) {
	        #ifdef DEBUG___
			cout << "end_candidate: M.pop front: " << M.front() <<endl;
			#endif
			M.pop_front();
			
		}

		//sample crosses t+a: add it to candidate list
		if (new_candidate) {
			while (!M.empty() && i->value >= M.back().value) {
				M.pop_back();
			}
			//detect if new maximum is found
			if (!M.empty()) {
				//if M non empty then t + a does not generate new maximum
				new_candidate = false;
			}
			//add candidate
			M.push_back(*i);

			//increment iterator
			i++;
		}

		//OUTPUT OF NEW MAXIMUM
	   	#ifdef DEBUG___
		   cout << "OUTPUT OF NEW MAXIMUM" << endl;
		#endif

		if (M.empty()) {                         
			#ifdef DEBUG___
			cout << "No candidate." << endl;
			#endif
			z->push_back(Sample(t, BOTTOM, 0));
		}
           //next best candidate
        else {
			//if ((!z->empty())&& (z->back().time==t)) 
				 //t = M.front().time;
			//	 z->push_back(Sample(M.front().time, M.front().value, 0));
			//else 
				z->push_back(Sample(t, M.front().value, 0));
		}
		//cout << "z: " << *z << endl;
		//NEXT EVENT DETECTION
	   	#ifdef DEBUG___
		   cout << "NEXT EVENT DETECTION" << endl;
		#endif
		if (!M.empty()) {

			if (i != y.end()) {
				if (i->time - a == M.front().time) {
					#ifdef DEBUG___
					cout << "M not empty and M.front+a == i->time" << endl << "=> new & end" << endl; 
					#endif
					t = M.front().time;
					new_candidate = true;
					end_candidate = true;
				} else if (i->time - a < M.front().time) {
					#ifdef DEBUG___
					cout << "M not empty and M.front+a > i->time" << endl << "=> new & !end" << endl; 
			    	#endif
					t = i->time - a;
					new_candidate = true;
					end_candidate = false;
				} else { //M.back().time < i->time - a
					#ifdef DEBUG___
					cout << "M not empty and M.front+a <= i->time" << endl << "=> !new & end" << endl; 
					#endif
					t = M.front().time;
					new_candidate = false;
					end_candidate = true;

				}
			} else {
				#ifdef DEBUG___
				cout << "i==y.end" << endl << "=> !new & end" << endl; 
				#endif
				t = M.front().time;
				new_candidate = false;
				end_candidate = true;
			}
		} else {
			if (i != y.end()) {
				#ifdef DEBUG___
				cout << "M empty & i != y.end()" << endl << " => !new & end " << endl; 
				#endif
				t = i->time - a;
				new_candidate = true;
				end_candidate = false;
			} else {
				#ifdef DEBUG___
				cout << "M empty and i==y.end" << endl << "!new & !end " << endl; 
				#endif
				new_candidate = false;
				end_candidate = false;
			}
		}
		cont = (new_candidate||end_candidate);
	}

    if(z->back().time==z->endTime) z->pop_back(); // from Breach ...

#ifdef DEBUG___
		cout << "OUT:" <<  *z  << endl; 
		cout << "<<<plateauMax:                            OUT.==================================================================================" << endl;
	#endif

	return z;
}


//copy of plateauMax, operator < switched with operator >, TOP replaces BOTTOM
Signal * plateauMin(Signal * x, double a) {
	bool new_candidate, end_candidate;
	double t, v;
	Signal::const_reverse_iterator j;
	Sequence M; //sorted in ascending times from front to back
	Sequence y; //minimum of x(t) and x(t-) at discontinuity points of x
	Sequence::iterator i;

	Signal * z = new Signal();
	z->beginTime = x->beginTime;
	z->endTime = x->endTime;

	//PRECOMPUTATION: list the local minimums of x
	t = x->endTime;
	v = TOP;
	for (j = x->rbegin(); j != x->rend(); j++) {
		y.push_front(Point(t, fmin(v, j->valueAt(t))));
		t = j->time;
		v = j->value;
	}
	y.push_front(Point(x->front().time, x->front().value));

	//INIT: read values in [0, a)
	i = y.begin();
	while (i->time < x->beginTime + a) {
		while (!M.empty() && i->value <= M.back().value) {
			M.pop_back();
		}
		M.push_back(*i);
		i++;
	}
	if (i->time == x->beginTime + a)
		new_candidate = true;
	else
		new_candidate = false;

	end_candidate = false;
	t = x->beginTime;

	// push first element in M
	//z->push_back(Sample(t, M.front().value, 0));

	//STEP
	bool cont = true;
	while (cont) {
		//UPDATE OF CANDIDATE LIST
		//candidate crosses t: remove it from list
		if (end_candidate) {
			M.pop_front();
		}

		//sample crosses t+a: add it to candidate list
		if (new_candidate) {
			while (!M.empty() && i->value <= M.back().value) {
				M.pop_back();
			}
			//detect if new minimum is found
			if (!M.empty()) {
				//if M non empty then t + a does not generate new minimum
				new_candidate = false;
			}
			//add candidate
			M.push_back(*i);

			//increment iterator
			i++;
		}

		//OUTPUT OF NEW MINIMUM
		//no candidate
        if (M.empty()) {
             z->push_back(Sample(t, TOP, 0)); //  => Why ? otherwise crashes :)
         }

        else //next best candidate
            {
		 	// To Be double checked
            //if ((!z->empty())&& (z->back().time==t)) 
			//	 z->push_back(Sample(M.front().time, M.front().value, 0));
			//else 
				z->push_back(Sample(t, M.front().value, 0));

			
			}

		//NEXT EVENT DETECTION
		if (!M.empty()) {
			if (i != y.end()) {
				if (i->time - a == M.front().time) {
					t = M.front().time;
					new_candidate = true;
					end_candidate = true;
				} else if (i->time - a < M.front().time) {
					t = i->time - a;
					new_candidate = true;
					end_candidate = false;
				} else { //M.back().time < i->time - a
					t = M.front().time;
					new_candidate = false;
					end_candidate = true;

				}
			} else {
				t = M.front().time;
				new_candidate = false;
				end_candidate = true;
			}
		} else {
			if (i != y.end()) {
				t = i->time - a;
				new_candidate = true;
				end_candidate = false;
			} else {
				new_candidate = false;
				end_candidate = false;
			}
		}
		cont = (new_candidate || end_candidate);
	}

    if(z->back().time==z->endTime) z->pop_back();

	return z;
}


}
