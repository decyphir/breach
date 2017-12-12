#include "signal.h"
#include "mex.h"

/* 
 * class Signal member functions
 */
Signal::Signal(double T, double V, int n) {
	double p_t, p_v; //previous times and value	
	
	beginTime=T;
	endTime=T;
	push_back(Sample(T, V, 0.));
	
}

Signal::Signal(double * T, double * V, int n) {
#ifdef DEBUG__
  cout << "Entering Signal::Signal" << endl;
#endif

  double p_t, p_v; //previous times and value	
	
  beginTime=T[0];
  endTime=T[n-1];
	
  if (n==1)
    push_back(Sample(T[0], V[0], 0.));
  else {
    for(int i=0; i < n-1; i++) {

#ifdef NO_LINEAR_INTERPOL
        push_back(Sample(T[i], V[i], 0.));
#else
        push_back(Sample(T[i], V[i], (V[i+1]-V[i]) / (T[i+1] - T[i])));
#endif
    }
    push_back(Sample(T[n-1], V[n-1], 0.));
  }
}

Signal::Signal(Sequence S) {
  double p_t, p_v; //previous times and value
  double d; // derivative
  
  if(S.size() == 0) {
    beginTime=0;
    endTime=0;
    Signal();
  }else{
    beginTime=S.front().time;
    p_t = S.front().time;
    p_v = S.front().value;
    S.pop_front();
    while(S.size() > 0){
      d = (S.front().value - p_v) / (S.front().time - p_t);
      push_back(Sample(p_t, p_v, d));
      p_t = S.front().time;
      p_v = S.front().value;
      S.pop_front();
    }
    push_back(Sample(p_t, p_v, 0.));
    endTime=p_t;
  }
}

int Signal::push_front(Sample P) {
// Specialize push_front - ensures time is properly increasing and not beyond endTime
    if (empty()) {
        if (P.time< endTime)
            std::deque<Sample>::push_front(P);
        else {
#ifdef DEBUG__
            mexPrintf("WARNING: time after endTime!\n");
#endif
            std::deque<Sample>::push_front(P);
            endTime=P.time;
        }
    }
    else
        if   (P.time<front().time) {
            std::deque<Sample>::push_front(P);
        }
        else {
#ifdef DEBUG__
            mexPrintf("WARNING: time not advancing!\n");
#endif
            return 1;
        }
    return 0;
}

//remove linear interpolations
void Signal::simplify() {
#ifdef DEBUG__
  mexPrintf("Entering Signal::simplify\n");
#endif
	push_back(front());
	pop_front();
	while (front().time != beginTime) {
		if( back().valueAt(front().time) != front().value || back().derivative != front().derivative ) {
			push_back(front());
		}
		pop_front();
	}
#ifdef DEBUG__
  mexPrintf("Leaving Signal::simplify\n");
#endif
}

void Signal::print() const {
    Signal::const_iterator i;
    
    if(begin() == end())
        mexPrintf("EMPTY\n");
    else {
        mexPrintf("beginTime: %g endTime:%g\n", beginTime , endTime);
        
        mexPrintf("time    value   derivative\n");
        for(i = begin(); i != end(); i++) {
            mexPrintf("%g %g %g\n",i->time, i->value, i->derivative);
        }
    }
}

void Signal::resize(double s, double t, double v) {

#ifdef DEBUG__
  cout << "Entering Signal::resize" << endl;
#endif
 
  Sample first;
  
  //trim or extend front of signal
  
  if(beginTime > s) {
    push_front(Sample(s, v, 0));
  }
  else {
    while(front().time < s) {
      first=front();
      pop_front();
    }
    if(front().time > s) {
      push_front(Sample(s, first.valueAt(s), 0));
    }
  }

  //trim or extend back of signal
  
  if(endTime < t) {
    push_back(Sample(endTime, v, 0));
  }
  else {
    while(back().time >=t) {
      pop_back();
    }
  }

  beginTime=s;
  endTime=t;
}

void Signal::shift(double a) {
  Signal::iterator i;
  
  beginTime=beginTime + a;
  endTime=endTime + a;
  
  for(i = begin(); i != end(); i++) {
    i->time=i->time + a; 
  }
}


Signal* Signal::reverse(){
  Signal::iterator i;
  Signal *z = new Signal();
  double prev_d = 0;
  Sample s;
  
  i = begin();
  z->beginTime = -endTime;
  z->endTime = -beginTime;

  while(i != end()){
    s = *i;
    z->push_front(Sample(-s.time, s.value, -prev_d));
    prev_d = s.derivative;
    i++;
  }
  
  return z;
}


/*
 * friend functions
 */
std::ostream & operator<<(std::ostream & out, const Point & point) {
	out << point.time << ";" << point.value ;
	return out;
}


std::ostream & operator<<(std::ostream & out, const Sample & sample) {
	out << sample.time << ";" << sample.value << ";" << sample.derivative;
	return out;
}


std::ostream & operator<<(std::ostream & out, const Sequence & M) {
	Sequence::const_iterator i;
	
	if(M.empty()) return out << "EMPTY" << std::endl;

	for(i = M.begin(); i != M.end(); i++) {
	        out << *i << std::endl;
	}
	return out;
}


std::ostream & operator<<(std::ostream & out, const Signal & y) {
	Signal::const_iterator i;
	
	if(y.begin() == y.end()) return out << "EMPTY" << std::endl;

	out << "def = [" << y.beginTime << ", " << y.endTime << ")" << std::endl;
	for(i = y.begin(); i != y.end(); i++) {
	        out << *i << std::endl;
	}
	return out;
}

