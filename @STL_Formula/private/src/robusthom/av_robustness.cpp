#include "robustness.h"
#include "window.h"
#include <list>

using namespace std;

/*---------------------------------------------------------------------------*
 *                           MAIN FUNCTIONS                                  *
 *---------------------------------------------------------------------------*/

Signal * constZero(double beginT, double endT) {
	double V[2] = {0, 0};
	double T[2] = {beginT, endT};
	return new Signal(T, V, 2);
}


Signal * computeLeftAvEventually(Signal * x, double a) {
#ifdef DEBUG__
  cout << "Entering computeBoundedEventually" << endl;
#endif

  double beginTime, endTime;
  Signal *zZero, *xPos,*z;
  Sample endpoint;

  beginTime = x->beginTime;
  endTime = x->endTime;

  zZero = constZero(beginTime, endTime);
  xPos = computeOr(x, zZero);

  //Add endpoint
  if(xPos->back().time < xPos->endTime){
    endpoint = Sample(xPos->endTime, xPos->back().valueAt(xPos->endTime), 0);
    xPos->push_back(endpoint);
  }

  z = plateauMaxAux(xPos,a);
  delete xPos;
  delete zZero;

  z->resize(beginTime, endTime, BOTTOM);
  z->simplify();

  return z;
}

Signal * computeRightAvEventually(Signal * x, double a) {
#ifdef DEBUG__
  cout << "Entering computeBoundedEventually" << endl;
#endif

  double beginTime, endTime;
  Signal *zZero, *xPos, *z;
  Sample endpoint;

  beginTime = x->beginTime;
  endTime = x->endTime;

  zZero = constZero(beginTime, endTime);
  xPos = computeOr(x, zZero);

  //Add endpoint
  if(!xPos->empty()){
    xPos->resize(beginTime-a, endTime, (*(x->begin())).value);
    xPos->resize(beginTime-a, endTime+a, (*(x->rbegin())).value);
  }
  if(xPos->back().time < xPos->endTime){
    endpoint = Sample(xPos->endTime, xPos->back().valueAt(xPos->endTime), 0);
    xPos->push_back(endpoint);
  }

  z = plateauMaxAux(xPos->reverse(), a)->reverse();
  delete xPos;
  delete zZero;

  z->shift(-a);
  z->resize(beginTime, endTime, BOTTOM);
  z->simplify();

  return z;

}


Signal * computeLeftAvGlobally(Signal * x, double a) {
#ifdef DEBUG__
  cout << "Entering computeBoundedEventually" << endl;
#endif

  double beginTime, endTime;
  Signal *zZero, *xPos,*z;
  Sample endpoint;

  beginTime = x->beginTime;
  endTime = x->endTime;

  zZero = constZero(beginTime, endTime);
  xPos = computeOr(x, zZero);

  //Add endpoint
  if(xPos->back().time < xPos->endTime){
    endpoint = Sample(xPos->endTime, xPos->back().valueAt(xPos->endTime), 0);
    xPos->push_back(endpoint);
  }

  z = plateauMinAux(xPos,a);
  delete xPos;
  delete zZero;

  z->resize(beginTime, endTime, BOTTOM);
  z->simplify();

  return z;

}

Signal * computeRightAvGlobally(Signal * x, double a) {
#ifdef DEBUG__
  cout << "Entering computeBoundedEventually" << endl;
#endif

  double beginTime, endTime;
  Signal *zZero, *xPos, *xPosRev, *z, *zRev;
  Sample endpoint;

  beginTime = x->beginTime;
  endTime = x->endTime;

  zZero = constZero(beginTime, endTime);
  xPos = computeOr(x, zZero);

  //Add endpoint
  if(!xPos->empty()){
    xPos->resize(beginTime-a, endTime, (*(x->begin())).value);
    xPos->resize(beginTime-a, endTime+a, (*(x->rbegin())).value);
  }
  if(xPos->back().time < xPos->endTime){
    endpoint = Sample(xPos->endTime, xPos->back().valueAt(xPos->endTime), 0);
    xPos->push_back(endpoint);
  }

  z = plateauMinAux(xPos->reverse(), a)->reverse();
  delete xPos;
  delete zZero;

  z->shift(-a);
  z->resize(beginTime, endTime, BOTTOM);
  z->simplify();

  return z;

}


// Computation for a signal s.t. t mapsto average of max x(t':t+a) over t' in [t,t+a].
// a is assumed to be smaller than length of x.
// This works like Lemire algorithm for "streaming min-max filter"
// and the window slides from left to right.
Signal * plateauMaxAux(Signal * x, double a) {
#ifdef DEBUG__
  cout << "Entering plateauMaxAux" << endl;
#endif
  bool cont = true;
  bool new_candidate, end_candidate;
  Sample newSample;

  double t; // current left endpoint of the window W
  Window W; // Monotonically decreasing window whose value at t' is max x in [t,t']
  Sequence zSeq;
  Signal * z;

  Signal::iterator i;

  new_candidate = true;
  end_candidate = false;

  //INIT: read values in [0, a)
  i=x->begin();

  while(i->time <= x->beginTime + a && i != x->end()) {
    newSample = *i;
    W.updateBackMax(newSample, new_candidate);
    i++;
  }

  //NEXT EVENT DETECTION
  new_candidate = (i->time - a == 0);
  end_candidate = false;
  t = x->beginTime;

  while(cont){
    // END CANDIDATE
    W.updateFront(t, end_candidate);
    // NEW CANDIDATE
    if(new_candidate){
      newSample = *i;
      i++;
    }else if(!W.empty()){
      if(i == x->end()){
        newSample = Sample(t+a, W.back().value, 0.);
      }else{
        newSample = Sample(t+a, W.back().valueAt(t+a), W.back().derivative);
      }
    }
    W.updateBackMax(newSample, new_candidate);

    // Update Sequence
    zSeq.push_back(Point(t, W.area / a));   

    //NEXT EVENT DETECTION
    // // Here the left most element of W is always in t+a.
    // // Hence we need to look up the second element of W 
    // // to estimate when the next end_candidate event happens.
    if(i == x->end()){
      // The case of no new_candidate exists
      new_candidate = false;
      end_candidate = (W.size() > 1);
      if (end_candidate) { t = W.at(1).time; }
    }else if(W.size() < 2){
      // The case of no end_candidate exists 
      // (i.e. W contains only its left endpoint)
      new_candidate = true;
      end_candidate = false;
      t = i->time - a;
    }else{
      new_candidate = (i->time - a) <= (W.at(1).time);
      end_candidate = (W.at(1).time) <= (i->time - a);
      t = fmin(i->time - a, W.at(1).time);
    }
    cont = (new_candidate || end_candidate);
  }
  
  if (t < x->endTime){
    zSeq.push_back(Point(x->endTime, W.area / a));
  }

  return new Signal(zSeq);
}


// Computation for a signal s.t. t mapsto average of min x(t':t+a) over t' in [t,t+a].
// a is assumed to be smaller than length of x.
// This works like Lemire algorithm for "streaming min-max filter"
// and the window slides from left to right.
Signal * plateauMinAux(Signal * x, double a) {
#ifdef DEBUG__
  cout << "Entering plateauMaxAux" << endl;
#endif
  bool cont = true;
  bool new_candidate, end_candidate;
  Sample newSample;
  double t; // current left endpoint of the window time
  Window W; // Monotonically increasing window whose value at t' is max x in [t,t']
  Sequence zSeq;
  Signal::iterator i;

  new_candidate = true;
  end_candidate = false;

  //INIT: read values in [0, a)
  i=x->begin();

  while(i->time <= x->beginTime + a && i != x->end()) {
    newSample = *i;
    W.updateBackMin(newSample, new_candidate);
    i++;
  }

  //NEXT EVENT DETECTION
  new_candidate = (i->time - a == 0);
  end_candidate = false;
  t = x->beginTime;

  while(cont){
    //END CANDIDATE
    W.updateFront(t, end_candidate);
    //NEW CANDIDATE
    if(new_candidate){
      newSample = *i;
      i++;
    }else if(!W.empty()){
      if(i == x->end()){
        newSample = Sample(t+a, W.back().value, 0.);
      }else{
        newSample = Sample(t+a, W.back().valueAt(t+a), W.back().derivative);
      }
    }
    W.updateBackMin(newSample, new_candidate);

    // Update Sequence
    zSeq.push_back(Point(t, W.area / a));   

    // NEXT EVENT DETECTION
    // // Here the left most element of W is always in t+a.
    // // Hence we need to look up the second element of W 
    // // to estimate when the next end_candidate event happens.
    if(i == x->end()){
      // The case of no new_candidate exists
      new_candidate = false;
      end_candidate = (W.size() > 1);
      if (end_candidate) { t = W.at(1).time; }
    }else if(W.size() < 2){
      // The case of no end_candidate exists 
      // (i.e. W contains only its left endpoint)
      new_candidate = true;
      end_candidate = false;
      t = i->time - a;
    }else{
      new_candidate = (i->time - a) <= (W.at(1).time);
      end_candidate = (W.at(1).time) <= (i->time - a);
      t = fmin(i->time - a, W.at(1).time);
    }
    cont = (new_candidate || end_candidate);
  }
  
  if (t < x->endTime){
    zSeq.push_back(Point(x->endTime, W.area / a));
  }
  return new Signal(zSeq);
}


