#include "robustness.h"
#include <list>
#include "mex.h"

using namespace std;

/*---------------------------------------------------------------------------*
 *                           MAIN FUNCTIONS                                  *
 *---------------------------------------------------------------------------*/

Signal * computeNot(Signal * y) {
    Signal * z = new Signal();
    Signal::iterator i;
    
    z->beginTime = y->beginTime;
    z->endTime = y->endTime;
    
    for(i = y->begin(); i != y->end(); i++) {
        z->push_back(Sample(i->time, -(i->value), -(i->derivative)));
    }
    return z;
}

Signal * computeAnd(Signal * x, Signal * y) {
    
#ifdef DEBUG__
    mexPrintf("Entering computeAnd\n");
#endif
    
    
    Signal::reverse_iterator i=x->rbegin();
    Signal::reverse_iterator j=y->rbegin();
    
    Signal * z=new Signal();
    
    z->beginTime=fmax(x->beginTime, y->beginTime);
    z->endTime=fmin(x->endTime, y->endTime);
    
    while(i->time >= z->endTime)
        i++;
    while(j->time >= z->endTime)
        j++;
    
    computePartialAnd(z, i, j, z->beginTime, z->endTime);
   
    z->simplify();
    return z;
}


Signal * computeOr(Signal * x, Signal * y) {
    Signal::reverse_iterator i=x->rbegin();
    Signal::reverse_iterator j=y->rbegin();
    
    Signal * z=new Signal();
    
    z->beginTime=fmax(x->beginTime, y->beginTime);
    z->endTime=fmin(x->endTime, y->endTime);
    
    while(i->time >= z->endTime)
        i++;
    while(j->time >= z->endTime)
        j++;
    
    computePartialOr(z, i, j, z->beginTime, z->endTime);
    z->simplify();
    
    return z;
}


Signal * computeEventually(Signal * x) {
#ifdef DEBUG__
    cout << "Entering computeEventually" << endl;
#endif
    
    Signal * z=new Signal();
    
    z->beginTime=x->beginTime;
    z->endTime=x->endTime;
    
    Signal::reverse_iterator i=x->rbegin();
    
    computePartialEventually(z,i,z->beginTime, z->endTime);
    
    return z;
}


Signal * computeUntil(Signal * x, Signal * y) {
#ifdef DEBUG__
    mexPrintf("Entering computeUntil\n");
#endif
    
    double s, t;
    double z_max=BOTTOM;
    Signal * z=new Signal();
    
    Signal::reverse_iterator i=x->rbegin();
    Signal::reverse_iterator j=y->rbegin();
    
    z->beginTime=fmax(x->beginTime, y->beginTime);
    z->endTime=fmin(x->endTime, y->endTime);
    
    while(i->time >= z->endTime) i++;
    while(j->time >= z->endTime) j++;
    
    s=z->beginTime;
    t=z->endTime;
    
    while(i->time > s) {
        
        computeSegmentUntil(z, *i, t, j, z_max);
        z_max=z->front().value;
        
        if(j->time == i->time) j++;
        t=i->time;
        i++;
    }
    
    if(i->time == s) computeSegmentUntil(z, *i, t, j, z_max);
    else computeSegmentUntil(z, Sample(s, i->valueAt(s), i->derivative), t, j, z_max);
    
    z->simplify();
    
    return z;
}


Signal * computeBoundedEventually(Signal * x, double a) {
#ifdef DEBUG__
    cout << "Entering computeBoundedEventually" << endl;
#endif
    //cout << ">  computeBoundedEventually: BREACH                 IN." << endl;
    //cout << "IN: " << *x << endl;
    
    Signal *z, *z1,
            *z2, *z3;
    
    z1=plateauMax(x, a);
    z2=new Signal(*x);
    z2->resize(x->beginTime + a, x->endTime + a, BOTTOM);
    z2->shift(-a);
    z3=computeOr(z2, z1);
    
    delete z1;
    delete z2;
    
    z=computeOr(x, z3);
    
    delete z3;
    z->simplify();
    
    //cout << "OUT: " << *z << endl;
    //cout << "<  computeBoundedEventually:    BREACH              OUT." << endl;
    return z;
}


Signal * computeBoundedGlobally(Signal * x, double a) {
    Signal *z, *z1, *z2, *z3;
    
    z1=plateauMin(x, a);
    
    z2=new Signal(*x);
    z2->resize(x->beginTime + a, x->endTime + a, TOP);
    z2->shift(-a);
    
    z3=computeAnd(z2, z1);
    delete z1;
    delete z2;
    
    z=computeAnd(x, z3);
    delete z3;
    
    z->simplify();
    
    return z;
}


Signal * computeTimedUntil(Signal * x, Signal * y, double a, double b) {
    Signal *z, *z1, *z2, *z3, *z4;
    
    if (a>0)
        z1=computeBoundedGlobally(x,a);
    
    z2=computeBoundedEventually(y,b-a);
    z3=computeUntil(x,y);
    z4=computeAnd(z2,z3);
    
    if (a>0) {
        z4->shift(-a);
        z=computeAnd(z1,z4);
        delete z1;
    }
    else
        z=computeAnd(x,z4);
    
    delete z2;
    delete z3;
    delete z4;
    
    return z;
}

/*---------------------------------------------------------------------------*
 *                         CONJUNCTION SUBROUTINES                           *
 *---------------------------------------------------------------------------*/

//PRECONDITIONS: j->time < t, i.time < t.
//POSTCONDITIONS: j->time <= i.time.
void computeSegmentAnd(Signal * z, const Sample & i, double t, Signal::reverse_iterator & j) {
    bool continued=false;
    double s=j->time;
    
    //mexPrintf("i.time: %g j->time:%g, t:%g\n" , i.time, j->time, t);
    
    // for every sample *j in (i.time, t)
    while(s > i.time) {
        if(i.valueAt(t) < j->valueAt(t)) {
            if (i.valueAt(s) > j->value) {
                t=i.timeIntersect(*j);
                z->push_front(Sample(t,i.valueAt(t),i.derivative));
                z->push_front(Sample(s,j->value,j->derivative));
                continued=false;
            }
            else continued=true;
        }
        else if (i.valueAt(t) == j->valueAt(t)) {
            if (i.valueAt(s) > j->value) {
                if(continued) {
                    z->push_front(Sample(t,i.valueAt(t), i.derivative));
                    continued=false;
                }
                z->push_front(Sample(s,j->value,j->derivative));
            }
            else continued=true;
        }
        else {
            if (i.valueAt(s) < j->value) {
                if(continued) {
                    z->push_front(Sample(t,i.valueAt(t),i.derivative));
                }
                t=i.timeIntersect(*j);
                
                z->push_front(Sample(t,j->valueAt(t),j->derivative));
                continued=true;
            }
            else {
                if(continued) {
                    z->push_front(Sample(t,i.valueAt(t),i.derivative));
                    continued=false;
                }
                z->push_front(Sample(s,j->value,j->derivative));
            }
        }
        
        //increment reverse iterator j
        t = s;
        j++;
        s = j->time;
    }
    
    
    //here we may have j->time < i.time
    // "i" values of z are no longer "continued"
    s=i.time;
    if(i.valueAt(t) < j->valueAt(t)) {
        if (i.value > j->valueAt(s)) {
            t=i.timeIntersect(*j);
            z->push_front(Sample(t,i.valueAt(t),i.derivative));
            z->push_front(Sample(s,j->valueAt(s),j->derivative));
        }
        else {
            z->push_front(i);
        }
    }
    else if (i.valueAt(t) == j->valueAt(t)) {
        if (i.value > j->valueAt(s)) {
            if(continued) {
                z->push_front(Sample(t,i.valueAt(t), i.derivative));
            }
            z->push_front(Sample(s,j->valueAt(s),j->derivative));
        }
        else {
            z->push_front(i);
        }
    }
    else {
        if (i.value < j->valueAt(s)) {
            if(continued) {
                z->push_front(Sample(t,i.valueAt(t),i.derivative));
            }
            t=i.timeIntersect(*j);
            z->push_front(Sample(t,j->valueAt(t),j->derivative));
            z->push_front(i);
            
        }
        else {
            if(continued) {
                z->push_front(Sample(t,i.valueAt(t),i.derivative));
            }
            z->push_front(Sample(s,j->valueAt(s),j->derivative));
        }
    }
    
}

void computePartialAnd(Signal * z, Signal::reverse_iterator & i, Signal::reverse_iterator & j, double s, double t) {
    
#ifdef DEBUG__
    mexPrintf("Entering computePartialAnd\n");
#endif
    
    while(i->time > s) {
        computeSegmentAnd(z,*i,t,j);
        if (j->time == i->time) j++;
        t=i->time;
        i++;
    }
    
    if(i->time == s) {
        computeSegmentAnd(z, *i, t, j);
    }
    else {
        computeSegmentAnd(z, Sample(s,i->valueAt(s),i->derivative), t, j);
    }
#ifdef DEBUG__
    mexPrintf("Leaving computePartialAnd\n");
#endif

}


/*---------------------------------------------------------------------------*
 *                         DISJUNCTION SUBROUTINES                           *
 *---------------------------------------------------------------------------*/

// copy of computeSegmentAnd, operator "<" switched with operator ">"
void computeSegmentOr(Signal * z, const Sample & i, double t, Signal::reverse_iterator & j) {
    bool continued=false;
    double s=j->time;
    
    // for every sample *j in (i.time, t)
    while(s > i.time) {
        if(i.valueAt(t) > j->valueAt(t)) {
            if (i.valueAt(s) < j->value) {
                t=i.timeIntersect(*j);
                z->push_front(Sample(t,i.valueAt(t),i.derivative));
                z->push_front(Sample(s,j->value,j->derivative));
                continued=false;
            }
            else continued=true;
        }
        else if (i.valueAt(t) == j->valueAt(t)) {
            if (i.valueAt(s) < j->value) {
                if(continued) {
                    z->push_front(Sample(t,i.valueAt(t), i.derivative));
                    continued=false;
                }
                z->push_front(Sample(s,j->value,j->derivative));
            }
            else continued=true;
        }
        else {
            if (i.valueAt(s) > j->value) {
                if(continued) {
                    z->push_front(Sample(t,i.valueAt(t),i.derivative));
                }
                t=i.timeIntersect(*j);
                z->push_front(Sample(t,j->valueAt(t),j->derivative));
                continued=true;
            }
            else {
                if(continued) {
                    z->push_front(Sample(t,i.valueAt(t),i.derivative));
                    continued=false;
                }
                z->push_front(Sample(s,j->value,j->derivative));
            }
        }
        
        //increment reverse iterator j
        t = s;
        j++;
        s = j->time;
    }
    
    
    //here we may have j->time < i.time
    // "i" values of z are no longer "continued"
    s=i.time;
    if(i.valueAt(t) > j->valueAt(t)) {
        if (i.value < j->valueAt(s)) {
            t=i.timeIntersect(*j);
            z->push_front(Sample(t,i.valueAt(t),i.derivative));
            z->push_front(Sample(s,j->valueAt(s),j->derivative));
        }
        else {
            z->push_front(i);
        }
    }
    else if (i.valueAt(t) == j->valueAt(t)) {
        if (i.value < j->valueAt(s)) {
            if(continued) {
                z->push_front(Sample(t,i.valueAt(t), i.derivative));
            }
            z->push_front(Sample(s,j->valueAt(s),j->derivative));
        }
        else {
            z->push_front(i);
        }
    }
    else {
        if (i.value > j->valueAt(s)) {
            if(continued) {
                z->push_front(Sample(t,i.valueAt(t),i.derivative));
            }
            t=i.timeIntersect(*j);
            z->push_front(Sample(t,j->valueAt(t),j->derivative));
            z->push_front(i);
        }
        else {
            if(continued) {
                z->push_front(Sample(t,i.valueAt(t),i.derivative));
            }
            z->push_front(Sample(s,j->valueAt(s),j->derivative));
        }
    }
    
}


void computePartialOr(Signal * z, Signal::reverse_iterator & i, Signal::reverse_iterator & j, double s, double t) {
    
    while(i->time > s) {
        computeSegmentOr(z,*i,t,j);
        if (j->time == i->time) j++;
        t=i->time;
        i++;
    }
    
    if(i->time == s) computeSegmentOr(z, *i, t, j);
    else computeSegmentOr(z, Sample(s,i->valueAt(s),i->derivative), t, j);
    
}



/*---------------------------------------------------------------------------*
 *                              UNTIL SUBROUTINES                            *
 *---------------------------------------------------------------------------*/

void computePartialEventually(Signal* z, Signal::reverse_iterator & i, double s, double t) {
    bool continued=false;
    double z_max=BOTTOM;
    while(i->time > s) {
        if(i->derivative >= 0) {
            if(z_max < i->valueAt(t)) {
                if(continued) {
                    z->push_front(Sample(t,z_max,0));
                }
                z_max=i->valueAt(t);
            }
            continued=true;
            //z->push_front(Sample(i->time, z_max, 0));
        }
        else if(i->valueAt(t) >= z_max) {
            if(continued) {
                z->push_front(Sample(t,z_max,0));
                continued=false;
            }
            z_max=i->value;
            z->push_front(*i);
        }
        else if(z_max >= i->value) {
            continued=true;
            //z->push_front(Sample(i->time, z_max, 0));
        }
        else {
            z->push_front(Sample(i->time + (z_max-i->value)/i->derivative, z_max, 0)); //time at which y reaches value next_z
            z->push_front(*i);
            z_max=i->value;
            continued=false;
        }
        
        t=i->time;
        i++;
    }
    
    //leftmost sample *i may not be on s
    //"z_max" values of z are not longer "continued".
    if(i->derivative >= 0) {
        if(z_max < i->valueAt(t)) {
            if(continued) {
                z->push_front(Sample(t,z_max,0));
            }
            z_max=i->valueAt(t);
        }
        z->push_front(Sample(s, z_max, 0));
    }
    else if(i->valueAt(t) >= z_max) {
        if(continued) {
            z->push_front(Sample(t,z_max,0));
        }
        z->push_front(Sample(s, i->valueAt(s), i->derivative));
    }
    else if(z_max >= i->valueAt(s)) {
        z->push_front(Sample(s, z_max, 0));
    }
    else {
        z->push_front(Sample(s + (z_max-i->value)/i->derivative, z_max, 0)); //time at which y reaches value next_z
        z->push_front(Sample(s, i->valueAt(s), i->derivative));
    }
    
}

void computeSegmentUntil(Signal * z, const Sample & i, double t, Signal::reverse_iterator & j, double z_max) {

#ifdef DEBUG__
    mexPrintf("Entering computeSegmentUntil\n");
#endif


    Signal *z1, *z2, *z3;
    Signal::reverse_iterator k, l;
    double s=i.time;
    
    if(i.derivative <= 0) {
        z1=new Signal();
        computeSegmentAnd(z1, i, t, j);
        
        z2=new Signal();
        k=z1->rbegin();
        computePartialEventually(z2, k, s, t);
        delete z1;
        
        l=z2->rbegin();
        computeSegmentOr(z, Sample(s, fmin(z_max, i.valueAt(t)), 0), t, l);
        delete z2;
    }
    else {
        z1=new Signal();
        computePartialEventually(z1, j, s, t);
        
        z2=new Signal();
        k=z1->rbegin();
        computeSegmentAnd(z2, i, t, k);
        delete z1;
        
        z1=new Signal();
        z3=new Signal();
        z3->push_front(Sample(s, z_max, 0));
        k=z3->rbegin();
        computeSegmentAnd(z1, i, t, k);
        delete z3;
        
        k=z1->rbegin();
        l=z2->rbegin();
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
#ifdef DEBUG__
    cout << "Entering plateauMax" << endl;
#endif
    bool new_candidate, end_candidate;
    double t, v;
    Signal::reverse_iterator j;
    Sequence M; //sorted in ascending times from front to back
    Sequence y; //maximum of x(t) and x(t-) at discontinuity points of x
    Sequence::iterator i;
    
    Signal * z=new Signal();
    z->beginTime=x->beginTime;
    z->endTime=x->endTime;
    
    
    //PRECOMPUTATION: list the local maximums of x
    t=x->endTime;
    v=BOTTOM;
    for(j=x->rbegin(); j!= x->rend(); j++) {
        y.push_front(Point(t, fmax(v,j->valueAt(t))));
        //    cout << "y.time:" << y.front().time << " y->value:" << y.front().value << endl;
        t=j->time;
        v=j->value;
    }
    y.push_front(Point(x->front().time, x->front().value));
    
    //INIT: read values in [0, a)
    i=y.begin();
    
    while(i->time < x->beginTime + a) {
        while(!M.empty() && i->value >= M.back().value) {
            M.pop_back();
        }
        //    cout << "i->time:" << i->time << " i->value:" << i->value << endl;
        
        M.push_back(*i);
        i++;
    }
    
    
    if(i->time == x->beginTime + a)
        new_candidate=true;
    else
        new_candidate=false;
    
    end_candidate=false;
    t=x->beginTime;
    
    //  cout << "M.front "  << M.front().value << endl;;
    
    //STEP
    bool cont=true;
    while(cont) {
        
        //UPDATE OF CANDIDATE LIST
        //candidate crosses t: remove it from list
        if(end_candidate) {
            //      cout << "end_candidate: M.pop front" << endl;
            M.pop_front();
        }
        
        //sample crosses t+a: add it to candidate list
        if(new_candidate) {
            //cout << "new_candidate: doing stuff" << endl;
            while(!M.empty() && i->value >= M.back().value) {
                M.pop_back();
            }
            //detect if new maximum is found
            if(!M.empty()) {
                //if M non empty then t + a does not generate new maximum
                new_candidate=false;
            }
            //add candidate
            M.push_back(*i);
            
            //increment iterator
            i++;
        }
        
        //OUTPUT OF NEW MAXIMUM
        //      cout << "Output new maximum" << endl;
        //no candidate
        if(M.empty()) {
            //cout << "No candidate" << endl;
            z->push_back(Sample(t, BOTTOM, 0));
        }
        //next best candidate
        else {
            z->push_back(Sample(t, M.front().value, 0));
            //	cout << "Pushed " << M.front().value << endl;
        }
        
        //NEXT EVENT DETECTION
        if(! M.empty()) {
            if(i != y.end()) {
                if(i->time - a == M.front().time) {
                    t=M.front().time;
                    new_candidate=true;
                    end_candidate=true;
                }
                else if(i->time - a < M.front().time) {
                    t=i->time - a;
                    new_candidate=true;
                    end_candidate=false;
                }
                else { //M.back().time < i->time - a
                    t=M.front().time;
                    new_candidate=false;
                    end_candidate=true;
                }
            }
            else {
                t=M.front().time;
                new_candidate=false;
                end_candidate=true;
            }
        }
        else {
            if(i != y.end()) {
                t=i->time - a;
                new_candidate=true;
                end_candidate=false;
            }
            else {
                new_candidate=false;
                end_candidate=false;
            }
        }
        cont = (new_candidate||end_candidate);
    }
    
    if(z->back().time==z->endTime) z->pop_back();
    
    return z;
}


//copy of plateauMax, operator < switched with operator >, TOP replaces BOTTOM
Signal * plateauMin(Signal * x, double a) {
    bool new_candidate, end_candidate;
    double t, v;
    Signal::reverse_iterator j;
    Sequence M; //sorted in ascending times from front to back
    Sequence y; //minimum of x(t) and x(t-) at discontinuity points of x
    Sequence::iterator i;
    
    Signal * z=new Signal();
    z->beginTime=x->beginTime;
    z->endTime=x->endTime;
    
    
    //PRECOMPUTATION: list the local minimums of x
    t=x->endTime;
    v=TOP;
    for(j=x->rbegin(); j!= x->rend(); j++) {
        y.push_front(Point(t, fmin(v,j->valueAt(t))));
        t=j->time;
        v=j->value;
    }
    y.push_front(Point(x->front().time, x->front().value));
    
    //INIT: read values in [0, a)
    i=y.begin();
    while(i->time < x->beginTime + a) {
        while(!M.empty() && i->value <= M.back().value) {
            M.pop_back();
        }
        M.push_back(*i);
        i++;
    }
    if(i->time == x->beginTime + a) new_candidate=true;
    else new_candidate = false;
    end_candidate=false;
    t=x->beginTime;
    
    //STEP
    bool cont=true;
    while(cont) {
        //UPDATE OF CANDIDATE LIST
        //candidate crosses t: remove it from list
        if(end_candidate) {
            M.pop_front();
        }
        
        //sample crosses t+a: add it to candidate list
        if(new_candidate) {
            while(!M.empty() && i->value <= M.back().value) {
                M.pop_back();
            }
            //detect if new minimum is found
            if(!M.empty()) {
                //if M non empty then t + a does not generate new minimum
                new_candidate=false;
            }
            //add candidate
            M.push_back(*i);
            
            //increment iterator
            i++;
        }
        
        //OUTPUT OF NEW MINIMUM
        //no candidate
        if(M.empty()) {
            z->push_back(Sample(t, TOP, 0));
        }
        //next best candidate
        else {
            z->push_back(Sample(t, M.front().value, 0));
        }
        
        //NEXT EVENT DETECTION
        if(! M.empty()) {
            if(i != y.end()) {
                if(i->time - a == M.front().time) {
                    t=M.front().time;
                    new_candidate=true;
                    end_candidate=true;
                }
                else if(i->time - a < M.front().time) {
                    t=i->time - a;
                    new_candidate=true;
                    end_candidate=false;
                }
                else { //M.back().time < i->time - a
                    t=M.front().time;
                    new_candidate=false;
                    end_candidate=true;
                    
                }
            }
            else {
                t=M.front().time;
                new_candidate=false;
                end_candidate=true;
            }
        }
        else {
            if(i != y.end()) {
                t=i->time - a;
                new_candidate=true;
                end_candidate=false;
            }
            else {
                new_candidate=false;
                end_candidate=false;
            }
        }
        cont = (new_candidate||end_candidate);
    }
    
    if(z->back().time==z->endTime) z->pop_back();
    
    return z;
}


