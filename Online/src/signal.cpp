#include "stdafx.h"
#include "signal.h"
#include "iomanip"

namespace CPSGrader {

    /* 
     * class Signal member functions
     */
    Signal::Signal(double T, double V, int n) {

        beginTime=T;
        endTime=T;
        push_back(Sample(T, V, 0.));

    }

    // TODO refactor piecewise constant vs piecewise linear...
    Signal::Signal(double * T, double * V, int n) {
#ifdef DEBUG__
        printf(">> Signal::Signal:                              IN." );
#endif

        beginTime=T[0];
        endTime=T[n-1];

        if (n==1)
            push_back(Sample(T[0], V[0], 0.));
        else {
            for(int i=0; i < n-1; i++) {
                push_back(Sample(T[i], V[i], (V[i+1]-V[i]) / (T[i+1] - T[i])));
                //push_back(Sample(T[i], V[i], 0.));
            }
            push_back(Sample(T[n-1], V[n-1], 0.));
        }

#ifdef DEBUG__
        printf("<< Signal::Signal:                            OUT." );
#endif
    }

    void Signal::addLastSample() {
        if (endTime> back().time) 
            push_back(Sample(endTime, back().valueAt(endTime),0.));
    }
    
    void Signal::appendSample(double t, double v) {

        if ((t<=endTime)&&size()>0)
            return;

        if (size()==0) {
            push_back(Sample(t,v,0.));
            beginTime=t;
        }
        else {
            back().derivative = (v-back().value) / (t - back().time);
            //back().derivative = 0.;
            push_back(Sample(t,v,0));
            endTime= t;
        }
    }

    void Signal::appendSample(double t, double v, double d) {

        //cout << "Appending t=" << t << " v=" << v << endl;

        if ((t<=endTime)&&size()>0)
            return;

        if (size()==0) {
            push_back(Sample(t,v,d));
            beginTime=t;
        }
        else {
            back().derivative = (v-back().value) / (t - back().time);
            //back().derivative = 0.;
            push_back(Sample(t,v,d));
            endTime= t;
        }
    }

    
    void Signal::appendSignal(Signal s) {

        Signal::const_iterator iter_s;

        for(iter_s = s.begin(); iter_s != s.end(); iter_s++) {
            appendSample((*iter_s).time,(*iter_s).value,(*iter_s).derivative);
        }

    }

    //remove linear interpolations
    void Signal::simplify() {
#ifdef DEBUG___
        printf(">>>Signal::simplify:                          IN." );
        cout << "IN: " << *this << endl;
#endif

        push_back(front());
        pop_front();
        while (front().time != beginTime) {
            if( back().valueAt(front().time) != front().value || back().derivative != front().derivative ) {
                push_back(front());
            }
            pop_front();
        }

        // check last sample
        if (back().time < endTime)
            push_back(Sample(endTime, back().valueAt(endTime), 0.));

#ifdef DEBUG___	
        cout << "OUT: " << *this << endl;
        printf("<<<Signal::simplify:                          OUT.\n";
#endif
    }

    void Signal::resize(double s, double t, double v) {
        // Resize signal to begin at time s and end at time t 

#ifdef DEBUG__
            printf(">>>Signal::resize:                            IN.\n");
        cout << "to start_time:" << s << " and end_time:" << t << endl;
        cout << "IN: " << *this << endl;
#endif

        if ( t<s-1e-14 ) {
            clear();
            beginTime=0.;
            endTime=0.;

#ifdef DEBUG__
            cout << "OUT(premature): " << *this << endl;
            printf("<<<Signal::resize:                            OUT.\n");
#endif
            return;
        }
        else 
            if (t < s)
                t = s;  // hope I don't regret this.
        Sample first;

        //trim or extend front of signal
        if(beginTime > s) {
            //double der = (front().value-v)/(front().time-s);
            push_front(Sample(s, front().value, 0));
        }
        else {
            while((!empty())&&(front().time < s)) {
                first=front();
                pop_front();
            }
            if (empty()) {
                //			cout << "push empty " << first << endl;
                push_front(Sample(s, first.valueAt(s), 0));
                if (endTime < s)
                    endTime = s;
            }
            else {
                if (front().time > s)  {
                    double val = first.valueAt(s);				
                    push_front(Sample(s,val,first.derivative));
                }
            }
        }
        //trim or extend back of signal
        if(endTime < t) {
            //		cout << "push_back here" << endl;
            if (back().value != v || back().derivative != 0.)
                push_back(Sample(endTime, v, 0));
        }
        else {
            while(!empty()&&back().time >t) {
                pop_back();
            }
        }
        if (empty()) {
            //		cout << "push_back empty" << endl;
            push_back(Sample(s, v, 0));
        }
        beginTime=s;
        endTime=t;
#ifdef DEBUG__
        cout << "OUT: " << *this << endl;
        printf("<<<Signal::resize:                            OUT.\n");
#endif
    }

    void Signal::shift(double a) {
        Signal::iterator i;

        beginTime=beginTime + a;
        endTime=endTime + a;

        for(i = begin(); i != end(); i++) {
            i->time=i->time + a;
        }
    }

    void Signal::removeInf() {

        while((!empty())&&
              (back().value==TOP ||
               back().derivative == TOP ||
               back().value==BOTTOM ||
               back().derivative == BOTTOM)) {
            pop_back();
        }
    }

    /*
     * friend functions
     */
    std::ostream & operator<<(std::ostream & out, const Point & point) {
        out << point.time << ";" << point.value ;
        return out;
    }

    std::ostream & operator<<(std::ostream & out, const Sample & sample) {

        out << std::setprecision(8) << std::setw(12) << sample.time << "  " << std::setprecision(8) << std::setw(12)  << sample.value << "  " << sample.derivative;
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

}
