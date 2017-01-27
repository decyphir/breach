#include "stdafx.h"
#include <transducer.h>
#include <algorithm>

namespace CPSGrader {

    /* init horizon */
    void transducer::init_horizon() {
    }

    void unary_transducer::init_horizon() {
        child->set_horizon(start_time, end_time);
        child->init_horizon();
    }

    void binary_transducer::init_horizon() {
        childL->set_horizon(start_time, end_time);
        childR->set_horizon(start_time, end_time);
        childL->init_horizon();
        childR->init_horizon();
    }

    void timed_unary_transducer::init_horizon() {

        //	cout << "timed_unary_transducer::init_horizon - start_time:" << start_time <<  endl;
        // checks whether a and b are given by parameters, and assign corresponding values
        double a,b;
        if (!get_param(I->begin_str,a)) a = I->begin;
        if (!get_param(I->end_str,b)) b = I->end;

        // update start_time and end_time of child
        child->set_horizon(start_time+a,end_time+b);
        child->init_horizon();

    }

    void timed_binary_transducer::init_horizon() {

        double a,b;
        if (!get_param(I->begin_str,a)) a = I->begin;
        if (!get_param(I->end_str,b)) b = I->end;

        // update start_time and end_time of children
        childL->set_horizon(start_time, end_time+b);
        childR->set_horizon(start_time+a, end_time+b);
        childL->init_horizon();
        childR->init_horizon();
    }

    // get end time complete for transducers
    double transducer::get_end_complete() {
        return z.endTime;
    }
    double transducer::get_end_complete_low() {
        return z_low.endTime;
    }
    double transducer::get_end_complete_up() {
        return z_up.endTime;
    }

    // get end time complete for timed transducers
    double timed_unary_transducer::get_end_complete() {
        double b;
        if (!get_param(I->end_str,b)) b = I->end;
        return child->z.endTime-b;
    }

    double timed_unary_transducer::get_end_complete_low() {
        double b;
        if (!get_param(I->end_str,b)) b = I->end;	
        return child->z_low.endTime-b;
    }

    double timed_unary_transducer::get_end_complete_up() {
        double b;
        if (!get_param(I->end_str,b)) b = I->end;
        return child->z_up.endTime-b;
    }


    /* Compute robustness */
    double not_transducer::compute_robustness() {
        //	cout << "Computing robustness of NOT"
        child->compute_robustness();
        z.compute_not(child->z);
        return z.front().value;
    }

    double ev_transducer::compute_robustness() {
#ifdef DEBUG__
        printf(">  ev_transducer::compute_robustness:         IN." );
        cout << "   I->a: " << I->begin << "   I->b: " << I->end << endl;
        cout << "   start_time:" << start_time << " end_time:" << end_time << endl;
#endif

        double a,b;
        if (!get_param(I->begin_str,a)) a = I->begin;
        if (!get_param(I->end_str,b)) b = I->end;

        child->compute_robustness();
        //if (z.empty()) // first run 

        z.compute_timed_eventually(child->z, a, b);

        //else // is there new stuff ? 
        //    {
        //       new z.endTime+b,  is anything beyond that 
        //       useless 
        //       child->z.endTime
        //    }
        

        double et = min(z.endTime,end_time);
        (child->z).resize(et-b, (child->z).endTime, 0.);

        z.resize(start_time,max(start_time,et), 0.);

#ifdef DEBUG__
        cout << "OUT:" << z << endl;
        cout << "<  ev_transducer::compute_robustness:         OUT." << endl;
#endif
        return z.front().value;
    }

    double alw_transducer::compute_robustness() {
#ifdef DEBUG__
        printf(">  alw_transducer::compute_robustness:        IN." );
        cout << "   I->a: " << I->begin << "   I->b: " << I->end << endl;
        cout << "   start_time:" << start_time << " end_time:" << end_time << endl;
#endif

        double a,b;
        if (!get_param(I->begin_str,a)) a = I->begin;
        if (!get_param(I->end_str,b)) b = I->end;
        child->compute_robustness();

        z.compute_timed_globally(child->z, a, b);
        double et =min(z.endTime,end_time);
        z.resize(start_time,max(start_time, et ),0.);

#ifdef DEBUG__
        cout << "OUT:" << z << endl;
        printf("<  alw_transducer::compute_robustness:        OUT." );
#endif
        return z.front().value;
    }

    double and_transducer::compute_robustness() {
#ifdef DEBUG__
        printf( ">  and_transducer::compute_robustness:        IN." );
#endif

        childL->compute_robustness();
        childR->compute_robustness();

        Signal z1= childL->z;
        Signal z2= childR->z;

        double et = max(z1.endTime, z2.endTime);
        z1.endTime = et;
        z2.endTime = et;
        z.compute_and(z1, z2);

        //	z.resize(start_time,max(start_time, z.endTime),0.);
#ifdef DEBUG__
        cout << "OUT:" << z << endl;
        printf("<  and_transducer::compute_robustness:        OUT." );
#endif
        return z.front().value;
    }

    double or_transducer::compute_robustness() {
#ifdef DEBUG__
        printf( ">  or_transducer::compute_robustness:         IN." );
#endif

        childL->compute_robustness();
        childR->compute_robustness();

        Signal z1= childL->z;
        Signal z2= childR->z;

        double et = max(z1.endTime, z2.endTime);
        z1.endTime = et;
        z2.endTime = et;

        z.compute_or(z1,z2);

#ifdef DEBUG__
        cout << "OUT:" << z << endl;
        printf("<  or_transducer::compute_robustness:         OUT." );
#endif
        return z.front().value;
    }

    double implies_transducer::compute_robustness() {

        childL->compute_robustness();
        childR->compute_robustness();

        Signal z1= childL->z;
        Signal z2= childR->z;
        double et = max(z1.endTime, z2.endTime);

        Signal *rob = computeImplies(&z1, &z2);
        z = *rob;
        delete rob;

        return z.front().value;
    }

    double until_transducer::compute_robustness() {

        double a,b;
        if (!get_param(I->begin_str,a)) a = I->begin;
        if (!get_param(I->end_str,b)) b = I->end;

        // update robustness of children
        childL->compute_robustness();
        childR->compute_robustness();

        z.compute_timed_until(childL->z, childR->z, a, b);
        z.resize(0.,z.endTime,0.);
        return z.front().value;
    }

    /* Utility functions */

    Signal transducer::get_signal() const {
        return z;
    }

    bool transducer::get_param(const string & param, double & val) {

        map<string, double>::iterator it;
        if ( (it= param_map.find(param)) != param_map.end()) {
            val = it->second;
            return true;
        }
        else{
            val = 0.;
            return false;
        }
    }

    std::ostream& operator<<(std::ostream& os, const transducer& T) {
        T.print(os);
        return os;
    }

    void transducer::print_trace() {
        for (auto ii = trace_data_ptr->begin(); ii != trace_data_ptr->end(); ii++){
            for (auto jj = (*ii).begin(); jj != (*ii).end(); jj++) {
                cout << *jj << " ";
            }
            cout << endl;
        }
    }
    
}

