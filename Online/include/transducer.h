#ifndef __TRANSDUCER_H
#define __TRANSDUCER_H
#include <map>
#include <stdlib.h>
#include <vector>
#include <iterator>
#include <string>
#include "robustness.h"
#include "interval.h"
#include "mex.h"

using namespace std;

#pragma warning (disable : 4482)

namespace CPSGrader {

    typedef vector<vector<double> > trace_data;

    /* Virtual classes */

    /** transducer virtual class. Provides a compute_robustness method. */
    class transducer {

    public:

        const trace_data *trace_data_ptr; // signal data to monitor vector of vector of double

        map<string, double> param_map;   //  parameter values
        map<string, int>  signal_map;    //  maps signal name to index in trace data 

        // interval of time for which the transducer needs to provide values
        double start_time, end_time;

        // z is neutral semantics, z_up upper bound, z_low lower bound
        Signal z, z_up, z_low;

        transducer(): start_time(0.), end_time(0.) {};
        
        virtual transducer * clone() const =0;

        virtual ~transducer() {};

        // Initializes horizons
        virtual void init_horizon()=0;

        inline void set_horizon(double s, double e) {
            start_time = s;
            end_time = e;
            init_horizon();
        }

        // set trace data 
        // TODO should be done at the constructor, parser and cloning level...
        virtual void set_trace_data_ptr(const trace_data &trace) {
            trace_data_ptr= &trace;
        }

        // compute quantitative semantics for current data
        virtual double compute_robustness()=0;
        virtual double compute_lower_rob();
        virtual double compute_upper_rob();

        // update quantitative semantics based on new data
        virtual double update_robustness();
        virtual double update_lower_rob(){return 0.;};
        virtual double update_upper_rob(){return 0.;};
                
        // get interval for which we have complete information to compute robustness
        virtual double get_end_complete();
        virtual double get_end_complete_low();
        virtual double get_end_complete_up();

        //TODO fix print mess
        virtual void print() const =0;
        virtual void print(ostream &os) const =0;

        void print_trace(); 

        // looks into param_map for a parameter value - returns success
        bool get_param(const string&, double &);
        Signal get_signal() const;

    };

    // unary transducers (one child)
    class unary_transducer: virtual public transducer {
    public:
        unary_transducer(): transducer() {
            child = nullptr;
        };

        explicit unary_transducer(transducer *_child) :
			transducer(), child(_child) {
            init_horizon();
        };

        void init_horizon();
        virtual void set_trace_data_ptr(const trace_data &trace) {
            trace_data_ptr= &trace;
            child->set_trace_data_ptr(trace);
        }

        // update quantitative semantics based on new data
        virtual double update_robustness();
        virtual double update_lower_rob(){return 0.;};
        virtual double update_upper_rob(){return 0.;};


        
        virtual ~unary_transducer() {
            delete child;
        };

        transducer *child;

    };

    // Binary transducers
    class binary_transducer: virtual public transducer {

    public:

        binary_transducer(): transducer() {
            childL = nullptr;
            childR = nullptr;
            init_horizon();
        };

        explicit binary_transducer(transducer *_childL,transducer *_childR) :
			transducer(), childL(_childL), childR(_childR) {
            init_horizon();
        };

        void init_horizon();

        virtual void set_trace_data_ptr(const trace_data &trace) {
            trace_data_ptr= &trace;
            childL->set_trace_data_ptr(trace);
            childR->set_trace_data_ptr(trace);
        }

        virtual ~binary_transducer() {
            delete childL;
            delete childR;
        };

        transducer *childL;
        transducer *childR;
    };



    // timed unary transducer (output depends on children output at the some time )
    class timed_unary_transducer: virtual public unary_transducer {
    public:

        interval *I;

        timed_unary_transducer():unary_transducer()  {
            I = new interval();
            init_horizon();
        };

        timed_unary_transducer(transducer *_child): unary_transducer(_child){
            I = new interval();
            init_horizon();
        };

        timed_unary_transducer(interval *_I, transducer *_child): unary_transducer(_child), I(_I){
            init_horizon();
        };

        void init_horizon();

        virtual double get_end_complete();
        virtual double get_end_complete_low();
        virtual double get_end_complete_up();

        virtual ~timed_unary_transducer() {
            delete I;
        };
    };


    // timed binary transducer (output depends on children output at some time )
    class timed_binary_transducer: virtual public binary_transducer {
    public:

        interval *I;

        timed_binary_transducer():binary_transducer()  {
            I = new interval();
            init_horizon();
        };

        timed_binary_transducer(transducer *_childL, transducer *_childR): binary_transducer(_childL, _childR){
            I = new interval();
            init_horizon();
        };

        timed_binary_transducer( transducer *_childL, interval *_I,transducer *_childR): binary_transducer(_childL, _childR), I(_I){
            init_horizon();
        };

        void init_horizon();

        // TODO
        //virtual double get_end_time_complete();
        //virtual double get_end_time_complete_low();
        //virtual double get_end_time_complete_up();

        virtual ~timed_binary_transducer() {
            delete I;
        };
    };

    /* Boolean transducers */
    class not_transducer: public unary_transducer {

    public:

        explicit not_transducer(transducer *_child):
            unary_transducer(_child) {
        };

        not_transducer* clone() const {
            transducer * child_clone= child->clone();
            return new not_transducer(child_clone);
        }

        // empty destructor: child is killed by unary_transducer destructor
        ~not_transducer() {};

        double compute_robustness();

        double compute_lower_rob();

        double compute_upper_rob();

        void print() const{
            print(cout);
        };

        void print(ostream &os) const {
            os << "not ";
            child->print(os);
        }
        ;

    };

    class and_transducer: public binary_transducer {
    public:

        explicit and_transducer(transducer *_childL, transducer *_childR) :
            transducer(), binary_transducer(_childL, _childR)
        {
        };

        virtual and_transducer* clone() const {
            transducer * childL_clone= childL->clone();
            transducer * childR_clone= childR->clone();
            return new and_transducer(childL_clone, childR_clone);
        }

        virtual ~and_transducer() {
        }

        double compute_robustness();

        double compute_lower_rob();

        double compute_upper_rob();


        void print() const{
            print(cout);
        };

        virtual void print(ostream &os) const {
            childL->print(os);
            os << " and ";
            childR->print(os);
        }
        ;

    };

    class or_transducer: public binary_transducer {

    public:

        explicit or_transducer(transducer *_childL, transducer *_childR) :
            transducer(), binary_transducer(_childL, _childR){
        };

        virtual or_transducer* clone() const {
            transducer * childL_clone= childL->clone();
            transducer * childR_clone= childR->clone();
            return new or_transducer(childL_clone, childR_clone);
        };

    
        virtual ~or_transducer() {
        };

        double compute_robustness();

        double compute_lower_rob();

        double compute_upper_rob();

        void print() const{
            print(cout);
        };

        virtual void print(ostream &os) const {
            childL->print(os);
            os << " or ";
            childR->print(os);
        }
        ;

    };

    class implies_transducer: public binary_transducer {

    public:

        explicit implies_transducer(transducer *_childL, transducer *_childR) :
            transducer(), binary_transducer(_childL, _childR){
        };

        virtual implies_transducer * clone() const {
            transducer * childL_clone= childL->clone();
            transducer * childR_clone= childR->clone();

            return new implies_transducer(childL_clone, childR_clone);
        }

        virtual ~implies_transducer() {
        };

        double compute_robustness();

        double compute_lower_rob();

        double compute_upper_rob();

        void print() const{
            print(cout);
        };

        virtual void print(ostream &os) const {
            childL->print(os);
            os << " => ";
            childR->print(os);
        }
        ;

    };

    /* Temporal operators */

    class ev_transducer: public timed_unary_transducer {
    public:

        explicit ev_transducer(interval *_I, transducer *_child) :
            unary_transducer(_child),timed_unary_transducer(_I, _child)   {
        };

        virtual ev_transducer* clone() const {
            transducer * child_clone= child->clone();
            interval *Iclone = new interval(*I);
            return new ev_transducer(Iclone, child_clone);
        }

        virtual ~ev_transducer() { // child is killed by unary transducer and interval by timed_unary_transducer
        }

        double compute_robustness();
        double compute_lower_rob();
        double compute_upper_rob();


        void print() const{
            print(cout);
        };

        virtual void print(ostream &os) const {
            os << "ev_";
            I->print(os);
            os << " (";
            child->print(os);
            os << ") ";
        }
        ;

    };

    class alw_transducer:  public timed_unary_transducer {
    public:

        explicit alw_transducer(interval *_I, transducer *_child) :
            unary_transducer(_child),timed_unary_transducer(_I, _child)   {
        };

        virtual alw_transducer* clone() const {
            transducer * child_clone= child->clone();
            interval *Iclone = new interval(*I);
            return new alw_transducer(Iclone, child_clone);
        }

        virtual ~alw_transducer() {
        }

        double compute_robustness();
        double compute_lower_rob();
        double compute_upper_rob();

        void print() const{
            print(cout);
        };

        virtual void print(ostream &os) const {
            os << "alw_";
            I->print(os);
            os << " (";
            child->print(os);
            os << ") ";
        }

    };

    class until_transducer: public timed_binary_transducer {

    public:
        explicit until_transducer(transducer *_childL, interval *_I, transducer *_childR) :
			transducer(), binary_transducer(_childL, _childR), timed_binary_transducer(_childL, _I, _childR) {
        }

        virtual until_transducer* clone() const {
            transducer * childL_clone= childL->clone();
            transducer * childR_clone= childR->clone();
            interval *Iclone = new interval(*I);
            return new until_transducer(childL_clone, Iclone, childR_clone);
        }

        virtual ~until_transducer() {
        }

        double compute_robustness();
        double compute_lower_rob();
        double compute_upper_rob();

        void print() const{
            print(cout);
        };

        virtual void print(ostream &os) const {
            os << " (";
            childL->print(os);
            os << ") until_";
            I->print(os);
            os << " (";
            childR->print(os);
            os << ") ";
        }
        ;

    };

    enum comparator {
        LESSTHAN, GREATERTHAN
    };

    /* Atoms and signal expressions */
    class stl_atom: public binary_transducer {

    public:
        comparator comp;

        stl_atom(): transducer(), binary_transducer(){
            comp = comparator::GREATERTHAN;
        };

        stl_atom(transducer *_childL, comparator cp, transducer *_childR) :
            transducer(),
            binary_transducer(_childL,_childR),
            comp(cp) {
        };

        stl_atom(transducer *_childL, string cp, transducer *_childR) :
            transducer(),
            binary_transducer(_childL,_childR)
        {
            if (cp.compare("<") == 0)
                comp = comparator::LESSTHAN;
            else
                comp = comparator::GREATERTHAN;
        }
        ;

        virtual stl_atom* clone() const {
            transducer * childL_clone= childL->clone();
            transducer * childR_clone= childR->clone();
            return new stl_atom(childL_clone, comp, childR_clone);
        };

        ~stl_atom() {
        };

        double compute_robustness();

        virtual void print(ostream &os) const {
            childL->print(os);
            switch (comp) {
            case (comparator::LESSTHAN):
                os << " < ";
                break;
            case (comparator::GREATERTHAN):
                os << " > ";
                break;
            }

            childR->print(os);

        }

        void print() const {
            print(cout);
        }
        ;
    };

    std::ostream& operator<<(std::ostream& os, const transducer& T);

}
#endif
