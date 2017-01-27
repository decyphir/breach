/*
 * signal_expr.h
 *
 *  Created on: Apr 11, 2014
 *      Author: alex
 */

#ifndef SIGNAL_EXPR_H_
#define SIGNAL_EXPR_H_

#include "transducer.h"

namespace CPSGrader {

    typedef vector<vector<double> >::const_iterator trace_data_iterator;

    /* signal_transducer: read data without transformation - has a variable name which must be
     * defined in the signal map   */
    class signal_transducer: public transducer {

    public:
        string variable;
        long int td_idx;

	signal_transducer() :
        variable ( "x" ) {
        };
        
        signal_transducer(string var) {
            variable = var;
        };
        

        void init_horizon() {};

        virtual signal_transducer* clone() const {

            signal_transducer * new_clone;
            new_clone = new signal_transducer(variable); 
            new_clone->td_idx = td_idx;
            new_clone->trace_data_ptr = trace_data_ptr;
            new_clone->signal_map = signal_map;
            new_clone->param_map = param_map;
            return new_clone;
        };

        virtual ~signal_transducer() {};
        
        virtual double compute_robustness();

        virtual void set_trace_data_ptr(const trace_data &trace) {
            trace_data_ptr = &trace;
            td_idx = 0;
        };
        
        void print(ostream &os) const {
            os << variable << "[" << "t" << "]";
        };
        
        void print() const {
            print(cout);
        };
        
    };

    class constant_transducer: public transducer {
    public:
        string param;
        double value;
        long int td_idx;

        constant_transducer() {
            value = 0.;
        }
        ;
        constant_transducer(double val) {
            value = val;
        }
        ;
        constant_transducer(const string &p);

        void init_horizon() {};

        virtual constant_transducer* clone() const {
            constant_transducer * T = new constant_transducer();
            T->td_idx = td_idx;
            T->trace_data_ptr = trace_data_ptr;
            T->param = param;
            T->value = value;
            return T;
        }

        virtual void set_trace_data_ptr(const trace_data &trace) {
            trace_data_ptr = &trace;
            td_idx = 0;
        };
        
        virtual ~constant_transducer() {
        };
        
        virtual double compute_robustness();

        void print(ostream &os) const {
            param.empty()? os << value : os << param;
        };

        void print() const {
            print(cout);
        };
    };


    
    /* Unary operators on signals */

    class abs_transducer: public unary_transducer {
    public:

        explicit abs_transducer(transducer *_child) :
        transducer(), unary_transducer(_child) {
        }

        virtual abs_transducer* clone() const {
            transducer * child_clone = child->clone();
            return new abs_transducer(child_clone);
        }

        virtual ~abs_transducer() {
        }

        double compute_robustness();

        void print() const {
            print(cout);
        }
        ;

        virtual void print(ostream &os) const {
            os << "abs( ";
            child->print(os);
            os << " )";
        }
        ;
    };

    /* binary signal operators - all operators assume that left and right operand have the
     * same number of samples.
     */

    class plus_transducer: public binary_transducer {
    public:

        explicit plus_transducer(transducer *_childL, transducer *_childR) :
        transducer(), binary_transducer(_childL, _childR) {
        };

        virtual plus_transducer* clone() const {
            transducer * childL_clone = childL->clone();
            transducer * childR_clone = childR->clone();
            return new plus_transducer(childL_clone, childR_clone);
        }

        virtual ~plus_transducer() {
        }

        double compute_robustness();
        void print() const {
            print(cout);
        }
        ;

        virtual void print(ostream &os) const {
            childL->print(os);
            os << "+";
            childR->print(os);
        }
        ;
    };

    class minus_transducer: public binary_transducer {
    public:

        explicit minus_transducer(transducer *_childL, transducer *_childR) :
		transducer(), binary_transducer(_childL, _childR) {
        };

        virtual minus_transducer* clone() const {
            transducer * childL_clone = childL->clone();
            transducer * childR_clone = childR->clone();
            return new minus_transducer(childL_clone, childR_clone);
        }

        virtual ~minus_transducer() {
        }

        double compute_robustness();
        void print() const {
            print(cout);
        }
        ;

        virtual void print(ostream &os) const {
            childL->print(os);
            os << "-";
            childR->print(os);
        }
        ;
    };

    class mult_transducer: public binary_transducer {
    public:

        explicit mult_transducer(transducer *_childL, transducer *_childR) :
        transducer(), binary_transducer(_childL, _childR) {
        };

        virtual mult_transducer* clone() const {
            transducer * childL_clone = childL->clone();
            transducer * childR_clone = childR->clone();
            return new mult_transducer(childL_clone, childR_clone);
        }

        virtual ~mult_transducer() {
        }

        double compute_robustness();

        void print() const {
            print(cout);
        }
        ;

        virtual void print(ostream &os) const {
            childL->print(os);
            os << "*";
            childR->print(os);
        }
        ;
    };


};
#endif /* SIGNAL_EXPR_H_ */
