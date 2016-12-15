#include "stdafx.h"
#include <transducer.h>
#include <algorithm>
#include <math.h>

namespace CPSGrader {

    double transducer::update_robustness(){

#ifdef DEBUG__
        cout << ">  transducer:update_robustness            IN."<< endl;
#endif

        compute_robustness();

        
#ifdef DEBUG__
        cout << "<  transducer:update_robustness            OUT."<< endl;
#endif
        return z.front().value;
    };


    double unary_transducer::update_robustness(){

#ifdef DEBUG__
        cout << ">  transducer:update_robustness            IN."<< endl;
#endif

        child->update_robustness();

        Sample last = child->z.back();
        child->z.clear();
        child->z.push_front(last);
        
#ifdef DEBUG__
        cout << "<  transducer:update_robustness            OUT."<< endl;
#endif
        return z.front().value;
    };

        
}
