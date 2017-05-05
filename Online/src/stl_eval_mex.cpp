// include stuff
#include "mex.h"
#include <string>
#include <iostream>
#include <sstream>
#include "stl_driver.h"

using namespace std;
using namespace CPSGrader;

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[] ) {


    stringstream fcout;
    
    if (nrhs<=1)
        mexErrMsgTxt("four inputs are expected: signal names, a formula (string), data (array), time interval (array of size 2).");
    
    /* read inputs: a string and data */
    char *signal_buf = mxArrayToString(prhs[0]);
    char *stl_buf = mxArrayToString(prhs[1]);   
    string phi_st = "signal "+ string(signal_buf) + "\n" + "phi:=" + string(stl_buf);
    
    size_t m = mxGetM(prhs[1]);
    size_t n=  mxGetN(prhs[1]);
    
    double *data_in = (double *) mxGetPr(prhs[1]);
    double *time_in = (double *) mxGetPr(prhs[2]);
    
    STLDriver stl_driver = STLDriver();
	
    vector<double> sample;
    for(int i = 0; i<n; i++){
        for(int j = 0; j<m; j++) {
            sample.push_back(data_in[j+ i*m]);
        }
        stl_driver.data.push_back(sample);
        sample.clear();
        //cout << endl;
    }
    
    transducer * phi; 
    double rob, rob_up, rob_low;
    rob = rob_up = rob_low = 0;
     
    Signal z, z_up, z_low;
	bool parse_success = stl_driver.parse_string(phi_st);

    if (parse_success) {
		phi = stl_driver.formula_map["phi"]->clone();
        phi->set_horizon(time_in[0], time_in[1]);
        phi->set_trace_data_ptr(stl_driver.data);
        
        cout << "Size of data:" <<        (phi->trace_data_ptr)->size() << endl;

           
        
#ifdef DEBUG__
        cout << "--------------------------------------------------------------" << endl;
        cout << "----------------------ROB-------------------------------------" << endl;
#endif
        rob    = phi->compute_robustness();
#ifdef DEBUG__
        cout << "--------------------------------------------------------------" << endl;
        cout << "---------------------ROB_UP-----------------------------------" << endl;
#endif
        rob_up = phi->compute_upper_rob();
#ifdef DEBUG__
        cout << "--------------------------------------------------------------" << endl;
        cout << "--------------------ROB_LOW-----------------------------------" << endl;
#endif
        rob_low= phi->compute_lower_rob();
#ifdef DEBUG__
        cout << "--------------------------------------------------------------" << endl;
        cout << "--------------------------------------------------------------" << endl;
#endif
        
        z =  phi->z;
        z_low = phi->z_low;
        z_up = phi->z_up;
    }
    else
    	mexErrMsgTxt("Problem parsing formula.");

    z.addLastSample();
    z_low.addLastSample();
    z_up.addLastSample();
    
    #ifdef DEBUG__
    cout << "rob:" << rob << endl;
    cout << "z:" << z  << endl;
    cout << "z_low:" << z_low << endl;
    cout << "z_up:" << z_up << endl;
    #endif

    /* compute and output robustness at 0 */
    int l =  z.size();
    plhs[0] = mxCreateDoubleMatrix(1,l, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,l, mxREAL);

    double *t_ptr = mxGetPr(plhs[0]);
    double *rob_ptr = mxGetPr(plhs[1]);
    Signal::const_iterator iter_z;
   
    int pos = 0;
    for(iter_z=z.begin(); iter_z!=z.end(); iter_z++) {
    	t_ptr[pos] = iter_z->time;
    	rob_ptr[pos] = iter_z->value;
    	pos++;
    }

    size_t l_low =  z_low.size();
    plhs[2] = mxCreateDoubleMatrix(1,l_low, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,l_low, mxREAL);

    double *t_low_ptr = mxGetPr(plhs[2]);
    double *rob_low_ptr = mxGetPr(plhs[3]);

    pos = 0;
    for(iter_z=z_low.begin(); iter_z!=z_low.end(); iter_z++) {
    	t_low_ptr[pos] = iter_z->time;
    	rob_low_ptr[pos] = iter_z->value;
    	pos++;
    }

    size_t l_up =  z_up.size();
    plhs[4] = mxCreateDoubleMatrix(1,l_up, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(1,l_up, mxREAL);

    double *t_up_ptr = mxGetPr(plhs[4]);
    double *rob_up_ptr = mxGetPr(plhs[5]);

    pos = 0;
    for(iter_z=z_up.begin(); iter_z!=z_up.end(); iter_z++) {
    	t_up_ptr[pos] = iter_z->time;
    	rob_up_ptr[pos] = iter_z->value;
    	pos++;
    }

    // clean
    mxFree(signal_buf);
    mxFree(stl_buf);
    delete phi;
}


//
