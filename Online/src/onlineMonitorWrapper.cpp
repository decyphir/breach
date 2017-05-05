#include <string>
#include <iostream>
#include <sstream>
#include "stl_driver.h"

#define S_FUNCTION_LEVEL 2
#define S_FUNCTION_NAME  onlineMonitorWrapper

// Need to include simstruc.h for the definition of the SimStruct and
// its associated macro definitions.
#include "simstruc.h"

#define SIG_IDX 0
#define STL_IDX 1
#define MAX_ROB_IDX 2
#define REFRESH_IDX 3

#define SIGNAL_STRING(S) ssGetSFcnParam(S,SIG_IDX)
#define STL_STRING(S)    ssGetSFcnParam(S,STL_IDX)
#define MAX_ROB(S)       ssGetSFcnParam(S,MAX_ROB_IDX)
#define REFRESH_RATE(S)  ssGetSFcnParam(S,REFRESH_IDX)

#define NPARAMS 4

#define NOUTPORT 2
#define UP_IDX 0
#define LOW_IDX 1


#if !defined(MATLAB_MEX_FILE)
/*
 * This file cannot be used directly with the Real-Time Workshop. However,
 * this S-function does work with the Real-Time Workshop via
 * the Target Language Compiler technology. See
 * matlabroot/toolbox/simulink/blocks/tlc_c/sfun_multiport.tlc
 * for the C version
 * matlabroot/toolbox/simulink/blocks/tlc_ada/sfun_multiport.tlc
 * for the Ada version
 */
# error This_file_can_be_used_only_during_simulation_inside_Simulink
#endif

using namespace std;
using namespace CPSGrader;

/*====================*
 * S-function methods *
 *====================*/

#define MDL_CHECK_PARAMETERS
#if defined(MDL_CHECK_PARAMETERS) && defined(MATLAB_MEX_FILE)
/* Function: mdlCheckParameters =============================================
 * Abstract:
 *    Validate our parameters to verify they are okay.
 */
static void mdlCheckParameters(SimStruct *S)
{
    
    if (mxIsEmpty(    ssGetSFcnParam(S,SIG_IDX)) ||
            mxIsSparse(   ssGetSFcnParam(S,SIG_IDX)) ||
            mxIsComplex(  ssGetSFcnParam(S,SIG_IDX)) ||
            mxIsLogical(  ssGetSFcnParam(S,SIG_IDX)) ||
            mxIsNumeric(  ssGetSFcnParam(S,SIG_IDX)) ||
            mxIsDouble(   ssGetSFcnParam(S,SIG_IDX)) ||
            !mxIsChar( ssGetSFcnParam(S,SIG_IDX))) {
        ssSetErrorStatus(S,"Parameter signal names must be a string.");
        return;
    }
    if (mxIsEmpty(    ssGetSFcnParam(S,STL_IDX)) ||
            mxIsSparse(   ssGetSFcnParam(S,STL_IDX)) ||
            mxIsComplex(  ssGetSFcnParam(S,STL_IDX)) ||
            mxIsLogical(  ssGetSFcnParam(S,STL_IDX)) ||
            mxIsNumeric(  ssGetSFcnParam(S,STL_IDX)) ||
            mxIsDouble(   ssGetSFcnParam(S,STL_IDX)) ||
            !mxIsChar( ssGetSFcnParam(S,STL_IDX))) {
        ssSetErrorStatus(S,"Parameter formula must be a string.");
        return;
    }
    
    if (mxIsEmpty(    ssGetSFcnParam(S,MAX_ROB_IDX)) ||
            mxIsSparse(   ssGetSFcnParam(S,MAX_ROB_IDX)) ||
            mxIsComplex(  ssGetSFcnParam(S,MAX_ROB_IDX)) ||
            mxIsLogical(  ssGetSFcnParam(S,MAX_ROB_IDX)) ||
            !mxIsNumeric(  ssGetSFcnParam(S,MAX_ROB_IDX)) ||
            !mxIsDouble(   ssGetSFcnParam(S,MAX_ROB_IDX)) ||
            mxIsChar( ssGetSFcnParam(S,MAX_ROB_IDX))) {
        ssSetErrorStatus(S,"Parameter max. robustness must be a positive double.");
        return;
    }

        if (mxIsEmpty(    ssGetSFcnParam(S,REFRESH_IDX)) ||
            mxIsSparse(   ssGetSFcnParam(S,REFRESH_IDX)) ||
            mxIsComplex(  ssGetSFcnParam(S,REFRESH_IDX)) ||
            mxIsLogical(  ssGetSFcnParam(S,REFRESH_IDX)) ||
            !mxIsNumeric(  ssGetSFcnParam(S,REFRESH_IDX)) ||
            !mxIsDouble(   ssGetSFcnParam(S,REFRESH_IDX)) ||
            mxIsChar( ssGetSFcnParam(S,REFRESH_IDX))) {
        ssSetErrorStatus(S,"Parameter refresh rate must be a positive double.");
        return;
    }

    
    char *signal_buf = mxArrayToString(SIGNAL_STRING(S));
    char *stl_buf = mxArrayToString(STL_STRING(S));
   
    string phi_st = "signal "+ string(signal_buf) + "\n" + "phi:=" + string(stl_buf);
    mxFree(signal_buf);
    mxFree(stl_buf);
    
    STLDriver *stl_driver = new STLDriver();
    bool parse_success = stl_driver->parse_string(phi_st);
    if (parse_success==0)
        ssSetErrorStatus(S,"Error while parsing STL Formula.");
    delete stl_driver;
    
}
#endif /* MDL_CHECK_PARAMETERS */

// Function: mdlInitializeSizes ===============================================
// Abstract:
//    The sizes information is used by Simulink to determine the S-function
//    block's characteristics (number of inputs, outputs, states, etc.).
static void mdlInitializeSizes(SimStruct *S)
{
    // Check the number of expected parameters
    ssSetNumSFcnParams(S, NPARAMS);
    // Parameter mismatch will be reported by Simulink
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return;
    }
  
    ssSetNumDiscStates(S, 2);
    
    mdlCheckParameters(S);
    if (ssGetErrorStatus(S) != NULL) {
        return;
    }
    ssSetSFcnParamTunable(S,STL_IDX,true);
    ssSetSFcnParamTunable(S,SIG_IDX,true);
    ssSetSFcnParamTunable(S,MAX_ROB_IDX,true);
    ssSetSFcnParamTunable(S,REFRESH_IDX,true);

    // Specify I/O
    if (!ssSetNumInputPorts(S, 1)) return;
    if(!ssSetInputPortDimensionInfo(S, 0, DYNAMIC_DIMENSION)) return;
    ssSetInputPortWidth(S, 0, DYNAMICALLY_SIZED);
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    
    if (!ssSetNumOutputPorts(S,1)) return;
    ssSetOutputPortWidth(S, 0, NOUTPORT);
    
    ssSetNumSampleTimes(S, 1);
    
    // Reserve place for C++ object
    ssSetNumPWork(S, 1);
    
    ssSetSimStateCompliance(S, USE_CUSTOM_SIM_STATE);
    
    ssSetOptions(S,
            SS_OPTION_WORKS_WITH_CODE_REUSE |
            SS_OPTION_EXCEPTION_FREE_CODE);
    
}
#if defined(MATLAB_MEX_FILE)
#define MDL_SET_INPUT_PORT_DIMENSION_INFO
static void mdlSetInputPortDimensionInfo(SimStruct *S,
        int_T port, const DimsInfo_T *dimsInfo)
{
    if(!ssSetInputPortDimensionInfo(S, port, dimsInfo)) return;
}
#endif

#if defined(MATLAB_MEX_FILE)
#define MDL_SET_OUTPUT_PORT_DIMENSION_INFO
static void mdlSetOutputPortDimensionInfo(SimStruct *S,
        int_T port, const DimsInfo_T *dimsInfo)
{
    if(!ssSetOutputPortDimensionInfo(S, port, dimsInfo)) return;
}
#endif

#if defined(MATLAB_MEX_FILE)
# define MDL_SET_DEFAULT_PORT_DIMENSION_INFO
/* Function: mdlSetDefaultPortDimensionInfo ===========================================
 * Abstract:
 *   In case no ports were specified, the default is an input port of width 2
 *   and an output port of width 1.
 */
static void mdlSetDefaultPortDimensionInfo(SimStruct *S)
{
    ssSetInputPortWidth(S, 0, 1);
    ssSetOutputPortWidth(S, 0, NOUTPORT);
}
#endif

#define MDL_INITIALIZE_CONDITIONS
/* Function: mdlInitializeConditions ========================================
 * Abstract:
 *    Initialize discrete states to initial min-max robustnesses
 */
static void mdlInitializeConditions(SimStruct *S)
{
    real_T *x0     = ssGetRealDiscStates(S);
    real_T max_rob = mxGetScalar(MAX_ROB(S));

    x0[UP_IDX]  =  max_rob;
    x0[LOW_IDX] = -max_rob;

}


// Function: mdlInitializeSampleTimes =========================================
// Abstract:
//   This function is used to specify the sample time(s) for your
//   S-function. You must register the same number of sample times as
//   specified in ssSetNumSampleTimes.
static void mdlInitializeSampleTimes(SimStruct *S)
{
    //ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    real_T refresh_rate = mxGetScalar(REFRESH_RATE(S));
    ssSetSampleTime(S, 0, refresh_rate);
    ssSetOffsetTime(S, 0, 0.0);
    ssSetModelReferenceSampleTimeDefaultInheritance(S);
    
}

// Function: mdlStart =======================================================
// Abstract:
//   This function is called once at start of model execution. If you
//   have states that should be initialized once, this is the place
//   to do it.
#define MDL_START
static void mdlStart(SimStruct *S)
{
    
    char *signal_buf = mxArrayToString(SIGNAL_STRING(S));
    char *stl_buf = mxArrayToString(STL_STRING(S));   
    string phi_st = "signal "+ string(signal_buf) + "\n" + "phi:=" + string(stl_buf);
    
    mxFree(signal_buf);
    mxFree(stl_buf);
    
    // Store new C++ object in the pointers vector
    STLDriver *stl_driver  = new STLDriver();
    
    stl_driver->parse_string(phi_st);    
    ssGetPWork(S)[0] = stl_driver;

}

#define MDL_UPDATE
static void mdlUpdate(SimStruct *S, int_T tid) {
  
    STLDriver *stl_driver = static_cast<STLDriver *>(ssGetPWork(S)[0]);
    double max_rob = mxGetScalar(MAX_ROB(S));
    Signal::BigM = max_rob;
    
    // Get data addresses of I/O
    InputRealPtrsType  u        =  ssGetInputPortRealSignalPtrs(S,0);
    int_T              u_width  =  ssGetInputPortWidth(S,0);
    real_T           *xd        =  ssGetRealDiscStates(S);
    real_T             T        =  ssGetT(S);


    vector<double> points;
    points.push_back(T);
    for(int_T i = 0; i< u_width; ++i ) {
        points.push_back(*u[i]);
    }
    
    stl_driver->data.push_back(points);
    
    double rob_up = max_rob;
    double rob_low= -max_rob;
    Signal z_up, z_low;
     
    if (stl_driver->data.size()>2) {
        transducer *phi = stl_driver->formula_map["phi"]->clone();
        phi->set_trace_data_ptr(stl_driver->data);
        
        rob_up  = phi->compute_upper_rob();
        rob_low = phi->compute_lower_rob();        
        delete phi;
    }  
   
    if (rob_up==TOP)
        rob_up= max_rob;
    if (rob_low==BOTTOM)
        rob_low= -max_rob;
    
    xd[UP_IDX] =  rob_up;
    xd[LOW_IDX] = rob_low;

}



// Function: mdlOutputs =======================================================
// Abstract:
//   compute the outputs of your S-function block.
static void mdlOutputs(SimStruct *S, int_T tid)
{
 
    real_T            *xd    = ssGetRealDiscStates(S);
    real_T            *y     = ssGetOutputPortRealSignal(S, 0);
       
    y[UP_IDX] =  xd[UP_IDX];
    y[LOW_IDX] = xd[LOW_IDX];
    
}


// Function: mdlTerminate =====================================================
// Abstract:
//   In this function, you should perform any actions that are necessary
//   at the termination of a simulation.  For example, if memory was
//   allocated in mdlStart, this is the place to free it.
static void mdlTerminate(SimStruct *S)
{
    // Retrieve and destroy C++ object
    STLDriver *stl_driver = static_cast<STLDriver *>(ssGetPWork(S)[0]);
    delete stl_driver;
}


// Required S-function trailer
#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
