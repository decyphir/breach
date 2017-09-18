#include "mex.h"
#include "signal.h"

void writeSignal(const Signal &y, mxArray * mx_times, mxArray* mx_values)  {
#ifdef DEBUG__
  mexPrintf("Entering wroteSignal\n");
#endif

  double * times = mxGetPr(mx_times);
  double * values = mxGetPr(mx_values);
  
  int k = 0;
  Signal::const_iterator i;  

  for(i = y.begin(); i != y.end(); i++) {
    times[k]= (*i).time; 
    values[k]=(*i).value; 
    k++;
  }
 
  i = y.end()-1;
  times[k] = y.endTime;
  values[k] = (*i).valueAt(y.endTime);

#ifdef DEBUG__
  mexPrintf("Leaving wroteSignal\n");
#endif 
}
