#include "mex.h"
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <ginac/ginac.h>
using namespace std;
using namespace GiNaC;


void mexFunction(int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ) 
{

  string fn(mxArrayToString(prhs[0]));
  
  char* param_st;
  ostringstream dfn_st;

  param_st = mxArrayToString(prhs[1]);

  parser reader;
  try {

    ex e = reader(fn);
    symtab table = reader.get_syms();
    symbol param;
    
    if (table.find(param_st) != table.end()) {
      param=ex_to<symbol>(table[param_st]);
      dfn_st << e.diff(param);
    }
    else 
      dfn_st << "0";
  } catch (exception &p) {
    cerr << p.what() << endl;
  }

  const char * st_out = dfn_st.str().c_str();
  plhs[0]=  mxCreateString(st_out);
  
  mxFree(param_st);

}
