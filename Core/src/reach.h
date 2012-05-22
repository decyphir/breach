#ifndef REACH_H
#define REACH_H

#include "traj.h"

class traj_set: public param_set {

  public:

  Array1D tspan; 
  vector<trajectory *> trajectories; 
  
  traj_set();
  
  traj_set(const mxArray * mxInput, int NumOfChildren) {
    const char **fnames;		// pointers to field names
    int ifield,nfields;
    mxArray *mx_epsi_input,*mx_pts_input,*mx_dim_input;
    
    nfields = mxGetNumberOfFields(mxInput);
    
    // allocate memory  for storing pointers
    fnames = (const char **)mxCalloc(nfields, sizeof(*fnames));
    
    // get field name pointers
    for (ifield=0; ifield< nfields; ifield++)
      {
	fnames[ifield] = mxGetFieldNameByNumber(mxInput,ifield);
	if(!strcmp(fnames[ifield],"pts"))
	  mx_pts_input = mxGetField(mxInput,0,fnames[ifield]);
	else if(!strcmp(fnames[ifield],"epsi"))
	  mx_epsi_input = mxGetField(mxInput,0,fnames[ifield]);
	else if(!strcmp(fnames[ifield],"dim"))
	  mx_dim_input = mxGetField(mxInput,0,fnames[ifield]);
      }
    
    /*assign pointers to input*/

    double * temp_pts = mxGetPr(mx_pts_input);		//data_pts_input is the pointer to input pts array
    M_pts = mxGetM(mx_pts_input);
    NumElems = mxGetN(mx_pts_input);
    
    double * temp_epsi = mxGetPr(mx_epsi_input);		// data_epsi_input is the pointer to inputepsi array
    M_epsi = mxGetM(mx_epsi_input);
    
    double * temp_dim = mxGetPr(mx_dim_input);		// data_dim_input is the pointer to input dim array
    
    pts = (double*)mxCalloc(M_pts*NumElems,sizeof(double));
    epsi = (double*)mxCalloc(M_epsi*NumElems,sizeof(double));
    dim = (double*)mxCalloc(M_epsi*1,sizeof(double));
    trajectories.resize(NumElems);
    NumRoots=NumElems;

    for(int i=0;i<NumElems;i++)
      {
	for(int j=0;j<M_pts;j++)
	  pts[i*M_pts+j] = temp_pts[i*M_pts+j]; 
	
	for(int j=0;j<M_epsi;j++)
	  epsi[i*M_epsi+j] = temp_epsi[i*M_epsi+j]; 
      }

    for(int j=0;j<M_epsi;j++)
		dim[j] = temp_dim[j];
    
    refine_history = (refine_tree *) mxMalloc(NumElems*sizeof(refine_tree));
    queue.resize(NumElems);

    for(int tree_ind=0;tree_ind<NumElems;tree_ind++) // this will create the roots of the trees
      {
	refine_tree root = new r_node(tree_ind,NumOfChildren);
	refine_history[tree_ind] = root;
	queue.at(tree_ind)=root;	
      }

  }; 

  bool refine_traj(int NumOfChildren, int (traj_set::*ptr_to_test)(refine_tree&, double[]),double tst_eps[]);
  int test_expansion(refine_tree&, double[]); 
  int test_nonlin_error(refine_tree&, double[]);
  void print_refine_history();

};


int CVM_Reach(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
int CVM_ComputeTrajs(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
int CVM_CmacComputeTraj(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

void TrajErr(const trajectory& global_traj, const Array1D& local_p, const Array1D& local_epsi, trajectory & local_traj, Array1D & Err);

void ComputeExpTraj(CmacVF *&C, Array2D pts, Array2D XS0, Array2D epsi, Array1D tspan, Array<trajectory*,1> &trajArray, Array2D &Xf, Array2D &XSf, Array2D &ExpaMax, Array1D &Keep );

void ComputeTrajs(CmacVF *&C, Array2D pts, Array1D tspan, Array<trajectory*,1> &trajArray,Array2D &Xf, Array1D &Keep);

#endif




