#ifndef PARAMSET_H
#define PARAMSET_H
#include "mex.h"
#include "string.h"
#include "matrix.h"
#include "math.h"
#include <vector>
#include <list>
#include <iostream>

using namespace std;

class r_node {

 public:
  
  int index;
  bool isleaf;
  int is_divided_in;
  r_node* *children;
  r_node* father;
  double val;
  int status;
  int NumChildren;
   
  /*constructors*/

  r_node();
  r_node(int x,int NumOfChildren)//constructor for root
  {
	index = x;
	isleaf = true;
	father = NULL;
	is_divided_in = -1;
	NumChildren = NumOfChildren;
	children = (r_node**)mxMalloc(NumOfChildren*sizeof(r_node*));
	for (int i = 0; i<NumOfChildren; i++)
	  children[i]=NULL;
   };
	
  r_node(r_node * x,int NumOfChildren) 
  {
	index = x->index;
	father = x;
	val = x->val;
	isleaf = true;
	is_divided_in = -1;
	NumChildren = NumOfChildren;
	children = (r_node**)mxMalloc(NumOfChildren*sizeof(r_node*));
	for (int i = 0; i<NumOfChildren; i++)
	  children[i]=NULL;

  };
  
  r_node(int a, r_node * x,double value,int NumOfChildren) 
  {
	index = a;
	father = x;
	val = value;
	isleaf = true;
	is_divided_in = -1;
	NumChildren = NumOfChildren;
	children = (r_node**)mxMalloc(NumOfChildren*sizeof(r_node*));
	for (int i = 0; i<NumOfChildren; i++)
	  children[i]=NULL;	
  };

};

typedef r_node* refine_tree;

class param_set {

 public:

  double *pts;     
  double *epsi;
  double *dim;
  int M_pts;
  int M_epsi;
  int NumElems;
  int NumRoots;
  vector<refine_tree> queue; // this maintains nodes that have to be tested 
  
  refine_tree * refine_history;
    
  /* constructors */ 

  param_set() {
    M_pts = M_epsi= NumElems = 1 ;
  };
  
  // create a param_set from a matlab structure 
  param_set(const mxArray * mxInput, int NumOfChildren)
  {
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
    NumRoots = NumElems;
    for(int tree_ind=0;tree_ind<NumElems;tree_ind++) // this will create the roots of the trees
      {
	refine_tree root = new r_node(tree_ind,NumOfChildren);
	refine_history[tree_ind] = root;
	queue.push_back(root);	
      }
  } 
  
  /* methods */

  int Test1(refine_tree,double[]);
  //int Test2(int,double[]);
  bool refine_all(int NumOfChildren, int (param_set::*ptr_to_test)(refine_tree, double[]),double[]); 
  void write_mxArray(mxArray * &mxOutput); 
  //refine_tree pop();
  double* Get_pts(int);
  double* Get_ep(int);


};

void CVM_Refine( int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] );
#endif


