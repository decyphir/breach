#include "mextools.h"

/**********************************************************************/
/* Blitz auxiliary function */

void BlitzArray_SetData(Array1D & A, int n, double * data) {

  Array1D Aux(data,shape(n),neverDeleteData,FortranArray<1>());
  A.resize(Aux.shape());
  A = Aux;

}

void BlitzArray_SetData(Array2D & A, int n, int m, double * data) {

  Array2D Aux(data,shape(n,m),neverDeleteData,FortranArray<2>());
  A.resize(Aux.shape());
  A = Aux;

}

void GetArray1D(Array1D & A, const mxArray * mxA){

  double * A_data;
  int n = mxGetNumberOfElements(mxA);
  A_data = mxGetPr(mxA);
  BlitzArray_SetData(A,n,A_data);

}


void GetArray2D(Array2D & A, const mxArray * mxA) {

  double * A_data;
  int n,m;

  n = mxGetM(mxA);
  m = mxGetN(mxA);

  A_data = mxGetPr(mxA);
  BlitzArray_SetData(A,n,m,A_data);
}

/***********************************************************************
/* Mex interface functions to extract fields in a Matlab structure    */

/* Read Functions */

char * GetString(const mxArray *input) {
  int buflen = (mxGetM(input) * mxGetN(input) * sizeof(mxChar)) + 1;
  char * res;
  char * err_msg;
  res= (char*) mxCalloc(buflen, sizeof(char));
  if (mxGetString(input,res,buflen)) {
    sprintf(err_msg,"Can't read String");
    mexErrMsgTxt(err_msg);
  }   
  return res;
}    

int GetIntField(const mxArray *PB, const char * FieldName) {
 mxArray *field = mxGetField(PB,0,FieldName);

 if (field == NULL) 
    return 0;

 int res = (int) mxGetScalar(field);
 return res;
}

int GetStringField(char* &str, const mxArray *PB,const char * FieldName) {

  mxArray *field = mxGetField(PB,0,FieldName);
 
  if (field == NULL) 
    return 0;
  else {

    int buflen = (mxGetM(field) * mxGetN(field) * sizeof(mxChar)) + 1;
    char * err_msg;
    str= (char*) mxCalloc(buflen, sizeof(char));
    if (mxGetString(field,str,buflen)) {
      sprintf(err_msg,"Can't read String field %s",FieldName);
      mexErrMsgTxt(err_msg);
    }   
    return 1;
  }
}    

double GetDoubleField(const mxArray *PB, const char * FieldName) {
  return mxGetScalar(mxGetField(PB,0,FieldName));
}

double * GetPrField(const mxArray *PB, const char * FieldName) {
  return mxGetPr(mxGetField(PB,0,FieldName));
}

/* Transforme un champ matrice d'une structure en Array1D (1 succes, 0 echec) */

int GetArrayField(Array1D &A, const mxArray *PB, const char * FieldName) {

  mxArray * mxA = mxGetField(PB,0,FieldName);
  
  double * A_data;
  int n;

  if (mxA == NULL) {
    return 0;
  }
  else {
    n = mxGetNumberOfElements(mxA);
    A_data = mxGetPr(mxA);
    BlitzArray_SetData(A,n,A_data);
    return 1;
  }
}

int GetArrayField(Array1D &A, const mxArray *PB, const char * FieldName, int idx) {

  mxArray * mxA = mxGetField(PB, idx,FieldName);
  
  double * A_data;
  int n;

  if (mxA == NULL) {
    return 0;
  }
  else {
    n = mxGetNumberOfElements(mxA);
    A_data = mxGetPr(mxA);
    BlitzArray_SetData(A,n,A_data);
    return 1;
  }
}

int GetArrayField(vector<int> &A, const mxArray *PB, const char * FieldName, int idx) {

  mxArray * mxA = mxGetField(PB,idx,FieldName);  
  double * A_data;
  int n;

  A.clear();
  if (mxA == NULL) {
    return 0;
  }
  else {
    n = mxGetNumberOfElements(mxA);
    A_data = mxGetPr(mxA);
    for (int i=0; i<n; i++) {
      A.push_back((int) A_data[i]);
    }


    return 1;
  }
}



/* La même pour les Array2D */

int GetArrayField(Array2D &A, const mxArray *PB, const char * FieldName) {

  mxArray * mxA = mxGetField(PB,0,FieldName);
  
  double * A_data;
  int n,m;

  if (mxA == NULL) {
    return 0;
  }

  else {
    n = mxGetM(mxA);
    m = mxGetN(mxA);
    A_data = mxGetPr(mxA);
    BlitzArray_SetData(A,n,m,A_data);
    return 1;
  }
}

int GetArrayField(Array2D &A, const mxArray *PB, const char * FieldName, int idx) {

  mxArray * mxA = mxGetField(PB,idx,FieldName);
  
  double * A_data;
  int n,m;

  if (mxA == NULL) {
    return 0;
  }

  else {
    n = mxGetM(mxA);
    m = mxGetN(mxA);
    A_data = mxGetPr(mxA);
    BlitzArray_SetData(A,n,m,A_data);
    return 1;
  }
}

int GetNField(const mxArray *PB, const char * FieldName) {
  return mxGetN(mxGetField(PB,0,FieldName));
}

/* Write Functions */

void SetArray(Array1D blitz_values, mxArray *&VF) {

  int n = blitz_values.extent(0);
  VF = mxCreateDoubleMatrix(1,n,mxREAL);
  double * values = mxGetPr(VF);
  memcpy(values,blitz_values.data(),n*sizeof(double));
}

void SetArray(Array2D blitz_values, mxArray *&VF) {

  int n = blitz_values.extent(0);
  int m = blitz_values.extent(1);
  Array2D matlab_ordered_blitz_values(n,m,FortranArray<2>());
  matlab_ordered_blitz_values = blitz_values;
  VF = mxCreateDoubleMatrix(n,m,mxREAL);
  double * values = mxGetPr(VF);
  memcpy(values,matlab_ordered_blitz_values.data(),n*m*sizeof(double));
 
}

void SetArray(Array<int,1> blitz_values, mxArray *&VF) {

  int n = blitz_values.extent(0);
  Array1D casted_blitz_values(n);
  int i;
  for(i=0;i<n;i++)
    casted_blitz_values(i) = (double) blitz_values(i);
  
  VF = mxCreateDoubleMatrix(1,n,mxREAL);
  double * values = mxGetPr(VF);
  memcpy(values,casted_blitz_values.data(),n*sizeof(double));
  
}



void SetArray(Array<int,2> blitz_values, mxArray *&VF) {

  int n = blitz_values.extent(0);
  int m = blitz_values.extent(1);
  Array2D matlab_ordered_blitz_values(n,m,FortranArray<2>());
  int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      matlab_ordered_blitz_values(i+1,j+1) = (double) blitz_values(i,j);
  
  VF = mxCreateDoubleMatrix(n,m,mxREAL);
  double * values = mxGetPr(VF);
  memcpy(values,matlab_ordered_blitz_values.data(),n*m*sizeof(double));
 
}

void SetField(Array1D blitz_values, mxArray *VF, const char * field_name) {

  int n = blitz_values.extent(0);
  mxArray * mxvalues = mxCreateDoubleMatrix(1,n,mxREAL);
  double * values = mxGetPr(mxvalues);
  memcpy(values,blitz_values.data(),n*sizeof(double));
  mxAddField(VF,field_name);
  int no_field = mxGetFieldNumber(VF,field_name);
  mxSetFieldByNumber(VF,0,no_field,mxvalues);
 
}

void SetField(Array2D blitz_values, mxArray *VF, const char * field_name) {

  int n = blitz_values.extent(0);
  int m = blitz_values.extent(1);

  Array2D matlab_ordered_blitz_values(n,m,FortranArray<2>());
  matlab_ordered_blitz_values = blitz_values;

  mxArray * mxvalues = mxCreateDoubleMatrix(n,m,mxREAL);

  double * values = mxGetPr(mxvalues);

  memcpy(values,matlab_ordered_blitz_values.data(),n*m*sizeof(double));

  mxAddField(VF,field_name);
  int no_field = mxGetFieldNumber(VF,field_name);
  mxSetFieldByNumber(VF,0,no_field,mxvalues);
 }

void SetField(Array<double,3> blitz_values, mxArray *VF, const char * field_name) {

  int n0 = blitz_values.extent(0);
  int n1 = blitz_values.extent(1);
  int n2 = blitz_values.extent(2);
  
  Array<double,3> matlab_ordered_blitz_values(n0,n1,n2,FortranArray<3>());
  matlab_ordered_blitz_values = blitz_values;

  int ndim =3; 
  int * dims = (int *) mxMalloc(ndim*sizeof(int));    

  dims[0] = n0;
  dims[1] = n1;
  dims[2] = n2;

  mxArray *mxvalues =mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);         
  double * values = mxGetPr(mxvalues);
  memcpy(values,matlab_ordered_blitz_values.data(),n0*n1*n2*sizeof(double));

  mxAddField(VF,field_name);
  int no_field = mxGetFieldNumber(VF,field_name);
  mxSetFieldByNumber(VF,0,no_field,mxvalues);

}

void SetField(Array<double,4> blitz_values, mxArray *VF, const char * field_name) {

  int n0 = blitz_values.extent(0);
  int n1 = blitz_values.extent(1);
  int n2 = blitz_values.extent(2);
  int n3 = blitz_values.extent(2);
  
  Array<double,4> matlab_ordered_blitz_values(n0,n1,n2,n3,FortranArray<4>());
  matlab_ordered_blitz_values = blitz_values;

  int ndim =4; 
int *  dims = (int *) mxMalloc(ndim*sizeof(int));    

  dims[0] = n0;
  dims[1] = n1;
  dims[2] = n2;
  dims[3] = n3;

  mxArray *mxvalues =mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);         
  double * values = mxGetPr(mxvalues);
  memcpy(values,matlab_ordered_blitz_values.data(),n0*n1*n2*sizeof(double));

  mxAddField(VF,field_name);
  int no_field = mxGetFieldNumber(VF,field_name);
  mxSetFieldByNumber(VF,0,no_field,mxvalues);
}

void CreateStruct(Array1D blitz_values, mxArray *& VF, const char * field_name){
		  
  int n = blitz_values.extent(0);
  mxArray * mxvalues = mxCreateDoubleMatrix(1,n,mxREAL);
  double * values = mxGetPr(mxvalues);
  memcpy(values,blitz_values.data(),n*sizeof(double));

  int dims_struct[2] = {1,1};
  const char * field_names[] = {field_name};
  int V_field;

  VF = mxCreateStructArray(2,dims_struct,1, field_names);
  V_field = mxGetFieldNumber(VF,field_name);
  mxSetFieldByNumber(VF,0,V_field,mxvalues);  

}

void CreateStruct(Array2D blitz_values, mxArray *& VF, const char * field_name){

  int n = blitz_values.extent(0);
  int m = blitz_values.extent(1);
  Array2D matlab_ordered_blitz_values(n,m,FortranArray<2>());
  matlab_ordered_blitz_values = blitz_values;
  mxArray * mxvalues = mxCreateDoubleMatrix(n,m,mxREAL);
  double * values = mxGetPr(mxvalues);
  memcpy(values,matlab_ordered_blitz_values.data(),n*m*sizeof(double));

  int dims_struct[2] = {1,1};
  const char * field_names[] = {field_name};
  int V_field;
  VF = mxCreateStructArray(2,dims_struct,1, field_names);
  V_field = mxGetFieldNumber(VF,field_name);
  mxSetFieldByNumber(VF,0,V_field,mxvalues);  

}

