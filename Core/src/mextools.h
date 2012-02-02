#ifndef MEXTOOLS_H
#define MEXTOOLS_H
#include <iostream>
#include <math.h>
#include <algorithm>
#include <time.h>
#include "mex.h"
#include <blitz/array.h>
#include <cstring>
#include <vector>

/* Misc */

using namespace blitz;

typedef Array<double,1> Array1D;
typedef Array<double,2> Array2D;

void BlitzArray_SetData(Array1D & A, int n, double * data);
void BlitzArray_SetData(Array2D & A, int n, int m, double * data);

/* Read */

void GetArray1D(Array1D & A, const mxArray * mxA);
void GetArray2D(Array2D & A, const mxArray * mxA);
char * GetString(const mxArray *input);
int GetIntField(const mxArray *PB,const char * FieldName);
int GetStringField(char * &str,const mxArray *PB,const char * FieldName);
double GetDoubleField(const mxArray *PB,const char * FieldName);
double * GetPrField(const mxArray *PB,const  char * FieldName); 
int GetArrayField(Array1D &A, const mxArray *PB, const char * FieldName);
int GetArrayField(Array2D &A, const mxArray *PB,const  char * FieldName);

int GetArrayField(Array1D &A, const mxArray *PB, const char * FieldName, int idx);
int GetArrayField(vector<int> &A, const mxArray *PB, const char * FieldName, int idx);
int GetArrayField(Array2D &A, const mxArray *PB, const char * FieldName, int idx);

int GetNField(const mxArray *PB, const char * FieldName);

/* Write */

void SetArray(Array1D blitz_values, mxArray *&VF);
void SetArray(Array2D blitz_values, mxArray *&VF);
void SetArray(Array<int,1> blitz_values, mxArray *&VF);
void SetArray(Array<int,2> blitz_values, mxArray *&VF);

void CreateStruct(Array1D blitz_values, mxArray *& VF, const char * field_name);
void CreateStruct(Array2D blitz_values, mxArray *& VF, const char * field_name);

void SetField(Array1D blitz_values, mxArray *VF, const char * field_name);
void SetField(Array2D blitz_values, mxArray *VF, const char * field_name);
void SetField(Array<double,3> blitz_values, mxArray *VF, const char * field_name);
void SetField(Array<double,4> blitz_values, mxArray *VF, const char * field_name);

#endif
