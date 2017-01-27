/*************************************************************************
 * MATLAB MEX ROUTINE LEMIRE_ND_ENGINE.C
 * [minval maxval] = LEMIRE_ND_ENGINE(a, window)
 *
 * PURPOSE: multiple 1D min/max running/filtering
 * Similar to LEMIRE_ENGINE but working on the second dimension, while
 * looping along the first and third dimensions. This MEX is used for
 * the engine for multidimensional min/max filtering
 *
 * INPUTS
 *  A: 3D arrays, logical and all numeric classes are supported
 *  window: scalar, size of the sliding window, must be >= 1
 * 
 * OUTPUTS
 *  minval, maxval: running min/max, vectors of dimension
 *  (length(A)-window+1), i.e.,
 *      minval(:,1,:) is min(A(:,1:win,:))
 *      minval(:,2,:) is min(A(:,2:win+1,:))
 *      ...
 *      minval(:,end,:) is min(A(:,end-window+1:end,:))
 *  The same indexing arrangement applies for maxval
 *
 * Note: if the data is complex, the imaginary part is ignored.
 *       window is limited to 2147483646 (2^31-2)
 *
 * Algorithm: Lemire's "STREAMING MAXIMUM-MINIMUM FILTER USING NO MORE THAN
 * THREE COMPARISONS PER ELEMENT" Nordic Journal of Computing, Volume 13,
 * Number 4, pages 328-339, 2006.
 *
 * Compilation:
 *  >> mex -O -v lemire_nd_engine.c % add -largeArrayDims on 64-bit computer
 *
 * see aldo: median filter, Kramer & Bruckner filter
 *
 * Author: Bruno Luong <brunoluong@yahoo.com>
 * History
 *  Original: 12/July/2009 
 *  13/July/2009: Correct bug for 64-bit machine (wedge index type)
 ************************************************************************/

#include "mex.h"
#include "matrix.h"

/* Uncomment this on older Matlab version where size_t has not been
 * defined */
/*
 * #define mwSize int
 * #define size_t int
 */

/* Define correct type depending on platform 
  You might have to modify here depending on your compiler */
#if defined(_MSC_VER) || defined(__BORLANDC__)
typedef __int64 int64;
typedef __int32 int32;
typedef __int16 int16;
typedef __int8 int08;
typedef unsigned __int64 uint64;
typedef unsigned __int32 uint32;
typedef unsigned __int16 uint16;
typedef unsigned __int8 uint08;
#else /* LINUX + LCC, CAUTION: not tested by the author */
typedef long long int int64;
typedef long int int32;
typedef short int16;
typedef char int08;
typedef unsigned long long int uint64;
typedef unsigned long int uint32;
typedef unsigned short uint16;
typedef unsigned char uint08;
#endif

/* Maximum (32-bit) integer */
#define MAXINT 0x7fffffff
            
/* This is the engine macro, used for different data type */
/* Note: j -> third dimension
 *       k -> first dimension
 *       i -> second (working) dimension */
#define ENGINE(a, minval, maxval, type) { \
for (j=0; j<q; j++) { \
    a = (type*)mxGetData(A) + j*stepA; \
    minval = (type*)mxGetData(MINVAL) + j*stepMinMax; \
    maxval = (type*)mxGetData(MAXVAL) + j*stepMinMax; \
    for (k=0; k<p; k++) { \
        nU = nL = 0; \
        Lfirst = Ufirst = 0; \
        Llast = Ulast = -1; \
        for (i=1; i<n; i++) { \
            left = (int)(i-window); \
            if (left >= 0) { \
                maxval[p*left] = a[p*(nU? U[Ufirst] : i-1)]; \
                minval[p*left] = a[p*(nL? L[Lfirst] : i-1)]; \
            } \
            if (a[p*i] > a[p*(i-1)]) \
            { \
                nL++; \
                if ((++Llast) == size) Llast = 0; \
                L[Llast] = i-1; \
                if (left == L[Lfirst]) { \
                    nL--; \
                    if ((++Lfirst) == size) Lfirst = 0; \
                } \
                while (nU) { \
                    if (a[p*i] <= a[p*U[Ulast]]) { \
                        if (left == U[Ufirst]) { \
                            nU--; \
                            if ((++Ufirst) == size) Ufirst = 0; \
                        } \
                        break; \
                    } \
                    nU--; \
                    if ((--Ulast) < 0) Ulast += size; \
                } \
            } \
            else { \
                nU++; \
                if ((++Ulast) == size) Ulast = 0; \
                U[Ulast] = i-1; \
                if (left == U[Ufirst]) { \
                    nU--; \
                    if ((++Ufirst) == size) Ufirst = 0; \
                } \
                while (nL) { \
                    if (a[p*i] >= a[p*L[Llast]]) { \
                        if (left == L[Lfirst]) { \
                            nL--; \
                            if ((++Lfirst) == size) Lfirst = 0; \
                        } \
                        break; \
                    } \
                    nL--; \
                    if ((--Llast) < 0) Llast += size; \
                } \
            } \
        } \
        left = (int)(n-window); \
        maxval[p*left] = a[p*(nU? U[Ufirst] : n-1)]; \
        minval[p*left] = a[p*(nL? L[Lfirst] : n-1)]; \
        a++; minval++; maxval++; \
    } \
} }
/* end ENGINE */

/* Define the name for Input/Output ARGUMENTS */
#define A prhs[0]
#define WINDOW prhs[1]
#define MINVAL plhs[0]
#define MAXVAL plhs[1]

/* Gateway of LEMIRE_ND_ENGINE */
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    
    mxClassID ClassID;
    
    /* Data pointers, which one are used depends on the class of A */
    double *adouble, *minvaldouble, *maxvaldouble;
    float *asingle, *minvalsingle, *maxvalsingle;
    int64 *aint64, *minvalint64, *maxvalint64;
    int32 *aint32, *minvalint32, *maxvalint32;
    int16 *aint16, *minvalint16, *maxvalint16;
    int08 *aint08, *minvalint08, *maxvalint08;
    uint64 *auint64, *minvaluint64, *maxvaluint64;
    uint32 *auint32, *minvaluint32, *maxvaluint32;
    uint16 *auint16, *minvaluint16, *maxvaluint16;
    uint08 *auint08, *minvaluint08, *maxvaluint08;

    mwSize i, n, window;
    int left, size;
    mwSize *U, *L; /* wedge */
    int nU, nL; /* wedge number of elements (0 is empty wedge) */
    int Ufirst, Lfirst, Ulast, Llast; /* Indices of two ends of the wedge */
    
    mwSize p, q, j, k;
    mwIndex stepA, stepMinMax;
    const mwSize* dimA;
    mwSize dimOut[3];
    
    /* Check number of arguments */
    if (nrhs!=2)
        mexErrMsgTxt("LEMIRE_ND_ENGINE: two arguments are required.");
    
    /* Get class of input matrix A */
    ClassID = mxGetClassID(A);
    
    /* Do not support on sparse */
    if (mxIsSparse(A))
        mexErrMsgTxt("LEMIRE_ND_ENGINE: First input A must be full.");       
    
    /* Get the number of elements of A */
    /* Get the size, MUST BE two or three, no check */
    dimA = mxGetDimensions(A);
    p = dimA[0];
    n = dimA[1];
    if (mxGetNumberOfDimensions(A)<3)
        q = 1; /* third dimension is singleton */
    else
        q = dimA[2];
   
    /* Window input must be double */
    if (mxGetClassID(WINDOW)!=mxDOUBLE_CLASS)
        mexErrMsgTxt("LEMIRE_ND_ENGINE: Second input WINDOW must be double.");
    
    /* Get the window size, cast it in mwSize */
    window = (mwSize)(*mxGetPr(WINDOW));
    
    if (window<1) /* Check if it's valid */
        mexErrMsgTxt("LEMIRE_ND_ENGINE: windows must be 1 or greater.");   
    if (window>n || window>MAXINT)
        mexErrMsgTxt("LEMIRE_ND_ENGINE: windows larger than data length.");
    
    /* Allocate wedges buffers for L and U, each is size (window+1) */
    size = (int)(window+1);
    L = mxMalloc((2*size)*sizeof(mwSize));
    if (L==NULL) mexErrMsgTxt("LEMIRE_ND_ENGINE: out of memory.");
    U = L + size;
    /* Initialize empty wedges L and U */
  
    /* Create output arrays */
    dimOut[0] = p; dimOut[1] = n-window+1; dimOut[2] = q;
    MINVAL = mxCreateNumericArray(3, dimOut, ClassID, mxREAL);
    MAXVAL = mxCreateNumericArray(3, dimOut, ClassID, mxREAL);
    
    /* Check if allocation is succeeded */
    if (MINVAL==NULL || MAXVAL==NULL)
         mexErrMsgTxt("LEMIRE_ND_ENGINE: out of memory.");
    
    /* Jump step of the third dimension */
    stepA = p*n; /* for A */
    stepMinMax = p*dimOut[1]; /* step for output */
    
    /* Call the engine depending on ClassID */
    switch (ClassID) {
        case mxDOUBLE_CLASS:
             ENGINE(adouble, minvaldouble, maxvaldouble, double);
             break;
        case mxSINGLE_CLASS:
             ENGINE(asingle, minvalsingle, maxvalsingle, float);
             break;
        case mxINT64_CLASS:
             ENGINE(aint64, minvalint64, maxvalint64, int64);
             break;
        case mxUINT64_CLASS:
             ENGINE(auint64, minvaluint64, maxvaluint64, uint64);
             break;
        case mxINT32_CLASS:
             ENGINE(aint32, minvalint32, maxvalint32, int32);
             break;
        case mxUINT32_CLASS:
             ENGINE(auint32, minvaluint32, maxvaluint32, uint32);
             break;
        case mxCHAR_CLASS:
             ENGINE(auint16, minvaluint16, maxvaluint16, uint16);
             break;
        case mxINT16_CLASS:
             ENGINE(aint16, minvalint16, maxvalint16, int16);
             break;
        case mxUINT16_CLASS:
             ENGINE(auint16, minvaluint16, maxvaluint16, uint16);
             break;
        case mxLOGICAL_CLASS:
             ENGINE(auint08, minvaluint08, maxvaluint08, uint08);
             break;
        case mxINT8_CLASS:
             ENGINE(aint08, minvalint08, maxvalint08, int08);
             break;
        case mxUINT8_CLASS:
             ENGINE(auint08, minvaluint08, maxvaluint08, uint08);
             break;
        default:
            mexErrMsgTxt("LEMIRE_ND_ENGINE: Class not supported.");
    } /* switch */
    
    /* Free the buffer */
    mxFree(L);
    
    return;
    
} /* Gateway LEMIRE_ND_ENGINE */
