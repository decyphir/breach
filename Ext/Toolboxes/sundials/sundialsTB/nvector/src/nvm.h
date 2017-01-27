/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2009-06-05 16:26:21 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Header file for the MNVECTOR Matlab interface.
 * -----------------------------------------------------------------
 */

#ifndef _NVM_H
#define _NVM_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <sundials/sundials_nvector.h>
#include "mex.h"

  /*
   * ------------------
   * Vector functions
   * ------------------
   */
  
  void InitVectors();
  N_Vector NewVector(int n);

  void PutData(N_Vector v, double *data, int n);
  void GetData(N_Vector v, double *data, int n);

#ifdef __cplusplus
}
#endif

#endif
