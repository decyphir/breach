#include "traj.h"

/*  Some Constants for DQ approximation of g*/

#define MIN_INC_MULT RCONST(1000.0)
#define ONE          RCONST(1.0)
#define ZERO         RCONST(0.0)
#define TWO          RCONST(2.0)

/*********** Subgate mex functions **********************/

int CVM_GetF(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
#if _DEBUG >= 1
  cout << "Entering CVM_GetF ..." << endl;
#endif

  /* Reading inputs */

  double t = mxGetScalar(prhs[0]);
  Array2D pts(N,1);
  GetArray2D(pts, prhs[1]);
  int nb_pts = pts.extent(1);
  
  cout << "N:" << N <<endl;

  /* Doing some stuff */

  Array2D fValues(N,nb_pts);
  InitFdata(f_data, cvm_Mdata->mx_data);
  
  Fdata* data = (Fdata*) f_data;
  N_Vector xv = N_VNew_Serial(N);
  N_Vector xdotv = N_VNew_Serial(N);

  realtype *xdata = N_VGetArrayPointer(xv);  
  realtype *xdotdata = N_VGetArrayPointer(xdotv);

  for (int j=0;j<nb_pts;j++) {
    
    for (int k=0; k<N; k++) {
      xdata[k] = pts(k,j);
      data->p[k] = xdata[k];
    }

    for (int k=N; k<data->dimp; k++)
      data->p[k] = pts(k,j);

    //    x = pts(Rx,j);
    //    p = pts(All,j);
    
    f(t,xv,xdotv,data);   

    for (int k=0; k<N; k++) {
      fValues(k,j) = xdotdata[k];
    }

    //    fValues(All,j)= xdot;
    
  }

  /* Writing outputs */ 
  
  SetArray(fValues,plhs[0]);

#if _DEBUG >= 1
  cout << "Leaving CVM_GetF ..." << endl;
#endif
  return 0;
}

int CVM_GetJf(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
#if _DEBUG >= 1
  cout << "Entering CVM_GetF ..." << endl;
#endif

  /* Reading inputs */

  double t = mxGetScalar(prhs[0]);
  Array2D pts(N,1);
  GetArray2D(pts, prhs[1]);
  int nb_pts = pts.extent(1);
  
  /* Doing some stuff */

  Fdata* data = (Fdata*) f_data;
  N_Vector xv = N_VNew_Serial(N);
  DenseMat Jv = DenseAllocMat(N,N);

  realtype *xdata = N_VGetArrayPointer(xv);  
  realtype **Jdata = Jv->data;
  
  Array1D x(xdata,shape(N),neverDeleteData);
  Array1D p(data->p,shape(data->dimp),neverDeleteData);
  Array2D J(*Jdata, shape(N,N), neverDeleteData);
  
  Range Rx(0,N-1);
  Range All = Range::all();   

  InitFdata(f_data, cvm_Mdata->mx_data);

  if (nb_pts>=1) {

    for (int j=0;j<nb_pts;j++) {
      J = 0;
      x = pts(Rx,j);
      p = pts(All,j);
      Jac(N,Jv,t,xv,NULL,f_data, NULL,NULL,NULL);             
    }

  }

  else {
    J = 0;
    x = pts(Rx,0);
    p = pts(All,0);
    Jac(N,Jv,t,xv,NULL,f_data, NULL,NULL,NULL);         
    SetArray(J,plhs[0]);
  }

  /* Writing outputs */ 
  

#if _DEBUG >= 1
  cout << "Leaving CVM_GetF ..." << endl;
#endif
  return 0;

}


int CVM_ComputeExpansion(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

#if _DEBUG >= 1
  cout << "Entering CVM_ComputeExpansion ..." << endl;
#endif

  /* Reading Inputs */

  Array2D pts(N,1);
  GetArrayField(pts, prhs[0],"pts");

  Array2D XS0(N*Ns,1);
  GetArrayField(XS0, prhs[0],"XS0");

  Array2D epsi(N,1);
  GetArrayField(epsi, prhs[0],"epsi");

  Array1D tspan(2);
  GetArray1D(tspan,prhs[1]);

  /* Doing some stuff */
  
  Array2D Xf(1,1);
  Array2D XSf(1,1);
  Array2D ExpaMax(1,1);

  ComputeExpansion(pts,XS0,epsi,tspan, Xf, XSf,ExpaMax);
 
  /* Writing outputs */ 

  plhs[0] = mxDuplicateArray(prhs[0]);
      
  SetField(Xf,plhs[0],"Xf");

  mxAddField(plhs[0],"XSf");
  SetField(XSf,plhs[0],"XSf"); 

  mxAddField(plhs[0],"ExpaMax");
  SetField(ExpaMax,plhs[0],"ExpaMax");


#if _DEBUG >= 1
   cout << "Leaving CVM_ComputeTraj ..." << endl;
#endif

   return 0;
}

 

int CVM_ComputeExpTraj(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

#if _DEBUG >= 1
  cout << "Entering CVM_ComputeExpTraj ..." << endl;
#endif
  
  /* Reading */

  Array2D pts(N,1);
  GetArrayField(pts, prhs[0],"pts");

  Array2D XS0(N*Ns,1);
  GetArrayField(XS0, prhs[0],"XS0");

  Array2D epsi(N,1);
  GetArrayField(epsi, prhs[0],"epsi");

  Array1D tspan(2);
  GetArray1D(tspan,prhs[1]);

  /* Doing some stuff */

  Array<trajectory*,1> trajArray(1);
  Array2D Xf(1,1);
  Array2D XSf(1,1);
  Array2D ExpaMax(1,1);

  ComputeExpTraj(pts, XS0, epsi, tspan, trajArray, Xf, XSf, ExpaMax);
 
  /* Writing outputs */ 

  plhs[0] = mxDuplicateArray(prhs[0]);
  mxArray * mxTraj;
  
  int nb_pts = trajArray.extent(0);  
  const char * fieldnames[4] = {"time", "X", "Expa", "U"};

  mxTraj = mxCreateStructMatrix(1,nb_pts,4, fieldnames);
  mxArray * mxX, *mxExpa, *mxU, *mxtime;

  for (int j = 0; j<nb_pts ; j++) {

    SetArray(*trajArray(j)->X,mxX);
    SetArray(*trajArray(j)->Expa,mxExpa);
    SetArray(*trajArray(j)->U,mxU);
    SetArray(*trajArray(j)->time,mxtime);

    mxSetField(mxTraj, j, "time",mxtime);
    mxSetField(mxTraj, j, "X",mxX);
    mxSetField(mxTraj, j, "Expa",mxExpa);
    mxSetField(mxTraj, j, "U",mxU);

  }
  
  SetField(Xf,plhs[0],"Xf");
  
  mxAddField(plhs[0],"ExpaMax");
  SetField(ExpaMax,plhs[0],"ExpaMax");

  mxAddField(plhs[0],"XSf");
  SetField(XSf,plhs[0],"XSf");

  mxAddField(plhs[0],"traj");
  mxSetField(plhs[0],0,"traj",mxTraj);

  /* Free memory */ 

  for (int j = 0; j<nb_pts ; j++) {
    delete trajArray(j);
  }

#if _DEBUG >= 1
  cout << "Leaving CVM_ComputeExpTraj ..." << endl;
#endif

  return 0;
}

int CVM_ComputeTrajSensi(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

#if _DEBUG >= 1
  cout << "Entering CVM_ComputeExpTraj ..." << endl;
#endif
  
  /* Reading */

  Array2D pts(N,1);
  GetArrayField(pts, prhs[0],"pts");

  /* set XS0 (same for all trajectories) */
  Array1D XS0(N*Ns);
  GetArrayField(XS0, prhs[0],"XS0");
  
  int i;
  for (int is=0; is<Ns; is++)
    for (i = 0; i<N; i++) 
      NV_Ith_S(yS0[is],i) = XS0(is*N+i);

  Array2D epsi(N,1);
  GetArrayField(epsi, prhs[0],"epsi");

  Array1D sensis(Ns);
  GetArrayField(sensis, prhs[0],"dim");
  
  Array1D tspan(2);
  GetArray1D(tspan,prhs[1]);

  /* Doing some stuff */

  Range All = Range::all();   
  int nb_pts = pts.extent(1);

  Array<trajectory*,1> trajArray(nb_pts);

  Array2D Xf(N,nb_pts);
  Array2D XSf(N*Ns,nb_pts);
  Array2D ExpaMax(N,nb_pts);
  Array1D jthExpaMax(N);  
  Array1D p0(pts.extent(0));
  Array1D Tspan;
  Array1D xS0(N*Ns);

  for (int j = 0; j<nb_pts ; j++) {

    p0 = pts(All, j);
    Tspan.resize(tspan.shape());  
    Tspan = tspan;
    trajArray(j) = new trajectory();
    trajArray(j)->p0 = &p0;
    trajArray(j)->ComputeTrajSensi(Tspan);
    trajArray(j)->ComputeExpa(epsi(Range::all(),j),jthExpaMax);
    ExpaMax(Range::all(), j)= jthExpaMax;
    int l = trajArray(j)->length; 
    Xf(All,j) = (*trajArray(j)->X)(All,l-1);
    XSf(All,j) = (*trajArray(j)->XS)(All,l-1);

  }
  
  /* Writing outputs */ 
 
  plhs[0] = mxDuplicateArray(prhs[0]);

  mxArray * mxTraj;
  
  const char * fieldnames[7] = {"param","sensis", "time", "X", "Expa","XS", "U"};

  mxTraj = mxCreateStructMatrix(1,nb_pts,7, fieldnames);
  mxArray * mxX, *mxExpa, *mxU, *mxtime, *mxXS, *mxParam, *mxsensis;

  //  SetArray(sensis,mxsensis);    

  for (int j = 0; j<nb_pts ; j++) {
    p0 = pts(All,j);
    SetArray(p0,mxParam);
    SetArray(sensis,mxsensis);    
    SetArray(*trajArray(j)->X,mxX);
    SetArray(*trajArray(j)->Expa,mxExpa);
    SetArray(*trajArray(j)->XS,mxXS);
    SetArray(*trajArray(j)->U,mxU);
    SetArray(*trajArray(j)->time,mxtime);

    mxSetField(mxTraj, j, "param",mxParam);
    mxSetField(mxTraj, j, "sensis",mxsensis);
    mxSetField(mxTraj, j, "X",mxX);
    mxSetField(mxTraj, j, "Expa",mxExpa);
    mxSetField(mxTraj, j, "XS",mxXS);
    mxSetField(mxTraj, j, "U",mxU);
    mxSetField(mxTraj, j, "time",mxtime);

  }
  
  SetField(Xf,plhs[0],"Xf");
  
  mxAddField(plhs[0],"ExpaMax");
  SetField(ExpaMax,plhs[0],"ExpaMax");

  mxAddField(plhs[0],"XSf");
  SetField(XSf,plhs[0],"XSf");

  mxAddField(plhs[0],"traj");
  mxSetField(plhs[0],0,"traj",mxTraj);

  /* Free memory */ 

  for (int j = 0; j<nb_pts ; j++) {
    delete trajArray(j);
  }

#if _DEBUG >= 1
  cout << "Leaving CVM_ComputeExpTraj ..." << endl;
#endif

  return 0;
}


int CVM_ComputeTraj(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

#if _DEBUG >= 1
  cout << "Entering CVM_ComputeTraj ..." << endl;
#endif

 /* Reading */

  Array2D pts(N,1);
  GetArrayField(pts, prhs[0],"pts");

  Array1D tspan(2);
  GetArray1D(tspan,prhs[1]);

  int test_choice;
  if (nrhs>=3)
    test_choice = mxGetScalar(prhs[2]);
  else
    test_choice = 0;
  
  /* Doing some stuff */

  Range All = Range::all();   
  int nb_pts = pts.extent(1);
  Array<trajectory*,1> trajArray(nb_pts);
  Array2D Xf(N,nb_pts);
  
  Array1D p0(pts.extent(0));
  Array1D Tspan;

  for (int j = 0; j<nb_pts ; j++) {

    p0 = pts(All, j);
    Tspan.resize(tspan.shape());  
    Tspan = tspan;
    trajArray(j) =  new trajectory();
    trajArray(j)->p0 = &p0;

    switch (test_choice) {
    case 0: 
      {
	(trajArray(j))->ComputeTraj(Tspan);
	break;
      }
    case 1: 
      {
	trajArray(j)->ComputeTraj(Tspan, &trajectory::test1);
	break;
      }
    }
  
    int l = trajArray(j)->length; 
    Xf(All,j) = (*trajArray(j)->X)(All,l-1);    
  
  }
 
  /* Writing outputs */ 

  plhs[0] = mxDuplicateArray(prhs[0]);
  mxArray * mxTraj;
  
  const char * fieldnames[4] = {"param", "time", "X", "U"};

  mxTraj = mxCreateStructMatrix(1,nb_pts,4, fieldnames);
  mxArray * mxX,*mxU, *mxtime, *mxParam;

  int j;
  for ( j = 0; j<nb_pts ; j++) {
    p0 = pts(All,j);
    SetArray(p0,mxParam);
    SetArray(*trajArray(j)->X,mxX);
    SetArray(*trajArray(j)->U,mxU);
    SetArray(*trajArray(j)->time,mxtime);
    mxSetField(mxTraj, j, "param",mxParam);
    mxSetField(mxTraj, j, "time",mxtime);
    mxSetField(mxTraj, j, "X",mxX);
    mxSetField(mxTraj, j, "U",mxU);
  }

  /* update status of trajectories if necessary */
  
  if (test_choice) {
    Array1D traj_status(nb_pts);
    mxArray * mxtraj_status;
    for ( j = 0; j<nb_pts ; j++) {
      traj_status(j) =  (double) (*trajArray(j)).status;      
     }
    SetArray(traj_status,mxtraj_status);
    mxAddField(plhs[0],"traj_status");
    mxSetField(plhs[0],0, "traj_status", mxtraj_status);
  }
 
  SetField(Xf,plhs[0],"Xf");
  mxAddField(plhs[0],"traj");
  mxSetField(plhs[0],0, "traj",mxTraj);
 
  /* Free memory */ 

  for (j = 0; j<nb_pts ; j++) {
    delete trajArray(j);
  }
  
#if _DEBUG >= 1
  cout << "Leaving CVM_ComputeTraj ..." << endl;
#endif
  return 0;
}


/*********** Implementations **********************/

void ComputeExpansion(Array2D pts, Array2D XS0, Array2D epsi, Array1D tspan, Array2D &Xf, Array2D &XSf,  Array2D &ExpaMax ) {

#if _DEBUG >= 1
  cout << "Entering ComputeExpansion ..." <<endl;
#endif

  int nb_traj = pts.extent(1);
  Xf.resize(N,nb_traj);
  Xf = 0;
  
  XSf.resize(N*Ns,nb_traj);
  XSf = 0;
  
  /* Sensitivity Reinitialisation */
  
  Array1D expa(N);
  Array2D xS(N,Ns,ColumnMajorArray<2>());
  xS = 0;
  double *xSdata = xS.data();
  ExpaMax.resize(N,nb_traj);
  ExpaMax = 0;
  
  /* CVodes Reinitialization */
  
  InitFdata(f_data, cvm_Mdata->mx_data);
  CVodeSetFdata(cvode_mem,f_data);

  int status, statusFSA;
  int itask = CV_ONE_STEP_TSTOP;
  int dimu;
  double tlast,tret,h;
  booleantype iret = FALSE;

  Fdata* data = (Fdata*) f_data;
  data->ns = Ns;

  realtype *xdata = N_VGetArrayPointer(y);
  Array1D x(xdata,shape(N),neverDeleteData);
  Array1D p(data->p,shape(data->dimp),neverDeleteData);

  double t0 = tspan(0);
  double tout = tspan(1);

  /* declarations for root (guard) functions  */

  data->ns = Ns;
  int *rootsfound = (int*) malloc((data->dimg)*sizeof(int));
 
  int ig;

  /* control */

  void* u_data;
  u_data = (Array1D*) new Array1D(1);  
  GetU(f_data, u_data);
  Array1D * u = (Array1D*) u_data; 
  dimu = (*u).extent(0);

  /* working variables */

  int i,is,j;

  Range Rx(0,N-1);
  Range All = Range::all();   
  
  for(j=0;j<nb_traj;j++) {

    x = pts(Rx,j); // new x0
    p = pts(All,j);

#ifdef _TD 
      x(N)=0.;
#endif

    UpdateFdata(t0,y,f_data,0,NULL);

    tret =t0;
    iret = FALSE;

    for (is=0; is<Ns; is++)
      for (i = 0; i<N; i++) 
	NV_Ith_S(yS0[is],i) = XS0(is*N+i,j);
    
    for (is=0; is<Ns; is++)
      GetData(yS0[is], &xSdata[is*N], N);
    
    for (i = 0; i< N; i++) {
      expa(i) = 0;
      for (is = 0; is <Ns; is++) {	
	expa(i) += abs(xS(i,is))*epsi(is,j);
      }
    }	  
    
    ExpaMax(Range::all(),j) = expa;

    /* Reinit solver */

    switch (itol) {
    case CV_SS:
      status = CVodeReInit(cvode_mem, f, t0, y, itol, reltol, &Sabstol);      
      statusFSA= CVodeSensReInit(cvode_mem,ism,yS0);  
      break;
    case CV_SV:
      status = CVodeReInit(cvode_mem, f, t0, y, itol, reltol, NV_abstol);
      statusFSA = CVodeSensReInit(cvode_mem,ism,yS0);
      break;
    }
    
    CVodeSetStopTime(cvode_mem,tout);
    
    /* Integrate system until tout */
    
    while (1) {
      tlast = tret;
      status = CVode(cvode_mem, tout, y, &tret, itask);
      
      /* break on CVode error */
      if (status < 0) {
	cout << "aie, problem." << endl;
	cout << "status (see cvodes guide):" << status << endl;
	cout << "tout:" << tout << endl;
	cout << "x:" << x << endl;
#if _DEBUG>=3
	cout << "time:" << tret << endl;
#endif
	CVodeGetDky(cvode_mem, tout, 1, y);
	cout << "dx:" << x << endl;

	cout << "CVODES has failed for some reason. Status= " << status << endl;
	mexErrMsgTxt("Dying...");	  

	return;
      }
      
      /* Test if tout was reached */
      CVodeGetCurrentStep(cvode_mem, &h);
      if ( (tret - tout)*h >= 0.0 ) {
	tret = tout;
	CVodeGetDky(cvode_mem, tout, 0, y);
	iret = TRUE;
      }      
   
      /* get sensitivity and compute expansion */

      CVodeGetSens(cvode_mem, tret, yS);
       
      /* Update fdata and resolve discontinuities */
      if (status == CV_ROOT_RETURN) { 
	
	CVodeGetRootInfo(cvode_mem, rootsfound);
	ig =0;
	
	for(i=0; i< data->dimg ; i++) { 
	  if (rootsfound[i]) {
	    ig = i;
	    break;
	  }
	}
	  
	UpdateFdata(tret,y,f_data,ig,yS);	  
	ComputeSensiJump(ig, tret, y, yS, f_data);      
	
	for (is=0;is < Ns; is++) {	  
	  
	  //	  cout << "xsJump[" << is << "]:" << endl;
	  //    N_VPrint_Serial(data->xsJump[is]);
	  
	  for(i=0; i< N; i++) 
	    NV_Ith_S(yS[is],i) += NV_Ith_S(data->xsJump[is],i);	
	}
	
	CVodeGetDky(cvode_mem, tret+h/10, 0, y);
	  
	switch (itol) {
	case CV_SS:
	  CVodeReInit(cvode_mem, f, tret,y, itol, reltol, &Sabstol);
	  CVodeSensReInit(cvode_mem,ism,yS);  
	  break;
	case CV_SV:
	  CVodeReInit(cvode_mem, f, tret,y, itol, reltol, NV_abstol);
	  CVodeSensReInit(cvode_mem,ism,yS);  
	  break;
	}
	CVodeSetInitStep(cvode_mem,h/2);
      }
      else {
	
	if (UpdateFdata(tret,y,f_data,0,NULL)) {
	  
	  CVodeGetDky(cvode_mem, tret, 0, y);
	  
	  switch (itol) {
	  case CV_SS:
	    CVodeReInit(cvode_mem, f, tret,y, itol, reltol, &Sabstol);
	    CVodeSensReInit(cvode_mem,ism,yS);  
	    break;
	  case CV_SV:
	    CVodeReInit(cvode_mem, f, tret,y, itol, reltol, NV_abstol);
	    CVodeSensReInit(cvode_mem,ism,yS);  
	    break;
	  }
	  CVodeSetInitStep(cvode_mem,h/2);
	}
      }
            
      for (is=0; is<Ns; is++)
	GetData(yS[is], &xSdata[is*N], N);
      
      for (i = 0; i< N; i++) {
	expa(i) = 0;
	for (is = 0; is <Ns; is++)
	  expa(i) += abs(xS(i,is))*epsi(is,j);
	ExpaMax(i,j) = max(expa(i), ExpaMax(i,j));
      }	  
      
      /* break if we need to */
      
      if(iret)  break;      
	  
    }

    /* last point */

    /* compute sensitivity */ 
    
    CVodeGetSens(cvode_mem, tret, yS);
    for (is=0; is<Ns; is++)
      GetData(yS[is], &xSdata[is*N], N);

    for (is = 0; is<Ns; is++) {
      for (i = 0; i< N; i++) 
	XSf(is*N+i,j) = xS(i,is);           
    }

    Xf(Range::all(),j) = x;
   
  }

  /* Free some memory  */ 

  free(rootsfound);

#if _DEBUG >= 1
  cout << "Leaving ComputeExpansion ..." << endl;  
#endif   
  
}

void ComputeExpTraj(Array2D pts, Array2D XS0, Array2D epsi, Array1D tspan, Array<trajectory*,1> &trajArray, Array2D &Xf, Array2D &XSf, Array2D &ExpaMax) {

#if _DEBUG >= 1
  cout << "Entering ComputeExpTraj ..." <<endl;
#endif

  int nb_traj = pts.extent(1);
  Xf.resize(N,nb_traj);
  Xf = 0;

  XSf.resize(N*Ns,nb_traj);
  XSf = 0;
  
  trajArray.resize(nb_traj);
  
  /* Sensitivity Reinitialisation */
  
  Array1D expa(N);
  Array2D xS(N,Ns,ColumnMajorArray<2>());
  xS = 0;
  double *xSdata = xS.data();
  ExpaMax.resize(N,nb_traj);
  ExpaMax = 0;

  /* CVodes Reinitialization */

  InitFdata(f_data, cvm_Mdata->mx_data);
  int status, statusFSA;
  int itask = CV_ONE_STEP_TSTOP;
  int dimu;
  double tret,h;
  booleantype iret = FALSE;

  Fdata* data = (Fdata*) f_data;
  data->ns = Ns;
  realtype *xdata = N_VGetArrayPointer(y);
  Array1D x(xdata,shape(N),neverDeleteData);
  Array1D p(data->p,shape(data->dimp),neverDeleteData);

  double t0 = tspan(0);
  double tout = tspan(1);

  /* declarations for root (guard) functions  */

  data->ns = Ns;
  int * rootsfound = (int*) malloc((data->dimg)*sizeof(int));
  int ig;

  /* control */

  void* u_data;
  u_data = (Array1D*) new Array1D(1);  
  GetU(f_data, u_data);
  Array1D * u = (Array1D*) u_data; 
  dimu = (*u).extent(0);

  /* working variables */

  Range Rx(0,N-1);
  Range All = Range::all();   


  int kmax = 5;
  int compteur = 0;
  int i,is,j;
  double tlast(t0);
  
  Array2D *nX,X;
  Array2D *nU,U;
  Array2D *nExpa,Expa;
  
  Array1D *ntime,time;

  for(j=0;j<nb_traj;j++) {

#if _DEBUG>=2
    cout << "----------------------------------------" <<endl;
    cout << " trajectory " << j << endl;
#endif

    x = pts(Rx,j); // new x0
    p = pts(All,j);

#ifdef _TD
	x(N)=0.;
#endif

    UpdateFdata(t0,y,f_data,0,NULL);

    tret = t0;
    compteur = 0; 
    iret = FALSE;

    /* Init output arrays  */
          
    time.resize(kmax);
    time(0) = t0;
    X.resize(N,kmax);
    X(All,0) = x;
    U.resize(dimu,kmax);
    U(All,0)=*u;
    Expa.resize(N,kmax);
    
    for (is=0; is<Ns; is++)
      for (i = 0; i<N; i++) 
	NV_Ith_S(yS0[is],i) = XS0(is*N+i,j);

    for (is=0; is<Ns; is++)
      GetData(yS0[is], &xSdata[is*N], N);
    
    for (i = 0; i< N; i++) {
      expa(i) = 0;
      for (is = 0; is <Ns; is++) {
	expa(i) += abs(xS(i,is))*epsi(is,j);
      }
    }	  

    Expa(All,0) = expa;
    ExpaMax(All,j) = expa;

    /* Reinit solver */

    switch (itol) {
    case CV_SS:
      status = CVodeReInit(cvode_mem, f, t0, y, itol, reltol, &Sabstol);
      statusFSA= CVodeSensReInit(cvode_mem,ism,yS0);
      break;
    case CV_SV:
      status = CVodeReInit(cvode_mem, f, t0, y, itol, reltol, NV_abstol);
      statusFSA = CVodeSensReInit(cvode_mem,ism,yS0);
      break;
    }

    CVodeSetStopTime(cvode_mem,tout);
  
    /* Integrate system until tout */
    
    while (1) {
      
      tlast = tret;

#if _DEBUG >=2
      cout << "Calling CVode .."<< " tout: "<< tout << endl;
      N_VPrint_Serial(y);
#endif
      status = CVode(cvode_mem, tout, y, &tret, itask);

      /* break and diagnose on CVode error */
      
      if (status < 0) {

	cout << "aie, problem." << endl;
	cout << "status (see cvodes guide):" << status << endl;
	cout << "tout:" << tout << endl;
	cout << "x:" << x << endl;

#if _DEBUG>=3
	cout << "time:" << time << endl;
	cout << "X:" << X << endl;
	cout << "U:" << U << endl;
#endif
	CVodeGetDky(cvode_mem, tout, 1, y);
	cout << "dx:" << x << endl;	

	cout << "CVODES complained for some reason. Status= " << status << endl;
	mexErrMsgTxt("Dying...");	      
	return;
      }

      /* Test if tout was reached */
      
      CVodeGetCurrentStep(cvode_mem, &h);
      if ( (tret - tout)*h >= 0.0 ) {
	tret = tout;
	CVodeGetDky(cvode_mem, tout, 0, y);
	iret = TRUE;
      }      

      /* Get sensitivities */

      CVodeGetSens(cvode_mem, tret, yS);
      
      /* Update fdata and resolve discontinuities */
      if (status == CV_ROOT_RETURN) { 
	
	CVodeGetRootInfo(cvode_mem, rootsfound);
	ig =0;
	
	for(i=0; i< data->dimg ; i++) { 
	  if (rootsfound[i]) {
	    ig = i;
	    break;
	  }
	}
	  
	UpdateFdata(tret,y,f_data,ig,yS);	  
	ComputeSensiJump(ig, tret, y, yS, f_data);      
	
	for (is=0;is < Ns; is++) {	  	 	  
	  for(i=0; i< N; i++) 
	    NV_Ith_S(yS[is],i) += NV_Ith_S(data->xsJump[is],i);	
	}
	
	CVodeGetDky(cvode_mem, tret+h/10, 0, y);
	  
	switch (itol) {
	case CV_SS:
	  CVodeReInit(cvode_mem, f, tret,y, itol, reltol, &Sabstol);
	  CVodeSensReInit(cvode_mem,ism,yS);  
	  break;
	case CV_SV:
	  CVodeReInit(cvode_mem, f, tret,y, itol, reltol, NV_abstol);
	  CVodeSensReInit(cvode_mem,ism,yS);  
	  break;
	}
	CVodeSetInitStep(cvode_mem,h/2);
      }
      else {
	
	if (UpdateFdata(tret,y,f_data,0,NULL)) {
	  
	  CVodeGetDky(cvode_mem, tret, 0, y);
	  
	  switch (itol) {
	  case CV_SS:
	    CVodeReInit(cvode_mem, f, tret,y, itol, reltol, &Sabstol);
	    CVodeSensReInit(cvode_mem,ism,yS);  
	    break;
	  case CV_SV:
	    CVodeReInit(cvode_mem, f, tret,y, itol, reltol, NV_abstol);
	    CVodeSensReInit(cvode_mem,ism,yS);  
	    break;
	  }
	  CVodeSetInitStep(cvode_mem,h/2);
	}
      }

      for (is=0; is<Ns; is++)
	GetData(yS[is], &xSdata[is*N], N);
      
      for (i = 0; i< N; i++) {
	expa(i) = 0;
	for (is = 0; is <Ns; is++)
	  expa(i) += abs(xS(i,is))*epsi(is,j);
	ExpaMax(i,j) = max(expa(i), ExpaMax(i,j));
      }	            


                
      /* get input value */
      
      GetU(f_data, u_data);  
      
      compteur++;
      if (compteur < kmax) {

	time(compteur) = tret;
	U(All,compteur)= *u;
	X(All,compteur)= x;
	Expa(All,compteur)= expa; 
	  
      }
      else  {// reached max number of point, must resize arrays. Try to guess how many points left.

	  kmax += max(1,min((int) ceil((tout-tret)/(tret-tlast)),kmax));

#if _DEBUG>=1
	  cout << "tout: "<< tout << " tret: "<< tret << " tlast: "<< tlast << endl;
	  cout << "kmax:" << kmax << endl;
#endif
	//	kmax *= (int) ceil((tout-t0)/(tret-t0))+1;
	
	X.resizeAndPreserve(N,kmax);
	U.resizeAndPreserve(dimu,kmax);
	Expa.resizeAndPreserve(N,kmax);
	time.resizeAndPreserve(kmax);
	  
	time(compteur) = tret;
	U(All,compteur)= *u;
	X(All,compteur)= x;     
	Expa(All,compteur)= expa; 
	
      }
            
      /* break if we need to */
      
      if(iret) { 
     
#if _DEBUG>=2 
	  cout << " End of trajectory." << endl;
#endif 
      
      break;      
      }
    }

    X.resizeAndPreserve(N,compteur+1);
    U.resizeAndPreserve(dimu,compteur+1);
    Expa.resizeAndPreserve(N,compteur+1);
    time.resizeAndPreserve(compteur+1);
        
    /* Mark new explored pts and only keep those in the trajectory */
 
    Xf(All, j) = X(All,compteur);
    for (is = 0; is<Ns; is++) {
      for (i = 0; i< N; i++) 
	XSf(is*N+i,j) = xS(i,is);
    }  
      
    trajArray(j) = new trajectory();
    nX = (trajArray(j)->X);
    nU = (trajArray(j)->U);
    nExpa = (trajArray(j)->Expa);
    ntime = (trajArray(j)->time);
    
    (*nX).resize(N,compteur+1);
    (*nU).resize(dimu,compteur+1);
    (*nExpa).resize(N,compteur+1);
    (*ntime).resize(compteur+1);
    
    (*nX) = X;
    (*nU) = U;
    (*nExpa) = Expa;
    (*ntime) = time;
  
  }

  free(rootsfound);

#if _DEBUG >= 1
  cout << "Leaving ComputeExpansion ..." << endl;  
#endif  

}

int giDQ(int ig,realtype t, N_Vector x, realtype* gx, realtype *dgi,void *g_data) {

  /* jac_data points to cvode_mem */
  Fdata* data= (Fdata*) g_data;

  realtype gnorm, minInc, inc, inc_inv, xjsaved, srur;

  realtype* gtemp = (realtype*) malloc( (gtemp[data->dimg])*sizeof(realtype));
  realtype  *x_data, *ewt_data;
  realtype h; 
  N_Vector nvgx;
  CVodeGetCurrentStep(cvode_mem,&h);
  				   
  int j;
  N_Vector ewt;
  ewt = N_VNew_Serial(N);
  CVodeGetErrWeights(cvode_mem, ewt);

  nvgx = N_VNew_Serial(data->dimg);
  N_VSetArrayPointer(gx,nvgx);

  
  /* Obtain pointers to the data for ewt, y */
  ewt_data = N_VGetArrayPointer(ewt);
  x_data   = N_VGetArrayPointer(x);
  
  /* Set minimum increment based on uround and norm of g */
  srur = RSqrt(DBL_EPSILON);

  gnorm = N_VWrmsNorm(nvgx, ewt);
  minInc = (gnorm != ZERO) ?
    (MIN_INC_MULT * ABS(h) * DBL_EPSILON * N * gnorm) : ONE;

  for (j = 0; j < N; j++) {
  
    xjsaved = x_data[j];
    inc = MAX(srur*ABS(xjsaved), minInc/ewt_data[j]);
    x_data[j] += inc;
    
    g(t, x, gtemp, f_data);
    x_data[j] = xjsaved;

    inc_inv = ONE/inc;
    dgi[j] = (gtemp[ig]-gx[ig])*inc_inv;
   }

  free(gtemp);
  return 0;

}

int ComputeSensiJump(int ig, realtype t, N_Vector x, N_Vector *xS, void* g_data) {

  //cout << "Entering ComputeSensiJump..." << endl;
  
  Fdata* data = (Fdata*) g_data;
  int i,is; 
  double dgs,dgf,dtau;

  N_Vector f1 = N_VNew_Serial(N);
  N_Vector f2 = N_VNew_Serial(N);

  N_Vector x1 = N_VNew_Serial(N);
  N_Vector x2 = N_VNew_Serial(N);
   
  realtype *dgi =  (realtype*) malloc(N*sizeof(realtype));
  realtype *gval=  (realtype*) malloc((data->dimg)*sizeof(realtype));
  
  for (i=0; i<data->dimg; i++)
    gval[i] =0;
  
  giDQ(ig,t,x,gval,dgi,g_data);
  realtype h;
  realtype tcur;
  
  CVodeGetCurrentStep(cvode_mem, &h);
  CVodeGetCurrentTime(cvode_mem,&tcur);

  realtype tinc = 200*DBL_EPSILON*(abs(tcur)+abs(h));
  realtype tm = max(t-tinc, tcur-h);
  realtype tp = min(t+tinc, tcur);
  
  CVodeGetDky(cvode_mem,tm, 0, x1);
  CVodeGetDky(cvode_mem, tp, 0, x2);
     
  f(tm,x1,f1,f_data);
  UpdateFdata(tp,x2,f_data,ig,yS);	  
  f(tp,x2,f2,f_data);

  for (is=0; is<Ns; is++) {
    
    dgs =0;
    dgf =0;
    
    for(i=0; i< N; i++) {
      dgs += dgi[i]*Ith(xS[is],i);  //  < dg , s_x0 >	  
      dgf += dgi[i]*Ith(f1,i);      //  < dg , f1 >
      Ith(data->xsJump[is],i)=0;    //   reset jump
    }
   
    dtau = dgs/dgf;      
    dtau = min(dtau,100.);
    
    for(i=0; i< N; i++) {	  
      //     cout <<  Ith(f2,i) - Ith(f1,i) << " ";
      Ith(data->xsJump[is],i) = dtau* (Ith(f2,i) - Ith(f1,i));	  
    }
  }          
  //  cout << "Leaving ComputeSensiJump..." << endl;

  free(dgi);
  free(gval);

  return 0;
}

void trajectory::ComputeTraj(Array1D& tspan) {

#if _DEBUG >= 1
  cout << "Entering traj::ComputeTraj ..." <<endl;
#endif
 
  /* CVodes Reinitialization */
 
#if _DEBUG >= 2
  cout << "Init CVodes..." << endl;
#endif  
 
  //  InitFdata(f_data, cvm_Mdata->mx_data);
  int status;
  int itask = CV_ONE_STEP_TSTOP;
  int dimu;
  double tlast,tret,h;
  booleantype iret = FALSE;

  Fdata* data = (Fdata*) f_data;
  realtype *xdata = N_VGetArrayPointer(y);
  Array1D x(xdata,shape(N),neverDeleteData);
  Array1D p(data->p,shape(data->dimp),neverDeleteData);

  current_x = &x;
  double t0 = tspan(0);
  double tout;

  /* control */ 

#if _DEBUG >= 2
  cout << "Control Variable ..." << endl;
#endif  

  void * u_data;
  u_data = (Array1D*) new Array1D(0);  
  GetU(f_data,  u_data);
  Array1D* u = (Array1D*) u_data;
  dimu = (*u).extent(0);

  /* working variables */

  Range Rx(0,N-1);
  Range All = Range::all();   

  int kmax=5;
  int compteur = 0;
  int i,j;
  int nb_points = tspan.extent(0);
  int statusgdk;

  if (nb_points == 2) {
    
    tout = tspan(1);          
    x = (*p0)(Rx); // new x0
    p = (*p0)(All);
#ifdef _TD
    x(N)=0.;
#endif	    
    UpdateFdata(t0,y,f_data,0,NULL);		
    compteur = 0;
    iret = FALSE;
    tret = t0;
	
    /* Reinit solver */
      
    switch (itol) {
    case CV_SS:
      CVodeReInit(cvode_mem, f, t0, y, itol, reltol, &Sabstol);      
      break;
    case CV_SV:
      CVodeReInit(cvode_mem, f, t0, y, itol, reltol, NV_abstol);
      break;
    }

    CVodeSetStopTime(cvode_mem,tout);
      
    /* Init root finding for hybrid */

    if (data->dimg>0){
      CVodeRootInit(cvode_mem, data->dimg, g, f_data);	  	  
    }
    
    /* Init Output Arrays */
    
    time->resize(kmax);
    (*time)(0) = t0;
    X->resize(N,kmax);
    (*X)(All,0) = x;
    U->resize(dimu,kmax);
    (*U)(All,0)=*u;
    
    /* Integrate system until tout */
	
    while (1) {
	  
      /* Integrate one step */
      tlast = tret;
      
#if _DEBUG >= 3
      cout << " CVodes Call ... " << endl;
#endif        
      status = CVode(cvode_mem, tout, y, &tret, itask);
	  
      /* break on CVode error */
      if (status < 0) {
	cout << "tout:" << tout << endl;
	cout << "x:" << x << endl;
#if _DEBUG >=3
	cout << "time:" << time << endl;
	cout << "X:" << *X << endl;
	cout << "U:" << *U << endl;
#endif
	CVodeGetDky(cvode_mem, tout, 1, y);
	cout << "dx:" << x << endl;
	cout << "CVODES failed miserably for some reason. Status= " << status << endl;
	mexErrMsgTxt("Dying...");	  
	return;
      }

      /* Test if tout was reached */
      CVodeGetCurrentStep(cvode_mem, &h);
      
      if ( (tret - tout)*h >= 0.0 ) {
	tret = tout;
	CVodeGetDky(cvode_mem, tout, 0, y);
	iret = TRUE;
      }      
	  
#if _DEBUG >= 3
      cout << " On avance, on avance.." << endl;
      cout << "tret:" << tret <<  endl;
#endif
	  
      /* Update f_data */
	  
      if (UpdateFdata(tret,y,f_data,0,NULL) || status == CV_ROOT_RETURN) {  //Update f_data. if changed at a discontinuity, reinit cvode 
	    
	CVodeGetDky(cvode_mem, tret, 0, y);
	
	switch (itol) {
	case CV_SS:
	  CVodeReInit(cvode_mem, f, tret,y, itol, reltol, &Sabstol);
	  break;
	case CV_SV:
	  CVodeReInit(cvode_mem, f, tret,y, itol, reltol, NV_abstol);
	  break;
	}
	CVodeSetInitStep(cvode_mem,h);
      }
      
      /* get input value */
      
#if _DEBUG >= 3
      cout << " Appel de GetU " << endl;
#endif       
      
      GetU(f_data, u_data);

      /* store step done, resizing arrays if necessary */

#if _DEBUG >= 3
      cout << " Resizing ... " << endl;
#endif       
      
      compteur++;
      if (compteur < kmax) {
	
#if _DEBUG >= 3
	cout << " no need. " << endl;
#endif       
	
	(*time)(compteur) = tret;
	(*U)(All,compteur)= *u;
	(*X)(All,compteur)= x;     
      }

      else  { // try to guess how many points left
	    
	kmax += max(1,min((int) ceil((tout-tret)/(tret-tlast)),kmax));
	
	(*X).resizeAndPreserve(N,kmax);
	(*U).resizeAndPreserve(dimu,kmax);
	(*time).resizeAndPreserve(kmax);
	
	(*time)(compteur) = tret;
	(*U)(All,compteur)= *u;
	(*X)(All,compteur)= x;     
      }
    
      /* break if we need to */
      
      if(iret)  break;      
      //      cout << "loop or not .. " << endl;
    }
    
#if _DEBUG >= 2
    cout << "Writing trajectory ..." << endl;
#endif
    
    /* resize */   
    length = compteur+1;
    (*X).resizeAndPreserve(N,compteur+1);
    (*U).resizeAndPreserve(dimu,compteur+1);
    (*time).resizeAndPreserve(compteur+1);    	
  }
  
  else { 
    
   if ((nb_points == 3)&&(tspan(1)==tspan(2))) {      
      tspan.resizeAndPreserve(2);
      nb_points = 2;
    }
   length = nb_points;
   (*X).resize(N,nb_points);
   (*U).resize(dimu,nb_points);
   (*time).resize(nb_points);
   (*time) = tspan;
    
#if _DEBUG >= 2
    cout << "Loop over trajectories..." << endl;
#endif
               
    x = (*p0)(Rx); // new x0
    p = (*p0)(All);
 
    UpdateFdata(t0,y,f_data,0,NULL);
            
    (*X)(All,0)= x;
    GetU(f_data, u_data);
    (*U)(All,0)=*u;
      
#ifdef _TD
    x(N)=0.;
#endif	
    iret = FALSE;
    
    /* Reinit solver */
    
    switch (itol) {
    case CV_SS:
      status = CVodeReInit(cvode_mem, f, t0, y, itol, reltol, &Sabstol);      
      break;
    case CV_SV:
      status = CVodeReInit(cvode_mem, f, t0, y, itol, reltol, NV_abstol);
      break;
    }
    
    /* Init root finding for hybrid */
    
    if (data->dimg>0) {
      CVodeRootInit(cvode_mem, data->dimg, g, f_data);
    }
    
    for(i=1; i<nb_points;i++) {
      
      iret = FALSE;
      tout=tspan(i);
      
      CVodeSetStopTime(cvode_mem,tout);

      /* Integrate system until tout */
      
      while (1) {
	
	/* Integrate one step */
	status = CVode(cvode_mem, tout, y, &tret, itask);
	
	/* break on CVode error */
	if (status < 0) return;   
	  
	/* Test if tout was reached */
	CVodeGetCurrentStep(cvode_mem, &h);
	
	if ( (tret - tout)*h >= 0.0 ) {
	  tret = tout;
	  CVodeGetDky(cvode_mem, tout, 0, y);
	  iret = TRUE;
	}      
	
	/* Update f_data */
	
	if (UpdateFdata(tret,y,f_data,0,NULL)||status == CV_ROOT_RETURN) { /* if f_data changed in a discontinuity, reinit cvode */
	  CVodeGetDky(cvode_mem, tret, 0, y);
	  
	  switch (itol) {
	  case CV_SS:
	    status = CVodeReInit(cvode_mem, f, tret,y, itol, reltol, &Sabstol);      
	    break;
	  case CV_SV:
	    status = CVodeReInit(cvode_mem, f, tret,y, itol, reltol, NV_abstol);
	  }
	  CVodeSetInitStep(cvode_mem,h);
	}
	  
	/* break if we need to */
	if(iret)  break;      
      }
	
      GetU(f_data, u_data);
      (*U)(All,i)= *u;
      (*X)(All,i)= x;	
    }       
  }

#if _DEBUG >= 1
  cout << "Leaving ComputeTraj." << endl;
#endif
}

void trajectory::ComputeTrajSensi(Array1D& tspan) {

 #if _DEBUG >= 1
  cout << "Entering trajectory::ComputeTrajSensi ..." <<endl;
 #endif

  /* Sensitivity Reinitialisation */
  
  Array2D xS(N,Ns,ColumnMajorArray<2>());
  xS = 0;
  double *xSdata = xS.data();

  /* CVodes Reinitialization */

  InitFdata(f_data, cvm_Mdata->mx_data);
  CVodeSetMaxConvFails(cvode_mem,20);

  int status, statusFSA;
  int itask = CV_ONE_STEP_TSTOP;
  int dimu;
  double tret,h;
  booleantype iret = FALSE;
  Fdata* data = (Fdata*) f_data;
  realtype *xdata = N_VGetArrayPointer(y);
  Array1D x(xdata,shape(N),neverDeleteData);
  Array1D p(data->p,shape(data->dimp),neverDeleteData);

  /* declarations for root (guard) functions  */

  data->ns = Ns;
  int * rootsfound  = (int*) malloc( sizeof(int)*(data->dimg));
  int ig;

  /* control */

  void* u_data;
  u_data = (Array1D*) new Array1D(0);  
  GetU(f_data, u_data);
  Array1D * u = (Array1D*) u_data; 
  dimu = (*u).extent(0);

  /* working variables */

  Range Rx(0,N-1);
  Range All = Range::all();   

  int kmax = 5;
  int compteur = 0;
  int i,is,j;
  double tlast,tcur;
  int nb_points = tspan.extent(0);

  double t0, tout;

  if (nb_points == 2) {

    t0 = tspan(0);
    tlast = tret = t0;
    tout = tspan(1);
    
    x = (*p0)(Rx); // new x0
    p = (*p0)(All);
    
#ifdef _TD
    x(N)=0.;
#endif

    //    InitFdata(f_data, cvm_Mdata->mx_data);
    UpdateFdata(t0,y,f_data,0,NULL);

    compteur = 0; 
    iret = FALSE;
    
    /* Init output arrays  */
    
    time->resize(kmax);
    (*time)(0) = t0;
    X->resize(N,kmax);
    (*X)(All,0) = x;
    U->resize(dimu,kmax);
    (*U)(All,0)=*u;
    XS->resize(N*Ns,kmax);
        
    for (is=0; is<Ns; is++)
      GetData(yS0[is], &xSdata[is*N], N);
      
    for (is=0; is<Ns; is++) 
      for (i=0; i <N; i++)
	(*XS)(i+is*N ,0) =xS(i,is) ;	    
    
    /* Reinit solver */

    switch (itol) {
    case CV_SS:
      status = CVodeReInit(cvode_mem, f, t0, y, itol, reltol, &Sabstol);
      statusFSA= CVodeSensReInit(cvode_mem,ism,yS0);
      break;
    case CV_SV:
      status = CVodeReInit(cvode_mem, f, t0, y, itol, reltol, NV_abstol);
      statusFSA = CVodeSensReInit(cvode_mem,ism,yS0);
      break;
    }

    CVodeSetStopTime(cvode_mem,tout);
  
    /* Init root finding for hybrid */
      
    if (data->dimg>0){
	CVodeRootInit(cvode_mem, data->dimg, g, f_data);	  	  
    }

    /* Integrate system until tout */
      
      while (1) {
	
	tlast = tret;
      
#if _DEBUG >=2
	cout << "Calling CVode .."<< " tout: "<< tout << endl;
	N_VPrint_Serial(y);
#endif
	status = CVode(cvode_mem, tout, y, &tret, itask);
        
	/* break and diagnose on CVode error */
	
	if (status < 0) {
	
	  cout << "aie, problem." << endl;
	  cout << "status (see cvodes guide):" << status << endl;
	  cout << "traj no: " << j << endl;
	  cout << "tout:" << tout << endl;	  
	  cout << "p:" << p << endl;
	  cout << "x:" << x << endl;
#if _DEBUG>=3
	  cout << "time:" << *time << endl;
	  cout << "X:" << *X << endl;
	  cout << "XS:" << *XS << endl;
	  cout << "U:" << *U << endl;
#endif
	  CVodeGetDky(cvode_mem, tout, 1, y);
	  cout << "dx:" << x << endl;	
	  cout << "CVODES failed for some reason. Yep. Status= " << status << endl;
	  mexErrMsgTxt("Dying...");	      
	
	  return;
	}
      
	/* Test if tout was reached */
	
	CVodeGetCurrentStep(cvode_mem, &h);
	if ( (tret - tout)*h >= 0.0 ) {
	  tret = tout;
	  CVodeGetDky(cvode_mem, tout, 0, y);
	  iret = TRUE;
	}      
      
	/* get sensitivity */
      
	CVodeGetSens(cvode_mem, tret, yS);
      
	/* Tests if a root is found  */
	
	if (status == CV_ROOT_RETURN) { 
	  
	  //cout << "Transitions !" << endl;

	  CVodeGetRootInfo(cvode_mem, rootsfound);
	  ig =0;
	
	  for(i=0; i< data->dimg ; i++) { 
	    if (rootsfound[i]) {
	      ig = i;
	      break;
	    }
	  }
	  
	  ComputeSensiJump(ig, tret, y, yS, f_data);      

	  for (is=0;is < Ns; is++) {	  
	    
	    //	  cout << "xsJump[" << is << "]:" << endl;
	    //    N_VPrint_Serial(data->xsJump[is]);
	  	   
	    for(i=0; i< N; i++) 
	      NV_Ith_S(yS[is],i) += NV_Ith_S(data->xsJump[is],i);	
	  }
	  
	  CVodeGetCurrentTime(cvode_mem,&tcur);
	  tret= min(tcur, tret+abs((tcur-tret)/1e6));
	  CVodeGetDky(cvode_mem, tret, 0, y);
	  
	  switch (itol) {
	  case CV_SS:
	    CVodeReInit(cvode_mem, f, tret,y, itol, reltol, &Sabstol);
	    CVodeSensReInit(cvode_mem,ism,yS);  
	    break;
	  case CV_SV:
	    CVodeReInit(cvode_mem, f, tret,y, itol, reltol, NV_abstol);
	    CVodeSensReInit(cvode_mem,ism,yS);  
	    break;
	  }
	  CVodeSetInitStep(cvode_mem,h/2);
	}
	else {
      
	  if (UpdateFdata(tret,y,f_data,0,NULL)) {
	    
	    CVodeGetDky(cvode_mem, tret, 0, y);
	    
	    switch (itol) {
	    case CV_SS:
	      CVodeReInit(cvode_mem, f, tret,y, itol, reltol, &Sabstol);
	      CVodeSensReInit(cvode_mem,ism,yS);  
	      break;
	    case CV_SV:
	      CVodeReInit(cvode_mem, f, tret,y, itol, reltol, NV_abstol);
	      CVodeSensReInit(cvode_mem,ism,yS);  
	      break;
	    }
	    CVodeSetInitStep(cvode_mem,h/2);
	  }
	}
    
	/* Store sensitivity matrix */
	
	for (is=0; is<Ns; is++)
	  GetData(yS[is], &xSdata[is*N], N);
		
         
      /* get input value */
	
	GetU(f_data, u_data);  
      
	compteur++;
	if (compteur < kmax) {
	  
	  (*time)(compteur) = tret;
	  (*U)(All,compteur)= *u;
	  (*X)(All,compteur)= x;	  
	  for (is=0; is<Ns; is++) 
	    for (i=0; i <N; i++)
	      (*XS)(i+is*N ,compteur) =xS(i,is) ;
	
	}
	else  {// reached max number of point, must resize arrays. Try to guess how many points left.
	
	  kmax += max(1,min((int) ceil((tout-tret)/(tret-tlast)),kmax));
	  X->resizeAndPreserve(N,kmax);
	  U->resizeAndPreserve(dimu,kmax);	  
	  XS->resizeAndPreserve(N*Ns,kmax);
	  time->resizeAndPreserve(kmax);
	  
	  (*time)(compteur) = tret;
	  (*U)(All,compteur)= *u;
	  (*X)(All,compteur)= x;     
	  for (is=0; is<Ns; is++) 
	    for (i=0; i <N; i++)
	      (*XS)(i+is*N ,compteur) =xS(i,is) ;	    
	}
            
	/* break if we need to */    
	
	if(iret) { 
	
#if _DEBUG>=2 
	  cout << " End of trajectory." << endl;
#endif 
      
	  break;      
	}
      }

      length = compteur+1;
      (*X).resizeAndPreserve(N,compteur+1);
      (*U).resizeAndPreserve(dimu,compteur+1);
      (*XS).resizeAndPreserve(N*Ns,compteur+1);
      (*time).resizeAndPreserve(compteur+1);
      
  }
  else {
    if ((nb_points == 3)&&(tspan(1)==tspan(2))) {      
      tspan.resizeAndPreserve(2);    
      nb_points = 2;
    }

    (*X).resize(N,nb_points);
    (*U).resize(dimu,nb_points);
    (*time).resize(nb_points);
    *time = tspan;    
    XS->resize(N*Ns,nb_points);
    length = nb_points;
    x = (*p0)(Rx); // new x0
    p = (*p0)(All);
    
#ifdef _TD
    x(N)=0.;
#endif
    
    //    InitFdata(f_data, cvm_Mdata->mx_data);
    UpdateFdata(t0,y,f_data,0,NULL);
    
    (*X)(All,0)= x;
    GetU(f_data, u_data);
    (*U)(All,0)=*u;
    
    iret = FALSE;
    
    /* Init output arrays  */
        
    for (is=0; is<Ns; is++)
      GetData(yS0[is], &xSdata[is*N], N);
    
    for (is=0; is<Ns; is++) 
      for (i=0; i <N; i++)
	(*XS)(i+is*N ,0) = xS(i,is) ;	    
    
    /* Reinit solver */

    switch (itol) {
    case CV_SS:
      status = CVodeReInit(cvode_mem, f, t0, y, itol, reltol, &Sabstol);
      statusFSA= CVodeSensReInit(cvode_mem,ism,yS0);
      break;
    case CV_SV:
      status = CVodeReInit(cvode_mem, f, t0, y, itol, reltol, NV_abstol);
      statusFSA = CVodeSensReInit(cvode_mem,ism,yS0);
      break;
    }
      
    /* Init root finding for hybrid */
      
    if (data->dimg>0){
      CVodeRootInit(cvode_mem, data->dimg, g, f_data);	  	  
    }

    /* Integrate system until tout */
    
    for(int k=1; k<nb_points;k++) {
      
      iret = FALSE;
      tout = tspan(k);
      
      CVodeSetStopTime(cvode_mem,tout);
      //	tlast = tret;
      
      while (1) {
	  
	/* Integrate one step */
	status = CVode(cvode_mem, tout, y, &tret, itask);
	
	/* break and diagnose on CVode error */
	  
	if (status < 0) {
	
	  cout << "aie, problem." << endl;
	  cout << "status (see cvodes guide):" << status << endl;
	  cout << "traj no: " << j << endl;
	  cout << "tout:" << tout << endl;
	  cout << "p:" << p << endl;
	  cout << "x:" << x << endl;
#if _DEBUG>=3
	  cout << "time:" << *time << endl;
	  cout << "X:" << *X << endl;
	  cout << "XS:" << *XS << endl;
	  cout << "U:" << *U << endl;
#endif
	  CVodeGetDky(cvode_mem, tout, 1, y);
	  cout << "dx:" << x << endl;	
	  cout << "CVODES failed for some reason. Yep. Status= " << status << endl;
	  mexErrMsgTxt("Dying...");	      
	
	  return;
	}
      
	/* Test if tout was reached */
	
	CVodeGetCurrentStep(cvode_mem, &h);
	if ( (tret - tout)*h >= 0.0 ) {
	  tret = tout;
	  CVodeGetDky(cvode_mem, tout, 0, y);
	  iret = TRUE;
	}      
	
	/* get sensitivity and compute expansion */
      
	CVodeGetSens(cvode_mem, tret, yS);
	  
	/* Update f_data. if changed at a discontinuity, reinit cvode */
	
	if (status == CV_ROOT_RETURN) { 

	  // cout << "Transitions !" << endl;
	  
	  CVodeGetRootInfo(cvode_mem, rootsfound);
	  ig =0;
	  
	  for(i=0; i< data->dimg ; i++) { 
	    if (rootsfound[i]) {
	      ig = i;
	      break;
	    }
	  }
	  
	  // UpdateFdata(tret,y,f_data,ig,yS);	  
	  ComputeSensiJump(ig, tret, y, yS, f_data);      
	  
	  for (is=0;is < Ns; is++) {	  
	    
	    //	  cout << "xsJump[" << is << "]:" << endl;
	    //    N_VPrint_Serial(data->xsJump[is]);
	    
	    for(i=0; i< N; i++) 
	      NV_Ith_S(yS[is],i) += NV_Ith_S(data->xsJump[is],i);	
	  }
	  
	  CVodeGetCurrentTime(cvode_mem,&tcur);
	  tret= min(tcur, tret+abs((tcur-tret)/1e6));
	  CVodeGetDky(cvode_mem, tret, 0, y);
	  
	  switch (itol) {
	  case CV_SS:
	    CVodeReInit(cvode_mem, f, tret,y, itol, reltol, &Sabstol);
	    CVodeSensReInit(cvode_mem,ism,yS);  
	    break;
	  case CV_SV:
	    CVodeReInit(cvode_mem, f, tret,y, itol, reltol, NV_abstol);
	    CVodeSensReInit(cvode_mem,ism,yS);  
	    break;
	  }
	  CVodeSetInitStep(cvode_mem,h/2);
	}
	  else {
	    
	    if (UpdateFdata(tret,y,f_data,0,NULL)) {
	    
	    CVodeGetDky(cvode_mem, tret, 0, y);
	    
	    switch (itol) {
	    case CV_SS:
	      CVodeReInit(cvode_mem, f, tret,y, itol, reltol, &Sabstol);
	      CVodeSensReInit(cvode_mem,ism,yS);  
	      break;
	    case CV_SV:
	      CVodeReInit(cvode_mem, f, tret,y, itol, reltol, NV_abstol);
	      CVodeSensReInit(cvode_mem,ism,yS);  
	      break;
	    }
	    CVodeSetInitStep(cvode_mem,h/2);
	    }
	  }
    	  	
	/* break if we need to */    
	if(iret) break;      
      }
      
      /* Store sensitivity matrix */
	
      for (is=0; is<Ns; is++)
	GetData(yS[is], &xSdata[is*N], N);
	               
      /* get input value and store things */
	
      GetU(f_data, u_data);  
	
      (*U)(All,k)= *u;
      (*X)(All,k)= x;	      

      for (is=0; is<Ns; is++) 
	for (i=0; i <N; i++)
	  (*XS)(i+is*N, k) =xS(i,is) ;       
    }	      
  }

  free(rootsfound);
#if _DEBUG >= 1
  cout << "Leaving trajectory::ComputeTrajSensi ..." << endl;  
#endif    
}

void trajectory::ComputeExpa(const Array1D& epsi, Array1D &ExpaMax) {

#if _DEBUG >= 1
  cout << "Entering trajectory::ComputeExpa ..." <<endl;
#endif

  if (length == 0) 
    return;

  Expa->resize(N, length);
  Array1D expa(N);

  ExpaMax = 0;
  int i, is;
  for(int k = 0; k<length; k++){    
    for (i = 0; i< N; i++) {
      expa(i) = 0;
      for (is = 0; is <Ns; is++)
	expa(i) += abs((*XS)(i+is*N,k))*epsi(is);

      ExpaMax(i) = max(expa(i), ExpaMax(i));
    }        
    (*Expa)(Range::all(),k)=expa;
  }
  
#if _DEBUG >= 1
  cout << "Leaving trajectory::ComputeExpa ..." << endl;  
#endif    
}

void trajectory::ComputeTraj(Array1D& tspan, int (trajectory::*test_function)()) {

#if _DEBUG >= 1
  cout << "Entering traj::ComputeTraj ..." <<endl;
#endif
   
  /* CVodes Reinitialization */
 
#if _DEBUG >= 2
  cout << "Init CVodes..." << endl;
#endif  
 
  InitFdata(f_data, cvm_Mdata->mx_data);
  int status;
  int itask = CV_ONE_STEP_TSTOP;
  int dimu;
  double tlast,tret,h;
  booleantype iret = FALSE;

  Fdata* data = (Fdata*) f_data;
  realtype *xdata = N_VGetArrayPointer(y);
  Array1D x(xdata,shape(N),neverDeleteData);
  Array1D p(data->p,shape(data->dimp),neverDeleteData);

  current_x = &x;
  
  double t0 = tspan(0);
  double tout;

  /* control */ 

#if _DEBUG >= 2
  cout << "Control Variable ..." << endl;
#endif  

  void * u_data;
  u_data = (Array1D*) new Array1D(0);  
  GetU(f_data,  u_data);
  Array1D* u = (Array1D*) u_data;
  dimu = (*u).extent(0);

  /* working variables */

  Range Rx(0,N-1);
  Range All = Range::all();   

  int kmax=5;
  int compteur = 0;
  int i,j;
  int nb_points = tspan.extent(0);
  int statusgdk;

  if (nb_points == 2) {
    
    tout = tspan(1);          
    x = (*p0)(Rx); // new x0
    p = (*p0)(All);

#ifdef _TD
    x(N)=0.;
#endif	    
    UpdateFdata(t0,y,f_data,0,NULL);		
    compteur = 0;
    iret = FALSE;
    tret = t0;
	
    /* Reinit solver */
      
    switch (itol) {
    case CV_SS:
      CVodeReInit(cvode_mem, f, t0, y, itol, reltol, &Sabstol);      
      break;
    case CV_SV:
      CVodeReInit(cvode_mem, f, t0, y, itol, reltol, NV_abstol);
      break;
    }

    CVodeSetStopTime(cvode_mem,tout);
      
    /* Init root finding for hybrid */

    if (data->dimg>0){
      CVodeRootInit(cvode_mem, data->dimg, g, f_data);	  	  
    }
    
    /* Init Output Arrays */
    
    time->resize(kmax);
    (*time)(0) = t0;
    X->resize(N,kmax);
    (*X)(All,0) = x;
    U->resize(dimu,kmax);
    (*U)(All,0)=*u;
    
    /* Integrate system until tout */
	
    while (1) {
	  
      /* Integrate one step */
      tlast = tret;
      
#if _DEBUG >= 3
      cout << " CVodes Call ... " << endl;
#endif        
      status = CVode(cvode_mem, tout, y, &tret, itask);

      /* break on CVode error */
      if (status < 0) {
	cout << "tout:" << tout << endl;
	cout << "x:" << x << endl;
#if _DEBUG >=3
	cout << "time:" << time << endl;
	cout << "X:" << *X << endl;
	cout << "U:" << *U << endl;
#endif
	CVodeGetDky(cvode_mem, tout, 1, y);
	cout << "dx:" << x << endl;
	cout << "CVODES failed miserably for some reason. Status= " << status << endl;
	mexErrMsgTxt("Dying...");	  
	return;
      }

      /* Test if tout was reached */
      CVodeGetCurrentStep(cvode_mem, &h);
      
      if ( (tret - tout)*h >= 0.0 ) {
	tret = tout;
	CVodeGetDky(cvode_mem, tout, 0, y);
	iret = TRUE;
      }      
	  
#if _DEBUG >= 3
      cout << " On avance, on avance.." << endl;
      cout << "tret:" << tret <<  endl;
#endif
	  
      /* Update f_data */
	  
      if (UpdateFdata(tret,y,f_data,0,NULL) || status == CV_ROOT_RETURN) {  //Update f_data. if changed at a discontinuity, reinit cvode 
	    
	CVodeGetDky(cvode_mem, tret, 0, y);
	
	switch (itol) {
	case CV_SS:
	  CVodeReInit(cvode_mem, f, tret,y, itol, reltol, &Sabstol);
	  break;
	case CV_SV:
	  CVodeReInit(cvode_mem, f, tret,y, itol, reltol, NV_abstol);
	  break;
	}
	CVodeSetInitStep(cvode_mem,h);
      }
      
      /* get input value */
      
#if _DEBUG >= 3
      cout << " Appel de GetU " << endl;
#endif       
      
      GetU(f_data, u_data);

      /* store step done, resizing arrays if necessary */

#if _DEBUG >= 3
      cout << " Resizing ... " << endl;
#endif       
      
      compteur++;
      if (compteur < kmax) {
	
#if _DEBUG >= 3
	cout << " no need. " << endl;
#endif       
	
	(*time)(compteur) = tret;
	(*U)(All,compteur)= *u;
	(*X)(All,compteur)= x;     
      }

      else  { // try to guess how many points left
	    
	kmax += max(1,min((int) ceil((tout-tret)/(tret-tlast)),kmax));
	
	(*X).resizeAndPreserve(N,kmax);
	(*U).resizeAndPreserve(dimu,kmax);
	(*time).resizeAndPreserve(kmax);
	
	(*time)(compteur) = tret;
	(*U)(All,compteur)= *u;
	(*X)(All,compteur)= x;     
      }
    
      /* break if we need to */    
  
      iret |= (this->*test_function)();
      if(iret)  break;      
    }
    
#if _DEBUG >= 2
    cout << "Writing trajectory ..." << endl;
#endif
    
    /* resize */   
    length = compteur+1;
    (*X).resizeAndPreserve(N,compteur+1);
    (*U).resizeAndPreserve(dimu,compteur+1);
    (*time).resizeAndPreserve(compteur+1);    	
  }
  
  else { 
    
   if ((nb_points == 3)&&(tspan(1)==tspan(2))) {      
      tspan.resizeAndPreserve(2);
      nb_points = 2;
    }
   length = nb_points;
   (*X).resize(N,nb_points);
   (*U).resize(dimu,nb_points);
   (*time).resize(nb_points);
   (*time) = tspan;
    
#if _DEBUG >= 2
    cout << "Loop over trajectories..." << endl;
#endif
               
    x = (*p0)(Rx); // new x0
    p = (*p0)(All);
 
    UpdateFdata(t0,y,f_data,0,NULL);
            
    (*X)(All,0)= x;
    GetU(f_data, u_data);
    (*U)(All,0)=*u;
      
#ifdef _TD
    x(N)=0.;
#endif	
    iret = FALSE;
    
    /* Reinit solver */
    
    switch (itol) {
    case CV_SS:
      status = CVodeReInit(cvode_mem, f, t0, y, itol, reltol, &Sabstol);      
      break;
    case CV_SV:
      status = CVodeReInit(cvode_mem, f, t0, y, itol, reltol, NV_abstol);
      break;
    }
    
    /* Init root finding for hybrid */
    
    if (data->dimg>0) {
      CVodeRootInit(cvode_mem, data->dimg, g, f_data);
    }
    
    for(i=1; i<nb_points;i++) {
      
      iret = FALSE;
      tout=tspan(i);
      
      CVodeSetStopTime(cvode_mem,tout);

      /* Integrate system until tout */
      
      while (1) {
	
	/* Integrate one step */
	status = CVode(cvode_mem, tout, y, &tret, itask);
	
	/* break on CVode error */
	if (status < 0) return;   
	  
	/* Test if tout was reached */
	CVodeGetCurrentStep(cvode_mem, &h);
	
	if ( (tret - tout)*h >= 0.0 ) {
	  tret = tout;
	  CVodeGetDky(cvode_mem, tout, 0, y);
	  iret = TRUE;
	}      
	
	/* Update f_data */
	
	if (UpdateFdata(tret,y,f_data,0,NULL)||status == CV_ROOT_RETURN) { /* if f_data changed in a discontinuity, reinit cvode */
	  CVodeGetDky(cvode_mem, tret, 0, y);
	  
	  switch (itol) {
	  case CV_SS:
	    status = CVodeReInit(cvode_mem, f, tret,y, itol, reltol, &Sabstol);      
	    break;
	  case CV_SV:
	    status = CVodeReInit(cvode_mem, f, tret,y, itol, reltol, NV_abstol);
	  }
	  CVodeSetInitStep(cvode_mem,h);
	}
	  
	/* break if we need to */
	if(iret)  break;      
      }
	
      GetU(f_data, u_data);
      (*U)(All,i)= *u;
      (*X)(All,i)= x;	
    }       
  }

#if _DEBUG >= 1
  cout << "Leaving ComputeTraj." << endl;
#endif
}

int trajectory::test1() {

  Fdata* data = (Fdata*) f_data;
  N_Vector xdotv = N_VNew_Serial(N);
  int res= 1;
  
  realtype tcur;
  CVodeGetCurrentTime(cvode_mem,&tcur);
  CVodeGetDky(cvode_mem, tcur, 1, xdotv);

  for(int i = 0; i<N; i++)  {
    switch (itol) {
    case CV_SS:
      res*=(abs(Ith(xdotv,i)) < Sabstol);
      break;
    case CV_SV:
      res*=(abs(Ith(xdotv,i)) < Ith(NV_abstol,i));
      break;
    }
    if (!res) {
      N_VDestroy(xdotv);
      return 0;
    }
  }
  //cout << "xdotv: "<< Ith(xdotv, 0) << "  " << Ith(xdotv, 1) << "  " << Ith(xdotv, 2) << "  res: " << res << endl;
  status = TRAJ_REACHED_SSTATE;
  N_VDestroy(xdotv);

  return res;
};
