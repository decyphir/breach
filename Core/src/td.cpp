#include "td.h"

static td_options * opt;
static td_stats * stats;

/******************* second mex subgate functions *********************/

int CVM_TD(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  opt = new td_options();
  stats = new td_stats();

  /* Reading inputs */

  Array2D pts(N,1);
  GetArrayField(pts, prhs[0],"pts");
  Array1D tspan(2);
  GetArray1D(tspan,prhs[1]);
  mxArray* mxDelta;

  if (nrhs == 3) {
    opt->lambda = GetDoubleField(prhs[2],"lambda");
    opt->eta = GetDoubleField(prhs[2],"eta");

    if (mxDelta = mxGetField(prhs[2],0,"delta"))
      GetArray1D(opt->delta,mxDelta);
    //cout << opt->delta<< endl; 
  }

  opt->gamma = exp(-S_G*tspan(1));
  
  /* Doing some stuff */
  
  InitFdata(f_data, cvm_Mdata->mx_data);
  TD( tspan,  pts);
  
  CmacVF* C = ((Fdata*) f_data)->C;
  double mem_used=C->cmac_memory_usage() ;

  mxArray * mxwmap;
  C->serialize_wmap(mxwmap);

  /* Writing outputs */ 

  plhs[0] = mxwmap;
  
  const char * fieldnames[4] = {"Av_Perf", "Av_BellRes", "Av_nb_new","mem_used"};
  mxArray* mxStats = mxCreateStructMatrix(1,1,4, fieldnames);

  mxSetField(mxStats, 0, "Av_Perf", mxCreateDoubleScalar(stats->Av_Perf));
  mxSetField(mxStats, 0, "Av_BellRes",mxCreateDoubleScalar(stats->Av_BellRes));
  mxSetField(mxStats, 0, "Av_nb_new",mxCreateDoubleScalar((double) stats->Av_nb_new));
  mxSetField(mxStats, 0, "mem_used",mxCreateDoubleScalar(mem_used));
  
  plhs[1] = mxStats;
  
}
  

int CVM_ComputeCosts(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  //  cout << "Entering CVM_ComputeCosts ..." << endl;

  opt = new td_options();
  stats = new td_stats();

  /* Reading inputs */
  
  Array2D pts(N,1);
  GetArrayField(pts, prhs[0],"pts");
  Array1D tspan(2);
  GetArray1D(tspan,prhs[1]);
  
  if (nrhs == 3) {
    opt->lambda = GetDoubleField(prhs[2],"lambda");
    opt->eta = GetDoubleField(prhs[2],"eta");
  }
  opt->gamma = exp(-S_G*tspan(1));

  /* Doing some stuff */
  
  Array2D Xf(N,1);
  Array<trajectory*,1> trajArray(1);
  Array1D Costs(1);
  ComputeCosts(tspan,pts,Xf,Costs,trajArray);

  /* Writing outputs */ 
  
  plhs[0] = mxDuplicateArray(prhs[0]);

  mxArray * mxTraj = Traj2mxStruct(trajArray);
  SetField(Xf,plhs[0],"Xf");
  SetField(Costs,plhs[0],"Costs");
  mxAddField(plhs[0],"traj");
  mxSetField(plhs[0],0,"traj",mxTraj);  

  const char * fieldnames[3] = {"Av_Perf", "Av_BellRes", "Av_nb_new"};
  mxArray* mxStats = mxCreateStructMatrix(1,1,3, fieldnames);

  mxSetField(mxStats, 0, "Av_Perf", mxCreateDoubleScalar(stats->Av_Perf));
  mxSetField(mxStats, 0, "Av_BellRes",mxCreateDoubleScalar(stats->Av_BellRes));
  mxSetField(mxStats, 0, "Av_nb_new",mxCreateDoubleScalar((double) stats->Av_nb_new));
  
  plhs[1] = mxStats;

  //  cout << "traj:" << *trajArray(0)->time << endl << *trajArray(0)->X << endl << *trajArray(0)->U << endl << *trajArray(0)->iV << endl;
  //  cout << "Leaving CVM_ComputeCosts ..." << endl;

}


/****************************** Implementation ************************************/

void ComputeCosts(Array1D tspan, Array2D pts,Array2D &Xf, Array1D &Costs,Array<trajectory*,1>  &trajArray) {
#if _DEBUG >= 1  
  cout << "Entering ComputeCosts ..." <<endl;
#endif
  
  InitFdata(f_data, cvm_Mdata->mx_data);
  int nb_traj = pts.extent(1);
  Xf.resize(N,nb_traj);
  Xf = 0;
  Costs.resize(nb_traj);
  trajArray.resize(nb_traj);

  /* CVodes Reinitialization */

  int status, statusgdk;
  int itask = CV_ONE_STEP_TSTOP;
  int dimu;
  double tlast,tret,h;
  booleantype iret = FALSE;

  Fdata* data = (Fdata*) f_data;
  CmacVF *C = (data)->C;

  realtype *xdata = N_VGetArrayPointer(y);
  Array1D x(xdata,shape(N),neverDeleteData);
  Array1D p(data->p,shape(data->dimp),neverDeleteData);

  double t0 = tspan(0);
  double tout;
  int  nb_points = tspan.extent(0);

  /* control */

  void* u_data;
  u_data = (Array1D*) new Array1D(1);  
  GetU(f_data, u_data);
  Array1D * u = (Array1D*) u_data; 
  dimu = (*u).extent(0);

  /* working variables */
  
  Range Rx(0,DIMX-1);
  Range All = Range::all();   

  int i,j;
  
  Array2D *X;
  Array2D *U;
  Array1D *time;
  Array1D *iV;
  int dimpts = pts.extent(0);


#ifdef _HYBRID

  double clok = 0.;
  double h_used;
  double h_min = data->h_min;
  double h_max = data->h_max;
  double backtrack;

#endif

  for (j=0;j<nb_traj;j++) {
    
    x = pts(Rx,j); // new x0
    p = pts(All,j);


#ifdef TD
      x(DIMX)=0.;
#endif
    iret = FALSE;
       
    tret = t0;
    //   UpdateFdata(tret,y,f_data,0,NULL);

#ifdef _HYBRID
    clok = 0.;
#endif

    /* Init output arrays  */
    
    trajArray(j) = new trajectory();
    X = (trajArray(j)->X);
    U = (trajArray(j)->U);
    time = (trajArray(j)->time);
    iV = (trajArray(j)->iV); 
    (*time).resize(nb_points);
    (*time) = tspan;
    (*X).resize(N,nb_points);
    (*X)(Range::all(),0) = x;
    (*U).resize(dimu,nb_points);
    (*U)(Range::all(),0)=*u;
    (*iV).resize(nb_points);
    (*iV)(0) = C->response(x(Range(0,DIMX-1)));
    
    /* CVodes Reinitialization */
    
    switch (itol) {
    case CV_SS:
      status = CVodeReInit(cvode_mem, f, t0, y, itol, reltol, &Sabstol);
      break;
    case CV_SV:
      status = CVodeReInit(cvode_mem, f, t0, y, itol, reltol, NV_abstol);
    }
    
    UpdateFdata(tret,y,f_data,0,NULL);
       
    if (status < 0) break;
    
    for(i=1; i<nb_points;i++) {
      
      iret = FALSE;
      tout=tspan(i);

      /* Integrate system until tout */
      
      while (1) {
	
	/* Integrate one step */
	tlast = tret;     
	CVodeSetStopTime(cvode_mem,tout);	

	status = CVode(cvode_mem, tout, y, &tret, itask);

	/* break on CVode error */
	if (status < 0) {
	  cout << "CVODES failed for some unknown reason. Status= " << status << endl;
	  mexErrMsgTxt("Dying...");
	  break;
	}

	/* Test if tout was reached */
	CVodeGetCurrentStep(cvode_mem, &h);

	if ( (tret - tout)*h >= 0.0 ) {
	  tret = tout;
	  CVodeGetDky(cvode_mem, tout, 0, y);
	  iret = TRUE;
	}
	  
	/* Update f_data */

#ifdef _HYBRID

	h_used = tret-tlast;
	clok += h_used;

	if (clok>= h_min) {  // check if spent enough in current mode 

	  if (UpdateFdata(tret,y,f_data,0,NULL)) { // Update f_data. if changed at a discontinuity, reinit cvode
	    
	    backtrack = max(h_used-h_max,0.); 
	    tret = tret - backtrack;

	    if (backtrack>0) { // if backtrack, need to retest if tout was reached
	      iret = FALSE;
	      if ( (tret - tout)*h >= 0.0 ) {
		tret = tout;
		CVodeGetDky(cvode_mem, tout, 0, y);
		iret = TRUE;
	      }      
	    }

	     statusgdk = CVodeGetDky(cvode_mem, tret, 0, y);

	     if (statusgdk) {

	       cout << "statusgdk: " << statusgdk << endl;
	       cout << "needed backtrack:" << backtrack << endl;
	       cout << "Pb backtracking"<< endl;	    
	       cout << "x:" << x << endl;
	       cout << "tret:" << tret << endl;
	       cout << "h_used:" << h_used << endl;
	       mexErrMsgTxt("CVODES failed for some reason.");

	     }
	     
	     switch (itol) {
	     case CV_SS:
	       CVodeReInit(cvode_mem, f, tret,y, itol, reltol, &Sabstol);
	       break;
	     case CV_SV:
	       CVodeReInit(cvode_mem, f, tret,y, itol, reltol, NV_abstol);
	       break;
	    }
	     CVodeSetInitStep(cvode_mem,h);
	     
	     clok = 0;
	  }
	  
	}

#else 	

	if (UpdateFdata(tret,y,f_data,0,NULL)) {
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

#endif
	
	/* break if we need to */

	if(iret)  break;

      }

      GetU(f_data, u_data);
      (*U)(Range::all(),i)= *u;
      (*X)(Range::all(),i)= x;
      (*iV)(i) = C->response(x(Range(0,N-2)));

    }

   Xf(Range::all(),j) = x;
   Costs(j)=x(N-1)+exp(-S_G*tout)*(*iV)(i);

    /* Eval Bellman residual for this trajectory */
    
    trajArray(j)->EvalBellRes();
    if (opt->lambda == -1.) 
      trajArray(j)->EvalTDUpdate();
    else
      trajArray(j)->EvalTDlambdaUpdate();

    stats->Av_Perf += trajArray(j)->perf;
    stats->Av_BellRes += abs((*trajArray(j)->BellRes)(0));

    //  trajArray(j)->affiche();
  }

  stats->Av_Perf = stats->Av_Perf/nb_traj;
  stats->Av_BellRes = stats->Av_BellRes/nb_traj;   
  cout << "Perf moyenne: " << stats->Av_Perf;
  cout << "  Bellman residuel moyen : " << stats->Av_BellRes << endl;
  
#if _DEBUG >= 1
  cout << "Leaving ComputeCosts." << endl;
#endif
}

void ComputeCosts(Array1D tspan, Array1D x0, trajectory & traj) { 
  
#if _DEBUG >= 1
  cout << "Entering ComputeCosts ..." <<endl;
#endif

  /* CVodes Reinitialization */
  
  int status;
  int itask = CV_ONE_STEP_TSTOP;
  int dimu;
  int dimx = N;
  double tret,h;
  booleantype iret = FALSE;


  Fdata* data = (Fdata*) f_data;
  CmacVF *C = (data)->C;

  realtype *xdata = N_VGetArrayPointer(y);
  Array1D x(xdata,shape(N),neverDeleteData);
  Array1D p(data->p,shape(data->dimp),neverDeleteData);


  x = x0;

  double t0 = tspan(0);
  tret = t0;
  double tlast,tout;
  int nb_points = tspan.extent(0);

 /* control */

  void* u_data;
  u_data = (Array1D*) new Array1D(1);  
  GetU(f_data, u_data);
  Array1D * u = (Array1D*) u_data; 

  /* working variables */

  int i,j;    
  iret = FALSE;
  int statusgdk;

#ifdef _HYBRID
  
  double clok = 0.;
  double h_used;
  double h_min = data->h_min;
  double h_max = data->h_max;
  double backtrack;

#endif

 /* Init output arrays  */

  (*(traj.X))(Range::all(),0) = x;
  (*(traj.U))(Range::all(),0)=*u;
  (*(traj.iV))(0) = C->response(x(Range(0,N-2)));
  
  //  cout << *(traj.time) << endl;
  //  cout << *(traj.X) << endl;
  //  cout << *(traj.U) << endl;
  //  cout << *(traj.iV) << endl;
  //  cout << "nb_points:" << nb_points << endl;

  /* CVodes Reinitialization */

  switch (itol) {
  case CV_SS:
    status = CVodeReInit(cvode_mem, f, t0, y, itol, reltol, &Sabstol);
    break;
  case CV_SV:
    status = CVodeReInit(cvode_mem, f, t0, y, itol, reltol, NV_abstol);
  }

  if (status < 0) {        
    cout << "CVODES failed for some unknown reason. Status= " << status << endl;
    mexErrMsgTxt("Dying...");
    return;
  }

  for(i=1; i<nb_points;i++) {
      
    iret = FALSE;
    tout=tspan(i);

    /* Integrate system until tout */

    while (1) {
      
      /* Integrate one step */
      CVodeSetStopTime(cvode_mem,tout);
      tlast = tret;
      status = CVode(cvode_mem, tout, y, &tret, itask);
	
      /* break on CVode error */
      if (status < 0) break;
      
      /* Test if tout was reached */
      CVodeGetCurrentStep(cvode_mem, &h);
      
      if ( (tret - tout)*h >= 0.0 ) {
	tret = tout;
	CVodeGetDky(cvode_mem, tout, 0, y);
	iret = TRUE;
      }
	  
	/* Update f_data */

#ifdef _HYBRID

	//	CVodeGetLastStep(cvode_mem,&h_used);

	h_used = tret-tlast;
	clok += h_used;

	if (clok>= h_min) {  // check if spent enough in current mode 

	  if (UpdateFdata(tret,y,f_data,0,NULL)) { // Update f_data. if changed at a discontinuity, reinit cvode
	    
	    //	    cout << "Switch !" << endl;
	     
	    backtrack = max(h_used-h_max,0.); 
	    //	    cout << "need backtrack:" << backtrack << endl;
	   
	    tret = tret - backtrack;

	    if (backtrack>0) { // if backtrack, need to retest if tout was reached
	      iret = FALSE;
	      if ( (tret - tout)*h >= 0.0 ) {
		tret = tout;
		CVodeGetDky(cvode_mem, tout, 0, y);
		iret = TRUE;
	      }      
	    }

	     statusgdk = CVodeGetDky(cvode_mem, tret, 0, y);
	     if (statusgdk) {

	       cout << "statusgdk: " << statusgdk << endl;
	       cout << "needed backtrack:" << backtrack << endl;
	       cout << "Pb backtracking"<< endl;	    
	       cout << "x:" << x << endl;
	       cout << "tret:" << tret << endl;
	       cout << "h_used:" << h_used << endl;

	       mexErrMsgTxt("CVODES failed for some reason.");
	     }

	    switch (itol) {
	    case CV_SS:
	      CVodeReInit(cvode_mem, f, tret,y, itol, reltol, &Sabstol);
	    break;
	    case CV_SV:
	      CVodeReInit(cvode_mem, f, tret,y, itol, reltol, NV_abstol);
	      break;
	    }
	    CVodeSetInitStep(cvode_mem,h);
	    clok = 0;

	  }	  	  
	}

#else

      /* Update f_data */
	
      if (UpdateFdata(tret,y,f_data,0,NULL)) {
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

#endif      


      /* break if we need to */

      if(iret)  break;
      
    }
 
    GetU(f_data, u_data);

    (*(traj.U))(Range::all(),i)= *u;
    (*(traj.X))(Range::all(),i)= x;
    (*(traj.iV))(i) = C->response(x(Range(0,N-2)));
  }

#if _DEBUG >= 1 
  cout << "Leaving ComputeCosts." << endl;  
#endif

}

void TD(Array1D tspan, Array2D pts) {

#if _DEBUG >= 1
  cout << "Entering TD ... " << endl;
#endif

  InitFdata(f_data, cvm_Mdata->mx_data);

  trajectory traj=trajectory();
  int nb_traj = pts.extent(1);
  int nb_points = tspan.extent(0);
  double dt;

  //  cout << "nb_points: "<< nb_points << endl;

  (*(traj.time)).resize(nb_points);
  (*(traj.time)) = tspan;  
  dt = tspan(1)-tspan(0);
  (*(traj.X)).resize(((Fdata*) f_data)->dimx,nb_points);
  (*(traj.U)).resize(((Fdata*) f_data)->dimu,nb_points);
  (*(traj.iV)).resize(nb_points);
  Array1D x(N);

  CmacVF* C = ((Fdata*) f_data)->C;

  Array2D X;
  Array1D V;
  Array<int,1> Xind;
  
  int i,j, nb_new, nb_new_moy(0);

  Range Rx(0,DIMX-1);

  if (opt->lambda == -1.) {
    //   cout << "TD(0)" << endl;
   for ( j =0 ; j<nb_traj; j++) {

     x = pts(Rx,j); // new x0

     ComputeCosts(tspan,x,traj);
     traj.EvalBellRes();
     traj.EvalTDUpdate();

     //     cout << "Vtilde:" << *(traj.Vtilde) << endl;
     //     cout << "X:" << (*traj.X)(Range(0,N-2),Range(0,nb_points-2)) << endl;
     //     C->trainAAt((*traj.X)(Range(0,N-2),Range(0,nb_points-2)),*(traj.Vtilde));

     C->train((*traj.X)(Range(0,N-2),Range(0,nb_points-2)),*(traj.Vtilde));
     stats->Av_Perf += traj.perf;
     stats->Av_BellRes += abs((*traj.BellRes)(0))/dt;

   }
  }
  
  else {

    //    cout << "TD(" << opt->lambda << ")" << endl;

    for (j =0 ; j<nb_traj; j++) {

      x = pts(Rx,j); // new x0

      ComputeCosts(tspan,x,traj);

      traj.EvalBellRes();
      traj.EvalTDlambdaUpdate();
      
      C->disperse_pts((*traj.X)(Range(0,N-2),Range(0,nb_points-2)), Xind);
      nb_new = Xind.extent(0);

      //    cout << "nb_new:" << nb_new<< endl;
      //      cout << Xind << endl;

      X.resize(N-1,nb_new);
      V.resize(nb_new);

      for (i=0;i<nb_new;i++) {
	X(Range(0,N-2),i) = (*traj.X)(Range(0,N-2),Xind(i));
	V(i) = (*(traj.Vtilde))(Xind(i));
      }
      //      cout << "X:"<< X << endl;

      C->train(X,V);

      stats->Av_nb_new += nb_new;
      stats->Av_Perf += traj.perf;
      stats->Av_BellRes += abs((*traj.BellRes)(0))/dt;      
    }
  }

  stats->Av_nb_new = stats->Av_nb_new/nb_traj;
  stats->Av_Perf = stats->Av_Perf/nb_traj;
  stats->Av_BellRes = stats->Av_BellRes/(nb_traj);

  //  cout << " Av. nb. new pts: " << stats->Av_nb_new;
  //  cout << " Av. Perf: " << stats->Av_Perf;
  //  cout << " Av. Bellman residual: " << stats->Av_BellRes << endl;  

}

mxArray* Traj2mxStruct(Array<trajectory*,1> trajArray) {

  //      cout << "Entering Traj2mxStruct ..." << endl;
  
  int nb_traj = trajArray.extent(0);  
  mxArray * mxTraj;
  const char * fieldnames[7] = {"perf","time", "X", "U","iV","BellRes","Vtilde"};
  
  mxTraj = mxCreateStructMatrix(1,nb_traj,7,fieldnames);
  mxArray *mxPerf, *mxX, *mxU, *mxtime, *mxiV, *mxBellRes, *mxVtilde;
  
  for (int j = 0; j<nb_traj ; j++) {
 
    mxPerf = mxCreateDoubleScalar(trajArray(j)->perf);
    SetArray(*trajArray(j)->X,mxX);
    SetArray(*trajArray(j)->U,mxU);
    SetArray(*trajArray(j)->time,mxtime);
    SetArray(*trajArray(j)->iV,mxiV);
    SetArray(*trajArray(j)->BellRes,mxBellRes);
    SetArray(*trajArray(j)->Vtilde,mxVtilde);
    
    mxSetField(mxTraj, j, "perf",mxPerf);
    mxSetField(mxTraj, j, "time",mxtime);
    mxSetField(mxTraj, j, "X",mxX);
    mxSetField(mxTraj, j, "U",mxU);
    mxSetField(mxTraj, j, "iV",mxiV);
    mxSetField(mxTraj, j, "BellRes",mxBellRes);
    mxSetField(mxTraj, j, "Vtilde",mxVtilde);

  }
  
  // cout << "Leaving Traj2mxStruct ... " << endl;

  return mxTraj;
  
 }

void trajectory::EvalBellRes() {

  // cout << "Entering eval_Bell_res ..." << endl;

  int nb_points = time->extent(0);
  (*BellRes).resize(nb_points);
  
  int i=0;
  double c,V,Vn,delta,t,tn;

  t = (*time)(0);
  tn = (*time)(1);

  double gamma = exp(-S_G*tn);
  double inv_g;

  c= (*X)(N-1,1);
  V = (*iV)(0);
  Vn = (*iV)(1);

  //  cout << "c: " << c << " gamma: " << gamma << " V:" << V << " Vn:" << Vn << endl; 
  
  (*BellRes)(0) = c+gamma*Vn-V;

  for(i=1;i<nb_points-1;i++) {

    gamma = exp(-S_G*(tn-t));
    inv_g = exp(S_G*t);
    t = tn;
    tn = (*time)(i+1);

    c = inv_g*((*X)(N-1,i+1)-(*X)(N-1,i));
    V = (*iV)(i);
    Vn = (*iV)(i+1);
   
    (*BellRes)(i) = c+gamma*Vn-V;

  }

  (*BellRes)(nb_points-1) = (*BellRes)(nb_points-2);
  //  cout << "Leaving EvalBellRes ..." << endl;

}

void trajectory::EvalTDUpdate() {

  //  cout << "Entering EvalTDUpdate ..." << endl;

 int nb_points = time->extent(0);

 // cout << "time:" << *time <<endl;
 // cout << "iV:" << *iV <<endl;
 // cout << "BellRes:" << *BellRes <<endl;

 (*Vtilde).resize(nb_points-1);  

 for (int i=0; i< nb_points-1 ; i++) {
   (*Vtilde)(i) = (*iV)(i)+.1*(*BellRes)(i);
 }

 // cout << "Vtilde:" << *iV <<endl;
 // cout << "Leaving EvalTDUpdate ..." << endl;
 
}

void trajectory::EvalTDlambdaUpdate() {

  //    cout << "Entering eval_TDlambda_update ..." << endl;
  
  int i;
  double xvTf,V;
  int nb_points = time->extent(0);
 (*Vtilde).resize(nb_points-1);  

  V = (*iV)(nb_points-1);
  xvTf = (*X)(N-1,nb_points-1);   
  double T = (*time)(nb_points-1);

  perf = xvTf+exp(-S_G*T)*(*iV)(nb_points-1); // Performance optimiste ...
          
  // cout << "************************************" << endl;
  
  //  TRACE("costs:", costs);
  //  TRACE("x0:", traj(Range(1,N),0));
  //  TRACE("xf:", traj(Range(1,N),iTf));
  //  TRACE("perf:", perf);
  //  cout << "perf: " << perf;
  //  cout << "  xvTf: " << xvTf;
  //  cout << "  exp(-S_G*T): "<< exp(-S_G*T); 
  //  cout << "  iV: " << (*iV)(nb_points-1) << endl;
 
  
  //  TRACE("iV(x0):", iV);
  //  cout << "Bellman_res:", Bellman_res);
  //  TRACE("iTf", iTf);
  
  //  cout << "Updates:" << endl;

  i = nb_points-2;
  double e = (*BellRes)(nb_points-2);
  
  while (i>=0) {      
    
    (*Vtilde)(i) = V+opt->eta*e;
    (*Vtilde)(i) = max(0.,(*Vtilde)(i));
    (*Vtilde)(i) = min(MAXCOST,(*Vtilde)(i));
    i--;      
    V = (*iV)(i);
    e *= (opt->lambda)*(opt->gamma);
    e += (*BellRes)(i);

  }

  //  cout << "Leaving EvalTDlambdaUpdate ..." << endl;
  
}

void trajectory::affiche()  {
    cout << "time:" << *time << endl;
    cout << "X:" << *X << endl;
    cout << "U:" << *U << endl;
    cout << "Expa:" << *Expa << endl;
    cout << "iV:" << *iV << endl;
    cout << "Vtilde:" << *Vtilde << endl;
  }

ostream& operator<<(ostream& os, const trajectory& traj){

  os << "time:" << *(traj.time) << endl;
  os << "X:" << *(traj.X) << endl;
  os << "U:" << *(traj.U)<< endl;
  os << "Expa:" << *(traj.Expa) << endl;
  os << "iV:" << *(traj.iV) << endl;
  os << "Vtilde:" << *(traj.Vtilde) << endl;
  return os;

};


void trajectory::write_mxTraj(mxArray* & mxTraj) {
  
  const char * fieldnames[3] = {"time", "X", "U"};
  mxTraj = mxCreateStructMatrix(1,1,3, fieldnames);
  mxArray * mxX,*mxU, *mxtime;

  SetArray(*X,mxX);
  SetArray(*U,mxU);
  SetArray(*time,mxtime);
  mxSetField(mxTraj, 0, "time",mxtime);
  mxSetField(mxTraj, 0, "X",mxX);
  mxSetField(mxTraj, 0, "U",mxU);

}

