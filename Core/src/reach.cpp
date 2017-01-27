#include "reach.h"


int CVM_CmacComputeExpTraj(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

#if _DEBUG >= 1
    cout << "Entering CVM_ComputeExpTraj ..." << endl;
#endif

  /* Reading */

  CmacVF *C;
  ReadCmac(prhs[0],C);

  Array2D pts(N,1);
  GetArrayField(pts, prhs[1],"pts");

  Array2D XS0(N*Ns,1);
  GetArrayField(XS0, prhs[1],"XS0");

  Array2D epsi(N,1);
  GetArrayField(epsi, prhs[1],"epsi");

  Array1D tspan(2);
  GetArray1D(tspan,prhs[2]);

  /* Doing some stuff */

  Array<trajectory*,1> trajArray(1);
  Array2D Xf(1,1);
  Array2D XSf(1,1);
  Array2D ExpaMax(1,1);
  Array1D Keep(1);

  ComputeExpTraj(C,pts, XS0, epsi, tspan, trajArray, Xf, XSf, ExpaMax, Keep);
 
  /* Writing outputs */ 

  int imem_used;
  imem_used = C->cmac_memory_usage();
  double mem_used= (double) imem_used;
  plhs[2] = mxCreateDoubleScalar(mem_used);

  mxArray * mxwmap;
  C->serialize_wmap(mxwmap);    
  plhs[0] = mxwmap;

  plhs[1] = mxDuplicateArray(prhs[1]);
  mxArray * mxTraj;
  
  int nb_pts = trajArray.extent(0);  
  const char * fieldnames[4] = {"time", "X", "Expa", "U"};

  mxTraj = mxCreateStructMatrix(1,nb_pts,4, fieldnames);
  mxArray * mxX, *mxExpa, *mxU, *mxtime, *mxKeep;

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
  
  SetField(Xf,plhs[1],"Xf");
  
  mxAddField(plhs[1],"ExpaMax");
  SetField(ExpaMax,plhs[1],"ExpaMax");

  mxAddField(plhs[1],"XSf");
  SetField(XSf,plhs[1],"XSf");

  mxAddField(plhs[1],"Keep");
  SetArray(Keep,mxKeep);
  mxSetField(plhs[1],0,"Keep",mxKeep);

  mxAddField(plhs[1],"traj");
  mxSetField(plhs[1],0,"traj",mxTraj);

  /* Free memory */ 

  for (int j = 0; j<nb_pts ; j++) {
    delete trajArray(j);
  }


#if DEBUG >= 1
  cout << "Leaving CVM_CmacComputeExpTraj ..." << endl;
#endif
}

int CVM_ComputeTrajs(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

#if _DEBUG >= 1
  cout << "Entering CVM_ComputeTrajs ..." << endl;
#endif

  /* Reading */
  CmacVF *C;
  ReadCmac(prhs[0],C);

  Array2D pts(N,1);
  GetArrayField(pts, prhs[1],"pts");

  Array1D tspan(2);
  GetArray1D(tspan,prhs[2]);

  /* Doing some stuff */
  Array<trajectory*,1> trajArray(1);
  Array2D Xf(1,1);
  Array1D Keep(1);

  ComputeTrajs(C, pts, tspan, trajArray, Xf, Keep);
  
  /* Writing outputs */ 

#if _DEBUG >= 1
  cout << "Writing outputs ..." << endl;
#endif

  double mem_used= (double) C->cmac_memory_usage();
  plhs[2] = mxCreateDoubleScalar(mem_used);

  mxArray * mxwmap;
  C->serialize_wmap(mxwmap);    
  plhs[0] = mxwmap;
  plhs[1] = mxDuplicateArray(prhs[1]);

  mxArray * mxTraj;
  
  int nb_pts = trajArray.extent(0);  
  const char * fieldnames[3] = {"time", "X", "U"};

  mxTraj = mxCreateStructMatrix(1,nb_pts,3, fieldnames);
  mxArray * mxX,*mxU, *mxtime, *mxKeep;

#if _DEBUG >= 1
  cout << "Writing trajectories ..." << endl;
#endif
  Array1D p0(pts.extent(0));
  for (int j = 0; j<nb_pts ; j++) {
    SetArray(*trajArray(j)->X,mxX);
    SetArray(*trajArray(j)->U,mxU);
    SetArray(*trajArray(j)->time,mxtime);

    mxSetField(mxTraj, j, "time",mxtime);
    mxSetField(mxTraj, j, "X",mxX);
    mxSetField(mxTraj, j, "U",mxU);

  }
  
  SetField(Xf,plhs[1],"Xf");
  mxAddField(plhs[1],"Keep");
  SetArray(Keep,mxKeep);
  mxSetField(plhs[1],0,"Keep",mxKeep);
  mxAddField(plhs[1],"traj");
  mxSetField(plhs[1],0,"traj",mxTraj);

  /* Free memory */ 

  for (int j = 0; j<nb_pts ; j++) {
    delete trajArray(j);
  }

  
#if _DEBUG >= 1
  cout << "Leaving CVM_ComputeTrajs ..." << endl;
#endif
}

int CVM_Reach(int nlhs,  mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  
#if _DEBUG >= 1
  cout << "Entering CVM_Reach ..." << endl;
#endif
   int NumOfChildren = 3;
   traj_set * P = new traj_set(prhs[0],NumOfChildren);
   GetArray1D(P->tspan, prhs[1]);
   
   const mxArray * mxOptions = prhs[2];

   double * expa_test = mxGetPr(mxGetField(mxOptions,0, "ExpaTest"));
   int MaxIter = (int) GetDoubleField(mxOptions, "MaxIter");
   int iter =0;
   bool done = false;

   P->print_refine_history();
   while(done == false) {
     cout << "iteration " << iter << endl;
     //     done = P->refine_traj(NumOfChildren,&traj_set::test_expansion,expa_test);
     done = P->refine_traj(NumOfChildren,&traj_set::test_nonlin_error,expa_test);
     // done = P->refine_all(NumOfChildren,&param_set::Test1,expa_test);
     if (++iter>MaxIter)
       break;
   }
   //  cout << " Refine history " << endl;
   //  P->print_refine_history();
   
   P->write_mxArray(plhs[0]);
   
   delete P; 

   if (!done)
     cout << "didn't converge" << endl;
#if _DEBUG >= 1
   cout << "Leaving CVM_Reach ..." << endl;
#endif
   
   return 1;

}

int CVM_CmacComputeTraj(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

#if _DEBUG >= 1
  cout << "Entering CVM_ComputeTraj ..." << endl;
#endif

 /* Reading */

  CmacVF *C;
  ReadCmac(prhs[0],C);
  
  Array2D pts(N,1);
  GetArrayField(pts, prhs[1],"pts");

  Array1D tspan(2);
  GetArray1D(tspan,prhs[2]);

  int stop_cycle=0;  
  if (nrhs>3) {
    stop_cycle = (int) mxGetScalar(prhs[3]);
  }
  
  /* Doing some stuff */

  Range All = Range::all();   
  int nb_pts = pts.extent(1);
  Array<trajectory*,1> trajArray(nb_pts);
  Array2D Xf(N,nb_pts);
  
  Array1D p0(pts.extent(0));
  Array1D Tspan;
  int l;
  double period;
  Array1D periods;
  if (stop_cycle) 
    periods.resize(nb_pts);

  for (int j = 0; j<nb_pts ; j++) {

    p0 = pts(All, j);
    Tspan.resize(tspan.shape());  
    Tspan = tspan;
    trajArray(j) =  new trajectory();
    trajArray(j)->p0 = &p0;
    period = (trajArray(j))->ComputeTraj(C,Tspan, stop_cycle);
    //    cout << "period found:" << period << endl;
    l = trajArray(j)->length; 
    Xf(All,j) = (*trajArray(j)->X)(All,l-1);    
    if (stop_cycle) {
      C->clear_cmac_weights();
      periods(j) = period;
    }
  }
  
  /* Writing trajectories */ 

  plhs[1] = mxDuplicateArray(prhs[1]);
  mxArray * mxTraj, *mxPeriod;
  
  const char * fieldnames[4] = {"param", "time", "X", "U"};

  mxTraj = mxCreateStructMatrix(1,nb_pts,4, fieldnames);
  mxArray * mxX,*mxU, *mxtime, *mxParam;

  int j;
  
  if (stop_cycle)
    mxAddField(mxTraj,"period");
  
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
    
    if (stop_cycle) {
      mxPeriod = mxCreateDoubleScalar(periods(j));    
      mxSetField(mxTraj,j, "period",mxPeriod);    
    }
  }

  SetField(Xf,plhs[1],"Xf");
  mxAddField(plhs[1],"traj");
  mxSetField(plhs[1],0, "traj",mxTraj);
 
  /* Writing cmac */
  
  double mem_used= (double) C->cmac_memory_usage();
  
  plhs[2] = mxCreateDoubleScalar(mem_used);
  
  mxArray * mxwmap;
  C->serialize_wmap(mxwmap);    
  plhs[0] = mxwmap;
  
  /* Free memory */ 
  
  for (j = 0; j<nb_pts ; j++) {
    delete trajArray(j);
  }
  
#if _DEBUG >= 1
  cout << "Leaving CVM_CmacComputeTraj ..." << endl;
#endif
  
}

double trajectory::ComputeTraj(CmacVF* &C, Array1D& tspan, int stop_cycle) {

#if _DEBUG >= 1
  cout << "Entering traj::ComputeTraj ..." <<endl;
#endif
 
  /* CVodes Reinitialization */
 
#if _DEBUG >= 2
  cout << "Init CVodes..." << endl;
#endif  
  
  InitFdata(f_data, cvm_Mdata->mx_data);
  int cvstatus;
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
  double tout, T;

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
      
 
  /* rootfinding */
  
  for (int i=0; i<N; i++) {
    data->delta[i]=1.1*C->qnt[i];    
  }
  CVodeRootInit(cvode_mem, data->dimg+N, g_delta, f_data);	  	  
  Array1D xp(data->xp, shape(N)); // x prec  
  xp = x(Range(0,N-1));
      
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
    cvstatus = CVode(cvode_mem, tout, y, &tret, itask);
      
    /* break on CVode error */
    if (cvstatus < 0) {
      cout << "tout:" << tout << endl;
      cout << "x:" << x << endl;
#if _DEBUG >=3
      cout << "time:" << time << endl;
      cout << "X:" << *X << endl;
      cout << "U:" << *U << endl;
#endif
      CVodeGetDky(cvode_mem, tout, 1, y);
      cout << "dx:" << x << endl;
      cout << "CVODES failed miserably for some reason. Status= " << cvstatus << endl;
      mexErrMsgTxt("Dying...");	  
      return -1.;
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
    
    if (UpdateFdata(tret,y,f_data,0,NULL) || cvstatus == CV_ROOT_RETURN) {  //Update f_data. if changed at a discontinuity, reinit cvode 
	    
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
    if ((cvstatus == CV_ROOT_RETURN)||(iret)) {   
      compteur++;
      if (compteur < kmax) {
	
#if _DEBUG >= 3
	cout << " no need. " << endl;
#endif       
	
	(*time)(compteur) = tret;
	(*U)(All,compteur)= *u;
	(*X)(All,compteur)= x;     
	xp = x(Range(0,N-1));
      }

      else  { // try to guess how many points left
	
	kmax += max(1,min((int) ceil((tout-tret)/(tret-tlast)),kmax));
	
	(*X).resizeAndPreserve(N,kmax);
	(*U).resizeAndPreserve(dimu,kmax);
	(*time).resizeAndPreserve(kmax);
	
	(*time)(compteur) = tret;
	(*U)(All,compteur)= *u;
	(*X)(All,compteur)= x;
	xp = x(Range(0,N-1));     
      }
      
      if (!iret) 
	if (stop_cycle &&(!iret)) {
	  T = C->mark_time_visited(x, tret);
	  if (T>0) {
	    iret= TRUE;
	    status = TRAJ_OSCILLATES;
	  }
	}     
    } // END if (status == CV_ROOT_RETURN)
    
    /* break if we need to */
    
    if(iret)  break;      
    //      cout << "loop or not .. " << endl;
  } // END while (1) 
    

#if _DEBUG >= 2
  cout << "Writing trajectory ..." << endl;
#endif
   
  /* resize */   
  length = compteur+1;
  (*X).resizeAndPreserve(N,compteur+1);
  (*U).resizeAndPreserve(dimu,compteur+1);
  (*time).resizeAndPreserve(compteur+1);    	
    
#if _DEBUG >= 1
  cout << "Leaving ComputeTraj." << endl;
#endif
  return T;
}

int traj_set::test_expansion(refine_tree &nod, double expa_test[] ){
    
  int nd = nod->index;
  double * pts_data = Get_pts(nd); 
  Array1D pts = Array1D(pts_data, shape(M_pts));
  //  cout << " pts: " << pts << endl;
  double * epsi_data = Get_ep(nd);
  Array1D eps = Array1D(epsi_data, shape(M_epsi));
  //  cout << "eps:" << eps << endl;

  nod->status=-1;
  trajectory traj;
  Array1D expa_max(N); 
  traj.p0 = &pts;
  traj.ComputeTrajSensi(tspan);
  traj.ComputeExpa(eps, expa_max);
  //  cout << "expa_max : " << expa_max << endl;

  for(int i=0 ; i < expa_max.extent(0) ; i++) {
    if (expa_max(i) > expa_test[i]) {
      int i_max =0;
      double epsi_max = eps(0);
      for(int j = 1; j< M_epsi; j++) {
	if (eps(j)>epsi_max) {
	  i_max = j;
	  epsi_max = eps(j);
	}
      }
      //   cout << "To be refined in " << i_max << ", caus' " << expa_max(i) << " greater than " << expa_test[i] << endl;  
      nod->status=i_max;
      break;
    }
  }
  
  return nod->status;
}

int traj_set::test_nonlin_error(refine_tree &node, double nonlin_err_test[] ) {
  
  //  cout << "Entering test_nonlin_error..." << endl;

  /* This version works for fixed time steps */ 

  /* firstly, compute trajectory */
  
  int nd_local = node->index;    
  double * pt_local_data = Get_pts(nd_local); 
  Array1D pt_local = Array1D(pt_local_data, shape(M_pts));

  Array1D expa_max(N); 
  double * eps_data = Get_ep(nd_local); 
  Array1D eps = Array1D(eps_data, shape(M_epsi));  
  double epsi_max;

  //  cout << "eps: " << eps << endl;

  trajectories[nd_local] = new trajectory(); 

  trajectory *local_traj= trajectories[nd_local];
    
  local_traj->p0 = &pt_local;
  local_traj->ComputeTrajSensi(tspan);
  local_traj->ComputeExpa(eps, expa_max);

  /* then check if the node is a root (return if yes) */

  r_node * global_node = node->father;
   
  if (global_node==NULL) {    
    return 0;//find_greatest_epsi(nd_local);
  }

  int nd_global = global_node->index;
  if (nd_global==nd_local) {
    global_node = global_node->children[2];
    nd_global = global_node->index;
  }

  trajectory* global_traj = trajectories[nd_global];  

  /* Compute local traj estim */
  Range All=Range::all();
  int traj_length= global_traj->length;
  //  cout << " global traj length:" << traj_length << endl;
  Array2D Xlocal_traj_estim(N,traj_length);

  Array1D dp = Array1D(N);
  double *pt_global = Get_pts(nd_global); 
  
  for(int i= 0; i<N; i++)
    dp(i)= pt_local_data[i]-pt_global[i];
  
  Xlocal_traj_estim = *global_traj->X;
  
  //cout << "Global traj: " << Xlocal_traj_estim << endl;

  //cout << "dp: " << dp << endl;
  int i, j, is, k, i_max; 
  //cout << "Global traj XS:"<< (*global_traj->XS) <<endl;

  for(k=0; k< traj_length ; k++) {
    for(is=0; is<Ns; is++)
      for(i=0; i<N; i++) {
	Xlocal_traj_estim(i, k) += ((*global_traj->XS)(i+is*N,k))*dp(is); 
      }
  }
  
  //  cout << "Xlocal_traj_estim: " << Xlocal_traj_estim << endl;
  // cout << "Xlocal_traj:" << *local_traj->X << endl;
  /* Compute Error */
 
  Array2D ErrX = Array2D(Xlocal_traj_estim.shape());
  Array2D ErrXs = Array2D((*global_traj->XS).shape());
  
  ErrX = Xlocal_traj_estim-(*local_traj->X);
  ErrXs = (*global_traj->XS)-(*local_traj->XS);

  //cout << "ErrX: " << ErrX << endl;
  //cout << "ErrXs: " << ErrXs << endl;

  Array2D ExpaErr(N, traj_length); 
  Array1D expa_err(N);
  node->status =-1;

  for(k = 0; k<traj_length; k++){    
    for (i = 0; i< N; i++) {
      expa_err(i) = abs(ErrX(i,k));
      for (is = 0; is <Ns; is++)
	expa_err(i) += abs(ErrXs(i+is*N,k))*eps(is);

      if ((expa_err(i)> nonlin_err_test[i])&&(node->status==-1)) {
	epsi_max =eps(0);
	i_max=0;
	for(j = 1; j< M_epsi; j++) {
	  if (eps(j)>epsi_max) {
	    i_max = j;
	    epsi_max = eps(j);
	  }
	}
	node->status = i_max;	 
	//	cout << "To be refined in " << i_max << ", caus' " << expa_err(i) << " greater than " << nonlin_err_test[i] << endl;  
      }
    }        
 
   ExpaErr(Range::all(),k)=expa_err;
    
  }

  //  cout << "ExpaErr: " << ExpaErr << endl;
  
  return node->status; 

}

bool traj_set::refine_traj(int NumOfChildren, int (traj_set::*ptr_to_test)(refine_tree&, double[]),double tst_eps[]) 
{

  //  cout << "Entering refine_traj" << endl;
  double *data_pts_out,*data_epsi_out,*data_dim_out,*data;

  int i,i_father,k=0,l=0,m,j,ind_pts,ind_ep,Ne=0;
  int tree_ind=NumElems-1; 
  r_node* father;
  double len,temp;  
  
  data = (double *)mxCalloc(NumOfChildren,sizeof(double));
  
  /* test all elements in the queue */

  vector<refine_tree>::reverse_iterator rit;
  vector<refine_tree>::iterator it;

  int stat;
  // cout << "queue size: " << queue.size() << endl; 
  for(rit= queue.rbegin();rit<queue.rend();++rit)
    {
      //      cout << "testing node:" << (*rit)->index << endl;
      (*rit)->status = (this->*ptr_to_test)(*rit,tst_eps);	  
      if((*rit)->status != -1) {
	//	cout << "ah, un pt à refiner: " << (*it)->status << endl; 
	Ne++;
      }
    }

  if (Ne==0) 
    return true;

  int NumElemsOld = NumElems;
  
  NumElems = (NumElemsOld + Ne*(NumOfChildren - NumOfChildren%2)); //New pts will contain these many pts 
  //cout << "New numelem: "  << NumElems << endl;
  
  data_pts_out = (double *)mxCalloc(M_pts*NumElems,sizeof(double));
  data_epsi_out = (double *)mxCalloc(M_epsi*NumElems,sizeof(double));
  data_dim_out = (double *)mxCalloc(M_epsi*1,sizeof(double));
  trajectories.resize(NumElems);

  /* copy previous pts at the beginning of new_pts */
  for (; k< M_pts*NumElemsOld ; k++) 
    data_pts_out[k] = pts[k];
  
  //cout << endl << " current k: " << k << endl;
  for (; l< M_epsi*NumElemsOld ; l++) 
    data_epsi_out[l] = epsi[l];
	  
  //cout << endl << " current l: " << l << endl;
  
  //cout << endl << " current tree_ind: " << tree_ind << endl;
	       
  /* check and refine nodes that need to be refined */
		
  vector<refine_tree> new_queue;

  for(it= queue.begin();it!=queue.end();++it)
    {
      father = *it;
      //      cout << "processing node:" << (*it)->index << endl;
      /* test if the node needs to be refined (check status) */
	    
      ind_ep = father->status;
      if(ind_ep != -1) // do *nothing* if status == -1
	{	      
	  father->isleaf = false;		
	  father->is_divided_in = (int)dim[ind_ep];
	  ind_pts = ((int)dim[ind_ep] - 1);
	  i_father = father->index;
	  father->val = *(pts+i_father*M_pts+ind_pts);
	  
	  temp = (double)(NumOfChildren/2.0);
	  len = *(epsi+i_father*M_epsi+ind_ep)/temp;
	  m = (NumOfChildren/2);		//m represents for how many times the loop should run
		
	  if(NumOfChildren%2 == 1)
	    {
	      data[0] = *(pts+i_father*M_pts+ind_pts);
	      for(j=1;j<=m;j++)		//data contains resp co_ordinates for new points
		{	
		  temp = j*len;
		  data[j] = data[0] - temp;
		  data[j+m] = data[0] + temp;
		}		
	    }//if ends
	  
	  else
	    {
	      data[0] = *(pts+i_father*M_pts+ind_pts) - len/2.0;
	      for(j=1;j<m;j++)		//data contains resp co_ordinates for new points
		{	
		  temp = j*len;
		  data[j] = data[0] - temp;
		  data[j+m] = data[0] + temp;
		}
	      data[j] = data[0] + j*len;
	    }//else ends
	  len = len/2.0;
		
	  for(m=0;m<NumOfChildren;m++)	//This loop is for populating pts and epsi
	    {
	      // test if the new pts is the same as its father
	      
	      if (father->val == data[m])
		{
		  refine_tree node = new r_node(father,NumOfChildren);
		  father->children[m] = node;
		  new_queue.push_back(node);			
		  /* pts is the same but we have to update epsi */
		  data_epsi_out[i_father*M_epsi+ind_ep] = len;			
		}
	      else
		{	       	
		  //cout << "increasing tree_ind:" << tree_ind << endl;
		  tree_ind++;
		  refine_tree node = new r_node(tree_ind,father,data[m],NumOfChildren);	    
		  father->children[m] = node;
		  new_queue.push_back(node);
			
		  for(k=0;k<M_pts; k++) 
		    data_pts_out[tree_ind*M_pts+k] = pts[i_father*M_pts+k];			  
				
		  data_pts_out[tree_ind*M_pts+ind_pts] = data[m];
			
		  for(l=0;l<M_epsi; l++) 
		    data_epsi_out[tree_ind*M_epsi+l]= epsi[i_father*M_epsi+l];

		  data_epsi_out[tree_ind*M_epsi+ind_ep]= len;
			  
		} //else for  if (father->new_val == data[m]) ends				
	    }// m-for loop ends       
	}	//if(ind_ep != -1) ends
    }//NumElems for loop ends	
		
  for(i=0;i<M_epsi;i++)
    data_dim_out[i] = dim[i];
	
  //updating the object(P) variables
  mxFree(pts);
  mxFree(epsi);
  mxFree(dim);

  queue = new_queue;
  pts = data_pts_out;
  epsi = data_epsi_out;
  dim = data_dim_out;
	
  return false;
}//method refine_all ends

void traj_set::print_refine_history(){
  cout << "Entering print refine history..." << endl;
  int i,j;
  r_node* nod;
  cout << "NumRoots: " << NumRoots << endl;
  for( i=0; i< NumRoots; i++) {
    cout << refine_history[i]->index << "  ";
  }
  cout << endl;
  
  for( i=0; i< NumRoots; i++) {
    
    nod = refine_history[i]->children[j];
    for (j = 0; j<3; j++){
      nod = refine_history[i]->children[j];
      if (nod != NULL)
	cout << nod->index << " ";	
      
      else 
	cout << "NaN "; 
    }
    cout << "    "<< endl;
  }
  cout << endl;
}

void ComputeExpTraj(CmacVF *&C, Array2D pts, Array2D XS0, Array2D epsi, Array1D tspan, Array<trajectory*,1> &trajArray, Array2D &Xf, Array2D &XSf, Array2D &ExpaMax, Array1D &Keep ) {

#if _DEBUG >= 1
   cout << "Entering ComputeExpTraj with CMAC..." <<endl;
#endif

  int nb_traj = pts.extent(1);
  Xf.resize(N,nb_traj);
  Xf = 0;

  XSf.resize(N*Ns,nb_traj);
  XSf = 0;
  
  Keep.resize(nb_traj);
  Keep = 0;

  trajArray.resize(nb_traj);
  
  /* Sensitivity Reinitialisation */
  
  Array1D expa(N);
  Array2D xS(N,Ns,ColumnMajorArray<2>());
  xS = 0;
  double *xSdata = xS.data();
  ExpaMax.resize(N,nb_traj);
  ExpaMax = 0;

  /* CVodes Reinitialization */

#if _DEBUG >= 2
  cout << "Init CVodes..." << endl;
#endif  

  InitFdata(f_data, cvm_Mdata->mx_data);
  int status, statusFSA;
  int itask = CV_ONE_STEP_TSTOP;
  int dimu;
  double tlast,tret(0.),h;
  booleantype iret = FALSE;

  Fdata* data = (Fdata*) f_data;
  data->ns = Ns;

  realtype *xdata = N_VGetArrayPointer(y);
  Array1D x(xdata,shape(N),neverDeleteData);
  Array1D p(data->p,shape(data->dimp),neverDeleteData);

  double t0 = tspan(0);
  double tout = tspan(1);

  /* rootfinding */

#if _DEBUG >= 2
  cout << "Root Finding..." << endl;
#endif  

  for (int i=0; i<DIMX; i++) {
    data->delta[i]=C->qnt[i]/3;
  }

#if _DEBUG >= 2
  cout << "Root Init ..." << endl;
#endif  

  CVodeRootInit(cvode_mem, DIMX, g_delta, f_data);
  Array1D xp(data->xp, shape(DIMX)); // x prec

  /* control */

#if _DEBUG >= 2
  cout << "Control Variable ..." << endl;
#endif  

  void* u_data;
  u_data = (Array1D*) new Array1D(1);  
  GetU(f_data, u_data);
  Array1D * u = (Array1D*) u_data; 
  dimu = (*u).extent(0);

  /* working variables */

  Range Rx(0,DIMX-1);
  Range All = Range::all();   

  int kmax = 5;
  int compteur = 0;
  int i,ir,is,j,nb_new_pts, nb_new_end_pts(0);
  int statusgdk;

  Array2D *nX,X;
  Array2D *nU,U;
  Array2D *nExpa,Expa;
  
  Array1D *ntime,time;
  Array<int,1> Xind;

#if _DEBUG >= 2
  cout << "Loop over trajectories..." << endl;
#endif

  for(j=0;j<nb_traj;j++) {

#if _DEBUG >= 2
    cout << "----------------------------------------" <<endl;
    cout << " tajectory " << j << endl;
#endif

    x = pts(Rx,j); // new x0
    p = pts(All,j);

#ifdef _TD
      x(DIMX)=0.;
#endif

    UpdateFdata(t0,y,f_data,0,NULL);
    xp = x(Range(0,DIMX-1));
    compteur = 0; 
    iret = FALSE;

    /* Init output arrays  */
          
    time.resize(kmax);
    time(0) = t0;
    tret = t0;
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
      CVodeReInit(cvode_mem, f, t0, y, itol, reltol, &Sabstol);
      CVodeSensReInit(cvode_mem,ism,yS0);
      break;
    case CV_SV:
      CVodeReInit(cvode_mem, f, t0, y, itol, reltol, NV_abstol);
      CVodeSensReInit(cvode_mem,ism,yS0);
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
	cout << "xp:" << xp << endl;
#if _DEBUG >= 3       
	cout << "time:" << time << endl;
	cout << "X:" << X << endl;
	cout << "U:" << U << endl;
#endif
	CVodeGetDky(cvode_mem, tout, 1, y);
	cout << "dx:" << x << endl;


	cout << "CVODES was unhappy for some reason. Status= " << status << endl;
	mexErrMsgTxt("Dying...");	      

	return;
      }
      
      /* Test if tout was reached */
      
      CVodeGetCurrentStep(cvode_mem, &h);
      if ( (tret - tout)*h >= 0.0 ) {
	tret = tout;
	statusgdk = CVodeGetDky(cvode_mem, tout, 0, y);
	if (statusgdk)
	  cout << "Pb testing tout 1"<< endl;
	iret = TRUE;
      }      

      /* get sensitivity and compute expansion */
      
      CVodeGetSens(cvode_mem, tret, yS);
      for (is=0; is<Ns; is++)
	GetData(yS[is], &xSdata[is*N], N);
      
      for (i = 0; i< N; i++) {
	expa(i) = 0;
	for (is = 0; is <Ns; is++)
	  expa(i) += abs(xS(i,is))*epsi(is,j);
	ExpaMax(i,j) = max(expa(i), ExpaMax(i,j));
      }	        
      
      /*  Update f_data. if changed at a discontinuity, reinit cvode */
      
      if (UpdateFdata(tret,y,f_data,0,NULL)) {
	
	CVodeGetDky(cvode_mem, tret, 0, y);
	
	switch (itol) {
	case CV_SS:
	  CVodeReInit(cvode_mem, f, tret,y, itol, reltol, &Sabstol);
	  break;
	case CV_SV:
	  CVodeReInit(cvode_mem, f, tret,y, itol, reltol, NV_abstol);
	  break;
	}
	CVodeSetInitStep(cvode_mem,h/2);
      }
      
      GetU(f_data, u_data);  
      
      if (status == CV_ROOT_RETURN) {
	
	compteur++;
	if (compteur < kmax) {
	  time(compteur) = tret;
	  U(All,compteur)= *u;
	  X(All,compteur)= x;
	  Expa(All,compteur)= expa; 
	  xp = x(Range(0,DIMX-1));

	}
	else  {// try to guess how many points left
	  
	  kmax += max(1,min((int) ceil((tout-tret)/(tret-tlast)),kmax));
	
	  X.resizeAndPreserve(N,kmax);
	  U.resizeAndPreserve(dimu,kmax);
	  Expa.resizeAndPreserve(N,kmax);
	  time.resizeAndPreserve(kmax);
	  
	  time(compteur) = tret;
	  U(All,compteur)= *u;
	  X(All,compteur)= x;     
	  Expa(All,compteur)= expa; 
	  xp = x(Range(0,DIMX-1));
	
	}
      
      }
    
      /* break if we need to */
      
      if(iret)  break;      
      
    }

    /* last point */

    if (status != CV_ROOT_RETURN) { 

      compteur++;
      X.resizeAndPreserve(N,compteur+1);
      U.resizeAndPreserve(dimu,compteur+1);
      Expa.resizeAndPreserve(N,compteur+1);
      time.resizeAndPreserve(compteur+1);

      time(compteur) = tret;
      U(All,compteur)= *u;
      X(All,compteur)= x;     
      Expa(All,compteur)= expa; 
      xp = x(Range(0,DIMX-1));

    }

    else {
      
      X.resizeAndPreserve(N,compteur+1);
      U.resizeAndPreserve(dimu,compteur+1);
      Expa.resizeAndPreserve(N,compteur+1);
      time.resizeAndPreserve(compteur+1);

    }
    
    /* Mark new explored pts and only keep those in the trajectory */
    C->mark_as_visited(X,Xind);
    nb_new_pts = Xind.extent(0);
 
    if (nb_new_pts) {
      //      cout << Xind(nb_new_pts-1) << "/" << compteur << endl;    
      
      if (Xind(nb_new_pts-1) >= compteur-2) {	
	Xf(All, nb_new_end_pts) = X(All,compteur);
	for (is = 0; is<Ns; is++) {
	  for (i = 0; i< N; i++) 
	    XSf(is*N+i,nb_new_end_pts) = xS(i,is);
	}  
	nb_new_end_pts++;
	Keep(j) =1.;
      }
    }
 
    trajArray(j) = new trajectory();
    nX = (trajArray(j)->X);
    nU = (trajArray(j)->U);
    nExpa = (trajArray(j)->Expa);
    ntime = (trajArray(j)->time);
    
    (*nX).resize(N,nb_new_pts);
    (*nU).resize(dimu,nb_new_pts);
    (*nExpa).resize(N,nb_new_pts);
    (*ntime).resize(nb_new_pts);
    
    for (i=0; i< nb_new_pts; i++) {
 
      ir = Xind(i);
      (*nX)(All,i) = X(All,ir);
      (*nU)(All,i) = U(All,ir);
      (*nExpa)(All,i) = Expa(All,ir);
      (*ntime)(i) = time(ir);

    }
    
    //  cout << "time:" << *time << endl;
    //  cout << "X:" << *X << endl;
    //  cout << "Expa:" << *Expa << endl;
    //  cout << "U:" << *U << endl;

  }

  Xf.resizeAndPreserve(N,nb_new_end_pts);
  XSf.resizeAndPreserve(N*Ns,nb_new_end_pts);  


#if _DEBUG >= 1
  cout << "Leaving ComputeExpTraj with CMAC ..." << endl;  
#endif
}

void ComputeTrajs(CmacVF *&C, Array2D pts, Array1D tspan, Array<trajectory*,1> &trajArray, Array2D &Xf, Array1D &Keep) {

#if _DEBUG >= 1
  cout << "Entering ComputeTrajs ..." <<endl;
#endif

  int nb_traj = pts.extent(1);
  trajArray.resize(nb_traj);  
  Keep.resize(nb_traj);
  Keep = 0;

  /* CVodes Reinitialization */

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

  double t0 = tspan(0);
  double tout = tspan(1);
  Xf.resize(N,nb_traj);

  /* rootfinding */
  
  for (int i=0; i<DIMX; i++)
    data->delta[i]=C->qnt[i]/3;
  
  CVodeRootInit(cvode_mem, DIMX, g_delta, f_data);
  Array1D xp(data->xp, shape(DIMX)); // x prec
  
  /* control */

  void* u_data;
  u_data = (Array1D*) new Array1D(1);  
  GetU(f_data, u_data);
  Array1D * u = (Array1D*) u_data; 
  dimu = (*u).extent(0);

  /* working variables */

  Range Rx(0,DIMX-1);
  Range All = Range::all();   

  int kmax=5;
  int compteur = 0;
  int i,ir,j,nb_new_pts, nb_new_end_pts(0);

  Array2D *nX,X;
  Array2D *nU,U;
  Array1D *ntime,time;
  Array<int,1> Xind;

#if _DEBUG >= 2
    cout << "Loop over trajectories..." << endl;
#endif
 
  for(j=0;j<nb_traj;j++) {

#if _DEBUG >= 2
    cout << "----------------------------------------" <<endl;
    cout << "trajectory #" << j << endl;
#endif

    x = pts(Rx,j); // new x0
    p = pts(All,j);
    
#ifdef _TD
	x(DIMX)=0.;
#endif

    UpdateFdata(t0,y,f_data,0,NULL);

    tret = t0;
    xp = x(Range(0,DIMX-1));
    compteur = 0;
    iret = FALSE;

    /* Init output arrays  */
          
    time.resize(kmax);
    time(0) = t0;
    X.resize(N,kmax);
    X(All,0) = x;
    U.resize(dimu,kmax);
    U(All,0)=*u;

    /* Reinit solver */

    switch (itol) {
    case CV_SS:
      status = CVodeReInit(cvode_mem, f, t0, y, itol, reltol, &Sabstol);      
       break;
    case CV_SV:
      status = CVodeReInit(cvode_mem, f, t0, y, itol, reltol, NV_abstol);
      break;
    }
    
    CVodeSetStopTime(cvode_mem,tout);

    /* Integrate system until tout */

    while (1) {
      tlast = tret;

#if _DEBUG >= 3
      cout << "CVode call ... " << endl;
#endif

      status = CVode(cvode_mem, tout, y, &tret, itask);

      /* break on CVode error */
      if (status < 0) {

	cout << "tout:" << tout << endl;
	cout << "traj no: " << j << endl;
	cout <<"p: " << p << endl; 
	cout << "x: " << x << endl;
	cout << "xp: " << xp << endl;
#if _DEBUG>=3
	cout << "time:" << time << endl;
	cout << "X:" << X << endl;
	cout << "U:" << U << endl;
#endif

	CVodeGetDky(cvode_mem, tout, 1, y);
	cout << "dx:" << x << endl;

	cout << "CVODES failed for some reason. Beuh. Status= " << status << endl;
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
 
      /* Update fdata and resolve discontinuities */
   
      if (UpdateFdata(tret,y,f_data,0,NULL)) { // Update f_data. if changed at a discontinuity, reinit cvode 
      
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
 
#if _DEBUG >= 3
  cout << "get new input ..." << endl;
#endif  
      GetU(f_data, u_data);     
            
      if (status == CV_ROOT_RETURN) {
	
	/* update compteur and resize vectors if need for more space */
	
	compteur++;
	if (compteur < kmax) {
	  time(compteur) = tret;
	  U(All,compteur)= *u;
	  X(All,compteur)= x;     
	  xp = x(Range(0,DIMX-1));
	}

	else  {// try to guess how many points left

	  kmax += max(1,min((int) ceil((tout-tret)/(tret-tlast)),kmax));
	  X.resizeAndPreserve(N,kmax);
	  U.resizeAndPreserve(dimu,kmax);
	  time.resizeAndPreserve(kmax);
	  
	  time(compteur) = tret;
	  U(All,compteur)= *u;
	  X(All,compteur)= x;     
	  xp = x(Range(0,DIMX-1));
	
	}
      
      }
    
      /* break if we need to */
      
      if(iret)  break;      
      
    }

    /* last point */

#if _DEBUG >= 2
      cout << "Last point ... " << endl;
#endif

      if (status != CV_ROOT_RETURN) { 

	compteur++;
	X.resizeAndPreserve(N,compteur+1);
	U.resizeAndPreserve(dimu,compteur+1);
	time.resizeAndPreserve(compteur+1);

	time(compteur) = tret;
	U(All,compteur)= *u;
	X(All,compteur)= x;     
	xp = x(Range(0,DIMX-1));
	
      }
      
      else {
	
     	X.resizeAndPreserve(N,compteur+1);
	U.resizeAndPreserve(dimu,compteur+1);
	time.resizeAndPreserve(compteur+1);
      }
      

    /* Mark new explored pts and only keep those in the trajectory */

#if _DEBUG >= 2
      cout << "Sparse points ... " <<  endl;
#endif
      
    C->mark_as_visited(X,Xind);
    nb_new_pts = Xind.extent(0);

    if (nb_new_pts) {
      //      cout << Xind(nb_new_pts-1) << "/" << compteur << endl;    

      if (Xind(nb_new_pts-1) >= compteur-2) {	
	Xf(All, nb_new_end_pts) = X(All,compteur);
	nb_new_end_pts++;
	Keep(j) =1.;
      }
    }
      
    /* Write Trajectory */

#if _DEBUG >= 2
      cout << "Write trajectory ... " << j << endl;
      cout << "time:" << time << endl;
      cout << "X:" << X << endl;
      cout << "U:" << U << endl;
#endif
      
      trajArray(j) = new trajectory();
      nX = (trajArray(j)->X);
      nU = (trajArray(j)->U);
      ntime = (trajArray(j)->time);
      
      (*nX).resize(N,nb_new_pts);
      (*nU).resize(dimu,nb_new_pts);
      (*ntime).resize(nb_new_pts);
      
      for (i=0; i< nb_new_pts; i++) {
	ir = Xind(i);
	(*nX)(All,i) = X(All,ir);
	(*nU)(All,i) = U(All,ir);
	(*ntime)(i) = time(ir);
      }

#if _DEBUG >= 2
      cout << "Xind:" << Xind << endl;
      cout << "ntime:" << *ntime << endl;
      cout << "nX:" << *nX << endl;
      cout << "nU:" << *nU << endl;
#endif
     
  }

  Xf.resizeAndPreserve(N,nb_new_end_pts);
  
#if _DEBUG >= 1
  cout << "Leaving ComputeTrajs ..." << endl;
#endif  
  
}


