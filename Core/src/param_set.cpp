#include "param_set.h"

double* param_set::Get_pts(int nd)
{
    double *point = (double *)mxCalloc(M_pts,sizeof(double));
    nd = nd*M_pts;
    for(int i=0;i<M_pts;i++)
        point[i] = pts[nd+i];
    return point;  
}   

double* param_set::Get_ep(int nd)
{
    double *epsilon = (double *)mxCalloc(M_epsi,sizeof(double));
    nd = nd*M_epsi;
    for(int i=0;i<M_epsi;i++)
        epsilon[i] = *(epsi+nd+i);
    return epsilon;  
}   

/* This method is for test strategy*/
int param_set::Test1(refine_tree nod,double test[])
{
	int i;
	int nd= nod->index;
	nod->status = -1;
	
	nd = nd*M_epsi;
	for(i=0;i<M_epsi;i++)
	  {	
	    //cout << "idx=" << nd+i << " :"<< epsi[nd+i] << " " << test[i] << endl;
	    if(epsi[nd+i] > test[i]) {
	      nod->status=i;
	      break;
	    }
	  }
	return (nod->status);
}

/*
refine_tree param_set :: pop()
{
	refine_tree node;
	vector<refine_tree>:: iterator itq;
	itq = queue.begin();
	node = *itq;
	queue.erase(itq);
	return node;
}*/

bool param_set::refine_all(int NumOfChildren, int (param_set::*ptr_to_test)(refine_tree, double[]),double tst_eps[]) 
{
	double *data_pts_out,*data_epsi_out,*data_dim_out,*data;
	int i,i_father,k=0,l=0,m,j,ind_pts,ind_ep,N=0;
	int tree_ind=NumElems-1; 
	refine_tree father;
	double len,temp;
		
	mexPrintf("\n-- NumElems is: %d\n",NumElems);	

	data = (double *)mxCalloc(NumOfChildren,sizeof(double));

	/* test all elements in the queue */
	
	vector<refine_tree>::iterator it;
	int stat;
	for(it= queue.begin();it!=queue.end();++it)
	{
	  //cout << "testing node:" << (*it)->index << endl;
	  (*it)->status = (this->*ptr_to_test)(*it,tst_eps);	  
	  if((*it)->status != -1) {
	    //cout << "ah, un pt à refiner: " << (*it)->status << endl; 
	    N++;
	  }
	}

	if (N==0) {
	  return true;
	}
	int NumElemsOld = NumElems;
	
	NumElems = (NumElemsOld + N*(NumOfChildren - NumOfChildren%2)); //New pts will contain these many pts 
	//cout << "New numelem: "  << NumElems << endl;

	data_pts_out = (double *)mxCalloc(M_pts*NumElems,sizeof(double));
	data_epsi_out = (double *)mxCalloc(M_epsi*NumElems,sizeof(double));
	data_dim_out = (double *)mxCalloc(M_epsi*1,sizeof(double));

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
	    //cout << "processing node:" << (*it)->index << endl;
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
					father->children[j] = node;
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

void param_set::write_mxArray(mxArray * &mxOutput)
{
	const char **fnames;
	mxArray *mx_pts_output,*mx_epsi_output,*mx_dim_output;
	
	fnames = (const char **)mxCalloc(3, sizeof(*fnames));
	fnames[0] = "pts";
	fnames[1] = "epsi";
	fnames[2] = "dim";
							// create 1X1 struct matrix for output
	mxOutput = mxCreateStructMatrix(1, 1, 3,fnames);
	
	mx_pts_output = mxCreateDoubleMatrix(M_pts,NumElems,mxREAL);
	mxFree(mxGetPr(mx_pts_output));
	mxSetPr(mx_pts_output,pts);
	mxSetField(mxOutput, 0, fnames[0],mx_pts_output);
	
	mx_epsi_output = mxCreateDoubleMatrix(M_epsi,NumElems,mxREAL);
	mxFree(mxGetPr(mx_epsi_output));
	mxSetPr(mx_epsi_output,epsi);
	mxSetField(mxOutput, 0, fnames[1],mx_epsi_output);

	mx_dim_output = mxCreateDoubleMatrix(M_epsi,1,mxREAL);
	mxFree(mxGetPr(mx_dim_output));
	mxSetPr(mx_dim_output,dim);
	mxSetField(mxOutput, 0, fnames[2],mx_dim_output);
}



void CVM_Refine( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ) 
{
  /*read the inputs*/

  double *dou_dim = (double *)mxGetPr(prhs[1]);
  int NumOfChildren = (int)*dou_dim;
  param_set * P = new param_set(prhs[0], NumOfChildren);
  double *test_epsi = (double *)mxGetPr(prhs[2]);
  
  bool done = false;
  
  while(done == false)
    done = P->refine_all(NumOfChildren,&param_set::Test1,test_epsi);
  
  P->write_mxArray(plhs[0]);
  
  delete P;
}

