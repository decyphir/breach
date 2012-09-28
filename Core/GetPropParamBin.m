function [p rob] = GetPropParamBin(Sys, phi, P, params, optim, p_interval, p_tol, traj) 
%  GETPROPPARAMBIN search values for parameters in a formula phi so that
%  phi is satisfied by a set of traces 
%     
%  Synopsis: [p rob] = GetPropParamBin(Sys, phi, Paramset, params_to_optim, optim, p_interval, p_tol, traces) 
%   
%  where : 
%  
%  - Sys             is a system structure for Breach
%  - phi             is an STL (QMITL) property
%  - ParamSet        is a Breach set of parameters with param values for phi 
%  - params_to_optim is a cell of property param names to find
%  - optim           is an array specifying whether each parameter should be maximized    
%                    (optim[i] = j) or minimized (optim[i] = -1); in the case of multiple  
%                    possible satisfying values, the order in which the parameter appears   
%                    determines the priority order in which they are
%                    optimized e.g. [1 2] means params{1} is maximized, then
%                    params{2} is maximized; [-2 1] means params{2} is minimized, then
%                    params{1} is maximized, etc                
%  - p_interval      is the search interval(s) for parameter values
%  - p_tol           Precision of the binary search in each parameter                  
%  - traces          is an array of trajectories
%      

  
  for i= 1:numel(params) 
    
    if (sign(optim( find(i==abs(optim))))< 0)
      pol(i) = 1;
      pb(i) = p_interval(i,2);
      pw(i) = p_interval(i,1);
    else
      pb(i) = p_interval(i,1);
      pw(i) = p_interval(i,2);
      pol(i) = -1;
    end
  end
      
  Pb = SetParam(P, params, pb');
  Pw = SetParam(P, params, pw');
       
  valb = QMITL_Eval(Sys, phi, Pb, traj, 0);
  valw = QMITL_Eval(Sys, phi, Pw, traj, 0);
    
  if  (min(valb.*valw) > 0)
    
    fprintf('\n Interval contains only sat or unsat params, try larger parameter region \n');
%    p_interval(:,1) = p_interval(:,1)/2;
%    p_interval(:,2) = p_interval(:,2)*2;
%
%    [p rob] = GetPropParamBin(Sys, phi,P, params, optim, p_interval, p_tol, traj);
    return
  end
      
  % Now we know that there are satisfiable values

  val = min(valb);
  
  for i_to_optim = abs(optim)      % optimize independently in the order
                                   % given in optim 
    
    fprintf('\nOptimizing %s ', params{i_to_optim});
    
    timeout=100;
    pimax = p_interval(i_to_optim,2);
    pimin = p_interval(i_to_optim,1);
    
    err = p_tol(i_to_optim);       
    
    rfprintf_reset();
    while (abs(pimax-pimin)>err)

      p_i = (pimax+pimin)/2;            
      Pb = SetParam(Pb, params(i_to_optim),p_i');            
      valb = pol(i)*QMITL_Eval(Sys, phi, Pb, traj , 0);
      val = min(valb);
      
      %fprintf('  pimin: %g  pimax: %g p_i: %g val %g \n', pimin, pimax, p_i, val);
      res = num2str(pimax);
      rfprintf(res);
      
      if (val>0)
        pimax = p_i;
        rob = val;
      else      
        pimin = p_i;
      end
                 
      timeout =timeout-1;
      if timeout==0
        fprintf('Time out !!');
        break;
      end               
    
    end % end while
    
    p(i_to_optim) = pimax;
    Pb = SetParam(Pb, params(i_to_optim),pimax');  
    
  end
  fprintf('\n');
  