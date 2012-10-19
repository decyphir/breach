function [p rob] = GetPropParamBin(Sys, phi, P, params, monotony, p_interval, p_tol, traj) 
%  GETPROPPARAMBIN search values for parameters in a formula phi so that
%  phi is satisfied by a set of traces - assumes monotonicity, to be
%  specified by the user
%     
%  Synopsis: [p rob] = GetPropParamBin(Sys, phi, Paramset, params_to_optim, monotony, p_interval, p_tol, traces) 
%   
%  where : 
%  
%  - Sys             is a system structure for Breach
%  - phi             is an STL (QMITL) property
%  - ParamSet        is a Breach set of parameters with param values for phi 
%  - params_to_optim is a cell of property param names to find
%  - monotony        is an array specifying the monotonicity of phi wrt each parameter. should be maximized    
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
    
    if (monotony(i)>0)
      pb(i) = p_interval(i,2);
      pw(i) = p_interval(i,1);
    else
      pb(i) = p_interval(i,1);
      pw(i) = p_interval(i,2);
    end
  
  end
      
  Pb = SetParam(P, params, pb');
  Pw = SetParam(P, params, pw');
       
  valb = QMITL_Eval(Sys, phi, Pb, traj, 0);
  valw = QMITL_Eval(Sys, phi, Pw, traj, 0);
  
  
%% Check if everybody is sat

  if (all(valw>=0))
    p = pw;
    rob = min(valw);      
    fprintf(['\n Warning: Interval contains only sat params, result may be not tight. Try larger parameter ' ...
               'region. \n']);      
    return
    
  end

  if (any(valb<0))
    p = pb;
    rob = max(valb);      
    fprintf(['\n Error: Interval contains only unsat params, result not tight. Try larger parameter ' ...
               'region. \n']);      
    return
  end

  traj= traj(valw<0);
  % Now we know that there are satisfiable values

  val = min(valb);
  rob = inf;
  for i= 1:numel(params)      % optimize independently in the order
                             % given in params 
    
    fprintf('\nOptimizing %s ', params{i});
    
    timeout=100;
    pimax = p_interval(i,2);
    pimin = p_interval(i,1);
    
    err = p_tol(i);       
    
    rfprintf_reset();
    while (abs(pimax-pimin)>err)

      p_i = (pimax+pimin)/2;            
      Pb = SetParam(Pb, params(i),p_i');            
      valb = QMITL_Eval(Sys, phi, Pb, traj , 0);
      val = min(valb);
      
      %fprintf('  pimin: %g  pimax: %g p_i: %g val %g \n', pimin, pimax, p_i, val);
      res = num2str(pimax);
      rfprintf(res);
      
      if (val>0)
        rob = min(val,rob);
        if monotony(i)<0
          pimin = p_i;
        else
          pimax = p_i;
        end
      else      
        if monotony(i)<0
          pimax = p_i;
        else
          pimin = p_i;
        end
      end
                 
      timeout =timeout-1;
      if timeout==0
        fprintf('Time out !!');
        break;
      end               
    
    end % end while
 
    
    if (val>0)
      p(i) = p_i;
    else
      if monotony(i)>0
        p(i) = pimax;
      else
        p(i) = pimin;
      end      
    end
    Pb = SetParam(Pb, params(i),p(i)');  
  end
  fprintf('\n');
  