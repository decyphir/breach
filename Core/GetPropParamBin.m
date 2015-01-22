function [p, rob] = GetPropParamBin(Sys, phi, P, params, monotony, p_interval, p_tol, traj,opt)
%GETPROPPARAMBIN search values for parameters in a formula phi so that phi
% is satisfied by a set of traces - assumes monotonicity, to be specified
% by the user
% 
% Synopsis: [p, rob] = GetPropParamBin(Sys, phi, P, params, monotony, p_interval, p_tol, traj)
% 
% Inputs:
%  - Sys        : is a system structure for Breach
%  - phi        : is an STL (QMITL) property
%  - P          : is a Breach set of parameters with param values for phi
%  - params     : is a cell of property param names to find
%  - monotony   : is an array specifying the monotonicity of phi wrt each
%                 parameter. should be maximized (optim[i] = j) or
%                 minimized (optim[i] = -1); in the case of multiple
%                 possible satisfying values, the order in which the
%                 parameter appears determines the priority order in which
%                 they are optimized e.g. [1 2] means params{1} is
%                 maximized, then params{2} is maximized; [-2 1] means
%                 params{2} is minimized, then params{1} is maximized,
%                 etc...
%  - p_interval : is the search interval(s) for parameter values
%  - p_tol      : Precision of the binary search in each parameter
%  - traj       : is an array of trajectories
% 
% Outputs:
%  - p   :  parameter values for phi
%  - rob :  corresponding robust satisfaction
%

if exist('opt','var')    
    if ~isfield(opt,'verbose')
        opt.verbose=1;
    end
else
  opt.verbose=1;
end

% NM: I assume that p_interval is a 2D array
num_params = numel(params);
pb = zeros(1,num_params);
pw = zeros(1,num_params);
for ii = 1:num_params
    if(monotony(ii)>0)
        pb(ii) = p_interval(ii,2);
        pw(ii) = p_interval(ii,1);
    else
        pb(ii) = p_interval(ii,1);
        pw(ii) = p_interval(ii,2);
    end
end

Pb = SetParam(P, params, pb');
Pw = SetParam(P, params, pw');

valb = QMITL_Eval(Sys, phi, Pb, traj, 0);
valw = QMITL_Eval(Sys, phi, Pw, traj, 0);

%% Check if everybody is sat

if all(valw>=0)
    p = pw;
    rob = min(valw);
    if opt.verbose
        fprintf(['Warning: Interval contains only sat params, result may be not tight. ' ...
                 'Try larger parameter region. \n']);
    end
    return;
    
end

%% Check if at least somebody is sat
if any(valb<0)
    p = pb;
    rob = max(valb);
    if opt.verbose
        fprintf(['Error: Interval contains only unsat params, result not tight. Try larger parameter ' ...
        'region. \n']);
    end
    return;
end

% Now we know that there are satisfying and unsatisfying values

val = min(valb);
rob = inf;
p = zeros(1,num_params); % initialize p
for ii = 1:num_params      % optimize independently in the order
                           % given in params
    if opt.verbose
        fprintf('Optimizing %s ', params{ii});
    end
    timeout = 100;
    pimax = p_interval(ii,2);
    pimin = p_interval(ii,1);
    
    err = p_tol(ii);
    
    if opt.verbose
        rfprintf_reset();
    end
    
    p_i = (pimax+pimin)/2; % initialize p_i in case we don't enter the loop
    while (abs(pimax-pimin)>err)
        
        p_i = (pimax+pimin)/2;
        Pb = SetParam(Pb, params(ii), p_i');
              
        valb = QMITL_Eval(Sys, phi, Pb, traj, 0);
        val = min(valb);
        
%       fprintf('  pimin: %g  pimax: %g p_i: %g val %g \n', pimin, pimax, p_i, val);
        res = num2str(p_i);
        if opt.verbose
            rfprintf(res);
        end
        if(val>0)
            rob = min(val,rob);
            if(monotony(ii)<0)
                pimin = p_i;
            else
                pimax = p_i;
            end
        else
            if(monotony(ii)<0)
                pimax = p_i;
            else
                pimin = p_i;
            end
        end
       
        timeout = timeout-1;
        if(timeout==0)
            fprintf('Time out !!');
            break;
        end
        
    end % end while
    
    if opt.verbose
        rfprintf(num2str(p_i));
    end    
    
    if(val>0)
        p(ii) = p_i;
    else
        if(monotony(ii)>0)
            p(ii) = pimax;
        else
            p(ii) = pimin;
        end
    end
    Pb = SetParam(Pb, params(ii),p(ii)');
    if opt.verbose
        fprintf('\n');
    end
end
end
