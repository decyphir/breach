function [p,rob,Pb] = GetPropParamBin(Sys, phi, P, prop_opt, traj, t_phi)
%GETPROPPARAMBIN search values for parameters in a formula phi so that phi
% is satisfied by a set of traces - assumes monotonicity, to be specified
% by the user
%
% Synopsis: [p, rob] = GetPropParamBin(Sys, phi, P, prop_opt [, traj])
%
% Inputs:
%  - Sys        : is a system structure for Breach
%  - phi        : is an STL property
%  - P          : is a Breach set of parameters with param values for phi
%  - prop_opt   : structure with the following fields
%     - params     : is a cell of property param names to find. In the case of multiple
%                    possible satisfying values, the order in which the
%                    parameter appears determines the priority order in which
%                    they are optimized
%     - monotony   : is an array specifying the monotonicity of phi wrt each
%                    parameter. should be maximized (monotony[i] = j) or
%                    minimized (monotony[i] = -1); note  monotony[i] is the
%                    monotonicity of param[i];
%     - ranges     : is the search interval(s) for parameter values
%     - p_tol      : Precision of the binary search in each parameter
%  - traj       : is an array of trajectories
%
% Outputs:
%  - p   :  parameter values for phi
%  - rob :  corresponding robust satisfaction
%

if ~exist('traj','var')
    traj = P.traj;
elseif isempty('traj')
    traj = P.traj;
end

if ~exist('t_phi','var')
    t_phi=0;
end


if exist('opt','var')
    if ~isfield(prop_opt,'verbose')
        opt.verbose=1;
    end
else
    opt.verbose=1;
end

if isfield(prop_opt, 'order')
    order = prop_opt.order;
else
    order = 1:numel(prop_opt.params);
end

params = prop_opt.params(order);
monotony = prop_opt.monotony(order);
ranges = prop_opt.ranges(order,:);

if ~isfield(prop_opt,'p_tol')
    has_p_tol = 0;
    p_tol = -1*ones(numel(params));
else
    has_p_tol = 1;
    p_tol = prop_opt.p_tol(order);
end

if ~isfield(prop_opt,'rel_tol')
    rel_tol = 1e-4;
else
    rel_tol = prop_opt.rel_tol;
end

num_params = numel(params);
pb = zeros(1,num_params);
pw = zeros(1,num_params);

% indices for best and worst corners in range
% note: these are scalar indices for 2D arrays
% recall that (i,j) == i+nb_row*j
i_best  = (1:num_params) + num_params*(monotony+1)/2;
i_worst = (1:num_params) + num_params*(1-monotony)/2;

pb = ranges(i_best);
pw = ranges(i_worst);

Pb = SetParam(P, params, pb');
Pw = SetParam(P, params, pw');
valb = STL_Eval(Sys, phi, Pb, traj, t_phi);
valw = STL_Eval(Sys, phi, Pw, traj, t_phi);

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
        fprintf(['Warning: Interval contains only unsat params, result not tight. Try larger parameter ' ...
            'region. \n']);
    end
    return;
end

%% Now we know that there are satisfying and unsatisfying values

val = min(valb);
rob = inf;
p = zeros(1,num_params);   % initialize p
pmax = ranges(:,2);
pmin = ranges(:,1);

rel_err = zeros(1,num_params);
for ip = 1:num_params
    rel_err(ip) = abs(pmax(ip) - pmin(ip))/abs(pmax(ip));
    p_i = (pmax(ip)+pmin(ip))/2; % initialize p_i in case we don't enter the loop
end

if opt.verbose
    rfprintf_reset();
end
    
for ip = 1:num_params
    if opt.verbose
        fprintf('Finding tight %s           ', params{ip});
    end
    while ((abs(pmax(ip)-pmin(ip))>p_tol(ip))&& rel_err(ip)>rel_tol)
        res = bisect(ip);
        if opt.verbose
            rfprintf(res);
        end
    end % end while
    
    if opt.verbose
        rfprintf(num2str(p_i));
    end
        
    Pb = SetParam(Pb, params(ip),p(ip)');
    if opt.verbose
        fprintf('\n');
    end
    
end

    function res =  bisect(ip)
        p_i = (pmax(ip)+pmin(ip))/2;
        Pb = SetParam(Pb, params(ip), p_i');
        valb = STL_Eval(Sys, phi, Pb, traj, t_phi);
        val = min(valb);
        
        res = num2str(p_i);
        if(val>=0)
            rob = min(val,rob);
            p(ip) = p_i;
            if(monotony(ip)<0)
                pmin(ip) = p_i;
            else
                pmax(ip) = p_i;
            end
        else
            if(monotony(ip)<0)
                pmax(ip) = p_i;
            else
                pmin(ip) = p_i;
            end
        end
        rel_err(ip) = abs(pmax(ip) - pmin(ip))/abs(pmax(ip));
    end
end

