function [p, val, vald] = STL_SEvalDiff(Sys, phi, P, tspan, params, tau, VERBOSE)
%STL_SEVALDIFF computes (d phi / d param) at time tau, for each parameter
% in params. FORMULA PARAMETERS ARE NOT MANAGED.
% 
% Synopsis: [p, val, vald] = STL_SEvalDiff(Sys, phi, P, tspan, params, tau[, VERBOSE])
% 
% Inputs:
%  - Sys     : the system
%  - phi     : the STL formula for which the derivative is computed
%  - P       : the parameter set. It may contain many parameter vectors.
%              The trajectories do not need to be computed, nor the
%              evaluation of phi.
%  - tspan   : time step for trajectories computation
%  - params  : indexes or names of the parameters and initial conditions
%              against which the derivative of phi is computed (may contain
%              many parameters). Parameters not in P.ParamList will not be
%              considerated. If params is empty or contains no valid
%              parameters, p, val and vald are set to [].
%  - tau     : time step at which phi is evaluated. It must be a scalar. It
%              is recommanded that tau is a time step of tspan.
%  - VERBOSE : (Optional, default=true) boolean indicating if the
%              computation progress bar is shown.
% 
% Outputs:
%  - p    : the values of the parameter described by params. It is an array
%           of dimension numel(params) x size(P.pts,2)
%  - val  : the truth value of phi at time tau. Dimension of val is
%           1 x size(P.pts, 2).
%  - vald : the derivative d phi/d p for all parameter p in params. Size of
%           vald is numel(params) x size(P.pts,2), except for formula
%           parameter, for which, vald is set to zero.
% 
% Example (Lorentz84):
%  format shortG
%  CreateSystem
%  P = CreateParamSet(Sys, {'a','b','x0'}, [0.1,0.4 ; 1,7 ; -1,1], 2);
%  P = SetParam(P, {'x1h', 'x1l', 'T'}, [.3; -.3; 5]);
%  STL_ReadFile('oscil_prop.stl');
%  [~, val, vald] = STL_SEvalDiff(Sys, x1_oscil, P, 0:0.1:10, P.dim, 0)
% 
% Other example (Lorentz84):
%  CreateSystem
%  phi_x0 = STL_Formula('phi_x0','x0[0]<-5');
%  P = CreateParamSet(Sys, {'x0','a'}, [-10,1;-0.25,0.25],2);
%  [~, val, vald] = STL_SEvalDiff(Sys, phi_x0, P, 0:0.1:10, {'x0','x1','a'}, 0)
%                          % here, d phi_x0 / d x0 is not nul, but both
%                          % d phi_x0 / d x1 and d phi_x0 / d a are zero.
% 
%See also SplotSensiBar STL_Formula CreateParamSet ComputeTrajSensi
%SPropSensi PPhiSensiLocal
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Check inputs

% dealing with params
if ~exist('VERBOSE','var')
    VERBOSE=1;
end

if isnumeric(params)
    i_params = params(params<=size(P.pts,1)); % keep only existing parameters
    i_params = i_params(i_params>0);
    params = P.ParamList(i_params);
else
    if ischar(params)
        params = {params};
    end
    [idx_valid,i_params] = ismember(params,P.ParamList);
    params = params(idx_valid); % keep only existing parameters
    i_params = i_params(idx_valid);
end
if isempty(i_params)
    p = [];
    val = [];
    vald = [];
    return;
end

%  Check P
if ~isfield(P,'traj_ref')
    P.traj_ref = 1:size(P.pts,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Doing stuff

p = P.pts(i_params, :);

% Compute trajectories if needed
P = ComputeTraj(Sys, P, tspan);

% Find the culpredicate and the right time
[val, mustar, tstar]  = STL_Culprit(Sys, phi, P, tau, [], VERBOSE);

nbParam = numel(params);
nbParamSet = size(P.pts,2);
vald = zeros(nbParam, nbParamSet);
for ii = 1:nbParamSet
    Pi = Sselect(P,ii);
    for jj = 1:nbParam
        param = params{jj};
        % diff the culpredicate
        dmustar = DiffFormula(mustar{ii}, param);
        
        if(i_params(jj)<=P.DimP) % REQUIRED UNTIL FORMULA PARAMETERS ARE MANAGED
            % Eval the derivative at the right time
            vald(jj,ii) = STL_Eval(Sys, dmustar, Pi, Pi.traj, tstar(ii));
        end
    end
end

end
