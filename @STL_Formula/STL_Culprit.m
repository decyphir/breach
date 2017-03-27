function [val0, mu, t, idx] = STL_Culprit(Sys, phi, P, tau, threshold, VERBOSE)
%STL_CULPRIT tries to find out who (which predicate) and when is
% responsible for the value of val0 = rho(phi, traj, tau), for all
% parameter vector of P. This is a brute force implementation: does not try
% to solve conflict (e.g. equality in min or max), and does not verify that
% the value is actually admissible (time might be wrong). It just
% decomposes phi into its predicates, evaluates them and looks for a time
% when value val is reached. It recomputes the evaluation of phi on P using
% the classic evaluation method.
% 
% Synopsis: [val0, mu, t, idx]  = STL_Culprit(Sys, phi, P, tau [, threshold[, VERBOSE]])
% 
% Inputs:
%  - Sys       : The system
%  - phi       : The formula for which we want to find the predicate. If
%                it is a predicate, mu is phi.
%  - P         : The parameter set containing trajectories against which
%                the formula is checked. It may contain many parameter
%                vectors. The trajectories must be computed. If the
%                trajectories are not computed, all output are equals to
%                [], and a warning is thrown. The evaluation of phi does
%                not need to be computed.
%  - tau       : Time step at which we ask for the predicate responsible
%                for the value of phi. It must be a scalar. It is strongly
%                recommended that tau belongs to traj.time, otherwise, the
%                time point t may be not found.
%  - threshold : (Optionnal, default or empty=1.0e-14) Indicates the
%                maximal difference between the evaluation of the formula
%                at a time point and its evalutation at time tau to say
%                that we have found the valid time point t.
%  - VERBOSE   : (Optional, default=true) boolean indicating if the
%                computation progress bar is shown.
% 
% Outputs:
%   - val0 Truth value of phi at time tau. Dimension of val0 is
%          1 x size(P.pts,2).
%   - mu   The predicate responsible for the value val0. It is a cell array
%          of dimension 1 x size(P.pts,2). If, for a parameter vector, the
%          predicate responsible for the evaluation of the formula is not
%          found, mu, at the corresponding index, is empty.
%   - t    time point such that rho(mu, traj, t) == val0. It is an array of
%          dimension 1 x size(P.pts,2). If, for a parameter vector, the
%          predicate responsible for the evaluation of the formula is not
%          found, t, at the corresponding index, is set to NaN.
%   - idx  Index of t in traj.time. If t does not belong to traj.time
%          (typically if phi is a predicate and tau doesn't belong to
%          traj.time), the index of the smallest time point greater that
%          tau is provided, or the last time point if there is no greater
%          time point greater than tau. If, for a parameter vector, the
%          predicate responsible for the evaluation of the formula is not
%          found, idx, at the corresponding index, is set to 0.
% 
%See also STL_SEvalDiff
%
% NEEDED IMPROVEMENTS ARE: FIND OUT WHICH PREDICATE IS RESPONSIBLE FOR THE
% VALUE OF PHI FOR EACH PARAMETER SET.
%

%%%%%%%%%%%%%%%%%
% check input
if ~exist('VERBOSE','var')
    VERBOSE = 1;
end

if(~exist('threshold','var')||isempty(threshold))
    threshold = 1e-14;
end

if ~isfield(P,'traj')
    warning('STL_Culprit:noTrajField','The parameter set P has no traj field.');
    val0 = [];
    mu = [];
    t = [];
    idx = [];
    return
end

%%%%%%%%%%%%%%%%%
% do stuff
[~,val0] = SEvalProp(Sys, P, phi, tau, [], 0, 'classic', VERBOSE);

numParamSet = size(P.pts,2);
mu = cell(1,numParamSet);
t = zeros(1,numParamSet);
idx = zeros(1,numParamSet);

if strcmp(phi.type, 'predicate')
    t = tau*ones(1,numParamSet);
    [mu{:}] = deal(phi); % fill mu with phi
    for ii = 1:numParamSet
        traj_time = P.traj{P.traj_ref(ii}).time;
        idx(ii) = find(traj_time>=tau,1,'first'); % index of tau if tau belong in traj_time
        if isempty(idx(ii))
            idx(ii) = traj_time(end); % there is no time point higher than tau
        end
    end
else
    %
    % TODO : IMPROVE THE IMPLEMENTATION OF THIS CODE
    %
    % Typically : do a recursive exploration by splitting at level 2 and
    %             looking for the subformula leading to this evaluation
    %             value (be carefull : how do we manage ev, alw and until?)
    %
    % LIMITATION : if tau does not belong to traj.time, we may not found a
    %               valid time instant
    %
    mus = STL_ExtractPredicates(phi);
    for ii = 1:numParamSet % for each param set
        traj = P.traj{P.traj_ref(ii});
        Pi = Sselect(P,ii);
        found = false;
        jj = 1;
        while ~found && jj<=numel(mus) %for each predicate, we look if the truth value, at one time
            mu_tmp = mus(jj);          %point, is equal to the initial evaluation for this param set
            val = STL_Eval(Sys, mu_tmp, Pi, traj, traj.time, 'classic');
            idx_tmp = find(abs(val-val0(ii))<threshold,1);
            if ~isempty(idx_tmp)
                mu{ii} = mu_tmp;
                t(ii) = traj.time(idx_tmp);
                idx(ii) = idx_tmp;
                found = true; % go for the next parameter vector
            end
            jj = jj + 1;
        end
        if(~found)
            t(ii) = NaN;
        end
    end
end

end
