function [val, tau] = STL_Eval_Gen(Sys, phi, P, trajs, partition, relabs, taus)
%STL_EVAL_Gen computes the satisfaction function of a property for one or
% many traces.
% 
% Synopsis: [val, tau] = STL_Eval_Gen(Sys, phi, P, partition, relabs, trajs[, taus])
% 
% Inputs:
%  - Sys    : the system
%  - phi    : a STL Formula
%  - P      : is a parameter set which contains one parameter vector only
%             used for properties parameters.
%  - trajs  : is a structure with fields X and time. It may contains many
%             trajectories. In this case, all will be checked wrt the
%             property parameter described in P.
%  - partition  : is the partition of signals given as an array of strings:
%             the robustness is computed in the signals of the partition.
%  - relabs : is a string indicating how to treat variables that are 
%             not in the partition: 'rel' for -inf/+inf or 'abs' for 0.
%  - taus   : (Optional, default=traj.time for each traj in trajs) is the
%             time, possibly an array, when to eval the satisfaction of the
%             property. All time points not belonging to traj.time will be
%             linearly interpolated.
%
% Outputs:
%  - val : is a cell array of dimension 1 x numel(trajs). Each cell
%          contains an line array describing the evaluation of phi at each
%          time of tau. If numel(trajs) is 1, val is the content of its
%          only cell to avoid a useless cell array (so, it is a line
%          array).
%  - tau : is a cell array of dimension 1 x numel(trajs). It contains time
%          points at which the formula is evaluated or interpolated. If the
%          parameter taus is provided, tau is equal to taus. If
%          numel(trajs) is 1, tau contains the content of the only cell
%          (this avoid a useless cell array), and thus becomes a line
%          array.
% 
%See also SEvalProp STL_Formula
%

%% formula is given directly as a string
if ischar(phi)
    STL_Formula('phi_tmp__', phi);
    switch nargin
        case 6
            [val, tau] = STL_EvalThom_Gen(Sys, phi_tmp__, P, trajs, partition, relabs);
        case 7
            [val, tau] = STL_EvalThom_Gen(Sys, phi_tmp__, P, trajs, partition, relabs, taus);
    end
    evalin('base', 'clear phi_tmp__');
    return
end

%% default parameters
if isfield(phi.params,'default_params')
    pnames = fieldnames(phi.params.default_params);
    for ip = 1:numel(pnames)
        % if P doesn't define the parameter, use default
        if FindParam(P,pnames{ip})> size(P.pts,1)
            pval = phi.params.default_params.(pnames{ip});
            if isscalar(pval)
                P = SetParam(P,pnames{ip},pval );
            else
                P.(pnames{ip}) = pval;
            end
        end
    end
end

switch nargin
    case 6
        [val, tau] = STL_EvalThom_Gen(Sys, phi, P, trajs, partition, relabs);
    case 7
        [val, tau] = STL_EvalThom_Gen(Sys, phi, P, trajs, partition, relabs, taus);
end

% Apply constant semantics, if needed
constSemanticsVal = 100; 
if strcmp(phi.semantics, 'constant')
    val(val >= 0) = constSemanticsVal;
    val(val < 0) = -constSemanticsVal;
end


end
