function [Pbest,val_best] = SOptimPropLog(Sys, phi, opts, verbose)
% SOPTIMPROPLOG answers either the better parameter set found if it is
% negative or all leading to a positive evaluation of phi. If
% opts.StopWhenFound is set to 1, the function stops as soon as it finds a
% parameter set positive (max) or negative (min).
% The value of the parameter not in opts.params in equal to Sys.p.
%
% Synopsis: [Pbest,val_best] = SOptimPropLog(Sys, phi, opts [, verbose])
%
% Inputs:
%   - phi
%   - opts : describes the options. It contains the following fields:
%      - tspan          See SOptimProp
%      - tau            See SOptimProp (default= first time point of tspan)
%      - params         See SOptimProp
%      - lbound         See SOptimProp
%      - ubound         See SOptimProp
%      - MaxIter        See SOptimProp
%      - OptimType      See SOptimProp (default= 'Max')
%      - StopWhenFound  See SOptimProp
%      - StopWhenFoundInit : more or less used (TO IMPLEMENT)
%      - nbSplit     (mandatory if OptimizeNInitPoints or OptimizeNboxes)
%              it defines the number of initial points or the number of
%              sub-space to explore.
%      - OptimizeGlobally    : if set to one, SOptimPropLog performs an
%              optimization in the parameter space defined by opts.lbound
%              and opts.ubound, starting the research of the optimum from
%              the center of this space 
%      - OptimizeNInitPoints : if set to one, SOptimPropLog performs an
%              optimization in the parameter space defined by opts.lbound
%              and opts.ubound, starting the research from nbSplit initial
%              points, uniformly randomly distributed over a log scale.
%      - OptimizeNBoxes      : if set to one, SOptimPropLog performs 
%              optimizations in opts.nbSplit sub-spaces quasi-randomly
%              defined on a logarithmic scale.
%   - verbose : set to one to print progress of the computation. Possible
%               values : 0, 1 or 2. (default= 0)
%
% Outputs:
%      - Pbest   : contains either the best parameter set found or the set of
%                    parameter sets verifying phi.
%      - val_best: contains the truth value of phi for the parameters set in
%                    Pbest
%See also SOptimProp

%% check inputs
if isfield(opts, 'tspan')
    tspan = opts.tspan;
elseif isfield(Sys, 'tspan')
    tspan = Sys.tspan;
else
    error('SOptimPropLog:noTspan','The opts field does not contain tspan field.');
end

if isfield(opts,'tau')
    tau = opts.tau;
    if tspan(1) > tau
        tspan = [tau tspan];
    end
else
    tau = tspan(1);
end

if isfield(opts,'OptimType')
    OptimType = lower(opts.OptimType); % to avoid case mistake, we convert to lower case
else
    OptimType = 'max';
end

if isfield(opts,'StopWhenFound')
    StopWhenFound = opts.StopWhenFound;
else
    StopWhenFound = 0;
end

if isfield(opts,'StopWhenFoundInit')
    StopWhenFound = StopWhenFound | opts.StopWhenFoundInit;
end

if isfield(opts,'OptimizeGlobally')
    OptimizeGlobally = opts.OptimizeGlobally;
else
    OptimizeGlobally = 0;
end

if isfield(opts,'OptimizeNInitPoints')
    OptimizeNInitPoints = opts.OptimizeNInitPoints;
else
    OptimizeNInitPoints = 0;
end

if isfield(opts,'OptimizeNBoxes')
    OptimizeNBoxes = opts.OptimizeNBoxes;
else
    OptimizeNBoxes = 0;
end

if OptimizeNInitPoints==1 || OptimizeNBoxes==1
    if ~isfield(opts,'nbSplit')
        error('SOptimPropLog:noNbSplit',...
                'There are no field nbSplit in opts, despite there is a field OptimizeNBoxes or OptimizeNInitPoints.');
    else
        nbSplit = opts.nbSplit;
    end
end

if ~exist('verbose','var')
    verbose = 0;
end

phi = QMITL_OptimizePredicates(Sys,phi); % optimization of the predicates
val_best = -inf;


%% OptimizeGlobally
if OptimizeGlobally == 1
    if(verbose>=1)
        fprintf('\n  Global optimization :\n');
    end
    Pall = CreateParamSet(Sys,opts.params,[opts.lbound;opts.ubound]');
    [val_best,Pbest] = SOptimProp(Sys,Pall,phi,opts);
    
    if StopWhenFound % check if we have found a solution
        if strcmp(OptimType,'max') && val_best>0
            return ;
        elseif strcmp(OptimType,'min') && val_best<0
            return ;
        end
    end
end

%% OptimizeNInitPoints
if OptimizeNInitPoints==1
    if(verbose>=1)
        fprintf('\n  Global optimization with %d initial values:\n',nbSplit);
    end
    Prl = CreateRandomLogParamSets(Sys,opts.params,[opts.lbound;opts.ubound]',nbSplit);
    [val_tmp,Ptmp] = SOptimProp(Sys,Prl,phi,opts);
    
    idx_pos = find(val_tmp>0); % keep the best or the positive ones
    if ~isempty(idx_pos)
        if(val_best>0)
            val_best = [val_best,val_tmp(idx_pos)];
            Pbest = SConcat(Pbest,Sselect(Ptmp,idx_pos));
        else
            val_best = val_tmp(idx_pos);
            Pbest = Sselect(Ptmp,idx_pos);
        end
    elseif(val_best<=0)% no positive in val_tmp, compare only if val_best is negative
        [val_tmp,idx_tmp] = max(val_tmp);
        if(val_tmp>val_best)
            val_best = val_tmp;
            Pbest = Sselect(Ptmp,idx_tmp);
        end
    end
    
    if StopWhenFound % check if we have found a solution
        if strcmp(OptimType,'max') && val_best(1)>0
            return ;
        elseif strcmp(OptimType,'min') && val_best(1)<0
            return ;
        end
    end
end

%% OptimizeNBoxes
if OptimizeNBoxes==1
    if(verbose>=1)
        fprintf('\n  Global optimization in %d sub-spaces.\n',nbSplit);
    end
    Prl = CreateRandomLogParamSets(Sys,opts.params,[opts.lbound;opts.ubound]',nbSplit);
    try
        %TODO : IF opt.StopWhenFoundInit IS SET TO 1, COMPUTE TRAJECTORIES
        % ONE BY ONE (see SOptimProp)
        Prl = ComputeTraj(Sys,Prl,tspan);
    catch
        error('SOptimPropLog:CVODESerror','Error during computation of initial trajectories. Try with opts.OptimizeNBoxes=0.');
    end
    [Prl,val_tmp] = SEvalProp(Sys, Prl, phi, tau);
    
    % we look for the maximum
    switch OptimType
        case 'max'
            [~,iv] = sort(-val_tmp);
        case 'min'
            [~,iv] = sort(val_tmp);
        case 'zero'
            [~, iv] = sort(abs(val_tmp));
    end
    
    k=0;
    for i=iv
        k=k+1;
        if(verbose>=2)
            fprintf('\n  Optimization of subspace %d/%d :\n', k, numel(iv));
        end
        Ptmp = Sselect(Prl,i);
        val = GetParam(Ptmp,params);
        epsi = GetEpsi(Ptmp,params);
        opts.lbound = (val-epsi)';
        opts.ubound = (val+epsi)';
        [val_tmp,Ptmp] = SOptimProp(Sys,Ptmp,phi,opts);
        switch OptimType
            case 'max'
                if(val_tmp > val_best)
                    Pbest = Ptmp;
                    val_best = val_tmp;
                end
                if(val_best>0 && StopWhenFound)
                    return;
                end
            case 'min'
                if(val_tmp < val_best)
                    Pbest = Ptmp;
                    val_best = val_tmp;
                end
                if(val_best<0 && StopWhenFound)
                    return;
                end
            case 'zero'
                if(abs(val_tmp) < abs(val_best))
                    Pbest = Ptmp;
                    val_best = val_tmp;
                end
        end
    end
end

end

