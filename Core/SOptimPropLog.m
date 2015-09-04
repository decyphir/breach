function [Pbest,val_best] = SOptimPropLog(Sys, phi, opts, verbose)
%SOPTIMPROPLOG computes either the best parameter vector found if it is
% negative or all parameter vector leading to a positive evaluation of phi.
% If opts.StopWhenFound is set to 1, the function stops as soon as it finds
% a parameter set positive (max) or negative (min). The value of the
% parameter not in opts.params are set equal to Sys.p.
%
% Synopsis: [Pbest, val_best] = SOptimPropLog(Sys, phi, opts [, verbose])
%
% Inputs:
%   - Sys  : the considered system
%   - phi  : the STL formula to verify
%   - opts : describes the options. It contains the following fields:
%      - tspan      The time domain computation of the trajectories. If
%                   not provided, either Sys must have a tspan field, or P
%                   must contains computed trajectories. Otherwise, an
%                   error is thrown.
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
%
%See also SOptimProp
%

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
    if(tspan(1) > tau)
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

if(OptimizeNInitPoints==1 || OptimizeNBoxes==1)
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

phi = STL_OptimizePredicates(Sys,phi); % optimization of the predicates
if strcmp(OptimType,'min')
    val_best = inf;
else
    val_best = -inf;
end


%% OptimizeGlobally
if(OptimizeGlobally)
    if(verbose>=1)
        fprintf('\n  Global optimization :\n');
    end
    Pall = CreateParamSet(Sys,opts.params,[opts.lbound;opts.ubound]');
    [val_best,Pbest] = SOptimProp(Sys,Pall,phi,opts);
    
    if(StopWhenFound) % check if we have found a solution
        if(strcmp(OptimType,'max') && val_best>0)
            return ;
        elseif(strcmp(OptimType,'min') && val_best<0)
            return ;
        end
    end
end

%% OptimizeNInitPoints
if(OptimizeNInitPoints)
    if(verbose>=1)
        fprintf('\n  Global optimization with %d initial values:\n',nbSplit);
    end
    Prl = CreateRandomLogParamSets(Sys,opts.params,[opts.lbound;opts.ubound]',nbSplit);
    [val_tmp,Ptmp] = SOptimProp(Sys,Prl,phi,opts);
    
    if strcmp(OptimType,'max') % keep the positive ones or the best if all negative
        idx_val = find(val_tmp>0);
        if ~isempty(idx_val)
            if(val_best>0)
                val_best = [val_best,val_tmp(idx_val)];
                Pbest = SConcat(Pbest,Sselect(Ptmp,idx_val));
            else
                val_best = val_tmp(idx_val);
                Pbest = Sselect(Ptmp,idx_val);
            end
        elseif(val_best<=0)% no positive in val_tmp, compare only if val_best is negative
            [val_tmp,idx_tmp] = max(val_tmp);
            if(val_tmp>val_best)
                val_best = val_tmp;
                Pbest = Sselect(Ptmp,idx_tmp);
            end
        end
    elseif strcmp(OptimType,'min') % keep the negative ones or the lowest if all positive
        idx_val = find(val_tmp<0);
        if ~isempty(idx_val)
            if(val_best<0)
                val_best = [val_best,val_tmp(idx_val)];
                Pbest = SConcat(Pbest,Sselect(Ptmp,idx_val));
            else
                val_best = val_tmp(idx_val);
                Pbest = Sselect(Ptmp,idx_val);
            end
        elseif(val_best>=0)% no negative in val_tmp, compare only if val_best is positive
            [val_tmp,idx_tmp] = max(val_tmp);
            if(val_tmp<val_best)
                val_best = val_tmp;
                Pbest = Sselect(Ptmp,idx_tmp);
            end
        end
    end
    
    if StopWhenFound % check if we have found a solution
        if(strcmp(OptimType,'max') && val_best(1)>0)
            return ;
        elseif(strcmp(OptimType,'min') && val_best(1)<0)
            return ;
        end
    end
end

%% OptimizeNBoxes
if(OptimizeNBoxes)
    if(verbose>=1)
        fprintf('\n  Global optimization in %d sub-spaces.\n',nbSplit);
    end
    Prl = CreateRandomLogParamSets(Sys,opts.params,[opts.lbound;opts.ubound]',nbSplit);
    try
        %TODO : IF opt.StopWhenFoundInit IS SET TO 1, COMPUTE TRAJECTORIES
        % ONE BY ONE (see SOptimProp)
        Prl = ComputeTraj(Sys,Prl,tspan);
    catch %#ok<CTCH>
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
    
    kk=0;
    for ii=iv
        kk=kk+1;
        if(verbose>=2)
            fprintf('\n  Optimization of subspace %d/%d :\n', kk, numel(iv));
        end
        Ptmp = Sselect(Prl,ii);
        val = GetParam(Ptmp,params);
        epsi = GetEpsi(Ptmp,params);
        opts.lbound = (val-epsi)';
        opts.ubound = (val+epsi)';
        [val_tmp,Ptmp] = SOptimProp(Sys,Ptmp,phi,opts);
        switch OptimType
            case 'max'
                if(val_tmp>0) % found a valid parameter set
                    if(val_best(1)>0) % val_best may be an array
                        Pbest = SConcat(Pbest,Ptmp);
                        val_best = [val_best,val_tmp]; %#ok<AGROW>
                    else
                        Pbest = Ptmp;
                        val_best = val_tmp;
                    end
                elseif(val_best<=0 && val_tmp>val_best)
                    Pbest = Ptmp; % Pbest and Ptmp are both negative, keep the best one
                    val_best = val_tmp;
                end
                if(val_best(1)>0 && StopWhenFound)
                    return;
                end
            case 'min'
                if(val_tmp<0) % found a valid parameter set
                    if(val_best(1)<0)
                        Pbest = SConcat(Pbest,Ptmp);
                        val_best = [val_best,val_tmp]; %#ok<AGROW>
                    else
                        Pbest = Ptmp;
                        val_best = val_tmp;
                    end
                elseif(val_best>=0 && val_tmp<val_best)
                    Pbest = Ptmp; % Pbest and Ptmp are not valid, keep the best one
                    val_best = val_tmp;
                end
                if(val_best(1)<0 && StopWhenFound)
                    return;
                end
            case 'zero'
                if(abs(val_tmp) < abs(val_best)) % keep the closest one to zero
                    Pbest = Ptmp;
                    val_best = val_tmp;
                end
        end
    end
end

end

