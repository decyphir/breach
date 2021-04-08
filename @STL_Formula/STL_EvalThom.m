function [val__, time_values__] = STL_EvalThom(Sys_, phi_, P_, trajs_, t)
%STL_EVALTHOM computes the satisfaction function of a property for one
% or many trajectory(ies). This function uses a variable time step robust
% monitoring algorithm.
%
% Synopsis: [val__, time_values__] = STL_EvalThom(Sys, phi, P, trajs[, t])
%
% Input:
%  - Sys   : is the system
%  - phi   : is a STL property
%  - P     : is a parameter set which contains one parameter vector only
%            used for properties parameters.
%  - trajs : is a structure with fields X and time. It may contains many
%            trajectories. In this case, all will be checked considering
%            the property parameters described in P.
%  - t     : (optional, default=traj.time) is the time point(s), so
%            possibly an array, when to eval the satisfaction of the
%            property. All time points not belonging to traj.time will be
%            linearly interpolated.
%
% Output:
%  - val__         : is a cell array of dimension 1 x numel(trajs). Each
%                    cell contains a line array describing the evaluation
%                    of phi at each time of time_values__. If numel(trajs)
%                    is 1, val__ is the content of its only cell to avoid a
%                    useless cell array (so, it is a line array).
%  - time_values__ : is a cell array of dimension 1 x numel(trajs). Each
%                    cell contains an array indicating time point(s) at
%                    which the formula is evaluated or interpolated. If the
%                    parameter t is provided, each cells of time_values__
%                    is equal to t. If numel(trajs) is 1, time_values__
%                    contains the content of the only cell (this avoid a
%                    useless cell array), and thus becomes a line array.
%
%See also STL_Eval SEvalProp
%

%
%  This works with piecewise affine signals.
%


%% defines the parameter as global variables so that they are available for
% all subsequent computations

global BreachGlobOpt;
if ~isempty(P_.ParamList)
    BreachGlobOpt.GlobVarsDeclare = ['global ', sprintf('%s ',P_.ParamList{:})]; % contains parameters and IC values (can remove IC if phi is optimized)
    eval(BreachGlobOpt.GlobVarsDeclare); % These values may be used in generic_predicate and GetValues
else
    BreachGlobOpt.GlobVarsDeclare = ''; % contains parameters and IC values (can remove IC if phi is optimized)
end

ii=1;
num_dim = size(P_.pts,1);
eval_str = [P_.ParamList(1:num_dim);num2cell(1:num_dim)];
eval_str = sprintf('%s=P_.pts(%d,ii);',eval_str{:});
eval(eval_str);

%% for each trajectory, compute values and times

numTrajs = numel(trajs_);
val__ = cell(1, numTrajs);
time_values__ = cell(1, numTrajs);

if isstruct(trajs_)||isa(trajs_, 'matlab.io.MatFile')
    trajs_ = {trajs_};
end

for ii=1:numTrajs % we loop on every traj in case we check more than one
    if (Psize_pts(P_)==1)
        Pii = P_;
    else
        Pii = Sselect(P_, ii);
        eval(eval_str); % needed, as parameters can change from one Pii to another
    end
    
    % Ensures that traj.X and traj.time are double precision
    traj.time = double(trajs_{ii}.time);
    traj.X = double(trajs_{ii}.X);
    
    % Robusthom doesn't like singular intervals - should be optimized one
    % of these days ...
    if exist('t','var')
        time_values__{ii} = t;
        if(numel(t)==1)
            tn = find(traj.time>t,1);
            if isempty(tn)
                interval = [t t+1]; % Maybe not the best choice, but we have to make one !
            else
                interval = [t traj.time(1,tn)];
            end
        else
            interval = [t(1) t(end)];
        end
        [val, time_values] = GetValues(Sys_, phi_, Pii, traj, interval);
        
        try
 %           if(numel(t)==1) % we handle singular times
 %               val__{ii} = val(1);
 %           else
                if isfield(BreachGlobOpt, 'disable_robust_linear_interpolation')&&BreachGlobOpt.disable_robust_linear_interpolation
                    val__{ii} = interp1(time_values, val, t, 'previous');
                else
                    val__{ii} = interp1(time_values, val, t);
                end
 %           end
        catch % if val is empty
            val__{ii} = NaN(1,numel(t));
        end
    else
        interval = [0 traj.time(1,end)];
        [val__ii, time_values__ii] = GetValues(Sys_, phi_, Pii, traj, interval);
        
        val__{ii} = val__ii(time_values__ii<=traj.time(1,end));
        time_values__{ii} = time_values__ii(time_values__ii<=traj.time(1,end));
        
        if (time_values__{ii}(end) < traj.time(end))
            vend__ = interp1(time_values__ii, val__ii,traj.time(1,end));
            time_values__{ii}(end+1) = traj.time(1,end);
            val__{ii}(end+1)= vend__;
        end
        
    end
end

if(numTrajs==1)
    val__ = val__{1};
    time_values__ = time_values__{1};
else
    if numel(t)==1
        val__ = cell2mat(val__);
        time_values__ = cell2mat(time_values__);
    end
end
end
%%

function [valarray, time_values] = GetValues(Sys, phi, P, traj, interval)
global BreachGlobOpt;
eval(BreachGlobOpt.GlobVarsDeclare);

switch(phi.type)
    
    case 'predicate'
        time_values = GetTimeValues(traj, interval);
        params = phi.params;
        params.Sys = Sys;
        params.P = P;
        evalfn = @(t) phi.evalfn(0, traj, t, params); % will call generic_predicate
        
        try   % works if predicate can handle multiple time values
            valarray = evalfn(time_values);
        catch %#ok<CTCH>
            valarray = arrayfun(evalfn, time_values);
        end
        
    case 'not'
        [valarray, time_values] = GetValues(Sys, phi.phi, P, traj, interval);
        valarray = - valarray;
        
    case 'or'
        [valarray1, time_values1] = GetValues(Sys, phi.phi1, P, traj, interval);
        [valarray2, time_values2] = GetValues(Sys, phi.phi2, P, traj, interval);
        
        switch phi.semantics
            case 'max'
                [time_values, valarray] = RobustOr(time_values1, valarray1, time_values2, valarray2);
            case 'add'
                % ||+
                [time_values, valarray] = robustAndPlus(time_values1, -valarray1, time_values2, -valarray2);
                valarray = -valarray;
            case 'vbool_v1'
                [time_values, valarray] = robustAndPlus_v1(time_values1, -valarray1, time_values2, -valarray2);
                valarray = -valarray;
            case 'MARV'
                % On this level, MARV is just standard robustness, since
                % MARV only applies to top-level "always"-operator.
                [time_values, valarray] = RobustOr(time_values1, valarray1, time_values2, valarray2);
            case 'constant'
                [time_values, valarray] = RobustOr(time_values1, valarray1, time_values2, valarray2);
            otherwise
                error('Unknown robust semantics (phi.semantics)');
        end
        
    case 'and'
        [valarray1, time_values1] = GetValues(Sys, phi.phi1, P, traj, interval);
        [valarray2, time_values2] = GetValues(Sys, phi.phi2, P, traj, interval);
        
        % JOHAN CHANGE
        switch phi.semantics
            case 'max'
                % Standard and
                [time_values, valarray] = RobustAnd(time_values1, valarray1, time_values2, valarray2);
            case 'add'
                % Koen's &+
                [time_values, valarray] = robustAndPlus(time_values1, valarray1, time_values2, valarray2);
            case 'vbool_v1'
                % Old additive semantics
                [time_values, valarray] = robustAndPlus_v1(time_values1, valarray1, time_values2, valarray2);
            case 'MARV'
                % On this level, MARV is just standard robustness, since
                % MARV only applies to top-level "always"-operator.
                [time_values, valarray] = RobustAnd(time_values1, valarray1, time_values2, valarray2);
            case 'constant'
                % Standard and
                [time_values, valarray] = RobustAnd(time_values1, valarray1, time_values2, valarray2);
            otherwise
                error('Unknown robust semantics (phi.semantics)');
        end
        
    case 'andn'
        n_phi = numel(phi.phin);
        valarray = cell(1,n_phi);
        time_values = cell(1,n_phi);
        for ii=1:n_phi
            [valarray{ii},time_values{ii}] = GetValues(Sys, phi.phin(ii), P, traj, interval);
        end
        [time_values, valarray] = RobustAndn(time_values,valarray);
        
    case '=>'
        [valarray1, time_values1] = GetValues(Sys, phi.phi1, P, traj, interval);
        [valarray2, time_values2] = GetValues(Sys, phi.phi2, P, traj, interval);
        valarray1 = -valarray1;
        
        switch phi.semantics
            case 'max'
                [time_values, valarray] = RobustOr(time_values1, valarray1, time_values2, valarray2);
            case 'add'
                % Standard implication, but with vbool andPlus
                [time_values, valarray] = robustAndPlus(time_values1, -valarray1, time_values2, -valarray2);
                valarray = -valarray;
            case 'vbool_v1'
                [time_values, valarray] = robustAndPlus_v1(time_values1, -valarray1, time_values2, -valarray2);
                valarray = -valarray;
            case 'MARV'
                % On this level, MARV is just standard robustness, since
                % MARV only applies to top-level "always"-operator.
                [time_values, valarray] = RobustOr(time_values1, valarray1, time_values2, valarray2);
            case 'constant'
                [time_values, valarray] = RobustOr(time_values1, valarray1, time_values2, valarray2);
            otherwise
                error('Unknown robust semantics (phi.semantics)');
        end
        
        
    case 'always'
        I___ = eval(phi.interval);
        I___ = max([I___; 0 0]);
        I___(1) = min(I___(1), I___(2));
        next_interval = I___+interval;
        [valarray, time_values] = GetValues(Sys, phi.phi, P, traj, next_interval);
        % JOHAN FIX
        % valarray is EMPTY if the formula is "true". The valarray is
        % assigned Inf at all time steps, which is then "removed" to
        % prevent unwanted behaviour.
        % Solution: If valarray is empty, set the valarray to be
        % true_value.
        if isempty(valarray)
            if isfield(phi.params.default_params,'true_value__')
                % true_value__ is defined for phi
                valarray = phi.params.default_params.true_value__;
            else
                % true_value__ is NOT defined for phi!!
                warning('true_value__ is not defined for phi! Using true_value__ = 100.')
                valarray = 100;
            end
            time_values = I___(1);
        end
        % END JOHAN FIX
        
        switch phi.semantics
            case 'max'
                if(I___(end)~=inf)
                    time_values = [time_values time_values(end)+I___(end)];
                    valarray = [valarray valarray(end)];
                end
                [time_values, valarray] = RobustEv(time_values, -valarray, I___);
                valarray = -valarray;
            case 'add'
                %[time_values, valarray] = RobustAvEvRight(time_values, -valarray, I___);
                %valarray = -valarray;
                [time_values, valarray] = RobustAlways(time_values, valarray, I___);
            case 'vbool_v1'
                %[time_values, valarray] = RobustAvEvRight(time_values, -valarray, I___);
                %valarray = -valarray;
                [time_values, valarray] = RobustAlways_v1(time_values, valarray, I___);
            case 'MARV'
                % On this level, MARV is just standard robustness, since
                % MARV only applies to top-level "always"-operator.
                if(I___(end)~=inf)
                    time_values = [time_values time_values(end)+I___(end)];
                    valarray = [valarray valarray(end)];
                end
                [time_values, valarray] = RobustEv(time_values, -valarray, I___);
                valarray = -valarray;
            case 'constant'
                if(I___(end)~=inf)
                    time_values = [time_values time_values(end)+I___(end)];
                    valarray = [valarray valarray(end)];
                end
                [time_values, valarray] = RobustEv(time_values, -valarray, I___);
                valarray = -valarray;
            otherwise
                error('Unknown objective function!');
        end
        
    case 'av_eventually'
        I___ = eval(phi.interval);
        I___ = max([I___; 0 0]);
        I___(1) = min(I___(1), I___(2));
        next_interval = I___+interval;
        [valarray1, time_values1] = GetValues(Sys, phi.phi, P, traj, next_interval);
        if(I___(end)~=inf)
            time_values1 = [time_values1 time_values1(end)+I___(end)];
            valarray1 = [valarray1 valarray1(end)];
        end
        [time_values, valarray] = RobustAvEvRight(time_values1, valarray1, I___);
        
    case 'eventually'
        I___ = eval(phi.interval);
        I___ = max([I___; 0 0]);
        I___(1) = min(I___(1), I___(2));
        next_interval = I___+interval;
        [valarray1, time_values1] = GetValues(Sys, phi.phi, P, traj, next_interval);
        
        switch phi.semantics
            case 'max'
                if(I___(end)~=inf)
                    time_values1 = [time_values1 time_values1(end)+I___(end)];
                    valarray1 = [valarray1 valarray1(end)];
                end
                [time_values, valarray] = RobustEv(time_values1, valarray1, I___);
            case 'add'
                [time_values, valarray] = RobustAlways(time_values1, -valarray1, I___);
                valarray = -valarray;
            case 'vbool_v1'
                [time_values, valarray] = RobustAlways_v1(time_values1, -valarray1, I___);
                valarray = -valarray;
            case 'MARV'
                % On this level, MARV is just standard robustness, since
                % MARV only applies to top-level "always"-operator.
                if(I___(end)~=inf)
                    time_values1 = [time_values1 time_values1(end)+I___(end)];
                    valarray1 = [valarray1 valarray1(end)];
                end
                [time_values, valarray] = RobustEv(time_values1, valarray1, I___);
            case 'constant'
                if(I___(end)~=inf)
                    time_values1 = [time_values1 time_values1(end)+I___(end)];
                    valarray1 = [valarray1 valarray1(end)];
                end
                [time_values, valarray] = RobustEv(time_values1, valarray1, I___);
            otherwise
                error('Unknown objective function!');
        end
    
    case 'once'
        I___ = eval(phi.interval);
        I___ = max([I___; 0 0]);
        I___(1) = min(I___(1), I___(2));
        
        
        next_interval = interval-[I___(1)+I___(2), I___(1)];
        next_interval(1) = max(0, next_interval(1));        
        [valarray1, time_values1] = GetValues(Sys, phi.phi, P, traj, next_interval);
        
        % Flipping time, taking into account constant interpolation with
        % previous 
        Tend__ =  time_values1(end)+1; 
        past_time_values1 = fliplr(Tend__-[time_values1 Tend__]);
        past_valarray1 =    fliplr([valarray1(1) valarray1]);               
        
        switch phi.semantics
            case 'add'
                [past_time_values, past_valarray] = RobustAlways(past_time_values1, -past_valarray1, I___);  
                past_valarray = -past_valarray;
            case 'max'
                [past_time_values, past_valarray] = RobustEv(past_time_values1, past_valarray1, I___);  
            case 'vbool_v1'
                [past_time_values, past_valarray] = RobustAlways_v1(past_time_values1, -past_valarray1, I___);  
                past_valarray = -past_valarray;
            case 'MARV'
                % On this level, MARV is just standard robustness, since
                % MARV only applies to top-level "always"-operator.
                [past_time_values, past_valarray] = RobustEv(past_time_values1, past_valarray1, I___);  
            case 'constant'
                [past_time_values, past_valarray] = RobustEv(past_time_values1, past_valarray1, I___);  
            otherwise
                error('Unknown objective function!');
        end

        % Flipping back
        time_values = fliplr(Tend__-[past_time_values Tend__]);
        valarray = fliplr([past_valarray(1) past_valarray]);                  

    case 'historically'
        I___ = eval(phi.interval);
        I___ = max([I___; 0 0]);
        I___(1) = min(I___(1), I___(2));
        
        next_interval = interval-[I___(1)+I___(2), I___(1)];
        next_interval(1) = max(0, next_interval(1));        
        [valarray1, time_values1] = GetValues(Sys, phi.phi, P, traj, next_interval);        
        
        % Flipping time, taking into account constant interpolation with
        % previous 
        Tend__ =  time_values1(end)+1; 
        past_time_values1 = fliplr(Tend__-[time_values1 Tend__]);
        past_valarray1 =    fliplr([valarray1(1) valarray1]);        
        
        switch phi.semantics
            case 'max'
                [past_time_values, past_valarray] = RobustEv(past_time_values1, -past_valarray1, I___);  
            case 'add'
                [past_time_values, past_valarray] = RobustAlways(past_time_values1, past_valarray1, I___);  
                past_valarray = -past_valarray;
            case 'vbool_v1'
                [past_time_values, past_valarray] = RobustAlways_v1(past_time_values1, past_valarray1, I___);  
                past_valarray = -past_valarray;
            case 'MARV'
                % On this level, MARV is just standard robustness, since
                % MARV only applies to top-level "always"-operator.
                [past_time_values, past_valarray] = RobustEv(past_time_values1, -past_valarray1, I___);  
            case 'constant'
                [past_time_values, past_valarray] = RobustEv(past_time_values1, -past_valarray1, I___);  
            otherwise
                error('Unknown objective function!');
        end
        
        % Flipping back
        time_values = fliplr(Tend__-[past_time_values Tend__]);
        valarray = -fliplr([past_valarray(1) past_valarray]);        
       
    case 'until'
        I___ = eval(phi.interval);
        I___ = max([I___; 0 0]);
        I___(1) = min(I___(1), I___(2));
        interval1 = [interval(1), I___(2)+interval(2)];
        interval2 = I___+interval;
        
        [valarray1, time_values1] = GetValues(Sys, phi.phi1, P, traj, interval1);
        [valarray2, time_values2] = GetValues(Sys, phi.phi2, P, traj, interval2);
        if(I___(end)~=inf)
            time_values1 = [time_values1 time_values1(end)+I___(end)];
            valarray1 = [valarray1 valarray1(end)];
            time_values2 = [time_values2 time_values2(end)+I___(end)];
            valarray2 = [valarray2 valarray2(end)];
        end
        [time_values, valarray] = RobustUntil(time_values1, valarray1, time_values2, valarray2, I___);
end

%%  Sanity checks

% time progress
isblocked = (diff(time_values)<=0);
if find(isblocked)
    %  warning('Some time values not strictly increasing for property %s', disp(phi));
    keep = [1 find(~isblocked)+1];
    valarray = valarray(keep);
    time_values = time_values(keep);
end

ibof = isnan(valarray)|isnan(time_values)|isinf(valarray)|isinf(time_values);
if ~isempty(find(ibof, 1))
    val_ok = valarray(~ibof);
    time_ok = time_values(~ibof);
    if ~isempty(val_ok)
        warning('STL_Eval:Inf_or_Nan', 'Some values are NaN or inf for property %s (use warning(''off'', ''STL_Eval:Inf_or_Nan'') to disable warning)', disp(phi));
        if numel(val_ok)==1
            valarray(1,:) = val_ok;
        else
            valarray = interp1(time_ok, val_ok, time_values, 'nearest');
        end
    else
        warning('STL_Eval:Inf_or_Nan', 'All values are NaN or inf for property %s', disp(phi));
    end
end


end

function time_values = GetTimeValues(traj,interval)
%GETTIMEVALUES provides time points belonging to traj.time strictly
% included in interval plus the bounds of interval.
%
% Note: for now it does not deal correctly with operations on time, e.g.
% x[t-a] of x[ g(t) ] where g is some function
%
% Also if we wanted to deal correctly with open or closed intervals that
% would be the place to look into. For now, the interpretation is rather
% that of closed intervals.
%

if(interval(1)==interval(end))
    time_values = interval(1);
    return ;
end

% first time instant
ind_ti = find(traj.time>=interval(1),1);
if isempty(ind_ti)
    if ~isempty(traj.time)
        time_values = [traj.time(1,end) traj.time(1,end)+1];
    else
        time_values = [];
    end
    return
end

if(traj.time(1,ind_ti)==interval(1))
    time_values = traj.time(1,ind_ti);
    ind_ti = ind_ti+1;
else
    time_values = interval(1);
end

% Last time instant
if(interval(end)==inf)
    time_values = [time_values traj.time(1,ind_ti:end)];
else
    ind_tf = find(traj.time >= interval(end),1);
    if isempty(ind_tf)
        time_values = [time_values traj.time(1,ind_ti:end) interval(end)];
    elseif(traj.time(1,ind_tf)==interval(end))
        time_values = [time_values traj.time(1,ind_ti:ind_tf)];
    else
        time_values = [time_values traj.time(1,ind_ti:ind_tf-1) interval(end)];
    end
end

end

