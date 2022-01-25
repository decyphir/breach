function [trajsSensitive, time_values__] = calculate_structural_sensitivity(Sys, phi, P, traj1, traj2, partition, relabs, t)


%% defines the parameter as global variables so that they are available for
% all subsequent computations

global BreachGlobOpt;
if ~isempty(P.ParamList)
    BreachGlobOpt.GlobVarsDeclare = ['global ', sprintf('%s ',P.ParamList{:})]; % contains parameters and IC values (can remove IC if phi is optimized)
    eval(BreachGlobOpt.GlobVarsDeclare); % These values may be used in generic_predicate and GetValues
else
    BreachGlobOpt.GlobVarsDeclare = ''; % contains parameters and IC values (can remove IC if phi is optimized)
end

ii=1;
num_dim = size(P.pts,1); 
eval_str = [P.ParamList(1:num_dim);num2cell(1:num_dim)];
eval_str = sprintf('%s=P.pts(%d,ii);',eval_str{:});
eval(eval_str);

%% for each trajectory, compute values and times

numTrajs = numel(traj1);
val__ = cell(1, numTrajs);
time_values__ = cell(1, numTrajs);

if isstruct(traj1)||isa(traj1, 'matlab.io.MatFile')
    traj1 = {traj1};
    traj2 = {traj2};
end

for ii=1:numTrajs % we loop on every traj in case we check more than one
    if (Psize_pts(P)==1)
        Pii = P;
    else
        Pii = Sselect(P, ii);
        eval(eval_str); % needed, as parameters can change from one Pii to another
    end
    
    % Ensures that traj.X and traj.time are double precision
    traj.time = double(traj1{ii}.time);
    traj.X = double(traj1{ii}.X);
    
    % Add information from the second traj aswell
    assert(all(traj1{ii}.time == traj2{ii}.time));
    traj.X2 = double(traj2{ii}.X);
    
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
        trajsSensitive = GetValues(Sys, phi, Pii, traj, partition, relabs, interval);
        
        try
 %           if(numel(t)==1) % we handle singular times
 %               val__{ii} = val(1);
 %           else
                if isfield(BreachGlobOpt, 'disable_robust_linear_interpolation')&&BreachGlobOpt.disable_robust_linear_interpolation
                    val__{ii} = interp1(time_values, val, t, 'previous');
                else
                    val__{ii} = interp1(time_values, val, t);
                end
%            end
        catch % if val is empty
            val__{ii} = NaN(1,numel(t));
        end
    else
        error('Not implemented yet');
        interval = [0 traj.time(1,end)];
        [val__ii, time_values__ii] = GetValues(Sys, phi, Pii, traj, partition, relabs, interval);
        
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

function [trajsSensitive, time_values] = GetValues(Sys, phi, P, traj, partition, relabs, interval)
% trajsSensitive shows, for each time value, where the subformulas do NOT
% have the same value. 1 - not same value. 0 - same value. 
global BreachGlobOpt;
eval(BreachGlobOpt.GlobVarsDeclare);

switch(phi.type)
    
    case 'predicate'
        time_values = GetTimeValues(traj, interval);
        params = phi.params;
        params.Sys = Sys;
        params.P = P;
        
        % Create the second trajectory
        traj2.time = traj.time;
        traj2.X = traj.X2;
        evalfn1 = @(t) phi.evalfn(0, traj, t, params); % will call generic_predicate
        evalfn2 = @(t) phi.evalfn(0, traj2, t, params); % will call generic_predicate
        
        try   % works if predicate can handle multiple time values
            valarrayTraj1 = evalfn1(time_values);
            valarrayTraj2 = evalfn2(time_values);
        catch err %#ok<CTCH>
            throw(err);
            valarray = arrayfun(evalfn, time_values);
        end
        
        trajsSensitive = (valarrayTraj1 ~= valarrayTraj2);
    
    case 'not'
        % Negating does not affect sensitivity
        [trajsSensitive, time_values] = GetValues(Sys, phi.phi, P, traj, partition, relabs, interval);
        
    case 'or'
        [trajsSensitive1, time_values1] = GetValues(Sys, phi.phi1, P, traj, partition, relabs, interval);
        [trajsSensitive2, time_values2] = GetValues(Sys, phi.phi2, P, traj, partition, relabs, interval);
        
        % Sensitivity either in first traj OR second traj gives sensitivity
        % of disjunction
        [trajsSensitive, time_values] = StaticSensitivity(time_values1, trajsSensitive1, time_values2, trajsSensitive2);
        
    case 'and'
        [trajsSensitive1, time_values1] = GetValues(Sys, phi.phi1, P, traj, partition, relabs, interval);
        [trajsSensitive2, time_values2] = GetValues(Sys, phi.phi2, P, traj, partition, relabs, interval);
        
        % Sensitivity either in first traj OR second traj gives sensitivity
        % of conjunction
        [trajsSensitive, time_values] = StaticSensitivity(time_values1, trajsSensitive1, time_values2, trajsSensitive2);
        
    case 'andn'
        error('not implemented');
        
    case '=>'
        error('not implemented');
        
    case 'always'
        I___ = eval(phi.interval);
        I___ = max([I___; 0 0]);
        I___(1) = min(I___(1), I___(2));
        next_interval = I___+interval;
        [trajsSensitive1, time_values] = GetValues(Sys, phi.phi, P, traj, partition, relabs, next_interval);
        
        [trajsSensitive, time_values] = TimedSensitivity(time_values, trajsSensitive1, I___);
        
    case 'av_eventually'
        error('Not implemented');
        
    case 'eventually'
        5;
        I___ = eval(phi.interval);
        I___ = max([I___; 0 0]);
        I___(1) = min(I___(1), I___(2));
        next_interval = I___+interval;
        [valarray1, time_values1] = GetValues(Sys, phi.phi, P, traj, partition, relabs, next_interval);
        
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
        
        [trajsSensitive1, time_values1] = GetValues(Sys, phi.phi, P, traj, partition, relabs, next_interval);  
        
        % Flipping time, taking into account constant interpolation with
        % previous 
        Tend__ =  time_values1(end)+1; 
        past_time_values1 = fliplr(Tend__-[time_values1 Tend__]);
        past_valarray1 =    fliplr([trajsSensitive1(1) trajsSensitive1]);  
        
        [past_valarray, past_time_values] = TimedSensitivity(past_time_values1, past_valarray1, I___); 
        
        % Flipping back
        time_values = fliplr(Tend__-[past_time_values Tend__]);
        trajsSensitive = fliplr([past_valarray(1) past_valarray]);  
        
    case 'historically'
        I___ = eval(phi.interval);
        I___ = max([I___; 0 0]);
        I___(1) = min(I___(1), I___(2));
        
        next_interval = interval-[I___(1)+I___(2), I___(1)];
        next_interval(1) = max(0, next_interval(1));        
                
        [trajsSensitive1, time_values1] = GetValues(Sys, phi.phi, P, traj, partition, relabs, next_interval);   
        
        % Flipping time, taking into account constant interpolation with
        % previous 
        Tend__ =  time_values1(end)+1; 
        past_time_values1 = fliplr(Tend__-[time_values1 Tend__]);
        past_valarray1 =    fliplr([trajsSensitive1(1) trajsSensitive1]);  
        
        [past_valarray, past_time_values] = TimedSensitivity(past_time_values1, past_valarray1, I___); 
        
        % Flipping back
        time_values = fliplr(Tend__-[past_time_values Tend__]);
        trajsSensitive = fliplr([past_valarray(1) past_valarray]);  

    case 'until'
        5;
        I___ = eval(phi.interval);
        I___ = max([I___; 0 0]);
        I___(1) = min(I___(1), I___(2));
        interval1 = [interval(1), I___(2)+interval(2)];
        interval2 = I___+interval;
        
        [valarray1, time_values1] = GetValues(Sys, phi.phi1, P, traj, partition, relabs, interval1);
        [valarray2, time_values2] = GetValues(Sys, phi.phi2, P, traj, partition, relabs, interval2);
        if(I___(end)~=inf)
            time_values1 = [time_values1 time_values1(end)+I___(end)];
            valarray1 = [valarray1 valarray1(end)];
            time_values2 = [time_values2 time_values2(end)+I___(end)];
            valarray2 = [valarray2 valarray2(end)];
        end
        [time_values, valarray] = RobustUntil(time_values1, valarray1, time_values2, valarray2, I___);
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

function [valarray, time_values] = StaticSensitivity(time_values1, valarray1, time_values2, valarray2)
% Like robustAndPlus, but for sensitivity instead of robustness

% Store all time instances that exist in the two time arrays time_values1
% and time_values2
time_values = union(time_values1, time_values2);

% Create the valarray that will contain sensitivity values
valarray = zeros(size(time_values));

% Handle cases if time values are nonexistent!
if isempty(time_values1)
    time_values = time_values2;
    valarray = valarray2;
    return
end

if isempty(time_values2)
    time_values = time_values1;
    valarray = valarray1;
    return
end

% Loop over all times
idx1new = 1;
idx2new = 1;
for k = 1:length(time_values)
    % Find the index we need in time_values1
    if idx1new == numel(time_values1)
        % Do nothing
    elseif time_values1(idx1new+1) <= time_values(k)
        idx1new = idx1new + 1;
    end
    
    % Find the index we need in time_values2
    if idx2new == numel(time_values2)
        % Do nothing
    elseif time_values2(idx2new+1) <= time_values(k)
        idx2new = idx2new + 1;
    end
    
    val1 = valarray1(idx1new);
    val2 = valarray2(idx2new);
    
    % If EITHER is sensitive, the output is sensitive
    valarray(k) = val1 || val2;
    
end

end

function [trajsSensitive, time_values] = TimedSensitivity(time_values, valarray, I___)
% Like RobustAlways, but for sensitivity values instead of robustness
% values

time_values = time_values - I___(1);
I___ = I___ - I___(1);

startIdxNew = 1;
endIdxNew = 1;

indicesToRemove = [];

for valIndex = 1:length(valarray)
    thisTime = time_values(valIndex);
    
    % Start time is currentTime + first element of I___
    startTime = thisTime + I___(1);
    
    % End time is startTime + second element of I___
    endTime = startTime + I___(2);
    
    % New getting partialTime and partialValarray
    for startCounter = startIdxNew:numel(time_values)
        if time_values(startCounter) >= startTime
            startIdxNew = startCounter;
            break;
        end
    end
    
    for endCounter = endIdxNew:numel(time_values)
        if time_values(endCounter) - endTime > 1e-10
            endIdxNew = max(endCounter-1, 1);
            break;
        end
    end
    
    if endCounter == numel(time_values)
        endIdxNew = endCounter;
    end
    
    % There is no point beyond the start point
    % => Set end index to same as start index
    if endIdxNew < startIdxNew
        endIdxNew = startIdxNew;
    end
    
    timeIntervalIndexNew = startIdxNew:endIdxNew;
    partialTimeNew = time_values(timeIntervalIndexNew);
    partialValarrayNew = valarray(timeIntervalIndexNew);
    
    partialRob = PartialTimedSensitivity(partialTimeNew, partialValarrayNew);
    trajsSensitive(valIndex) = partialRob;
    
    if valIndex > 2
        if time_values(valIndex) == time_values(valIndex - 1) ...
                && valarray(valIndex) == valarray(valIndex - 1)
            indicesToRemove = [indicesToRemove valIndex]; %#ok<*AGROW>
        end
    end
    
end
% 
% % Remove subsequent indices that are identical in both time and value
% valarray(indicesToRemove) = [];
% time_values(indicesToRemove) = [];

end

function partialSensitivity = PartialTimedSensitivity(time_values, valarray)

partialSensitivity = any(valarray);

end

