function [val__, tau__] = STL_EvalClassic(Sys, phi, P, trajs, t,dt__)
%STL_EVALCLASSIC computes the satisfaction function of a property for one
% or many trajectory(ies). This function uses a fixed time step robust
% monitoring algorithm. NOTE: if the horizon of phi is outside of traj.time, 
% extrapolates the trajectories with constant values to ]-inf, +inf[. 
%
% Synopsis: [val__, tau__] = STL_Eval(Sys, phi, P, trajs[, t,dt])
% 
% Inputs:
%  - Sys   : the system
%  - phi   : STL formula
%  - P     : Parameter set (might include property parameters) 
%  - trajs : array of trajectories 
%  - t     : (Optional, default=trajs.time for each traj in trajs) if t
%            doesn't belong to traj.time, the truth value is interpolated.
%            If t is a single time point and does not belong to traj.time,
%            the truth value answered is the truth value of the formula at
%            the closest time point belonging to traj.time, greater than t.
% 
% Outputs:
%  - val__ : The truth value of the trajs at time points tau__.
%  - tau__ : 
%
%See also STL_Eval SEvalProp
%


global BreachGlobOpt;

try %#ok<TRYNC>
    if BreachGlobOpt.AltEval
        warning('STL_EvalClassic:otherMethod',...
            'STL_EvalClassic is called, but as the system global variable AltEval is true, STL_EvalThom is used.');
        [val__, tau__, status__] = STL_EvalThom(Sys, phi, P, trajs, t);
        return;
    end
end

BreachGlobOpt.GlobVarsDeclare = ['global ', sprintf('%s ',P.ParamList{:})]; % contains parameters and IC values (can remove IC if phi is optimized)
eval(BreachGlobOpt.GlobVarsDeclare); % These values may be used in generic_predicate and GetValues

ii =1;
eval_str = [P.ParamList;num2cell(1:numel(P.ParamList))];
eval_str = sprintf('%s=P.pts(%d,ii);',eval_str{:});
eval(eval_str);

numTrajs = numel(trajs);
val__ = cell(1, numTrajs);
tau__ = cell(1, numTrajs);
exist_t = exist('t','var'); % boolean indicating if t is provided as arg of STL_EvalClassic
exist_dt = exist('dt__','var');


for ii=1:numTrajs
    
    if (Psize_pts(P)==1)
        Pii = P;
    else
        Pii = Sselect(P, ii);
        eval(eval_str); % needed, as parameters can change from one Pii to another        
    end
    
    % Ensures that traj.X and traj.time are double precision
    trajs(ii).time = double(trajs(ii).time);
    trajs(ii).X = double(trajs(ii).X);

    
    if ~exist_t
        t = trajs(ii).time;
    elseif isempty(t)
        t = trajs(ii).time;
    end
    
    if exist_dt % we resample - assumes t exists 
        new_time = trajs(ii).time(1):dt__:trajs(ii).time(end);
        new_X = interp1(trajs(ii).time, trajs(ii).X', new_time, 'linear',inf)';
        trajs(ii).time=new_time;
        trajs(ii).X = new_X;
        if (size(new_X, 2)==1)
            trajs(ii).X = trajs(ii).X';
        end
    else
        dt__ = trajs(ii).time(2)-trajs(ii).time(1);
    end

    interval = [0 t(end)];
    [valarray, time_values] = GetValues(Sys, phi, P, trajs(ii), interval,dt__);
    
    if isscalar(valarray)
        val__{ii} = valarray*ones(1, numel(time_values)); % Why valarray could be a scalar and not time_values ?
        tau__{ii} = time_values;
    elseif(numel(time_values) == numel(t) && all(time_values==t))
        val__{ii} = valarray;
        tau__{ii} = time_values;
    else
        val__{ii} = interp1(time_values, valarray, t, 'nearest', valarray(end));
        tau__{ii} = t;
    end
end

if(numTrajs==1)
    val__ = val__{1};
    tau__ = tau__{1};
end

end


function [valarray__, time_values__] = GetValues(Sys, phi, P, traj, interval,dt__)

global BreachGlobOpt;
eval(BreachGlobOpt.GlobVarsDeclare);

switch (phi.type)
    
   case 'predicate'
        
       params = phi.params;
       params.Sys = Sys;
       params.P = P;
       evalfn = @(t) phi.evalfn(0, traj, t, params);

       
       
       time_values__ = GetTimeValues(traj, interval);
        
      
        
        
        
        
        
        valarray__ = evalfn(time_values__);    
        
        
        
    case 'not'
        [valarray__, time_values__] = GetValues(Sys, phi.phi, P, traj, interval,dt__);
        valarray__ = - valarray__;
        
    case 'or'
        [valarray1, time_values1] = GetValues(Sys, phi.phi1, P, traj, interval,dt__);
        [valarray2, ~] = GetValues(Sys, phi.phi2, P, traj, interval,dt__);
        [valarray1, valarray2] = SyncValues(valarray1, valarray2);
        time_values__ = time_values1(1:numel(valarray1));
        valarray__ = max(valarray1 ,valarray2);
        
    case 'and'
        [valarray1, time_values1] = GetValues(Sys, phi.phi1, P, traj, interval,dt__);
        [valarray2, ~] = GetValues(Sys, phi.phi2, P, traj, interval,dt__);
        [valarray1, valarray2] = SyncValues(valarray1, valarray2);
        time_values__ = time_values1(1:numel(valarray1));
        valarray__ = min(valarray1, valarray2);
        
    case 'andn'
        [valarray__, time_values__] = GetValues(Sys, phi.phin(1), P, traj, interval,dt__);
        for ii=2:numel(phi.phin)
            valarray_tmp = GetValues(Sys, phi.phin(ii), P, traj, interval,dt__);
            [valarray__, valarray_tmp] = SyncValues(valarray__, valarray_tmp);
            time_values__ = time_values__(1:numel(valarray__));
            valarray__ = min(valarray__, valarray_tmp);
        end
        
    case '=>'
        [valarray1, time_values1] = GetValues(Sys ,phi.phi1, P, traj, interval,dt__);
        [valarray2, ~] = GetValues(Sys, phi.phi2, P, traj, interval,dt__);
        [valarray1, valarray2] = SyncValues(valarray1, valarray2);
        
        time_values__ = time_values1(1:numel(valarray1));
        valarray__ = max(-valarray1, valarray2);
        
    case 'always'
        I__ = eval(phi.interval);
        
        I__ = max([I__; 0 0]);
        I__(1) = min(I__(1), I__(2));
        
        next_interval = I__+interval;
        [valarray__, time_values__] = GetValues(Sys, phi.phi, P, traj, next_interval,dt__);
        
        win = max(ceil((I__(2)-I__(1))/dt__)+1, 1);
        
        if(win>=numel(valarray__))
            valarray__ = lim_inf(valarray__);
        else
            valarray__ = [valarray__ valarray__(end)*ones(1,win-1)];
            valarray__ = minmaxfilt1(valarray__, win, 'min');
        end
        time_values__ = time_values__-I__(1);
        if(time_values__(1)<0)
            time_values__ = time_values__-time_values__(1);
        end
        time_values__ = time_values__(1:numel(valarray__));
        
    case 'eventually'
        I__ = eval(phi.interval);
        I__ = max([I__; 0 0]);
        I__(1) = min(I__(1), I__(2));
        
        
        next_interval = I__+interval;
        [valarray__, time_values__] = GetValues(Sys, phi.phi, P, traj, next_interval,dt__);
        
        win = max(ceil((I__(2)-I__(1))/dt__)+1,1);
        
        if (win>=numel(valarray__))
            valarray__ = -lim_inf(-valarray__);
        else
            valarray__ = [valarray__ valarray__(end)*ones(1,win-1)];
            valarray__ = minmaxfilt1(valarray__, win, 'max');
        end
        
        time_values__ = time_values__-I__(1);
        if (time_values__(1)<0)
            time_values__ = time_values__-time_values__(1);
        end
        time_values__ = time_values__(1:numel(valarray__));
        
    case 'until'
        I__ = eval(phi.interval);
        I__ = max([I__; 0 0]);
        I__(1) = min(I__(1), I__(2));
        interval1 = [interval(1), I__(2)+interval(2)];
        interval2 = I__+interval;
        
        [valarray1, time_values1]= GetValues(Sys, phi.phi1, P, traj, interval1,dt__);
        [valarray2, ~]= GetValues(Sys, phi.phi2, P, traj, interval2,dt__);
        [valarray1, valarray2] = SyncValues(valarray1, valarray2);
        time_values__ = time_values1(1:numel(valarray1));
        
        if (I__(2) == inf)
            if (I__(1)==0)
                valarray__ = until_inf(valarray1, valarray2, 0, -1);
            else
                i1_interval= max(ceil(I__(1)/dt__)+1,1);
                valarray1 = [valarray1 valarray1(end)*ones(1, i1_interval-1)];
                
                minphi1_win = minmaxfilt1(valarray1, i1_interval);
                valarray__ = until_inf(valarray1, valarray2, i1_interval, -1, minphi1_win);
            end
        else
            i1_interval= max(ceil(I__(1)/dt__)+1, 1);
            i2_interval= ceil(I__(2)/dt__)+1;
            
            add_N = max(i1_interval, i2_interval)-1;
            valarray1 = [valarray1 valarray1(end)*ones(1,add_N)];
            valarray2 = [valarray2 valarray2(end)*ones(1,add_N)];
            
            minphi1_win = minmaxfilt1(valarray1,i1_interval);
            
            valarray__ = until_inf(valarray1, valarray2, i1_interval, i2_interval, minphi1_win);
            valarray__ = valarray__(1:numel(time_values__));
        end
        
end

end

function time_values = GetTimeValues(traj, interval)
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

if(traj.time(ind_ti)==interval(1))
    time_values = traj.time(ind_ti);
    ind_ti = ind_ti+1;
else
    time_values = interval(1);
end

% Last time instant
if(interval(end)==inf)
    time_values = [time_values traj.time(ind_ti:end)];
else
    ind_tf = find(traj.time >= interval(end),1);
    if isempty(ind_tf)
        time_values = [time_values traj.time(ind_ti:end) interval(end)];
    elseif(traj.time(ind_tf)==interval(end))
        time_values = [time_values traj.time(ind_ti:ind_tf)];
    else
        time_values = [time_values traj.time(ind_ti:ind_tf-1) interval(end)];
    end
end

end
    
    
    
function [v1, v2] = SyncValues(v1, v2)
l1 = numel(v1);
l2 = numel(v2);
if l1==l2
    return;
end
if(l1>l2)
    % v2 = [v2 zeros(1,l1-l2)];
    v1 = v1(1:l2);
else
    %   v1 = [v1 zeros(1,l2-l1)];
    v2 = v2(1:l1);
end

end

