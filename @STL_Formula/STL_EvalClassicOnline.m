function [val__, tau__] = STL_EvalClassicOnline(Sys, phi, P, trajs, t,dt__)
%STL_EVALCLASSICONLINE computes the satisfaction function interval (lower bound, upper bound) of a property for one
% or several trajectory(ies). This function uses a fixed time step robust monitoring algorithm.
%
% Synopsis: [val__, tau__] = STL_EvalClassicOnline(Sys, phi, P, trajs[, t])
% 
% Inputs:
%  - Sys   : the system
%  - phi   : 
%  - P     : 
%  - trajs : 
%  - t     : (Optional, default=trajs.time for each traj in trajs) if t
%            doesn't belong to traj.time, the truth value is interpolated.
%            If t is a single time point and does not belong to traj.time,
%            the truth value answered is the truth value of the formula at
%            the closest time point belonging to traj.time, greater than t.
% 
% Outputs:
%  - val__ : The truth value interval of the trajs at time points tau__.
%  - tau__ : Time
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

eval_str = [P.ParamList;num2cell(1:numel(P.ParamList))];
eval_str = sprintf('%s=P.pts(%d,1);',eval_str{:});
eval(eval_str);

numTrajs = numel(trajs);
val__ = cell(1, numTrajs);
tau__ = cell(1, numTrajs);
exist_t = exist('t','var'); % boolean indicating if t is provided as arg of STL_EvalClassic
exist_dt = exist('dt__','var');

for ii=1:numTrajs
    if ~exist_t
        t = trajs{ii}.time;
    elseif isempty(t)
        t = trajs{ii}.time;
    end
    
    if exist_dt % we resample - assumes t exists 
        new_time = trajs{ii}.time(1):dt__:trajs{ii}.time(end);
        new_X = interp1(trajs{ii}.time, trajs{ii}.X', new_time, 'linear',inf)';
        trajs{ii}.time=new_time;
        trajs{ii}.X = new_X;
        if (size(new_X, 2)==1)
            trajs{ii}.X = trajs{ii}.X';
        end
    else
        dt__ = trajs{ii}.time(2)-trajs{ii}.time(1);
    end
    
    interval = [t(1), t(end)];
    [valarray, time_values] = GetValues(Sys, phi, P, trajs{ii}, interval,dt__);
    
    if isscalar(valarray)
        val__{ii} = valarray*ones(1, numel(time_values)); 
        tau__{ii} = time_values;
    elseif(numel(time_values) == numel(t) && all(time_values==t))
        val__{ii} = valarray;
        tau__{ii} = time_values;
    else
        val__{ii} = interp1(time_values, valarray', t, 'nearest', inf)'; % nearest interpolation, otherwise we get NaN when interpolating at t when value at the right of t is inf
        tau__{ii} = t;
    end
end

if(numTrajs==1)
    val__ = val__{1};
    tau__ = tau__{1};
end

end

function [valarray__, time_values__] = GetValues(Sys, phi, P, traj, interval, dt__)

global BreachGlobOpt;
eval(BreachGlobOpt.GlobVarsDeclare);

switch (phi.type)
    
    case 'predicate'
        time_values__ = GetTimeValues(traj, interval);
        params = phi.params;
        params.Sys = Sys;
        params.P = P;
        evalfn = @(t) generic_predicate_online(0, traj, t, params);

        valarray__ = evalfn(time_values__);
        valarray__ = [valarray__; valarray__ ];
        
        if isempty(time_values__)
            time_values__ = [interval(1) interval(1)+dt__];
            valarray__ = [inf inf;-inf -inf];        
        elseif (interval(2)>time_values__(end)+eps) 
            time_out__ = time_values__(end)+dt__;
            val_out__ = [inf; -inf];  
            time_values__ = [time_values__ time_out__ ];
            valarray__ = [valarray__ val_out__];           
        end
        
    case 'not'
        [valarray__, time_values__] = GetValues(Sys, phi.phi, P, traj, interval,dt__);
        valarray__ = - valarray__([2 1],:);
        
    case 'or'
        [valarray1, time_values1] = GetValues(Sys, phi.phi1, P, traj, interval,dt__);
        [valarray2, ~] = GetValues(Sys, phi.phi2, P, traj, interval, dt__);
        [valarray1, valarray2] = SyncValues(valarray1, valarray2);
        time_values__ = time_values1(1:size(valarray1,2));
        valarray__ = max(valarray1 ,valarray2);
        
    case 'and'
        [valarray1, time_values1] = GetValues(Sys, phi.phi1, P, traj, interval,dt__);
        [valarray2, ~] = GetValues(Sys, phi.phi2, P, traj, interval,dt__);
        [valarray1, valarray2] = SyncValues(valarray1, valarray2);
        time_values__ = time_values1(1:size(valarray1,2));
        valarray__ = min(valarray1, valarray2);
        
    case 'andn'
        [valarray__, time_values__] = GetValues(Sys, phi.phin(1), P, traj, interval,dt__);
        for ii=2:numel(phi.phin)
            valarray_tmp = GetValues(Sys, phi.phin(ii), P, traj, interval,dt__);
            [valarray__, valarray_tmp] = SyncValues(valarray__, valarray_tmp);
            time_values__ = time_values__(1:size(valarray__,2));
            valarray__ = min(valarray__, valarray_tmp);
        end
        
    case '=>'
        [valarray1, time_values1] = GetValues(Sys ,phi.phi1, P, traj, interval,dt__);
        [valarray2, ~] = GetValues(Sys, phi.phi2, P, traj, interval,dt__);
        [valarray1, valarray2] = SyncValues(valarray1, valarray2);
        
        time_values__ = time_values1(1:size(valarray1,2));
        valarray__ = max(-valarray1, valarray2);
        
    case 'always'
        I__ = eval(phi.interval);        
        I__ = max([I__; 0 0]);            % set negative values to 0
        I__(1) = min(I__(1), I__(2));     % ensures I(2)>= I(1)
        
        next_interval = I__ + interval;   % potentially needed horizon for subformula 
        
        [valarray__, time_values__] = GetValues(Sys, phi.phi, P, traj, next_interval,dt__);
        
        win = max(ceil((I__(2)-I__(1))/dt__)+1, 1);
        
        if (win>size(valarray__,2))
            val_sup  =  lim_inf(valarray__(1,:));
            val_inf =   zeros(1,numel(valarray__(2,:)))-inf;
            valarray__ = [  val_sup ; ...
                            val_inf];
        elseif (win==size(valarray__,2))
            val_sup  =  lim_inf(valarray__(1,:));
            val_inf =   lim_inf(valarray__(2,:));
            valarray__ = [  val_sup ; ...
                            val_inf];
        else
            valarray__ = [valarray__ [zeros(1,win-1)+inf; ...
                                      zeros(1,win-1)-inf] ];
            valarray__ = [minmaxfilt1(valarray__(1,:), win, 'min') ; ...
                          minmaxfilt1(valarray__(2,:), win, 'min') ];
        end
        
        
        time_values__ = time_values__-I__(1);
        if(time_values__(1)<0)
            time_values__ = time_values__-time_values__(1);
        end
        valarray__ = valarray__(:, 1:numel(time_values__));

    case 'eventually'
        I__ = eval(phi.interval);
        I__ = max([I__; 0 0]);
        I__(1) = min(I__(1), I__(2));
        
        next_interval = I__+interval;
        [valarray__, time_values__] = GetValues(Sys, phi.phi, P, traj, next_interval,dt__);
        
        win = max(ceil((I__(2)-I__(1))/dt__)+1,1);

        if (win>size(valarray__,2))
            val_sup =   zeros(1,numel(valarray__(1,:)))+inf;
            val_inf  =  -lim_inf(-valarray__(2,:));
            valarray__ = [  val_sup ; ...
                            val_inf];
        elseif (win==size(valarray__,2))
            val_sup  =  -lim_inf(-valarray__(1,:));
            val_inf =   -lim_inf(-valarray__(2,:));
            valarray__ = [  val_sup ; ...
                            val_inf];
        else
            valarray__ = [valarray__ [zeros(1,win-1)+inf; ...
                                      zeros(1,win-1)-inf] ];
            valarray__ = [minmaxfilt1(valarray__(1,:), win, 'max') ; ...
                          minmaxfilt1(valarray__(2,:), win, 'max') ];
        end
        
        time_values__ = time_values__-I__(1);
        if (time_values__(1)<0)
            time_values__ = time_values__-time_values__(1);
        end
        valarray__ = valarray__(:, 1:numel(time_values__));

    case 'until' % TODO (bounded time)
        I__ = eval(phi.interval);
        I__ = max([I__; 0 0]);
        I__(1) = min(I__(1), I__(2));
        interval1 = [interval(1), I__(2)+interval(2)];
        interval2 = I__+interval;
        
        [valarray1, time_values1]= GetValues(Sys, phi.phi1, P, traj, interval1,dt__);
        [valarray2, ~]= GetValues(Sys, phi.phi2, P, traj, interval2,dt__);

%        [valarray1, valarray2] = SyncValues(valarray1, valarray2);
        time_values__ = time_values1(1:size(valarray1,2));
        
        val_sup = get_until(I__,time_values__, dt__, valarray1(1,:), valarray2(1,:));
        val_inf = get_until(I__,time_values__, dt__, valarray1(2,:), valarray2(2,:));
        
        valarray__ = [  val_sup ; ...
                        val_inf];
        
end

end

function time_values = GetTimeValues(traj, interval)
ind_ti = find(traj.time>= interval(1),1);
ind_tf = find(traj.time> interval(end),1);

if (ind_ti==1)
    ind_ti = 2;
end

if isempty(ind_tf)
    time_values = traj.time(ind_ti-1:end);
else
    time_values = traj.time(ind_ti-1:ind_tf); 
end

if (interval(1)<traj.time(1))
    time_values = [interval(1) time_values];
end

end

function [v1, v2] = SyncValues(v1, v2)

l1 = size(v1, 2);
l2 = size(v2, 2);
if l1==l2
    return;
end
if(l1>l2)
    % v2 = [v2 zeros(1,l1-l2)];
    v1 = v1(:,1:l2);
else
    %   v1 = [v1 zeros(1,l2-l1)];
    v2 = v2(:,1:l1);
end
end

function valarray__ = get_until(I__,time_values__,dt__,valarray1, valarray2)

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

