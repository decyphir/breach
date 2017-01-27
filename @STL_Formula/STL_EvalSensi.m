function [val, vald] = STL_EvalSensi(phi, traj, t)
%STL_EVALSENSI function probably deprecated, or not finished.
% 
% NOT RECOMMANDED TO USE
% 
% Synopsis: [val, vald] = STL_EvalSensi(phi, traj, t)
% 
% Inputs:
%   - phi  
%   - traj 
%   - t    time step when we ask for. If it is a scalar, we consider that 
%          the time step start from this time step and finish to the end of
%          the trajectory.
%
% Output:
%    - val
%    - vald
%

if (numel(t)==1)
    
    ind_t = find(traj.time>t,1);
    interval = [t traj.time(ind_t)];
else
    interval = [t(1) t(end)];
end
[valarray, valdarray, time_values] = GetValues(phi,traj,interval);
val = interp1(time_values, valarray, t,'linear',valarray(end));
vald = interp1(time_values, valdarray', t,'linear',0)';

end


function [valarray, valdarray, time_values] = GetValues(phi,traj,interval)
%
% GETVALUES
%
%   Synopsis : [valarray, valdarray, time_values] = STL_EvalSensi(phi, traj, interval)
%
%   Input :
%    - phi
%    - traj
%    - interval is [tmin,tmax]. we want the ??? between tmin and tmax
%
%   Output :
%     - valarray
%     - valdarray
%     - time_values
%

dt = traj.time(2)-traj.time(1);
valdarray = [];

switch (phi.type)
    
    case 'predicate'
        time_values = GetTimeValues(traj,interval);
        evalfn = @(t) phi.evalfn(1,traj,t);
        
        try  % works if predicate can handle multiple time values
            [valarray, valdarray] = evalfn(time_values);
        catch
            [valarray, valdarray] = arrayfun(evalfn, time_values);
        end
        
    case 'not'
        [valarray, valdarray, time_values] = GetValues(phi.phi,traj,interval);
        valarray = - valarray;
        valdarray = - valdarray;
        
    case 'or'
        [valarray1, valdarray1, time_values1] = GetValues(phi.phi1,traj,interval);
        [valarray2, valdarray2, time_values2] = GetValues(phi.phi2,traj,interval);
        [valarray1, valarray2] = SyncValues(valarray1,valarray2);
        [valdarray1, valdarray2] = SyncValues(valdarray1,valdarray2);
        
        time_values = time_values1(1:numel(valarray1));
        [valarray, indmax]  = max([valarray1 ; valarray2]);
        
        valdarray = [valdarray1;  valdarray2];
        NS = size(valdarray1,1);
        indN = repmat([1:size(valdarray1,1)]',[1 numel(indmax)]);
        indN = indN+repmat(2*NS*(0:size(valdarray,2)-1)+(indmax-1)*NS,[2 1]);
        
        valdarray = valdarray(indN);
        
    case 'and'
        [valarray1, valdarray1, time_values1] = GetValues(phi.phi1,traj,interval);
        [valarray2, valdarray2, time_values2] = GetValues(phi.phi2,traj,interval);
        [valarray1, valarray2] = SyncValues(valarray1,valarray2);
        [valdarray1, valdarray2] = SyncValues(valdarray1,valdarray2);
        
        [valarray, indmin]  = min([valarray1 ; valarray2]);
        
        valdarray = [valdarray1;  valdarray2];
        NS = size(valdarray1,1);
        indN = repmat([1:size(valdarray1,1)]',[1 numel(indmin)]);
        indN = indN+repmat(2*NS*(0:size(valdarray,2)-1)+(indmin-1)*NS,[2 1]);
        
        valdarray = valdarray(indN);
        
    case '=>'
        [valarray1, valdarray1, time_values1] = GetValues(phi.phi1,traj,interval);
        [valarray2, valdarray2, time_values2] = GetValues(phi.phi2,traj,interval);
        [valarray1, valarray2] = SyncValues(valarray1,valarray2);
        [valdarray1, valdarray2] = SyncValues(valdarray1,valdarray2);
        
        time_values = time_values1(1:numel(valarray1));
        [valarray, indmax]  = max([-valarray1 ; valarray2]);
        
        valdarray = [valdarray1;  valdarray2];
        NS = size(valdarray1,1);
        indN = repmat([1:size(valdarray1,1)]',[1 numel(indmax)]);
        indN = indN+repmat(2*NS*(0:size(valdarray,2)-1)+(indmax-1)*NS,[2 1]);
        valdarray = valdarray(indN);
        
    case 'always'
        
        next_interval = phi.interval+interval;
        [valarray, valdarray, time_values] = GetValues(phi.phi, traj, next_interval);
        if (phi.interval(2) == inf)
            [valarray, indx] = lim_inf_indx(valarray);
            valdarray = valdarray(indx);
        else
            % Warning: NOT IMPLEMENTED
            win = max(floor((phi.interval(2)-phi.interval(1))/dt),1);
            valarray = minmaxfilt1(valarray,win, 'min');
        end
        time_values = time_values-phi.interval(1);
        if (time_values(1)<0)
            time_values = time_values-time_values(1);
        end
        time_values = time_values(1:numel(valarray));
        
    case 'eventually'
        
        next_interval =  phi.interval+interval;
        [valarray, valdarray, time_values] = GetValues(phi.phi, traj, next_interval);
        
        if (phi.interval(2) == inf)
            [valarray, indx] = lim_inf_indx(-valarray);
            valarray= -valarray;
            valdarray = valdarray(indx);
        else
            
            % Warning: NOT IMPLEMENTED
            win = max(floor((phi.interval(2)-phi.interval(1))/dt),1);
            valarray = minmaxfilt1(valarray,win, 'max');
        end
        time_values = time_values-phi.interval(1);
        if (time_values(1)<0)
            time_values = time_values-time_values(1);
        end
        time_values = time_values(1:numel(valarray));
        
    case 'until'        % Warning: NOT IMPLEMENTED
        interval1  =  [interval(1), phi.interval(2)+interval(2)];
        interval2 =  phi.interval+interval;
        
        [valarray1, time_values1]= GetValues(phi.phi1, traj, interval1);
        [valarray2, time_values2]= GetValues(phi.phi2, traj, interval2);
        
        time_values = time_values1(1:numel(valarray1));
        
        if (phi.interval(2) == inf)
            if (phi.interval(1)==0)
                valarray = until_inf(valarray1,valarray2,0, -1);
            else
                i1_interval= floor(phi.interval(1)/dt);
                minphi1_win = minmaxfilt1(valarray1,i1_interval);
                valarray = until_inf(valarray1,valarray2, i1_interval,-1, minphi1_win);
            end
        else
            i1_interval= floor(phi.interval(1)/dt);
            i2_interval= floor(phi.interval(2)/dt);
            minphi1_win = minmaxfilt1(valarray1,i1_interval);
            valarray = until_inf(valarray1,valarray2, i1_interval,i2_interval, minphi1_win);
        end
        
end

end

function time_values = GetTimeValues(traj,interval)

ind_ti = find(traj.time>=interval(1),1);
ind_tf = find(traj.time>interval(end),1);

if (ind_ti==1)
    ind_ti = 2;
end

if isempty(ind_tf)
    time_values = traj.time(ind_ti-1:end);
else
    time_values = traj.time(ind_ti-1:ind_tf-1);
end

end

function [v1, v2] = SyncValues(v1, v2)
l1 = numel(v1);
l2 = numel(v2);
if l1==l2
    return;
end
if (l1>l2)
    % v2 = [v2 zeros(1,l1-l2)];
    v1 = v1(1:l2);
else
    %   v1 = [v1 zeros(1,l2-l1)];
    v2 = v2(1:l1);
end

end

