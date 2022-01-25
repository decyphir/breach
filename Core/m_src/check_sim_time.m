function tspan = check_sim_time(tspan)
%   tspan = check_sim_time(tspan)
%
%   Checks/enforces the following requirements:
%      - It must be a one dimensional array
%      - tspan(end) should be positive
%      - tspan can be a scalar, in which case it is treated as [0 t_in]
%      - tspan must be strictly increasing. If not, try to sort it and
%        remove duplicates  (with a warning)
%      - tspan(1) should be 0. Adds it with a warning if needed

if  ~all(diff(tspan)>0)
    warning('check_sim_time:non_monotonic_time', 'Simulation time should be strictly monotonic, increasing. Will sort and continue.' )
    tspan = sort(unique(tspan));
end

if tspan(end)<0
    error('check_time:neg_time', 'Cannot set negative time.')
end

if isscalar(tspan)  % standard case
    tspan = [0 tspan];
end

if size(tspan, 1)~=1
    tspan= tspan';
end

if tspan(1) ~= 0
    warning('check_sim_time:time_starts_at_non_zero', 'Simulation time should start at time 0. Adding 0 and trying to continue.')
    tspan = [0 tspan];
end
end
