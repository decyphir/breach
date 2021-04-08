function [time_values, valarray] = alw_monitor_MARV(time_values, valarray, I___)
% Calculate standard rho
rho = min(valarray);

% Calculate rho1
rho1 = 0;

for k = 1:length(time_values)-1
    rho1 = rho1 + valarray(k)*(time_values(k+1) - time_values(k));
end
rho1 = rho1 / (time_values(end) - time_values(1));

if rho < 0
    % Do nothing, since we want to keep the valarray as it is
else
    valarray = rho1*ones(size(valarray));
end

end