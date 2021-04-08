function [time_values, valarray] = robustAndPlus(time_values1, valarray1, time_values2, valarray2)
% This is the implementation of robustness value according to Koen's "&+"
% modified robustness for "AND"
% Somewhat inspired by function computeAnd in robustness.cpp
% (found in breach/@STL_Formula/private/src/robusthom)

% Store all time instances that exist in the two time arrays time_values1
% and time_values2
time_values = union(time_values1, time_values2);

% Create the valarray that will contain robustness values
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

% Loop over all time and calculate robustness
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
    
    % Store robustness values existing at these times 
    val1 = valarray1(idx1new);
    val2 = valarray2(idx2new);
    
    % Apply Koen's &+
    if val1 < 0 && val2 < 0
        valarray(k) = val1 + val2;
    elseif val1 > 0 && val2 > 0
        valarray(k) = 1/((1/val1) + (1/val2));
    else
        valarray(k) = min(val1, val2);
    end
    
end

end