function [time_values, valarray] = robustOrPlus(time_values1, valarray1, time_values2, valarray2)
% This is the implementation of robustness value according to Koen's "or+"
% modified robustness for "OR"
% Somewhat inspired by function computeAnd in robustness.cpp
% (found in breach/@STL_Formula/private/src/robusthom)

% Store all time instances that exist in the two time arrays time_values1
% and time_values2
time_values = union(time_values1, time_values2);

% Create the valarray that will contain robustness values
valarray = zeros(size(time_values));

idxToRemove = [];

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
for k = 1:length(time_values)
    % Find the index we need in both time_values1 and time_values2
    idx1 = find(time_values1 >= time_values(k), 1);
    idx2 = find(time_values2 >= time_values(k), 1);
    
    % Can't find the index if it is after the whole time range - fix this
    % manually
    if isempty(idx1)
        idx1 = length(time_values1);
    end
    if isempty(idx2)
        idx2 = length(time_values2);
    end
    
    % Store robustness values existing at these times 
    val1 = valarray1(idx1);
    val2 = valarray2(idx2);
    
    % Apply Koen's "or+"
    if (val1 < 0 && val2 < 0) || (val1 > 0 && val2 > 0)
        valarray(k) = val1 + val2;
    else
        valarray(k) = max(val1, val2);
    end
    
    % If the robustness value is the same as the last time instant, it is
    % not necessary and can be removed
    if k > 1
        if valarray(k) == valarray(k-1)
            idxToRemove = [idxToRemove k];
        end
    end
end

% Remove the values and times not needed
% This is currently commented becuase: It can result into having only ONE
% element in valarray and time_values, this yields an error when Breach
% uses RobustOr! (mexw64)
%valarray(idxToRemove) = [];
%time_values(idxToRemove) = [];

end