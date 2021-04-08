function [time_values, valarray] = alw_monitor_RobustAlways(time_values, valarray, I___)

indicesToRemove = [];
time_values = time_values - I___(1);
I___ = I___ - I___(1);

startIdxNew = 1;
endIdxNew = 1;

for valIndex = 1:length(valarray)
    thisTime = time_values(valIndex);
    
    % Start time is currentTime + first element of I___
    startTime = thisTime + I___(1);
    
    % End time is startTime + second element of I___
    endTime = startTime + I___(2);
    
%     timeIntervalIndex = find(time_values >= startTime & time_values <= endTime);
    
%     partialTime = time_values(timeIntervalIndex);
%     partialValarray = valarray(timeIntervalIndex);
    
    % New getting partialTime and partialValarray
    for startCounter = startIdxNew:numel(time_values)
        if time_values(startCounter) >= startTime
            startIdxNew = startCounter;
            break;
        end
    end
    
    for endCounter = endIdxNew:numel(time_values)
        if time_values(endCounter) > endTime
            endIdxNew = endCounter - 1;
            break;
        elseif endCounter == numel(time_values) && time_values(endCounter) <= endTime
            endIdxNew = endCounter;
        end
    end
    
    timeIntervalIndexNew = startIdxNew:endIdxNew;
    partialTimeNew = time_values(timeIntervalIndexNew);
    partialValarrayNew = valarray(timeIntervalIndexNew);
%     assert(all(timeIntervalIndex == timeIntervalIndexNew));
    
%     if isempty(partialTime)
%         % There is no signal in the time we are looking at
%         % We have to find the signal value from the point defined before
%         % this one
%         timeIntervalIndex = find(time_values <= endTime, 1, 'last');
%         partialTime = time_values(timeIntervalIndex);
%         partialValarray = valarray(timeIntervalIndex);
%     end
    
    partialRob = PartialRobustAlways(partialTimeNew, partialValarrayNew);
    
    % Check if we can reduce the size of time_values and valarray by
    % removing this element (if it is the same as the element before
    if valIndex == 1 || valIndex == length(valarray)
        % We cannot remove the first or the last element
        valarray(valIndex) = partialRob;
    elseif valarray(valIndex - 1) == partialRob
        % The value is the same as the previous one
        % We don't need to add a new value - instead, we can remove this
        % time step from the time_values and valarray vectors later on
        indicesToRemove(end+1) = valIndex; %#ok<*AGROW>
    else
        valarray(valIndex) = partialRob;
    end
    
end

% Remove the indices that are just "duplicates" of the value at the
% previous time stamp
time_values(indicesToRemove) = [];
valarray(indicesToRemove) = [];

end

function partialRob = PartialRobustAlways(time_values, valarray)
% Handle the case with only ONE value
% (We have no time to integrate over -> division by zero)
% Fix this by just returning the single value
if numel(valarray) == 1
    partialRob = valarray;
    return
end


rho = min(valarray(1:end));

if rho < 0
    % We are trying to speed up the process using vectorized calculations!
    timeDiff = diff(time_values);
    timeDiff = [timeDiff timeDiff(end)]; % Add last timeDiff element again
    timeDiffLessThanZero = timeDiff(valarray < 0);
    valLessThanZero = valarray(valarray < 0);
    partialRobNew = timeDiffLessThanZero.*valLessThanZero;
    partialRob = sum(partialRobNew);
    
    % Assert that the partialRob is negative - otherwise our additive
    % semantics are not sound with regards to the standard semantics
    assert(partialRob < 0);
else
    % The spec does not fail
    % Take the inverse of the integral of the inverse
    if rho == 0
        % Avoid division by zero: Just set the robustness to zero
        partialRob = 0;
    else
        % We are trying to speed up the process using vectorized calculations!
        timeDiff = diff(time_values);
        timeDiff = [timeDiff timeDiff(end)]; % Add last timeDiff element again
        reciprocalValues = 1./valarray;
        partialRobNew = sum(timeDiff.*reciprocalValues);
        partialRob = 1/partialRobNew;
    end
end
end