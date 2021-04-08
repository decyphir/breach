function [time_values, valarray] = RobustAlways(time_values, valarray, I___)

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
    
    partialRob = PartialRobustAlways(partialTimeNew, partialValarrayNew);
    valarray(valIndex) = partialRob;
    
    if valIndex > 2
        if time_values(valIndex) == time_values(valIndex - 1) ...
                && valarray(valIndex) == valarray(valIndex - 1)
            indicesToRemove = [indicesToRemove valIndex]; %#ok<*AGROW>
        end
    end
    
end

% Remove subsequent indices that are identical in both time and value
valarray(indicesToRemove) = [];
time_values(indicesToRemove) = [];

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
    % Note: We only need to assert this if the spec is actually positive
    % for time greater than zero (otherwise it is instantaneously positive
    % and we do not care)
    if partialRob > 0
        assert(sum(timeDiffLessThanZero) == 0);
    end
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