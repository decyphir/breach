function [true_intrvls, false_intrvls] = get_truth_intervals(val)
% get_truth_intervals assumes piecewise constant

true_vals = (val > 0|isnan(val));  % NaN true by default
true_intrvls =[];
false_intrvls =[];

if all(true_vals)
  true_intrvls = [1 numel(val)];
elseif any(true_vals)
    idx_start = find(diff([0.5 true_vals]));
    idx_intrvls =  [idx_start(1:end-1)' idx_start(2:end)'];
    if true_vals(1) == true  % first interval is true
        true_intrvls = idx_intrvls(1:2:end,:);
        false_intrvls = idx_intrvls(2:2:end,:);
    else
        true_intrvls = idx_intrvls(2:2:end,:);
        false_intrvls = idx_intrvls(1:2:end,:);
    end
    
    if true_vals(end)>0&&~all(true_vals>0)
        true_intrvls(end+1,:) = [false_intrvls(end,2) numel(val)];
    else
        false_intrvls(end+1,:) = [true_intrvls(end,2) numel(val)];
    end
else
    false_intrvls = [1 numel(val)];
end

% remove NaNs intervals from false
remove_nan = [];
for ifalse = 1:size(false_intrvls,1)
    if all(isnan(val(false_intrvls(ifalse,1):false_intrvls(ifalse,2)-1)))
        remove_nan(end+1) = ifalse;
    end
end

if ~isempty(remove_nan)
    false_intrvls(remove_nan,:) = [];
end

% remove NaNs intervals from true
remove_nan = [];
for itrue = 1:size(true_intrvls,1)
    if all(isnan(val(true_intrvls(itrue,1):true_intrvls(itrue,2)-1)))
        remove_nan(end+1) = itrue;
    end
end

if ~isempty(remove_nan)
    true_intrvls(remove_nan,:) = [];
end


end
