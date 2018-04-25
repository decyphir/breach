function highlight_truth_intervals(tau, val, true_color, true_alpha, false_color, false_alpha)

if ~exist('true_color', 'var')||isempty(true_color)
    true_color = 'g';
end
if ~exist('true_alpha', 'var')||isempty(true_alpha)
    true_alpha = 0.3;
end
if ~exist('false_color', 'var')||isempty(false_color)
    false_color = 'r';
end
if ~exist('false_alpha', 'var')||isempty(false_alpha)
    false_alpha = 0.3;
end

[itrue, ifalse] = get_truth_intervals(val);

if true_alpha ~=0
    for ii = 1:size(itrue,1)
        highlight_interval(gca, [tau(itrue(ii,1)), tau(itrue(ii,2))],true_color, true_alpha);
    end
end
if false_alpha ~= 0
    for ii = 1:size(ifalse,1)
        highlight_interval(gca, [tau(ifalse(ii,1)), tau(ifalse(ii,2))],false_color,false_alpha);
    end
end

end