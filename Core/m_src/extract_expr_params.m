function  params = extract_expr_params(expr)
%extract_params_expr extracts names of  parameters involved in expr
%
% Synopsis: params = extract_expr_params(expr)
%
% Input:
%  - expr : string
%
% Output:
%  - params: list of parameter names


params = {};
[~,~, ~, matches, tokens] = regexp(expr, '(\<\w+\>)');
for im=1:numel(matches)
    varname = tokens{im}{1};
    if isvarname(varname)&& ...
            all(exist(varname) ~= [2 3 5 6])  % checks for m-files, mex files, builtin, p-files
        params{end+1} = tokens{im}{1};
    end
end

reserved = { 'alw', 'ev', 'and', 'or', 'not', 'until', 't', ...
    'abs', 'sin', 'cos', 'exp','tan', 'norm','sqrt'};
params = setdiff(params, reserved);
params = unique(params);

end
