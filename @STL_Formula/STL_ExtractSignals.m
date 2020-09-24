function [signals, params,p0] = STL_ExtractSignals(phi)
%STL_ExtractSignals extracts names of signals and parameters involved in phi
%
% Synopsis: signals = STL_ExtractSignals(phi)
%
% Input:
%  - phi : an STL formula
%
% Output:
%  - signals: list of signals involved in phi (note: currently only detec-
%    ted from patterns of the form signal_id[t])
%  - params: list of parameters in the formula
%. - p0: default values for parameters

st_phi = disp(phi,0);
signals = {};
[~,~, ~, matches, tokens] = regexp(st_phi, '(\<\w+\>)[.+?\]');
for im=1:numel(matches)
    signals{end+1} = tokens{im}{1};
end
sreserved = {'alw_', 'ev_','until_', 'once_', 'hist_'};
signals = setdiff(signals, sreserved);


if nargout>1
    %% params
    params = {};
    [~,~, ~, matches, tokens] = regexp(st_phi, '(\<\w+\>)');
    for im=1:numel(matches)
        varname = tokens{im}{1};
        if isvarname(varname)
            if all(exist(varname) ~= [2 3 5 6])  % checks for m-files, mex files, builtin, p-files
                params{end+1} = tokens{im}{1};
            end
        end
    end
    
    reserved = [ sreserved signals  {'alw', 'ev','once', 'hist','and', 'or', '=>', 'not', 'until', 't', ...
        'abs', 'sin', 'cos', 'exp','tan', 'norm','sqrt'}];
    params = setdiff(params, reserved);
    params = unique(params);
    
    %% p0
    if nargout>=3
        p0 = zeros(1,numel(params));
        param_struct= get_params(phi);
        for ip = 1:numel(params)
            if isfield(param_struct, params{ip});
                p0(ip) = param_struct.(params{ip});
            else
                p0(ip) =0;
            end
        end
    end
end
end