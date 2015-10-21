function [iX, iP, t, args] = GetArgSensiBar(DimX, ParamList, args, props)
%
% Init default parameters for sensitivity analysis.
% Might be obsolete
%
if ~exist('args','var')
    
    default_var = [];
    
    for ip = 1:DimX
        
        sip = size(ParamList{ip},2);
        sdv = size(default_var,2);
        
        if(sip > sdv)
            default_var = [default_var repmat(' ', size(default_var,1), sip-sdv)];
            param = ParamList{ip};
        elseif ( sdv > sip )
            param = [ParamList{ip} repmat(' ', 1 , sdv-sip)];
        else
            param = ParamList{ip};
        end
        
        default_var = [default_var; param];
        
    end
    
    default_param = [];
    
    for ip = 1:numel(ParamList)
        
        sip = size(ParamList{ip},2);
        sdv = size(default_param,2);
        
        if(sip > sdv)
            default_param = [default_param repmat(' ', size(default_param,1), sip-sdv)];
            param = ParamList{ip};
        elseif ( sdv > sip )
            param = [ParamList{ip} repmat(' ', 1 , sdv-sip)];
        else
            param = ParamList{ip};
        end
        
        default_param = [default_param; param];
        
    end
    defaults = {default_var, default_param, '' };
else
    defaults = args;
end

% Call input dialog box

prompt = {'Choose variables to get sensitivities from ', ...
    'Choose parameters', 'tspan for computing trajectories and sensitivities ([] uses tspan of pre-computed trajectories)'};
name = 'Choose variables and parameters to plot sensitivities';

numlines = [10; 10; 1];

answer = inputdlg(prompt,name,numlines,defaults);

args = answer;
vars = answer{1};
params = answer{2};

iX = {};
for ii = 1:size(vars,1)
    iX = {iX{:}, strtrim(vars(ii,:))};
end

iP = {};
for ii = 1:size(params,1)
    iP = {iP{:}, strtrim(params(ii,:))};
end

t = [];
try
    t = eval(answer{3});
catch
    return
end

end
