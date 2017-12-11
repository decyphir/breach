function [params_u, idx_u] = FindParamsInput(Sys, InputName)
% FINDPARAMSINPUT returns the name and indices of control point parameters

params_u = {};    
idx_u = [];     

if isfield(Sys, 'InputList')
    InputNames = Sys.InputList;
else
    return;
end

if isfield(Sys, 'InputOpt')&&~isempty(Sys.InputOpt)
    
    InputOpt  = Sys.InputOpt;
     
    switch InputOpt.type
        
        case 'UniStep'
            for ku = 1:numel(InputNames)
                for k = 1:InputOpt.cp(ku)
                    params_u = {params_u{:} [InputNames{ku} '_u' num2str(k-1)]};
                end
            end
            
        case 'VarStep'
            for ku = 1:numel(InputNames)
                for k = 1:InputOpt.cp(ku)
                    params_u = {params_u{:} [InputNames{ku} '_dt' num2str(k-1)] [InputNames{ku} '_u' num2str(k-1)]};
                end
            end          
    end
end

if nargin==2
    res = regexp(params_u, ['^' InputName]);
    res = ~cellfun(@isempty,res);
    params_u = params_u(res);
end

if nargout ==2 
    idx_u = FindParam(Sys, params_u);
end

end