function res = updateSensiRes(res, res_new)

% checks that res_new has the same parameters and dimensions
if numel(res_new)~=numel(res)
    error('The two arguments must have same dimensions.');
end

for idx_res = 1:numel(res)
        
    if ~isequal(res{idx_res}.params,res_new{idx_res}.params)
        error('parameters are different in the new results.');                
    end
    res{idx_res}.rob = [res{idx_res}.rob res_new{idx_res}.rob]; % concat robustness results
    res{idx_res}.EE = [res{idx_res}.EE res_new{idx_res}.EE];   % and elementary effects
    
    num_params = numel(res{idx_res}.params);
    num_paths = size(res{idx_res}.EE,2)/num_params; % new number of paths then updates sensitivities    
    
    res{idx_res}.mu = 1/num_paths*sum(res{idx_res}.EE, 2);
    res{idx_res}.mustar = 1/num_paths* sum(abs(res{idx_res}.EE), 2);
    
    MU = repmat(res{idx_res}.mu,[1 size(res{idx_res}.EE,2)]);
    MUSTAR = repmat(res{idx_res}.mustar,[1 size(res{idx_res}.EE,2)]);
    res{idx_res}.sigma = sqrt(1/(num_paths-1)* sum((res{idx_res}.EE-MU).^2,2));
    res{idx_res}.sigmastar = sqrt(1/(num_paths-1)* sum((abs(res{idx_res}.EE)-MUSTAR).^2,2));
            
end