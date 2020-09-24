function display_morris_result(res, varargin)

opt.SortBy = 'mu';
opt = varargin2struct_breach(opt,varargin{:});

switch opt.SortBy
    case 'mu' 
        [~ , order] = sort(res.mu);
    case 'mustar' 
        [~ , order] = sort(res.mustar);
    case 'sigma' 
        [~ , order] = sort(res.sigma);
    case 'sigmastar' 
        [~ , order] = sort(res.sigmastar);        
end

fprintf('mu        mustar    sigma     sigmastar\n') 

for ip = 1:numel(res.params)   
    idx = order(ip);    
    %format = '%+5.5e';
    format = '%+5.5G';
    
    st_format = [format '   ' format '   ' format '   ' format '   ---> %s\n' ]; 
    fprintf(st_format, ...
       res.mu(idx), res.mustar(idx), res.sigma(idx), res.sigmastar(idx), res.params{idx}) 
end



end

