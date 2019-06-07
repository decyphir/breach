function [params, p0, domains] = read_uitable_params(h_uitable)
% read_uitable_params Takes a uitable element and reads parameter values and domains

data = get(h_uitable, 'Data');

for irow = 1:size(data,1)
    
    % param name
    params{irow} = data{irow,1};
    
    
    if isnumeric(data{irow,2})
        p0(irow,1)  =data{irow,2};
    else
        p0(irow,1)  =0;
    end
    
    if ~isempty(data{irow,3}) 
        domain = str2num(data{irow,3});
    else
        domain=[];
    end
    
    type = data{irow,5};
    enum  = data{irow, 4};
        
    switch type
        case 'enum'
            if isempty(enum)
                enum = data{irow,2};
            else
                enum = str2num(enum);
            end
            domains(irow) = BreachDomain(type, domain, enum);
        otherwise 
            domains(irow) = BreachDomain(type, domain);
    end
    
     p0(irow,1) = domains(irow).checkin( p0(irow,1));
     data{irow,2} = p0(irow,1);
    
    %
end

set(h_uitable, 'Data',data);

end