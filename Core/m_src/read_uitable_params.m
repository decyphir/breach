function [params, p0, domains] = read_uitable_params(h_uitable)
% read_uitable_params Takes a uitable element and reads parameter values and domains

data = get(h_uitable, 'Data');

for irow = 1:size(data,1)
    params{irow} = data{irow,1};
    if isnumeric(data{irow,2})
        p0(irow,1)  =data{irow,2};
    else
        p0(irow,1)  =0;
    end
    
    domain = str2num(data{irow,3});

    type = data{irow,5};
    enum  = data{irow, 4};
    
    domains(irow) = BreachDomain(type, domain);
    
%     if ~isempty(domain)
%         if isempty(enum)
%             if domain(2)-domain(1) >0
%                 domains(irow) = BreachDomain(type, [doml domu]);
%             else
%                 domains(irow) = BreachDomain(type);
%             end
%         else
%             enum = str2num(enum);
%             domains(irow) = BreachDomain(type, enum);
%             p0(irow,1) = domains(irow).checkin(p0(irow,1));
%         end
%     end
    
end
end