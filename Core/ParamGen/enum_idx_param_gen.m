classdef enum_idx_param_gen < param_gen
    % equal_param_gen enforces one param_out to be always equal to param_in
    
    properties
        domain_out
    end
    
    methods
        function this = enum_idx_param_gen(param, domain)
            if isnumeric(domain)
               domain = BreachDomain('enum', domain) ;
            end
            this.params = {[param '_enum_idx']};
            this.domain = {BreachDomain('int', [1 length(domain.enum)])};
            this.domain_out = {domain};
            this.domain_out{1}.domain = []; 
            this.params_out = {param};
            this.p0 = 1;
        end
        
        function p_out = computeParams(this, p_in)
            p_out = this.domain_out{1}.enum(p_in); 
        end
    end
     
end