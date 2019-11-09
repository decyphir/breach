classdef param_gen < handle
    properties
        params
        domain
        domain_out
        params_out
        p0  % default values for params
    end
    methods (Abstract)
       p_out = computeParams(this, p_in)
    end
    methods
       function assign_params(this, p)
            % assign_params fetch parameters and assign them in the current context
            for ip = 1:numel(this.params)
                assignin('caller', this.params{ip},p(ip));
            end
        end
    end
      
end