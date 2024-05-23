classdef equal_param_gen < param_gen
    % equal_param_gen enforces a set of param_out to be always equal to one param_in

    properties
        n_out
    end

    methods
        function this = equal_param_gen(param_in, param_out, p0)
            this.params = {param_in};
            
            
            if ~iscell(param_out)
                param_out={param_out};
            end
            this.params_out = param_out;
            this.n_out= numel(this.params_out);

            if nargin<3
                p0=0;
            end
            this.p0 =p0;
        end

        function p_out = computeParams(this, p_in)
            p_out = repmat(p_in,this.n_out,1);            
        end
    end

end