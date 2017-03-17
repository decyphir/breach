classdef signal_gen <handle
    % signal_gen An abstract simple class to generate single or multi-dimensional signals of some type.
    %
    % signal_gen Properties
    %    signals - cell array of names for the generated signals
    %    params  - cell array of parameter names
    %    p0      - default values for paramters, array with same dimension as params
    %
    % signal_gen Methods
    %   computeSignals (Abstract) - Takes parameter values and time and returns a signal
    %
    % See also constant_signal_gen, fixed_cp_signal_gen,
    % var_cp_signal_gen, pulse_signal_gen, step_signal_gen
    
    properties
        signals % names of the signals it generates
        params  % parameters such as control points, etc
        p0      % default values
        domain
    end
    
    
    methods
        function params = getParamNames(this)
            params = this.params;
        end
        
        function args= getSignalGenArgs(this)
            args = {};
        end
        
        function pval = get_param(this, pname)
           pidx = find(strcmp(this.params, pname));
           if ~isempty(pidx)
              pval = this.p0(pidx);  
           end
        end
        
        function set_param(this, pname, pval)
           pidx = find(strcmp(this.params, pname));
           if ~isempty(pidx)
              this.p0(pidx)=pval;  
           end
        end
        
    end
    
    methods (Abstract)
        computeSignals(this,p,time) % compute the signals
    end
    
end

