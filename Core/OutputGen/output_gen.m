classdef output_gen < handle

    properties % all are cells of strings, except for domains
        in_signals                   % input signals needed to compute outputs 
        out_signals                 % output signals 
        in_params                   % parameters of the output gen 
        p0                               % default values for parameters
        out_values                  % names says it 
        domains = containers.Map()  % maps all of the above to their respective domains 
   end

    methods (Abstract)
        [val, tout, sig_out] = eval(this, time, X, p) % returns output values, (output) signals, output times
    end

    methods
        function new = copy(this)
            % copy operator, works with R2010b or newer.
            objByteArray = getByteStreamFromArray(this);
            new = getArrayFromByteStream(objByteArray);
        end
    end

end