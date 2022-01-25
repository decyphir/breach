classdef from_file_pattern_signal_gen < signal_gen
    properties
        pattern                                   % filename pattern, e.g., 'data_firstparam%d_secondparam%d.mat'
        data_fmt = 'timed_variables' % other could be 'timed_array', 'struct_with_time', etc
    end
    methods
        function this = from_file_pattern_signal_gen(signals,fname)
            this.signals =signals;
            fname = regexprep(fname, '\\','/');
            [this.params,this.p0, this.pattern] = guess_filename_params(fname);                         
        end
        
        function [X, time] = computeSignals(this, p, time)
            fname = sprintf(this.pattern,p);
            X = zeros(numel(this.signals), numel(time));
            try 
                st = load(fname);
            catch
                warning(['Could not read file ' fname '. Returning NaN trace.'] )
                X(:,:) = NaN;
                return;
            end
            for isig = 1:numel(this.signals)
                % assumes that signals are in variable with [time  values]
                sig = st.(this.signals{isig});
                t_sig =  sig(:,1);
                v_sig = sig(:,2);
                x = interp1(t_sig, v_sig, time', 'linear', 'extrap'); 
                X(isig, :) = x';
            end
            
        end
    end
end