classdef from_file_signal_gen < signal_gen
    properties
        file_list                                    % filename pattern, e.g., 'data_firstparam%d_secondparam%d.mat'
        data_fmt = 'timed_variables' % other could be 'timed_array', 'struct_with_time', etc
        domain  
    end
    methods
        function this = from_file_signal_gen(fname, signals)
            this.signals = signals;
            this.file_list = dir(fname);        
            this.params = {'file_idx'};
            this.p0 = 1;
            this.domain = BreachDomain('int', [1 numel(this.file_list)]);
        end
        
        function X= computeSignals(this, p, time)
            X = zeros(numel(this.signals), numel(time));
            fname = this.file_list(p).name;
            try 
                st = load(fname);
            catch
                warning(['Could not read file ' fname '. Returning NaN trace.'] )
                X(:,:) = NaN;
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