classdef from_workspace_signal_gen < signal_gen
% from_workspace_signal_gen meant to be used with a From Workspace block.
    properties
        data_fmt = 'timed_variables' 
        org_signals 
    end
    methods
        function this = from_workspace_signal_gen(signals)
            this.signals =signals;          
            this.params = {};
            for isig = 1:numel(signals)
                    this.org_signals.(signals{isig}) = evalin('base',[signals{isig} ';']);
            end
        end
        
        function [X, time]= computeSignals(this, p, time)

            X = zeros(numel(this.signals), numel(time));
            for isig = 1:numel(this.signals)
                % assumes that signals are in variable with [time  values]
                sig = this.org_signals.(this.signals{isig});
                t_sig = sig(:,1);
                v_sig = sig(:,2);
                x = interp1(t_sig, v_sig, time', 'linear', 'extrap'); 
                X(isig, :) = x';
            end
            
        end
    end
end