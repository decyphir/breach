classdef constant_signal_gen < signal_gen
    % constant_signal_gen a class to generate simple multidimensional constant signals
    %
    % constant_signal_gen Methods
    %    constant_signal_gen - constructor takes a cell array of signal
    %    names and optional values. Parameters have same names as signals with '_u0' suffix.
    %
    % See also signal_gen.

    methods         
        function this = constant_signal_gen(signals, p0)
            
            if ischar(signals)
                signals = {signals};
            end
            
            nb_signals = numel(signals);
            this.signals = signals;
            
            this.params = {};
            for i_s = 1:nb_signals
                this.params = { this.params{:} [this.signals{i_s} '_u0']};
            end
            
            this.params_domain = repmat(BreachDomain(), 1, nb_signals); 
            this.signals_domain = repmat(BreachDomain(), 1, nb_signals); 
            
            if nargin == 2
                this.p0 = p0;
            else
                this.p0 = zeros(nb_signals,1);
            end
            
        end
            
        function [X, time] = computeSignals(this,p, time) % compute the signals
            if numel(p) ~= numel(this.signals)
                error('Wrong number of parameters for computing constant signal.' )
            end
            if size(p,1) ==1
                p = p';
            end
            X = repmat(p, 1, numel(time)); 
        end
        
        function type = getType(this)
            type = 'constant';
        end
        
    end
            
end


