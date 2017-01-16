classdef my_signal_generator < signal_gen
    properties
       lambda  % some parameter for the signal generator - this won't be visible from Breach API
    end
    methods
        % The constructor must name the signals and the parameters needed to construct them.
        function this = my_signal_generator(lambda)
            this.lambda = lambda;
            this.signals = {'Engine_Speed', 'Pedal_Angle'};
            this.params = {'my_param_for_Engine', 'my_param_for_Pedal'};
            this.p0 = [1000 50]; % default values     
        end
        
        % The class must implement a method with the signature below
        function [X, time] = computeSignals(this, p, time)
            % p contains values for the parameters in the declared order
            my_param_for_Engine = p(1);
            my_param_for_Pedal = p(2); 
                         
            % Constructs signals as some function of time and parameters
            Engine_Speed = this.lambda*cos(time)+ my_param_for_Engine; 
            Pedal_Angle = this.lambda*sin(time)+ my_param_for_Pedal; 
            
            % The signals must be returned as rows of X, in the declared order 
            X = [ Engine_Speed; Pedal_Angle ];
        end    
    end
end