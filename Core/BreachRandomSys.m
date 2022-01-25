classdef BreachRandomSys < BreachSignalGen
%BreachRandomSys Creates a random system from control point signal
%generator

    methods
        function this = BreachRandomSys(sigs, methods, ntraces)
            if isnumeric(sigs)
                nsigs = sigs;
                sigs = {};
                for isig=1:nsigs
                   sigs{isig}= ['x' num2str(isig)];
                end         
            elseif ischar(sigs)
                sigs=  {sigs};      
            end            
            nsigs=numel(sigs);
            
            if ~exist('methods','var')||isempty(methods)
               methods= 'spline';
            end
            
            sg = random_signal_gen(sigs, methods);
            this.InitSignalGen({sg});
            
            if exist('ntraces','var')
                %TODO                
                %this.SetParam();                
            end
            
        end
    
    
        function GenSim(this, n, seed, time)
            % Computes n signals with seed
            if nargin<2
                n=1;
            end
            if nargin<3
                seed = 0;
            end
            if nargin<4
                time = 0:.01:5;
            end
            
            seed_params = this.expand_param_name('seed');            
            num_seed_val = numel(seed_params)*n; 
            values = reshape(seed+1:seed+num_seed_val, numel(seed_params), n);
            this.SetParam(seed_params, values);
            this.Sim(time);            
        end
    
    
    end
end