classdef BreachSignalGen < BreachSystem
    
    % Class to generate signals of different types 
    properties
        signalGenerators
        dt_default=1e-3 % in case no time step is provided
    end
    
    methods
        %% Constructor
        % SignalGen(type, name, varargin)
        function this = BreachSignalGen(signalGenerators)
            
            this.signalGenerators= signalGenerators;
            % we need to declare parameters, signals, p0, and simfn
            
            signals ={}; 
            params = {};
            p0=[];
            for isg = 1:numel(signalGenerators)
                signals = {signals{:}, signalGenerators{isg}.signals{:}};
                params = {params{:}, signalGenerators{isg}.params{:}}; 
                p0sg = signalGenerators{isg}.p0;
                if size(p0sg,2) >1
                    p0sg = p0sg';
                end
                p0 = [p0; p0sg ];
            end
            p0 = [zeros(numel(signals),1) ; p0 ];
            this.Sys = CreateExternSystem('BreachSignalGen', signals, params, p0, @(Sys, tspan, p)breachSimWrapper(this, Sys, tspan, p));
            this.Sys.tspan =0:.01:10;
            this.P = CreateParamSet(this.Sys);
            this.P.epsi(:,:)=0 ;
            
            if isaSys(this.Sys) % Note: we ignore initial conditions for now in ParamRanges
                                % OK for Simulink, less so for ODEs...
                this.ParamRanges = [this.Sys.p(this.Sys.DimX+1:end) this.Sys.p(this.Sys.DimX+1:end)];
                this.SignalRanges = [];
            end
        end
        
        
        function [tspan, X] = breachSimWrapper(this, Sys, tspan, p)
            
            if numel(tspan)==1
               tspan = 0:this.dt_default:tspan; 
            elseif numel(tspan)==2
               tspan = tspan(1):this.dt_default:tspan(2); 
            end
            
            p = p(this.Sys.DimX+1:end);
            cur_ip =1;
            cur_is =1;
            for isg = 1:numel(this.signalGenerators)
               np = numel(this.signalGenerators{isg}.params);
               p_isg = p(cur_ip:cur_ip+np-1);
               ns = numel(this.signalGenerators{isg}.signals);
               X(cur_is:cur_is+ns-1, :) = this.signalGenerators{isg}.computeSignals(p_isg, tspan); 
               cur_ip = cur_ip+ np;
               cur_is = cur_is+ ns;
            end
                
        end
        
         
    end
    
end
