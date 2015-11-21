classdef BreachSimulinkSystem < BreachOpenSystem
    %BreachSimulinkSystem specializes BreachObject to Simulink systems
    
    methods
        
        % See CreateSimulinkSystem for arguments help
        function this = BreachSimulinkSystem(varargin)
            
            switch nargin
                case 0 % do nothing
                case 1 % Should be a Sys structure - check that type is right
                    
                    inSys = varargin{1};
                    if isaSys(inSys) && strcmp(Sys.type, 'Simulink')
                        this.Sys = inSys;
                        this.ResetParamSet();
                    elseif exist(inSys)==4  %  create Simulink system with default options
                        this.Sys = CreateSimulinkSystem(inSys);
                    else
                        error('BreachSimulinkSystem with one argument assumes that the argument is a system structure of type Simulink or the name of a Simulink model.')
                    end
                otherwise % creates a Simulink system
                    this.Sys = CreateSimulinkSystem(varargin{:});
            end
            
            if isaSys(this.Sys)
                this.P = CreateParamSet(this.Sys);
                this.ParamRanges = [this.Sys.p(this.Sys.DimX+1:end) this.Sys.p(this.Sys.DimX+1:end)];
                this.SignalRanges = [];
            end
            
        end
        
    end
end
