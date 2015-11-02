classdef BreachSimulinkSystem < BreachObject
    %BreachSimulinkSystem specializes BreachObject to Simulink systems 
       
    methods
        
        % Constructor
        function BrObj = BreachSimulinkSystem(varargin)
            BrObj = BrObj@BreachObject();

            try % try first to import a simulink model
                BrObj.Sys = CreateSimulinkSystem(varargin{:});
            catch
                switch nargin
                    case 1
                        inSys = varargin{1};
                        simulink_error = lasterr;
                        try
                            BrObj.Sys = inSys;
                            BrObj.P = CreateParamSet(BrObj.Sys);
                        catch
                            error(simulink_error);
                        end
                end
            end
            BrObj.P = CreateParamSet(BrObj.Sys);
            BrObj.Sys.ParamRanges = [BrObj.Sys.p BrObj.Sys.p];
            BrObj.Sys.SignalRanges = [];
        end
                 
    end
    
    
end
