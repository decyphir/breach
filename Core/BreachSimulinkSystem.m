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
        end
 
        function disp(BrObj)
            disp(['BreachObject interfacing model ' BrObj.Sys.mdl '.']);
            disp(' ');
            disp( 'Signals:')
            disp( '-------')
            cellfun(@disp, BrObj.Sys.ParamList(1:BrObj.Sys.DimX));
            disp(' ')
            disp('Inputs:')
            disp('------')
            cellfun(@disp, BrObj.Sys.InputList);
            disp(' ')
            disp('Parameters:')
            disp('----------')
            
            plist = BrObj.Sys.ParamList(BrObj.Sys.DimX+1:end);
            for ip = 1:numel(plist)
                fprintf('%s: %g\n', plist{ip}, GetParam(BrObj.P,plist{ip}));
            end
        end
                
    end
    
    
end
