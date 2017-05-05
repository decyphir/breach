classdef BreachFileSystem < BreachSystem
    % BreachFileSystem specializes BreachObject to monitoring traces
    % defined in files.
    % TODO Documentation 
    
    
    properties
        nb_traces
    end
    
    methods
        
        % Constructor
        function BrObj = BreachFileSystem(name, vars)
            
            BrObj = BrObj@BreachObject();
            BrObj.nb_traces= 0;
            params = {'no_file'};
            BrObj.Sys = CreateExternSystem(name, vars, params, 1);
            BrObj.Sys.tspan = 0:.1:1; % Default - should never be used
            BrObj.P = CreateParamSet(BrObj.Sys);
            
        end
        
        % TODO some error handling
        function ReadTrace(BrObj, trace)
            % if input is a string, assumes it's a file
            
            Pnew = CreateParamSet(BrObj.Sys);
            Pnew = SetParam(Pnew, 'no_file', BrObj.nb_traces+1);
            
            if ischar(trace)
                traj = load_traj(trace);
                traj.param(end+1) = BrObj.nb_traces+1;
                Pnew.traj=traj;
                Pnew.traj_ref = 1;
                Pnew.traj_to_compute =  []; 
                if BrObj.nb_traces == 0
                    BrObj.P = Pnew;
                else
                    BrObj.P = SConcat(BrObj.P, Pnew);
                end
                BrObj.nb_traces = BrObj.nb_traces+1;
                BrObj.P.traj_ref = 1:BrObj.nb_traces;
                BrObj.P.traj_to_compute =  []; 
            end
            
        end
        
        % TODO for CPSGrader style test plans
        function ReadSTLFile(BrObj,varargin)
            % TODO
        end
        
        
        % Plot signals
        function PlotSignals(BrObj, varargin)
            SplotVar(BrObj.P, varargin{:});
        end
        
        function disp(BrObj)
            disp(['BreachFileSystem ' BrObj.Sys.name '.']);
        end
        
    end
    
    
end
