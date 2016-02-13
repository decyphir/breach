classdef BreachTraceSystem < BreachSystem
    
    methods
        % constructor - takes signal names and an optional trace
        function this = BreachTraceSystem(signals, trace)
            
            if (nargin==0)
                return;
            end
            
            if ischar(signals)          
                if exist(signals, 'file')
                    [~, ~, ext] = fileparts(signals);
                    switch (ext)
                        case '.csv'
                            % signals is a CSV files
                            fid = fopen(signals,'r');
                            if(fid==-1)
                                error(['Couldn''t open file ' signal_names]);
                            end
                            tline = strtrim(fgetl(fid));
                            signal_names = strsplit(tline,',');
                            signal_names= signal_names(2:end);
                            fclose(fid);
                            trace = csvread(signals,1);
                    end
                end
            % simout data
            elseif isa(signals,'Simulink.SimulationOutput')
                [time, X, signal_names] = simout2X(signals);
                trace = [time' X'];
            else
                signal_names = signals;
            end
            
            % assumes now that we have signal names
            this.Sys = CreateExternSystem('TraceObject', signal_names, {'no_trace'},1);
            if exist('trace', 'var')
                this.AddTrace(trace);
            end
            
        end
        
        % counts number of traces
        function  nb_traces= CountTraces(this)
            if isfield(this.P,'traj')
                nb_traces = numel(this.P.traj);
            else
                nb_traces=0;
            end
        end
        
        % Add a trace, either from file or from array
        % TODO checks dimensions of signals and data
        function AddTrace(this, trace)
            
            if ischar(trace)
                traj = load_traj(trace);
            elseif isa(trace,'Simulink.SimulationOutput')
                [time, X] = simout2X(signals);
                traj.X = X;
                traj.time = time;
                traj.param = trace(1,2:end);
            elseif isnumeric(trace)
                traj.X = trace(:, 2:end)';
                traj.time = trace(:,1)';
                traj.param = trace(1,2:end);
            end
            
            Pnew = CreateParamSet(this.Sys);
            Pnew.epsi(:,:) = 0;
            
            nb_traces =this.CountTraces();
            
            traj.param(end+1) = nb_traces+1;
            Pnew.Xf = traj.X(:,end);
            Pnew.traj=traj;
            Pnew.traj_ref = 1;
            Pnew.traj_to_compute =  [];
            Pnew.pts(1:Pnew.DimP,1) = traj.param';
            if nb_traces == 0
                this.P = Pnew;
            else
                this.P = SConcat(this.P, Pnew);
            end
            this.P.traj_ref = 1:nb_traces+1;
            this.P.traj_to_compute =  [];
            this.P.pts(this.P.DimX+1,:) = 1:nb_traces+1; % index traces
            this.UpdateParamRanges();
            
        end
        
        function disp(this)
            disp(['BreachTraceSystem with ' num2str(this.CountTraces()) ' traces.']);
        end

        
    end
    
end
