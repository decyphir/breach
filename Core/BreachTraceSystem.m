classdef BreachTraceSystem < BreachSystem
    % BreachTraceSystem  a BreachSystem class to handle traces with no
    % simulator
    
    methods
        % constructor - takes signal names and an optional trace
        function this = BreachTraceSystem(signals, trace)
            InitBreach;
            if (nargin==0)
                return;
            end
            
            if isscalar(signals) && isnumeric(signals)
                ndim =  signals;
                signal_names = cell(1,ndim);
                for is = 1:ndim
                    signal_names{is} = ['x' num2str(is)];
                end
            elseif ischar(signals)
                if exist(signals, 'file')
                    [~, ~, ext] = fileparts(signals);
                    switch (ext)
                        case '.csv'
                            % signals is a CSV files
                            fid = fopen(signals,'r');
                            if(fid==-1)
                                error(['Couldn''t open file ' signals]);
                            end
                            tline = strtrim(fgetl(fid));
                            signal_names = strsplit(tline,',');
                            signal_names= signal_names(2:end);
                            for i_sig = 1:numel(signal_names)
                                sig =strtrim(signal_names{i_sig});
                                signal_names{i_sig} = regexprep(sig,'\W','_');
                            end
                            
                            fclose(fid);
                            trace = csvread(signals,1);
                    end
                end
                
                % simout data
            elseif isa(signals,'Simulink.SimulationOutput')
                [time, X, signal_names] = simout2X(signals);
                trace = [time' X'];
                
                % default signals should be a cell array of strings
            else
                signal_names = signals;
            end
            
            % assumes now that we have signal names
            this.Sys = CreateExternSystem('TraceObject', signal_names, {'trace_id'},1);
            this.P = CreateParamSet(this.Sys);
            
            if exist('trace', 'var')
                this.AddTrace(trace);
            end
            
            %  Default domains
            for ip = this.P.DimP
                this.Domains(ip) = BreachDomain();
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
                if size(trace,1) >=1
                    traj.param = trace(1,2:end);
                else
                    traj.param = zeros(1, size(trace,2)+1);
                end
                
                
            elseif isnumeric(trace)
                traj.X = trace(:, 2:end)';
                traj.time = trace(:,1)';
                if size(trace,1) >=1
                    traj.param = trace(1,2:end);
                else
                    traj.param = zeros(1, size(trace,2)-1);
                end
            elseif isstruct(trace)
                traj = trace;
                traj.param=traj.param(1:end-1);
            end
            
            Pnew = CreateParamSet(this.Sys);
            Pnew.epsi(:,:) = 0;
            
            nb_traces =this.CountTraces();
            
            traj.param(end+1) = nb_traces+1;
            Pnew.traj={traj};
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
            this.Sys.tspan = traj.time;
        end
        
        function AddRandomTraces(this,n_traces, n_samples, amp, end_time)
            
            if ~exist('n_traces', 'var')
                n_traces= 1;
            end
            if ~exist('n_samples', 'var')
                n_samples = 100;
            end
            if ~exist('amp', 'var')
                amp = 4;
            end
            if ~exist('end_time', 'var')
                end_time = 100;
            end
            
            dimx = this.Sys.DimX;
            dimp = this.Sys.DimP;
            for it = 1:n_traces
                traj.time = linspace(0,end_time,n_samples);
                traj.X = amp*rand([dimx n_samples])-amp*rand();
                traj.param = zeros(1,dimp);
                this.AddTrace(traj);
            end
            
        end
        
        
        function st = disp(this)
            st = ['BreachTraceSystem with ' num2str(this.CountTraces()) ' traces.'];
            if nargout<1
                disp(st);
            end
        end
        
    end
    
end
