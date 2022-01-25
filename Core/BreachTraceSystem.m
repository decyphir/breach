classdef BreachTraceSystem < BreachSystem
    % BreachTraceSystem  a BreachSystem class to handle traces with no
    % simulator
    
    methods
        % constructor - takes signal names and an optional trace
        function this = BreachTraceSystem(signals, trace, params)
            InitBreach;
            if (nargin==0)
                return;
            end
            
            if ischar(signals)
                signal_names = {signals};
            elseif isscalar(signals) && isnumeric(signals)
                ndim =  signals;
                signal_names = cell(1,ndim);
                for is = 1:ndim
                    signal_names{is} = ['x' num2str(is)];
                end
                
                
            elseif isstruct(signals)&&all(isfield(signals, {'time', 'outputs', 'inputs'}))
                trace1 = signals;
                signal_names = [trace1.outputs.names, trace1.inputs.names];
            elseif isstruct(signals)
                if all(isfield(signals, {'signals', 'time'})) % traj with signal names
                    trace1 = signals;
                    signal_names = trace1.signals.names;
                elseif all(isfield(signals,  {'time', 'X'})) % traj, but no signal names
                    trace1 = signals;
                    ndim = size(trace1.X,1);
                    signal_names = cell(1,ndim);
                    for is = 1:ndim
                        signal_names{is} = ['x' num2str(is)];
                    end
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
            if exist('params', 'var')&&~isempty(params)
                this.Sys = CreateSystem(signal_names, params,zeros(numel(signal_names)+numel(params),1));
            else
                this.Sys = CreateSystem(signal_names, {},zeros(numel(signal_names),1));
            end
            
            this.P = CreateParamSet(this.Sys);
            
            if exist('trace1', 'var')&&~isempty(trace1)
                this.AddTrace(trace1);
            end
            if exist('trace', 'var')&&~isempty(trace)
                this.AddTrace(trace);
            end
            
            %  Default domains
            for ip = 1:this.P.DimP
                this.Domains(ip) = BreachDomain();
            end
            
        end
        
        function nb_traces= CountTraces(this)
            % CountTraces counts number of traces
            
            if isfield(this.P,'traj')
                nb_traces = numel(this.P.traj);
            else
                nb_traces=0;
            end
        end
        
        function AddTrace(this, trace, params)
            % Add a trace, either from file or from array
            % TODO checks dimensions of signals and data
            
            signals = this.GetSignalNames();
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
            elseif isstruct(trace)
                if all(isfield(trace,{'time', 'X'}))
                    traj= trace;
                elseif all(isfield(trace,{'time', 'outputs', 'inputs'})) %  reading one struct obtained from a SaveResult
                    traj.time= trace.time;
                    traj.X = zeros(numel(signals),numel(traj.time));
                    for isig = 1:numel(signals)
                        idx_sig = find(strcmp(trace.inputs.names, signals{isig}),1);
                        if isempty(idx_sig)
                            idx_sig = find(strcmp(trace.outputs.names, signals{isig}),1);
                            if isempty(idx_sig)
                                error('BreachTraceSystem:signal_not_found', 'Signal %s not found', signals{isig});
                            else
                                traj.X(isig,:) = trace.outputs.values(idx_sig,:);
                            end
                        else
                            traj.X(isig,:) = trace.inputs.values(idx_sig,:);
                        end
                    end
                    
                elseif all(isfield(trace,{'time', 'signals'})) %  reading one struct obtained from a SaveResult
                    traj.time= trace.time;
                    signals = this.GetSignalNames();
                    traj.X = zeros(numel(signals),numel(traj.time));
                    for isig = 1:numel(signals)
                        idx_sig = find(strcmp(trace.signals.names, signals{isig}),1);
                        if isempty(idx_sig)
                            error('BreachTraceSystem:signal_not_found', 'Signal %s not found', signals{isig});
                        end
                        traj.X(isig,:) = trace.signals.values(idx_sig,:);
                    end
                end
            elseif isnumeric(trace)
                try
                    traj.X = trace(:, 2:end)';
                    traj.time = trace(:,1)';
                    traj.param = this.Sys.p';
                    if size(trace,1) >=1
                        traj.param(1:this.Sys.DimX) = trace(1,2:end);
                    else
                        traj.param(1:this.Sys.DimX) = zeros(1, size(trace,2)+1);
                    end
                catch
                    error('BreachTraceSystem:wrong_trace', 'Problem with trace format. Should be of size (N,%d)', this.P.DimX+1);
                end
                
            elseif isa(trace, 'matlab.io.MatFile')
                traj = trace;
            end
            
            nb_traces =this.CountTraces();
            
            if ~isfield(traj, 'param')&&~isa(traj, 'matlab.io.MatFile')
                traj.param = this.Sys.p';
            end
            if exist('params', 'var')
                num_p =this.P.DimP-this.P.DimX;
                if numel(params)~= num_p
                    error('BreachTraceSystem:wrong_param', 'Wrong number of parameters (last argument), should be %d', this.P.DimP-this.P.DimX);
                else
                    traj.param(this.DimX+1:end) = reshape(params, 1, num_p);
                end
                
            end
            new_pts  = traj.param';
            
            if nb_traces == 0
                Pnew = CreateParamSet(this.Sys);
                Pnew.epsi(:,:) = 0;
                Pnew.traj={traj};
                Pnew.traj_ref = 1;
                Pnew.traj_to_compute =  [];
                Pnew.pts(1:Pnew.DimP,1) =  new_pts;
                this.P = Pnew;
                this.Sys.tspan = traj.time;
            else
                % add pts and epsi
                this.P.pts(:,end+1) = new_pts;
                this.P.epsi(:,end+1) = this.P.epsi(:,end);
                this.P.traj{end+1} = traj;
                this.P.traj_ref(end+1) = numel(this.P.traj);
            end
            
            % FIXME still need to get rid of legacy useless Xf field...
            if isa(traj, 'matlab.io.MatFile')  % avoid loading matfile for Xf
                Xf= zeros(1,this.P.DimX);
            else
                Xf = traj.X(:,end);
            end
            if isfield(this.P, 'Xf')
                this.P.Xf(:,end+1)= Xf;
            else
                this.P.Xf= Xf;
            end
        end
        
        function GenTraceIdx(this, name)
            % GenTraceIdx generate unique index for traces
            
            if ~exist('name','var')||isempty(name)
                name = 'trace_idx';
            end
            
            [~, name_exists] = FindParam(this.Sys, name);
            num_traces = this.CountTraces();
            
            if num_traces>=0&&~name_exists
                for it = 1:num_traces
                    pts = this.P.pts(:,it);  % WARNING: assumes no traj_ref trick...
                    traj = this.P.traj{it};
                    traj.param= [traj.param it];
                    if it == 1 % create new Sys (add trace_idx param)
                        this.Sys.ParamList = [this.Sys.ParamList name];
                        this.Sys.DimP = this.Sys.DimP+1;
                        this.Sys.p = [this.Sys.p; 1];
                        Pnew = CreateParamSet(this.Sys);
                        Pnew.epsi = [this.P.epsi(:,1) ; 0 ];
                        Pnew.traj={traj};
                        Pnew.traj_ref = 1;
                        Pnew.traj_to_compute =  [];
                        Pnew.pts(1:Pnew.DimP,1) =  [pts; 1];
                    else
                        % add pts and epsi
                        Pnew.pts(:,end+1) = [pts; it];
                        Pnew.epsi(:,end+1) = [this.P.epsi(:,end); 0];
                        Pnew.traj{end+1} = traj;
                        Pnew.traj_ref(end+1) = numel(this.P.traj);
                    end
                end
                this.P = Pnew;
                this.Domains(end+1)= BreachDomain('int');
            end
        end
        
        function AddRandomTraces(this,n_traces, n_samples, amp, end_time)
            % AddRandomTraces Initially used to test monitoring algo
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
        
        function Sim(varargin)
            % BreachTraceSystem.Sim(varargin) does nothing - traces are added
            % using AddTrace method
            
        end
        
        function varargout = disp(this)
            st = ['BreachTraceSystem with ' num2str(this.CountTraces()) ' traces.\n'];
            if nargout == 0
                varargout = {};
                fprintf(st);
            else
                varargout{1} = st;
            end
            
        end
        
    end
    
end
