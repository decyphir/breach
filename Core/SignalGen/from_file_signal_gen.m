classdef from_file_signal_gen < signal_gen
    properties
        file_list
        file_name
        var_name
        time_var_name
        ignore_time = false   % if true, will ignore time given as parameter of Sim command and use time data instead
        data_fmt = 'timed_variables' % other could be 'timed_array', 'struct_with_time', etc
        pts  % stores possible values of parameters
    end
    methods
        function this = from_file_signal_gen(signals, fname, varname,params)
            
            if ~exist('fname', 'var')
                fname=  '*.mat';
            end
            if ~iscell(fname)
                fname= {fname};
            end
            
            % First thing is detecting files
            for ifn=1:numel(fname)
                dir_file_list = dir(fname{ifn});
                pathstr = fileparts(fname{ifn}); % dir does not keep the path...
                for ifnl = 1:numel(dir_file_list)
                    if ~isempty(pathstr)
                        dir_file_list(ifnl).name = [pathstr filesep dir_file_list(ifnl).name];
                    end
                end
                if isempty(this.file_list)
                    this.file_list = dir_file_list;
                else
                    this.file_list = [this.file_list dir_file_list];
                end
            end
            if isempty(this.file_list)
                error('from_file_signal_gen:no_trace_file', 'No trace file.')
            end
            
            this.params = {'file_idx'};
            this.p0 = 1;
            this.params_domain = BreachDomain('enum', [],  [1:numel(this.file_list)]);
            this.pts = [1:numel(this.file_list)];
            
            % open first file
            fname = this.file_list(1).name;
            st = load(fname);
            vars = fieldnames(st);
            
            % look for a time array - let's make that mandatory to begin
            % with.
            
            itime =find(strcmpi(vars, 'time'), 1);
            if isempty(itime)
                error('The file must contain a 1D array named ''time'' (case insensitive) with properly ordered time samples.');
            else
                time = st.(vars{itime});
                this.time_var_name = vars{itime};
                if any(diff(time)<0)
                    error('The file must contain a 1D array named ''time'' (case insensitive) with properly ordered time samples.');
                end
            end
            
            % Next is detecting all parameters- we're looking for
            % constant scalar parameters defined in all files - and all
            % signals - they have to be defined  in all files.
            if exist('params','var')
                if ischar(params)
                    params = {params};
                end
            else
                params= {}; 
            end
            
            
            signals_all = {};
            for iv = 1:numel(vars)
                v = vars{iv};
                if isempty(params)||(ismember(v,params)) % if params is specified, make sure v is in
                    if (iv ~= itime)&&isnumeric(st.(v))   % ignore time and everything not numeric
                        if isscalar(st.(v)) % this is a pararmeter!
                            this.params = [this.params v];
                            this.p0(end+1) = st.(v);
                        elseif length(st.(v))==length(time) % looks like  a signal
                            signals_all = [signals_all  v];
                        end
                    end
                end
            end
            % go over all files to fetch values for the parameters and
            % create enum domains
            for ifile = 1:numel(this.file_list)
                fname = this.file_list(ifile).name;
                st = load(fname);
                vars = fieldnames(st);
                for ip = 2:numel(this.params)
                    this.pts(ip, ifile) = st.(this.params{ip});   % TODO some try catch for parameter defined in first file, but not in other(s)
                end
            end
            for ip = 2:numel(this.params)
                this.params_domain(ip)= BreachDomain('enum', [], this.pts(ip, :));
            end
            
            % Match signals with signals_all
            if ischar(signals)
                signals = {signals};
            end
            
            if isempty(signals)
                signals = signals_all;
            end
            
            this.signals = signals;
            
            if ~exist('varname', 'var')||isempty(varname)
                varname = signals{1};
            end
            
            this.file_name = fname;
            this.var_name= varname;
            
            if nargin==3
                this.data_fmt = 'timed_array';
                this.var_name = varname;
            end
        end
        
        function [X, time] = computeSignals(this, p, time) % returns a p in pts
            
            
            fname = this.file_list(p(1)).name;
            
            try
                st = load(fname);
            catch
                warning(['Could not read file ' fname '. Returning NaN trace.'] )
                X(:,:) = NaN;
            end
            % fetch time data
            if isfield(st, this.time_var_name)
                t_sig = st.(this.time_var_name);
            else
                error('from_file_signal_gen:no_time_variable' ,'File %s does not contain time data as variable  %s.', fname, this.time_var_time);
            end
            
            
            if this.ignore_time
                if size(t_sig,1)==1
                    time =t_sig;
                else
                    time = t_sig';
                end
                
                X = zeros(numel(this.signals), numel(t_sig));
            else
                X = zeros(numel(this.signals), numel(time));
            end

            
            % fetch signals data
            for isig = 1:numel(this.signals)
                % assumes that signals are in variable with [time  values]
                if isfield(st, this.signals{isig})
                    sig = st.(this.signals{isig});
                else
                    error('from_file_signal_gen:missing_signal_variable' ,'File %s does not contain signal %s data.', fname, this.signals{isig});
                end
                
                if length(sig)~=length(t_sig)
                    error('from_file_signal_gen:inconsistent_size_time_signal', 'In file %s, signal %s size is inconsistent with size of time variable %s', fname, this.signal{isig}, this.time_var_name);
                end
                
                if size(sig, 2) > size(sig, 1) % make sure we have a column vector
                    sig = sig';
                end
                
                if size(sig,2) == 2
                    sig = sig(:,2);    % robust to timed variable - maye should issue a warning
                else
                    sig = sig(:,1);
                end
                
                if this.ignore_time
                    X(isig, :) = sig';
                else
                    x = interp1(t_sig, sig, time', 'linear', 'extrap');
                    X(isig, :) = x';
                end
            end
        end
        
               
        function args = getSignalGenArgs(this)
            args = {'file_name','var_name'};
        end
        
        
        
    end
end