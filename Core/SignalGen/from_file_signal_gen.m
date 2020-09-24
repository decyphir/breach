classdef from_file_signal_gen < signal_gen
    properties
        file_list
        file_name
        var_name
        sg_name
        time_var_name
        params_from_file
        ignore_time = false   % if true, will ignore time given as parameter of Sim command and use time data instead
        init_data_script
        pts  % stores possible values of parameters
    end
    methods
        function this = from_file_signal_gen(signals, fname, varname,params_from_file,sg_name,init_data_script, p0)
            
            if ~exist('fname', 'var')
                [filenames, paths] = uigetfile( ...
                    {  '*.mat','MAT-files (*.mat)'}, ...
                    'Pick one or more files', ...
                    'MultiSelect', 'on');
                
                if isequal(filenames,0) % cancel
                    fname = {};
                else
                    if ~iscell(filenames)
                        filenames= {filenames};
                    end
                    fname = cellfun( @(c)([ paths c  ] ), filenames,'UniformOutput',false);
                end
            end
            if ~iscell(fname)
                fname= {fname};
            end
            
            if exist('init_data_script', 'var')
                if isnumeric(init_data_script)&&~exist('p0','var') % Backward compatibility trick..
                    p0 = init_data_script;
                else
                    this.init_data_script = init_data_script;
                end
            end
                
            % First thing is detecting files
            for ifn=1:numel(fname)
                files_list = {};
                dir_file_list = dir(fname{ifn});   % try using dir (if directory or wildcard like *.mat)
                if ~isempty(dir_file_list)  
                    pathstr = fileparts(fname{ifn}); % dir does not keep the path...
                    for id = 1:numel(dir_file_list)
                        if ~isempty(pathstr)
                            files_list{id} = [pathstr filesep dir_file_list(id).name];
                        elseif isdir(fname{ifn})
                            files_list{id} = [fname{ifn} filesep dir_file_list(id).name];
                        else
                            files_list{id} = fname{ifn};
                        end
                    end
                else % try finding file in the path 
                    pfe = which(fname{ifn});
                    if ~isempty(pfe)
                        files_list = {pfe};
                    elseif  exist(fname{ifn}, 'file')
                        files_list = fname(ifn);
                    else
                        files_list = {};
                    end
                end
                
                
                if isempty(this.file_list)
                    this.file_list = files_list;
                else
                    this.file_list = [this.file_list files_list];
                end
            end
            
            if isempty(this.file_list)
                error('from_file_signal_gen:no_trace_file', 'No trace file.')
            end
            
            
            this.params = {'file_idx'};
            if exist('p0', 'var')&&~isempty(p0)
                this.p0 = p0;
                read_p0 = false;
            else
                this.p0 = 1;
                read_p0 = true;
            end
            
            this.params_domain = BreachDomain('enum', [],  [1:numel(this.file_list)]);
            this.pts = [1:numel(this.file_list)];
            
            % open first file
            fname = this.file_list{1};
            st = this.load_data(fname);
            vars = fieldnames(st);
            
            % look for a time array - let's make that mandatory to begin
            % with.
            
            itime =find(strcmpi(vars, 'time'), 1);
            if isempty(itime)
                itime = find(strcmpi(vars, 'ecutime'),1);
                if isempty(itime)
                    error('The file must contain a 1D array named ''time'' (case insensitive) with properly ordered time samples.');
                end
            end
            time = st.(vars{itime});
            this.time_var_name = vars{itime};
            if any(diff(time)<0)
                error('The file must contain a 1D array named ''time'' (case insensitive) with properly ordered time samples.');
            end
            
            % Next is detecting all parameters- we're looking for
            % constant scalar parameters defined in the files - and all
            % signals - they have to be defined  in all files.
            if exist('params_from_file','var')
                if ischar(params_from_file)&&~strcmp(params_from_file,'all')
                    params_from_file = {params_from_file};
                end
            else
                params_from_file= {};
            end
            this.params_from_file = params_from_file;
            
            signals_all = {};
            
            
            for iv = 1:numel(vars)
                v = vars{iv};
                if (iv ~= itime)&&isnumeric(st.(v))   % ignore time and everything not numeric
                    if isscalar(st.(v))&&...
                            (ischar(params_from_file)&&strcmp(params_from_file, 'all'))||...
                            isequal(v,params_from_file)||...
                            iscell(params_from_file)&&ismember(v, params_from_file) % this is a parameter
                        
                        this.params = [this.params v];
                        if read_p0
                            this.p0(end+1) = st.(v);
                        end
                    elseif length(st.(v))==length(time) % looks like  a signal
                        signals_all = [signals_all  v];
                    end
                end
            end
            
            % go over all files to fetch values for the parameters and
            % create enum domains
            for ifile = 1:numel(this.file_list)
                fname = this.file_list{ifile};
                st = this.load_data(fname);
                vars = fieldnames(st);
                for ip = 2:numel(this.params)
                    this.pts(ip, ifile) = st.(this.params{ip});   % TODO some try catch for parameter defined in first file, but not in other(s)
                end
            end
            for ip = 2:numel(this.params)
                this.params_domain(ip)= BreachDomain('enum', [], unique(this.pts(ip, :)));
            end
            
            if nargin==0||isempty(signals)
                signals = signals_all;
            end
            
            if ischar(signals)
                signals = {signals};
            end
            
            this.signals = signals;
            if ~exist('varname', 'var')||isempty(varname)
                varname = signals{1};
            end
            
            this.file_name = fname;
            this.var_name= varname;
            
            if ~exist('sg_name', 'var')||isempty(sg_name)
                sg_name =  varname;
            end
            
            this.params{1} = [sg_name '_file_idx'];
            this.sg_name = sg_name;
        end
        
        function st = load_data(this,fname)
            % loads, apply script if needed
            if isempty(this.init_data_script)
                st= load(fname);            
            else
                load(fname);
                try 
                    run(this.init_data_script);
                catch
                    error('from_file_signal_gen:load_data:pb_init_script', 'Error while executing script %s', this.init_data_script);
                end
                
                vars = who;
                st = struct();
                for iv= 1:numel(vars)
                    st.(vars{iv}) = eval(vars{iv});
                end
                
            end
        end
        
        function [X, time] = computeSignals(this, p, time) % returns a p in pts
            
            fname = this.file_list{p(1)};
            X = nan(numel(this.signals), length(time));
            
            try
                st = this.load_data(fname);
            catch
                warning(['Could not read file ' fname '. Returning NaN trace.'] )
                return;
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
            args = {'file_name','var_name', 'params_from_file', 'sg_name', 'init_data_script'};
        end
        
        function sgs = split(this)
        % split splits one signal_gen into one signal_gen per signal
             for i = 1:numel(this.signals)
                 sgs{i} = this.copy();
                 sgs{i}.signals = {this.signals{i}};
                 sgs{i}.var_name = this.signals{i};
                 sgs{i}.sg_name = this.signals{i};
             end
        end
        
    end
end