classdef from_file_signal_gen < signal_gen
    properties
        file_list
        file_name
        var_name
        data_fmt = 'timed_variables' % other could be 'timed_array', 'struct_with_time', etc
    end
    methods
        function this = from_file_signal_gen(signals, fname, varname, data_fmt)
            
            if ~exist('fname', 'var')
                fname=  '*.mat';
            end
            
            if ischar(signals)
                signals = {signals};
            end
            
            this.signals = signals;
            if ~iscell(fname)
                fname= {fname};
            end
            
            if ~exist('varname', 'var' )
                varname = signals{1};
            end
            
            if (nargin >= 4)  
                this.data_fmt = data_fmt;
            end
            
            this.file_name = fname;
            this.var_name= varname;
            
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
            
            this.params = {'file_idx'};
            this.p0 = 1;
            this.params_domain = BreachDomain('int', [1 numel(this.file_list)]);
            if nargin==3 
                this.data_fmt = 'timed_array';
                this.var_name = varname;
            end
        end
        
        function X= computeSignals(this, p, time)
            X = zeros(numel(this.signals), numel(time));
            fname = this.file_list(p).name;
            try
                st = load(fname);
            catch
                warning(['Could not read file ' fname '. Returning NaN trace.'] )
                X(:,:) = NaN;
            end
            
            switch (this.data_fmt)
                
                case 'timed_variables'
                    
                    for isig = 1:numel(this.signals)
                        % assumes that signals are in variable with [time  values]
                        try
                            sig = st.(this.var_name{isig});
                            t_sig =  sig(:,1);
                            v_sig = sig(:,2);
                            x = interp1(t_sig, v_sig, time', 'linear', 'extrap');
                            X(isig, :) = x';
                        catch
                            warning(['Could not read variable ' this.signals{isig} '. Returning NaN'] )
                            X(isig, :) = NaN;
                        end
                    end
                    
                case 'timed_array'
                    
                    for isig = 1:numel(this.signals)
                        % assumes that signals are in variable with [time  values]
                        sig = st.(this.var_name);
                        t_sig =  sig(:,1);
                        v_sig = sig(:,isig+1);
                        x = interp1(t_sig, v_sig, time', 'linear', 'extrap');
                        X(isig, :) = x';
                    end
                    
                    
            end
            
            
        end
        
        
        function args = getSignalGenArgs(this)
            args = {'file_name','var_name'};
        end
        
        
        
    end
end