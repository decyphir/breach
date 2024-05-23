classdef TA_signal_gen2 < var_cp_signal_gen
    properties
        TA_file
        num_evts
        labels
        wordgen_exe = 'wordgen'
        poly = 5
        expected_duration = 0
        template_in
        template_out
        params_cp                
        verbose=0 
        min_dt = 1e-9
    end
    
    methods
        function this= TA_signal_gen2(sigs, TA_file, labels, num_evt, method)
            
            if ~exist('method', 'var')
                method = 'previous';
            end
            
            if ~iscell(sigs)
                sigs= {sigs};
            end
            if ~iscell(labels)
                labels= {labels};
            end
                        
            this = this@var_cp_signal_gen(sigs, repmat(num_evt+2,1, numel(sigs)), method);
            
            this.signals = [this.signals 'timeword'];
            this.signals_domain(end+1) = BreachDomain('enum',[], 1:numel(labels));
            
            this.params_cp = this.params; 
            this.TA_file = TA_file;
            
            [this.wordgen_exe found] = find_wordgen_exe();
            
            if ~found
                warning('wordgen executable not found.');
            end

            
            this.labels=labels;
            
            this.num_evts = num_evt;
            this.params = {'time_scale'};
                        
            for ie = 0:num_evt-1
                this.params{end+1} = ['e' num2str(ie) '_dt'];
                this.params{end+1} = ['e' num2str(ie) '_branching' ];
                for i_l=1:numel(labels)
                    for i_sig=1:numel(sigs)
                        this.params{end+1} = ['e' num2str(ie) '_' sigs{i_sig} '_' labels{i_l} '_val'];
                    end
                end                            
            end
            
            this.params = this.params;            
            this.params_domain = repmat(BreachDomain(),1,numel(this.params));            
            this.p0  = zeros(numel(this.params),1);
            this.p0(1) = 1; % timescale
            idx_dt_branching = this.get_idx_dt_branching();
            this.p0(idx_dt_branching) = 0.5;
            
            this.template_in = repmat('%g[%g]', 1, num_evt);
            this.template_out = [repmat('%g[%s] ', 1, num_evt-1) '%g[%s]'];
            
        end
                    
        function [X,time] = computeSignals(this,p, time) % compute the signals
            [p_cp, dts, labels_idx] = computeParams(this, p);           
            X = computeSignals@var_cp_signal_gen(this, p_cp, time);            
            dts(dts<this.min_dt) = this.min_dt;
            ts = cumsum(dts);
            fprintf('%g ', sum(dts));
            X(end+1, :) = interp1([0 ts(1:end-1)], labels_idx, time, 'previous', 'extrap');        
        end
                
        function [p_cp, dts, labels_idx] = computeParams(this, p)
            
            exe = this.wordgen_exe;
            in_file = this.TA_file;
            num_pts = size(p,2);
            
            p_cp = zeros(numel(this.params_cp),1);
                        
            
            cmd_template = [strrep(exe,'\','\\') ' ' strrep(in_file,'\','\\') ' --template "%s"'  ...
                     ' --poly ' num2str(this.poly)];

            if this.expected_duration > 0
                cmd_template = [cmd_template ...
                    '--expected-duration' num2str(this.expected_duration)];
            end

            cmd_template = [cmd_template ...                     
                     ' -v ' num2str(this.verbose)];
            
            if ispc % for some reason caml cannot write cache file on windows
            cmd_template = [cmd_template ...
                     ' --no-cache ']; 
            end

            cmd_template = [cmd_template ...
                     ' --traj 1 --exact-rational --output-format timeword res.txt'];


            for ipt = 1:num_pts
                
                time_scale = p(1,ipt);
                idx_dt_branching = this.get_idx_dt_branching();
                p_wordgen = p(idx_dt_branching,ipt);
                
                word_template = sprintf(this.template_in, p_wordgen); % form wordgen template for this run 
                    
                cmd = sprintf(cmd_template, word_template);
                
                
%                 cmd = [exe ' ' in_file ' --template "' sprintf(this.template_in, p_wordgen) ...
%                     '" --poly ' num2str(this.poly)  ...
%                     ' --traj 1 --exact-rational --output-format timeword res.txt'];
                
                if this.verbose
                    disp(cmd);
                end
                stat = system(cmd);
                
                if stat
                    warning('wordgen returned error, or incorrect command: %s', cmd);
                    p_cp(2:2:end, ipt) = 1;
                else
                    fid= fopen('res.txt');
                    tline = fgetl(fid);
                    fclose(fid);
                    tok = regexp(tline, '(.*?)\[(\w+)\] ', 'tokens');
                                                            
                    for itok=1:numel(tok)
                        idx = find(strcmp(tok{itok}{2},this.labels), 1);
                        if isempty(idx)
                            error('Unknown label %s',tok{itok}{2});                            
                        else
                            dts(itok) = str2num(tok{itok}{1})*time_scale;
                            labels_idx(itok) = idx;
                        end                        
                    end
                    %fprintf('%g \n', sum(dts)/time_scale);
                    
                    idx = 1; 
                    for isig = 1:numel(this.signals)-1   % assign values from label to signals and control points                      
                        for i_tok = 1:this.num_cp(isig)-1
                            if i_tok<= numel(dts)
                                ilabel = labels_idx(i_tok);
                                iv = this.get_cp_idx(i_tok, isig, ilabel);
                                p_cp(idx, ipt) = p(iv);
                                p_cp(idx+1, ipt) = max(dts(i_tok), this.min_dt);
                                time_durations(idx) = max(dts(i_tok), this.min_dt);
                            end
                            idx = idx+2; % jump over dt..
                        end
                        if i_tok<= numel(dts)
                            p_cp(idx, ipt) = p(this.num_evts+1+(numel(this.signals)*(labels_idx(i_tok)-1)+isig),1) ;
                        end
                        idx = idx+1;                        
                    end
                                  
                end
                
            end
        end
        
        function type = getType(this)
            type = 'varstep';
        end
        
        function args = getSignalGenArgs(this)
            args = {'num_cp','TA_file','labels','method'};
        end
        
    end
    
    methods (Access=private)
        
        function idx = get_idx_dt_branching(this)
                num_param_per_evt = 2+(numel(this.signals)-1)*numel(this.labels);                
                idx_dt = [0 rem((1:numel(this.params)-1)-1,num_param_per_evt)==0];
                idx_branching = [0 rem((1:numel(this.params)-1)-1,num_param_per_evt)==1];
                idx= idx_dt|idx_branching;                                           
        end
        
        function iv = get_cp_idx(this, ie, isig, ilabel)
            num_param_per_evt = 2+(numel(this.signals)-1)*numel(this.labels);                
            iv = 1+...                             % timescale
                 (ie-1)*num_param_per_evt+...      % num events
                 2+...                             % dt,branching 
                 (ilabel-1)*(numel(this.signals)-1)+... % skips labels 
                 isig;
        end
        
        function i_labels = get_labels_idx(this)
           i_labels= this.num_evts+1:this.num_evts+numel(this.signals)*numel(this.labels);
        end
                
        
    end
    
    
end
