classdef TA_signal_gen < var_cp_signal_gen
    properties
        TA_file
        num_evts
        labels
        wordgen_exe = './wordgen'
        poly = 5
        template_in
        template_out
        params_cp                
        verbose=0
        min_dt = 1e-9
    end
    
    methods
        function this= TA_signal_gen(sigs, TA_file, labels, num_evt, method)
            
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
            
            this.params_cp = this.params; 
            this.TA_file = TA_file;
            
            
            this.labels=labels;
            
            this.num_evts = num_evt;
            this.params =  repmat({'params_e'}, 1, num_evt);
            
            for ie = 0:num_evt-1
                this.params{ie+1} = ['d_e' num2str(ie)];
            end
            this.params{end+1} = 'time_scale';
            
            for i_l=1:numel(labels)
                for i_sig=1:numel(sigs)
                    this.params{end+1} = [sigs{i_sig} '_' labels{i_l} '_val'];
                end
            end
          
            this.params = this.params;
            
            this.params_domain = repmat(BreachDomain(),1,numel(this.params));            
            this.p0  = zeros(numel(this.params),1);
            this.p0(1:num_evt,1) = 0.5;
            
            this.template_in = repmat('%g[0.5]', 1, num_evt);
            this.template_out = [repmat('%g[%s] ', 1, num_evt-1) '%g[%s]'];
            
        end
                    
        function X = computeSignals(this,p, time) % compute the signals
            p_cp = computeParams(this, p);
            X = computeSignals@var_cp_signal_gen(this, p_cp, time);            
        end
                
        function p_cp = computeParams(this, p_evt)
            
            exe = this.wordgen_exe;
            in_file = this.TA_file;
            num_pts = size(p_evt,2);
            
            p_cp = zeros(numel(this.params_cp),1);
            
            for ipt = 1:num_pts
                
                time_scale = p_evt(this.num_evts+1,ipt);
                p = p_evt(1:this.num_evts, ipt);                                
                               
                cmd = [exe ' ' in_file ' --template "' sprintf(this.template_in, p) ...
                    '" --poly ' num2str(this.poly) ' --traj 1 -v 0 --output-format timeword res.txt'];
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
                    tok = regexp(tline, '(.*?)\[(\w+)\] ', 'tokens');
                                                            
                    for itok=1:numel(tok)
                        idx = find(strcmp(tok{itok}{2},this.labels), 1);
                        if isempty(idx)
                            error('Unknown label %s',tok{itok}{1});                            
                        else
                            dt(itok) = str2num(tok{itok}{1})*time_scale;
                            idx_label(itok) = idx;
                        end                        
                    end
                    
                    idx = 1; 
                    for isig = 1:numel(this.signals)   % assign values from label to signals and control points                      
                        for i_tok = 1:this.num_cp(isig)-1
                            if i_tok<= numel(dt)
                                p_cp(idx, ipt) = p_evt(this.num_evts... % events
                                                       +1 ...           % time_scale
                                                       + (numel(this.signals)*(idx_label(i_tok)-1)+isig) ...  
                                                       ,1) ;
                                p_cp(idx+1, ipt) = max(dt(i_tok), this.min_dt);
                            end
                            idx = idx+2; % jump over dt..
                        end
                        if i_tok<= numel(dt)
                            p_cp(idx, ipt) = p_evt(this.num_evts+1+(numel(this.signals)*(idx_label(i_tok)-1)+isig),1) ;
                        end
                        idx = idx+1;
                        
                    end
                    
                    fclose(fid);
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
        
        function i_evt = get_evt_idx(this)
           i_evt = 1:this.num_evts;                      
        end
        
        function i_labels = get_labels_idx(this)
            i_labels= this.num_evts+1:this.num_evts+numel(this.signals)*numel(this.labels);
        end
                
        
    end
    
    
end
