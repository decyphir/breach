classdef TA_param_gen < param_gen
    properties
       TA_file 
       num_events  
       wordgen_exe = './wordgen'       
       poly = 5
       template_in
       template_out       
       val_map       
       p0_out
       verbose=0
       min_dt = 1e-9
    end
    
    methods
        function this= TA_param_gen(sig, TA_file, num_events)
            this.TA_file = TA_file;
            this.num_events = num_events;
            this.params = repmat({'params_e'}, 1, num_events);          
            this.params_out = repmat({'params_dt'}, 1, num_events);
            
            this.val_map= containers.Map();
            
            this.val_map('a') = 0;
            this.val_map('b') = 1;
            this.val_map('c') = 2;
            this.val_map('d') = 3;
            this.val_map('e') = 4;
                        
            for ie = 0:num_events-1
                this.params{ie+1} = [sig '_e' num2str(ie)];
                this.params_out{2*ie+1} = [sig '_u' num2str(ie)];
                this.params_out{2*ie+2} = [sig '_dt' num2str(ie)];
            end
            this.params{end+1} = 'time_scale';
            this.params{end+1} = 'alpha'; % linear transform of mapping
            this.params{end+1} = 'beta';
            
            this.params_out{2*num_events+1} = [sig '_u' num2str(num_events)];                        
            
            this.domain = repmat(BreachDomain(),1,numel(this.params));
            this.domain_out = repmat(BreachDomain(),1,numel(this.params_out));
            
            this.p0  = [.5*ones(num_events,1); ... 
                         [1 .5 0]']; % time_scale, alpha, beta
            this.template_out = [repmat('%g ', 1, num_events-1) '%g'];
            
        end
        
        function set_template_in(this, template_in)                       
           [~,~,~, matches, tokens] = regexp(template_in, '\[(\w+)\]');
           this.p0_out = zeros(2*this.num_events+1,1);
           this.template_in = template_in;
           for im = 1:numel(matches)
              this.p0_out(2*im-1,1) = this.val_map(tokens{im}{1});
           end
        end
                       
        function p_out = computeParams(this, p_in)
            
            exe = this.wordgen_exe;
            in_file = this.TA_file;
            num_pts = size(p_in,2);        
                        
            p_out = repmat(this.p0_out, 1, num_pts);
            
            for ipt = 1:num_pts
                
                time_scale = p_in(end-2,ipt);
                alpha= p_in(end-1, ipt);
                beta = p_in(end, ipt);
                p = p_in(1:end-3, ipt);
                p_out(1:2:end,ipt) = alpha*(p_out(1:2:end,ipt) +beta);                
                
                cmd = [exe ' ' in_file ' --template "' sprintf(this.template_in, p) ...
                      '" --poly ' num2str(this.poly) ' --traj 1 -v 0 res.txt'];                
                if this.verbose
                    disp(cmd);
                end
                stat = system(cmd);
                                
                if stat
                    warning('wordgen returned error, or incorrect command: %s', cmd);
                    p_out(2:2:end, ipt) = 1;
                else
                    f= fopen('res.txt');
                    dts = fscanf(f,this.template_out)*time_scale;
                    dts(dts==0) = this.min_dt;
                    p_out(2:2:end, ipt) = dts;
                    fclose(f);
                end
                
            end
        end
        
    end

    
end
