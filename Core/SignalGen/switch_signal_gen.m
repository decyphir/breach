classdef switch_signal_gen < signal_gen
    properties
        sgs
        sg_1st_param_idx
        p_1st_idx
    end
    
    methods
        function this= switch_signal_gen(sgs)
            this.sgs = sgs;
            this.params = {'sg_num'};
            this.params_domain = BreachDomain('int', [1 numel(sgs)]); 
            this.p0 = 1;
            this.signals = sgs{1}.signals;
            this.signals_domain = sgs{1}.signals_domain; % discutable ...   
            
            for isg = 1:numel(sgs)

                sg = sgs{isg};
                this.sg_1st_param_idx(isg) = numel(this.params)+1;
                if isempty(sg.params_domain)
                    sg.params_domain =  repmat(BreachDomain(),1, numel(sg.params));
                end
                this.params_domain = [this.params_domain sg.params_domain];
                for ip = 1:numel(sg.params)
                    this.params = [this.params ['sg' num2str(isg) '_' sg.params{ip}]];
                end
                
                if size(sg.p0,1)>1
                    this.p0 = [this.p0 sg.p0'];
                else
                    this.p0 = [this.p0 sg.p0];
                end
                
            end
            
        end
        
        function  [X, time] = computeSignals(this, p, time)
            idx_sg = p(1);
            sg = this.sgs{idx_sg};
            idx_p=this.sg_1st_param_idx(idx_sg);
            num_p= numel(sg.params);
            X= sg.computeSignals(p(idx_p:idx_p+num_p-1),time);
        end
        
    end
end