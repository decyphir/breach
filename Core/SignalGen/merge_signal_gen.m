classdef merge_signal_gen < signal_gen
    
    properties
        sgs
    end
    
    methods
        function this=  merge_signal_gen(sgs)
            this.sgs = sgs;
                
            for isg = 1:numel(sgs)
                sg= sgs{isg};
                this.signals = [this.signals  sgs{isg}.signals];
                this.params = [this.params sgs{isg}.params]; 
                
                % domains 
                num_sig = numel(sg.signals);
                if isempty(sg.signals_domain)
                    this.signals_domain = [this.signals_domain repmat(BreachDomain(),1, num_sig)];
                else
                    this.signals_domain = [this.signals_domain sg.signals_domain];
                end
                
                num_par = numel(sg.params);
                if isempty(sg.params_domain)
                    this.params_domain = [this.params_domain repmat(BreachDomain(),1, num_par)];
                else
                    this.params_domain = [this.params_domain sg.params_domain];
                end
                
                % default values
                p0sg = sg.p0;
                if size(p0sg,1) >1
                    p0sg = p0sg';
                end
                this.p0 = [this.p0 p0sg ];
            end
            
        end
        
        function [X, time] = computeSignals(this,p, time)
            X = zeros(0, numel(time));
            for isg = 1:numel(this.sgs)
                sg = this.sgs{isg};
                psg = p(1:numel(sg.params));
                p = p(numel(sg.params)+1:end);
                X(end+1:end+numel(sg.signals),:) = sg.computeSignals(psg,time);
            end
        end
    end
    
    
end
