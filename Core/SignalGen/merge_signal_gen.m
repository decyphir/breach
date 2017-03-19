classdef merge_signal_gen < signal_gen
    
    properties
        sgs
    end
    methods
        function this=  merge_signal_gen(sgs)
            this.sgs = sgs;
            this.signals = cellfun(@(c) (c.signals),sgs, 'UniformOutput',false);
            this.signals = [this.signals{:}];
            this.params = cellfun(@(c) (c.params),sgs, 'UniformOutput',false);
            this.params = [this.params{:}];
            this.p0 =  cell2mat(cellfun(@(c) (c.p0) ,sgs, 'UniformOutput',false ));
        end
        
        function X = computeSignals(this,p, time)
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
