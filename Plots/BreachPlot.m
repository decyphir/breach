classdef BreachPlot < handle
    
    properties
        BrSet
        fig_handle
        params
        listener
    end
    
    methods
        
        function this = BreachPlot(BrSet, params)
            if nargin == 0
                return;
            end
            this.BrSet = BrSet;
            this.fig_handle = figure;
            this.params = params; 
            this.addBrPListener
            this.updatePlot();
        end
     
        function addBrPListener(this)
            mco = ?BreachSet;
            prop_list = mco.PropertyList;
            for ip = 1:numel(prop_list)
                if strcmp(prop_list(ip).Name, 'P')
                    prop  = prop_list(ip);
                    break
                end
            end
            this.listener = event.proplistener(this.BrSet, prop, 'PostSet', @(o,e)updatePlot(this,o,e));
            set(this.fig_handle, 'DeleteFcn', @(o,e)(this.delBrPListener));
        end
        
        function delBrPListener(this)
            delete(this);
        end
        
        function updatePlot(this, o, e)
            figure(this.fig_handle);
            cla;
            this.BrSet.PlotParams(this.params);
        end
            
    end
end