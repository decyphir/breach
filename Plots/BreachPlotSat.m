classdef BreachPlotSat < BreachPlot
    % BreachPlotSat Maps PlotRobustSat to a listening plot
    
    properties
        spec
        depth
        tau
        ipts
    end
    
    methods
        
        function this = BreachPlotSat(BrSet, spec, depth, tau, ipts)
            this.BrSet = BrSet;
            this.fig_handle = figure;
            
            this.spec = spec;
            this.depth = depth;
            this.tau = tau;
            this.ipts = ipts;
            
            this.addBrPListener()
            this.updatePlot();
        end
        
        function updatePlot(this, o, e)
            figure(this.fig_handle);
            clf;
            this.BrSet.PlotRobustSat(this.spec, this.depth, this.tau, this.ipts);
        end
        
    end
end