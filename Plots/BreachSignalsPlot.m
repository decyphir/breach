classdef BreachSignalsPlot < handle
    
    properties
        BrSet
        Fig
        Axes
        itraj
    end
    
    methods
        
        function this = BreachSignalsPlot(BrSet, signals, itraj)
            
            switch nargin
                case 0
                return;
                case 1
                    signals = BrSet.P.ParamList{1};
                    itraj = 1;
                case 2
                    itraj = 1;
            end
            
            this.BrSet = BrSet;
            this.Fig = figure;
            this.itraj = itraj;
            
            if exist('signals','var')
                if ischar(signals)
                    signals = {signals};
                end
            else
                signals = {BrSet.P.ParamList{1}};
            end
            
            for is = 1:numel(signals)
                this.AddSignals(signals{is}, is);
            end
            
            
        end
               
        function AddAxes(this, pos)
        % AddAxes add new axe at specified position
            if nargin==1
                pos = numel(this.Axes)+1;
            end
            num_ax_old = numel(this.Axes);
            if numel(this.Axes)>0
                Xlim = get(this.Axes(1),'XLim');
            end 
            for ia = 1:num_ax_old;
              if ia < pos 
                  subplot(num_ax_old+1, 1, ia, this.Axes(ia))
              else
                  subplot(num_ax_old+1, 1, ia+1, this.Axes(ia))
              end
            end
            
            ax = subplot(num_ax_old+1, 1, pos);
            axes(ax);
            grid on;
                  
            if exist('Xlim', 'var')
                set(ax, 'XLim', Xlim);
            end
            this.Axes = [this.Axes(1:pos-1) ax this.Axes(pos:end)];
            figure(this.Fig);
            if numel(this.Axes)>1
                linkaxes(this.Axes, 'x');
            end
            
            cm = uicontextmenu;
            m_top_sigs = uimenu(cm, 'Label', 'Plot signal');
   
            if isa(this.BrSet, 'BreachRequirement')
                m_top_phis = uimenu(cm, 'Label', 'Highlight false intervals');
            end
            
            uimenu(cm, 'Label', 'Add Above', 'Callback', @(o,e)ctxtfn_add_axes_above(ax,o,e));
            uimenu(cm, 'Label', 'Add Below', 'Callback', @(o,e)ctxtfn_add_axes_below(ax, o,e));
            uimenu(cm, 'Label', 'Delete', 'Callback', @(o,e)ctxtfn_delete_axes(ax, o,e));
           
            for is = 1:this.BrSet.P.DimX
                sig= this.BrSet.P.ParamList{is};
                uimenu(m_top_sigs, 'Label', sig, 'Callback', @(o,e)ctxtfn_add_signal(ax,sig,o,e));
            end
            
            if isa(this.BrSet, 'BreachRequirement')
            
                for is = 1:numel(this.BrSet.formulas)
                    sig= get_id(this.BrSet.formulas{is}.formula);
                    uimenu(m_top_phis, 'Label', sig, 'Callback', @(o,e)ctxtfn_highlight_false(ax,sig,o,e));
                end
            end
            
            set(this.Axes(pos), 'UIContextMenu', cm);
            
            function ctxtfn_add_axes_above(ax, ~,~)
                for ia = 1:numel(this.Axes)
                    if isequal(ax,this.Axes(ia))
                        break;
                    end
                end
                this.AddAxes(ia);
            end
            
            function ctxtfn_add_axes_below(ax, ~,~)
                for ia = 1:numel(this.Axes)
                    if isequal(ax,this.Axes(ia))
                        break;
                    end
                end
                this.AddAxes(ia+1);
            end
            
            function ctxtfn_delete_axes(ax, ~,~)
                for ia = 1:numel(this.Axes)
                    if isequal(ax,this.Axes(ia))
                        break;
                    end
                end
                this.DeleteAxes(ia);
            end
           
            function ctxtfn_add_signal(ax, sig, ~,~)
                this.plot_signal(sig, ax);
            end
            
            function ctxtfn_highlight_false(ax, sig, ~,~)
                this.HighlightFalse(sig, ax);
            end
           
            
            
            
            
        end
        
        function DeleteAxes(this, pos)
        % DeleteAxe Remove axe  at specified position
            
            num_ax_old = numel(this.Axes);
            this.Axes(pos).delete;
            this.Axes = [this.Axes(1:pos-1) this.Axes(pos+1:end)];
            
            for ia = 1:num_ax_old-1;
                    subplot(num_ax_old-1, 1, ia, this.Axes(ia))
            end
            
            figure(this.Fig);
            linkaxes(this.Axes, 'x');
            
        end
    
        function AddSignals(this,sigs, ax)
            if ischar(sigs)
                sigs = {sigs};
            end
            
            if ~exist('ax', 'var')||isempty(ax)
                ax = numel(this.Axes)+1;
            end
            
            if isnumeric(ax)
                if ax==0 
                    this.AddAxes(1);
                    ax = this.Axes(1);
                elseif ax > numel(this.Axes) 
                    this.AddAxes(); 
                    ax =this.Axes(end);
                else
                    ax = this.Axes(ax);
                end
            end    
            if ~isa(ax, 'matlab.graphics.axis.Axes')
                error('Argument should be an Axes object or an index.')
            end
            
            for is = 1:numel(sigs)
                sig = sigs{is};
                this.plot_signal(sig, ax);
            end
        end
        
        function HighlightFalse(this, sig, ax)
            if nargin<3
                ax = this.Axes(1);
            end
            axes(ax);
            hold on;
            traj = this.BrSet.P.traj{this.itraj};
            idx = FindParam(this.BrSet.P, sig);
            tau = traj.time;
            val = traj.X(idx,:);
            highlight_truth_intervals(tau,val, 'g', 0, 'r', 0.3);
            l = legend('-DynamicLegend');
            set(l, 'Interpreter', 'None');
            
        end
        
    end
     
    methods (Access=protected)
        
        function plot_signal(this, sig, ax)
            axes(ax);
            hold on;
            traj = this.BrSet.P.traj{this.itraj}; 
            idx = FindParam(this.BrSet.P, sig);
            plot(traj.time , traj.X(idx,:), 'DisplayName', sig);
            l = legend('-DynamicLegend');
            set(l, 'Interpreter', 'None');
        end
        
    end
end