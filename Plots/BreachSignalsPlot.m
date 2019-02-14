classdef BreachSignalsPlot < handle
    
    properties
        BrSet
        Fig
        Axes
        ipts
    end
    
    methods
        
        function this = BreachSignalsPlot(BrSet, signals, ipts)
            
            switch nargin
                case 0
                    return;
                case 1
                    signals = BrSet.P.ParamList{1};
                    ipts = 1;
                case 2
                    ipts = 1;
            end
            
            this.BrSet = BrSet;
            this.Fig = figure;
            this.ipts = ipts;
            
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
        
        function ax = AddAxes(this, pos)
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
            ax = this.AddAxesContextMenu(ax);
            % default to horizontal zoom mode
            h = zoom;
            set(h,'Motion','horizontal','Enable','off');
        
        end
        
        function ax = AddAxesContextMenu(this, ax)
            cm = uicontextmenu(this.Fig);
            set(ax, 'UIContextMenu', cm);
            
            % all signals
            m_top_sigs = uimenu(cm, 'Label', 'Plot');
            m_sig = uimenu(m_top_sigs, 'Label', 'signal');
            for is = 1:this.BrSet.P.DimX
                sig= this.BrSet.P.ParamList{is};
                uimenu(m_sig, 'Label', sig, 'Callback', @(o,e)ctxtfn_add_signal(ax,sig,o,e));
            end
            
            %  signals by attributes
            try
                signature = this.BrSet.P.traj{this.ipts}.signature;
            catch
                signature  = this.BrSet.GetSignature();
            end
            
            % Aliases
            if ~isempty(this.BrSet.sigMap)
                att = 'aliases';
                keys = this.BrSet.sigMap.keys();
                m = uimenu(m_top_sigs, 'Label',  att);
                for is = 1:numel(keys)
                    sig= keys{is};
                    uimenu(m, 'Label', sig, 'Callback', @(o,e)ctxtfn_add_signal(ax,sig,o,e));
                end
            end
            
            for iatt = 1:numel(signature.signal_attributes) % build uimenu for attributes
                att =signature.signal_attributes{iatt};
                f = [att 's_idx'];
                signals_att = unique(signature.signals(signature.(f)));
                m = uimenu(m_top_sigs, 'Label',  att);
                for is = 1:numel(signals_att)
                    sig= signals_att{is};
                    uimenu(m, 'Label', sig, 'Callback', @(o,e)ctxtfn_add_signal(ax,sig,o,e));
                end
            end
            
            if ismember('quant_sat_signal',signature.signal_attributes)
                m = uimenu(cm, 'Label', 'Highlight false intervals');
                att ='quant_sat_signal';
                f = [att 's_idx'];
                signals_att = unique(signature.signals(signature.(f)));
                if ismember('predicate',signature.signal_attributes)
                    att = 'predicate';
                    f = [att 's_idx'];
                    signals_att = unique([signals_att signature.signals(signature.(f))]);
                end
                for is = 1:numel(signals_att)
                    sig= signals_att{is};
                    uimenu(m, 'Label', sig, 'Callback', @(o,e)ctxtfn_highlight_false(ax,sig,o,e));
                end
            end
            uimenu(cm, 'Label', 'Add axes above','Separator', 'on', 'Callback', @(o,e)ctxtfn_add_axes_above(ax,o,e));
            uimenu(cm, 'Label', 'Add axes below', 'Callback', @(o,e)ctxtfn_add_axes_below(ax, o,e));
            uimenu(cm, 'Label', 'Reset axes','Separator', 'on', 'Callback', @(o,e)ctxtfn_reset_axes(ax, o,e));
            uimenu(cm, 'Label', 'Delete axes', 'Callback', @(o,e)ctxtfn_delete_axes(ax, o,e));
                
            
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
            
            function ctxtfn_reset_axes(ax, ~,~)
                cla(ax);
                title('');
                legend off;
            end
            
            function ctxtfn_add_signal(ax, sig, ~,~)
                this.plot_signal(sig, ax);
                this.update_legend(ax);
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
            if ~isempty(this.Axes)
                linkaxes(this.Axes, 'x');
            end
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
            this.update_legend(ax);
            
        end
        
        function int_false= HighlightFalse(this, sig, ax)
            if ~exist('ax', 'var')||isempty(ax)
                ax = this.Axes(end);
            end
            
            axes(ax);
            hold on;
            traj = this.BrSet.P.traj{this.ipts};
            idx = FindParam(this.BrSet.P, sig);
            tau = traj.time;
            val = traj.X(idx,:);
            int_false = highlight_truth_intervals(tau,val, 'g', 0, 'r', 0.3);
            set(ax,'UserData', int_false);
            this.update_legend(ax);
            
        end
        
    end
    
    methods (Access=protected)
        
        function update_legend(this, ax)
            l = legend('-DynamicLegend');
            c = flipud(get(ax, 'Children'));
            num_patch = 0;
            lc = [];
            st = {};
            for idx = 1:numel(c)
                if isa(c(idx),'matlab.graphics.primitive.Patch')
                    if num_patch==0
                        lc(end+1) = c(idx);
                        patch_idx = numel(lc);
                        st{end+1} = 'F';
                    end
                    num_patch = num_patch+1;
                else
                    lc(end+1) = c(idx);
                    st{end+1} = get(c(idx), 'DisplayName');
                end
                    
            end
            
            if num_patch>0
                set(lc(patch_idx),'DisplayName' ,['Violations (' num2str(num_patch) ')']);
                st{patch_idx} = get(lc(patch_idx),'DisplayName'); 
            end
            l= legend(lc,st);
            set(l, 'Interpreter', 'None');
        end
        
        function plot_signal(this, sig, ax)
            axes(ax);
            hold on;
            itraj = this.BrSet.P.traj_ref(this.ipts);
            time = this.BrSet.P.traj{itraj}.time;
            sig_values = this.BrSet.GetSignalValues(sig, itraj);
            if ~isempty(sig_values)
                l = plot(time , sig_values, 'DisplayName', sig);
            end
        end
    end
end