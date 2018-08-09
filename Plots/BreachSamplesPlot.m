classdef BreachSamplesPlot < handle
    
    properties
        BrSet
        Fig
        params
        data
        summary
        signature
        ax
    end
        
    
    properties (Access=protected)
        x_axis = 'idx'
        y_axis = 'auto'
        z_axis
        idx_tipped
        pos_plot
        neg_plot
        all_plot
    end
    
    methods
        
        function this = BreachSamplesPlot(BrSet, params)
        % BreachSamplesPlot Initial implementation meant to navigate the summary of a BreachRequirement evaluation    
            switch nargin
                case 0
                    return;
            end
            
            this.BrSet = BrSet;
            this.Fig = figure;
            
            if exist('params','var')
                if ischar(params)
                    this.params = {params};
                else
                    this.params =params;
                end
            else
                this.params = {};
            end
       
             this.summary = BrSet.GetSummary(); 
             this.signature = this.summary.signature;
             this.update_plot();
        
        end
        
        function update_data(this)
            
            num_samples = size(this.BrSet.P.pts,2);
            
            % variables
            if isfield(this.signature,'variables_idx')
                this.data.variables = this.signature.params(this.signature.variables_idx);
            else 
                this.data.variables = {};
            end
            
            all_pts = 1:num_samples;
            if isa(this.BrSet, 'BreachRequirement')
                vals = this.summary.requirements.rob;
            else
                vals = 0*all_pts'+1;
            end
            
            this.data.all_pts.idx = all_pts;
            
            % satisfied requirements
            vals_pos = vals';
            vals_pos(vals'<=0) = 0;
            
            % falsified requirements
            vals_neg = vals';
            vals_neg(vals'>=0) = 0;
            
            % idx pos and neg
            num_vals_pos = sum(vals_pos>=0&vals_neg==0,1);
            num_vals_neg = sum(vals_neg<0,1);
            idx_pos = num_vals_pos  >0;
            idx_neg = num_vals_neg >0;
            
            if any(idx_pos)
                this.data.pos_pts.idx= all_pts(idx_pos);
                this.data.pos_pts.idx_traj = this.BrSet.P.traj_ref(this.data.pos_pts.idx);
                this.data.pos_pts.v_sum_pos = sum(vals_pos(:,idx_pos),1);
                this.data.pos_pts.v_num_pos = num_vals_pos(idx_pos);
            end
            
            if any(idx_neg)
                this.data.neg_pts.idx= all_pts(idx_neg);
                this.data.neg_pts.idx_traj = this.BrSet.P.traj_ref(this.data.neg_pts.idx);
                this.data.neg_pts.v_sum_neg = sum(vals_neg(:, idx_neg),1);
                this.data.neg_pts.v_num_neg = -num_vals_neg(idx_neg);
            end
            
        end
        
        function update_plot(this)
            
            this.update_data();
            
            if isa(this.BrSet, 'BreachRequirement')
                B = this.BrSet.BrSet;
            else
                B = this.BrSet;
            end       
            
            figure(this.Fig);
            if isempty(this.ax)
                this.ax = axes();
            else
                axes(this.ax)
            end
            cla;
            grid on;
    
            has_pos = isfield(this.data, 'pos_pts');
            has_neg = isfield(this.data, 'neg_pts');
            if has_pos
                pos_idx = this.data.pos_pts.idx;
                switch this.x_axis
                    case 'idx'
                        xdata_pos = pos_idx;
                    otherwise  % assumes x_axis is a parameter name
                        xdata_pos = B.GetParam(this.x_axis, pos_idx);
                end
                
                switch this.y_axis
                    case 'auto'
                        if ~isa(this.BrSet,'BreachRequirement')||(isa(this.BrSet,'BreachRequirement')&&numel(this.BrSet.req_monitors)==1)||...
                                has_neg&&~has_pos||...
                                has_pos&&~has_neg
                            ydata_pos = this.data.pos_pts.v_sum_pos;
                            plot_this = @plot_sum;
                        else
                            ydata_pos = this.data.neg_pts.v_num_neg;
                            plot_this = @plot_num;
                        end
                        
                    case 'sum'
                        ydata_pos = this.data.pos_pts.v_sum_pos;
                        plot_this = @plot_sum;
                    case 'num'
                        ydata_pos = this.data.pos_pts.v_num_pos;
                        plot_this = @plot_num;
                    otherwise
                        ydata_pos = B.GetParam(this.y_axis, pos_idx);
                        plot_this = @plot_param;
                end
            end
            
            if has_neg
                neg_idx = this.data.neg_pts.idx;
                
                switch this.x_axis
                    case 'idx'
                        xdata_neg = neg_idx;
                    otherwise  % assumes parameter name
                        xdata_neg = B.GetParam(this.x_axis, neg_idx);
                end
                
                switch this.y_axis
                    case 'auto'
                        if ~isa(this.BrSet,'BreachRequirement')||(isa(this.BrSet,'BreachRequirement')&&numel(this.BrSet.req_monitors)==1)||...
                                has_neg&&~has_pos||...
                                has_pos&&~has_neg
                            ydata_neg = this.data.neg_pts.v_sum_neg;
                            plot_this = @plot_sum;
                        else
                            ydata_neg = this.data.neg_pts.v_num_neg;
                            plot_this=@plot_num;
                        end
                    case 'sum'
                        ydata_neg = this.data.neg_pts.v_sum_neg;
                        plot_this = @plot_sum;
                    case 'num'
                        ydata_neg = this.data.neg_pts.v_num_neg;
                        plot_this = @plot_num;
                    otherwise
                        ydata_neg = B.GetParam(this.y_axis, neg_idx);
                        plot_this =  @plot_param;
                end
            end
            
            plot_this();
            
            function plot_param()
                if has_pos
                    this.pos_plot = plot(xdata_pos,ydata_pos,'.g', 'MarkerSize', 20);
                end
                hold on;
                if has_neg
                    this.neg_plot = plot(xdata_neg,ydata_neg,'.r', 'MarkerSize', 20);
                end
                grid on;
                xlabel(this.x_axis, 'Interpreter', 'None');
                ylabel(this.y_axis, 'Interpreter', 'None');
            end
                    
            function plot_sum()
                if has_pos
                    this.pos_plot = plot(xdata_pos,ydata_pos,'.g', 'MarkerSize', 20);
                end
                hold on;
                if has_neg
                     this.neg_plot = plot(xdata_neg,ydata_neg,'.r', 'MarkerSize', 20);
                end
                grid on;
                xlabel(this.x_axis, 'Interpreter', 'None');
                ylabel('Cumulative satisfactions/violations');
            end
     
            function plot_num()
                if has_pos
                    ydata_pos = this.data.pos_pts.v_num_pos;
                    this.pos_plot = bar(xdata_pos, ydata_pos ,0.5,'g');
                end
                hold on;
                grid on;
                if has_neg
                    ydata_neg = this.data.neg_pts.v_num_neg;
                    this.neg_plot = bar(xdata_neg, ydata_neg ,0.5,'r');
                    set(gca, 'YLim', [min(ydata_neg)-.1, max(ydata_pos)+.1],  'Ytick', ceil(min(ydata_neg)-.1):1:floor(max(ydata_pos)+.1));
                end
                xlabel(this.x_axis, 'Interpreter', 'None');
                ylabel('Num. requirement falsified/satisfied');
            end   
            h = title('Left click on data to get details, right click to plot signals', 'FontWeight', 'normal', 'FontSize', 10);
            
            
            %% Datacursor mode customization
            cursor_h = datacursormode(this.Fig);
            cursor_h.UpdateFcn = @myupdatefcn;
            cursor_h.SnapToDataVertex = 'on';
            datacursormode on
            
            function [txt] = myupdatefcn(obj,event_obj)
                pos = event_obj.Position;
                ipos = find(event_obj.Target.XData==pos(1)&event_obj.Target.YData==pos(2),1); 
                if isequal(this.neg_plot, event_obj.Target)
                    i_pts_req = neg_idx(ipos);
                elseif isequal(this.pos_plot, event_obj.Target)
                    i_pts_req = pos_idx(ipos);
                end
                
                this.idx_tipped = i_pts_req;
                
                if isa(this.BrSet,'BreachRequirement')
                    data_trace_idx_ = this.BrSet.GetParam('data_trace_idx_',i_pts_req);
                    txt{1} = ['req_val_idx_:' num2str(i_pts_req)] ;
                    txt{2} = ['data_trace_idx_:' num2str(data_trace_idx_)] ;
                    txt{3} = '--------------';
                    for irr = 1:numel(this.summary.requirements.names)
                        txt{end+1} = [this.summary.requirements.names{irr} ':' num2str(this.summary.requirements.rob(i_pts_req, irr))];
                    end
                else
                    txt{1} = ['Sample idx:' num2str(i_pts_req)] ;
                end
                
                if isfield(this.signature, 'variables_idx')
                    txt{end+1} = '--------------';
                    for irr = 1:numel(this.signature.variables_idx)
                        var_name = this.signature.params{this.signature.variables_idx(irr)};
                        var_value =B.GetParam(var_name, i_pts_req);
                        txt{end+1} = [var_name ': ' num2str(var_value)];
                    end
                end
            end
            
            %% Context menu
            cm = uicontextmenu;
            uimenu(cm, 'Label', 'Open signals plot','Callback', @ctxtfn_signals_plot)
            
            if isa(this.BrSet, 'BreachRequirement')
                top_diag = uimenu(cm, 'Label', ['Plot diagnosis']);
                for ir = 1:numel(this.summary.requirements.names)
                    uimenu(top_diag,'Label', this.summary.requirements.names{ir},'Callback', @(o,e)ctxtfn_plot_diagnosis(ir, o, e));
                end
            end
            top_x = uimenu(cm, 'Label', ['Change x-axis']);
            uimenu(top_x, 'Label', 'idx','Callback',@(o,e)(this.set_x_axis('idx')));
            for ip = 1:numel(this.data.variables)
                uimenu(top_x, 'Label', this.data.variables{ip},'Callback',@(o,e)(this.set_x_axis(this.data.variables{ip})));
            end
         
            top_y = uimenu(cm, 'Label', ['Change y-axis']);
            uimenu(top_y, 'Label', 'auto','Callback',@(o,e)(this.set_y_axis('auto')));
            uimenu(top_y, 'Label', 'sum','Callback',@(o,e)(this.set_y_axis('sum')));
            uimenu(top_y, 'Label', 'num','Callback',@(o,e)(this.set_y_axis('num')));
            for ip = 1:numel(this.data.variables)
                uimenu(top_y, 'Label', this.data.variables{ip},'Callback',@(o,e)(this.set_y_axis(this.data.variables{ip})));
            end
         
            set(cursor_h, 'UIContextMenu', cm);
         
            function ctxtfn_plot_diagnosis(ir, o,e)
                if isempty(this.idx_tipped)
                    it = 1;
                else
                    it = this.idx_tipped(1);
                end
                F = this.BrSet.PlotDiagnosis(ir, it);
                set(F.Fig,'Name', ['Trace idx= ' num2str(it)]);
            end
            
            function ctxtfn_signals_plot(o,e)
                if isempty(this.idx_tipped)
                    it = 1;
                else
                    it = this.idx_tipped(1);
                end
                sig = this.signature.signals{1};
                F = BreachSignalsPlot(this.BrSet,sig, it);
                set(F.Fig,'Name', ['Trace idx= ' num2str(it)]);
            end
            
        end
  
        function set_x_axis(this, param)
            current_axis = this.x_axis;
            try
                this.x_axis = param;
                cla(this.ax,'reset');
                this.update_plot();
            catch ME
                warning('BreachSamplesPlot:set_axis_fail', 'set_x_axis failed with error %s',   ME.message );
                this.set_x_axis(current_axis)
            end
        end
            
        function set_y_axis(this, param)
            current_axis = this.y_axis;
            try
                this.y_axis = param;
                cla(this.ax,'reset');
                this.update_plot();
            catch ME
                warning('BreachSamplesPlot:set_axis_fail', 'set_y_axis failed with error %s',   ME.message );
                this.set_y_axis(current_axis)
            end
        
        end
        
    end
end

