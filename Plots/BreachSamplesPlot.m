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
        z_axis = 'none (2D plot)'
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
            
            this.summary = BrSet.GetSummary();
            this.signature = this.summary.signature;
            
            if exist('params','var')
                if ischar(params)
                    this.params = {params};
                else
                    this.params =params;
                end
            else
                this.params = {};
            end
            
            this.update_data();
            
            %% default axis
            if ~isa(this.BrSet, 'BreachRequirement')
                var = this.data.variables;
                if numel(var)>=1
                    this.x_axis = var{1};
                end
                if numel(var)>=2
                    this.y_axis = var{2};
                end
                if numel(var)>=3
                    this.z_axis = var{3};
                end
            end
            
            this.update_plot();
            
        end
        
        function update_data(this)
            % update_data computes plottable quantities
            
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
            
            %% robustnesses of each requirement
            
            if isa(this.BrSet, 'BreachRequirement')
                this.data.req_names = {};
                for ir = 1:numel(this.BrSet.req_monitors)
                    name = this.BrSet.req_monitors{ir}.name;
                    this.data.req_names{end+1} = name;
                    vals_req = vals(:,ir);
                    vals_pos = vals_req';
                    vals_pos(vals_req'<=0) = 0;
                    
                    % falsified requirements
                    vals_neg = vals_req';
                    vals_neg(vals_req'>=0) = 0;
                    
                    % idx pos and neg
                    num_vals_pos = sum(vals_pos>=0&vals_neg==0,1);
                    num_vals_neg = sum(vals_neg<0,1);
                    idx_pos = num_vals_pos >0;
                    idx_neg = num_vals_neg >0;
                    
                    this.data.(name).pos_pts.idx= [];
                    this.data.(name).pos_pts.idx_traj = [];
                    this.data.(name).pos_pts.rob = [];
                    
                    this.data.(name).neg_pts.idx= [];
                    this.data.(name).neg_pts.idx_traj = [];
                    this.data.(name).neg_pts.rob = [];
                    
                    if any(idx_pos)
                        this.data.(name).pos_pts.idx= all_pts(idx_pos);
                        this.data.(name).pos_pts.idx_traj = this.BrSet.P.traj_ref(this.data.(name).pos_pts.idx);
                        this.data.(name).pos_pts.rob = vals_pos(idx_pos);
                    end
                    
                    if any(idx_neg)
                        this.data.(name).neg_pts.idx= all_pts(idx_neg);
                        this.data.(name).neg_pts.idx_traj = this.BrSet.P.traj_ref(this.data.(name).neg_pts.idx);
                        this.data.(name).neg_pts.rob = vals_neg(idx_neg);
                    end
                end
            end
            %% sum and num
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
            
            figure(this.Fig);
            if isempty(this.ax)
                this.ax = axes();
            else
                axes(this.ax)
            end
            cla;
            grid on;
            
            if isa(this.BrSet, 'BreachRequirement')
                this.update_req_plot();
            else
                this.update_brset_plot();
            end
            
        end
        
        function update_brset_plot(this)
            B = this.BrSet;
            idx = this.data.pos_pts.idx;
            
            % xdata
            switch this.x_axis
                case 'idx'
                    xdata = idx;
                    x_label = 'idx sample';
                otherwise  % assumes x_axis is a parameter name
                    xdata = B.GetParam(this.x_axis, idx);
                    x_label = this.x_axis;
            end
            
            % ydata
            switch this.y_axis
                case 'idx'
                    ydata = idx;
                    y_label = 'idx sample';
                case 'auto'
                    if isempty(this.data.variables)
                        y_label = B.P.ParamList(B.P.DimX+1);
                    else
                        y_label = this.data.variables{1};
                    end
                    ydata = B.GetParam(y_label,idx);
                otherwise  % assumes y_axis is a parameter name
                    ydata = B.GetParam(this.y_axis, idx);
                    y_label = this.y_axis;
            end
            
            % zdata
            none_z = 'none (2D plot)';
            switch this.z_axis
                case none_z
                    this.pos_plot = plot(xdata,ydata,'.b', 'MarkerSize', 20);
                    hold on;
                    grid on;
                    xlabel(x_label, 'Interpreter', 'None');
                    ylabel(y_label, 'Interpreter', 'None');
                otherwise  % 3d plot - assumes z_axis is a parameter name
                    zdata = B.GetParam(this.z_axis, idx);
                    this.pos_plot = plot3(xdata,ydata,zdata,'.b', 'MarkerSize', 20);
                    hold on;
                    grid on;
                    xlabel(x_label, 'Interpreter', 'None');
                    ylabel(y_label, 'Interpreter', 'None');
                    zlabel(this.z_axis, 'Interpreter', 'None');
            end
            
            
            %% Datacursor mode customization
            cursor_h = datacursormode(this.Fig);
            cursor_h.UpdateFcn = @myupdatefcn;
            cursor_h.SnapToDataVertex = 'on';
            
            function [txt] = myupdatefcn(obj,event_obj)
                pos = event_obj.Position;
                ipos = find(event_obj.Target.XData==pos(1)&event_obj.Target.YData==pos(2),1);
                i_pts = idx(ipos);
                
                this.idx_tipped = i_pts;
                txt{1} = ['Sample idx:' num2str(i_pts)] ;
                if isfield(this.summary, 'file_names')
                    [p, f ,e ] = fileparts(this.summary.file_names{i_pts});
                    txt{end+1} = ['File Name: ' f e];
                end
                
                if isfield(this.signature, 'variables_idx')
                    txt{end+1} = '--------------';
                    for irr = 1:numel(this.signature.variables_idx)
                        var_name = this.signature.params{this.signature.variables_idx(irr)};
                        var_value =B.GetParam(var_name, i_pts);
                        txt{end+1} = [var_name ': ' num2str(var_value)];
                    end
                end
            end
            
            %% Context menu
            
            cm_proj = uicontextmenu;
            
            top_x = uimenu(cm_proj, 'Label', ['Change x-axis']);
            uimenu(top_x, 'Label', 'idx','Callback',@(o,e)(this.set_x_axis('idx')));
            for ip = 1:numel(this.data.variables)
                uimenu(top_x, 'Label', this.data.variables{ip},'Callback',@(o,e)(this.set_x_axis(this.data.variables{ip})));
            end
            
            top_y = uimenu(cm_proj, 'Label', ['Change y-axis']);
            uimenu(top_y, 'Label', 'idx','Callback',@(o,e)(this.set_y_axis('idx')));
            for ip = 1:numel(this.data.variables)
                uimenu(top_y, 'Label', this.data.variables{ip},'Callback',@(o,e)(this.set_y_axis(this.data.variables{ip})));
            end
            
            top_z = uimenu(cm_proj, 'Label', ['Change z-axis']);
            uimenu(top_z, 'Label', none_z,'Callback',@(o,e)(this.set_z_axis(none_z)));
            for ip = 1:numel(this.data.variables)
                uimenu(top_z, 'Label', this.data.variables{ip},'Callback',@(o,e)(this.set_z_axis(this.data.variables{ip})));
            end
            
            set(this.ax, 'UIContextMenu', cm_proj);
            set(this.Fig, 'UIContextMenu', cm_proj);
            
            cm = uicontextmenu;
            uimenu(cm, 'Label', 'Open signals plot','Callback', @ctxtfn_signals_plot)
            set(cursor_h, 'UIContextMenu', cm);
            
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
            
            top_x = uimenu(cm, 'Label', ['Change x-axis']);
            uimenu(top_x, 'Label', 'idx','Callback',@(o,e)(this.set_x_axis('idx')));
            for ip = 1:numel(this.data.variables)
                uimenu(top_x, 'Label', this.data.variables{ip},'Callback',@(o,e)(this.set_x_axis(this.data.variables{ip})));
            end
            
            top_y = uimenu(cm, 'Label', ['Change y-axis']);
            uimenu(top_y, 'Label', 'idx','Callback',@(o,e)(this.set_y_axis('idx')));
            for ip = 1:numel(this.data.variables)
                uimenu(top_y, 'Label', this.data.variables{ip},'Callback',@(o,e)(this.set_y_axis(this.data.variables{ip})));
            end
            
            top_z = uimenu(cm, 'Label', ['Change z-axis']);
            uimenu(top_z, 'Label', none_z,'Callback',@(o,e)(this.set_z_axis(none_z)));
            for ip = 1:numel(this.data.variables)
                uimenu(top_z, 'Label', this.data.variables{ip},'Callback',@(o,e)(this.set_z_axis(this.data.variables{ip})));
            end
            
            
        end
        
        
        function update_req_plot(this)
            
            B = this.BrSet.BrSet;
            none_z = 'none (2D plot)';
            
            has_pos = isfield(this.data, 'pos_pts');
            has_neg = isfield(this.data, 'neg_pts');
            if has_pos
                pos_idx = this.data.pos_pts.idx;
                
                switch this.z_axis
                    case none_z
                    case 'idx'
                        zdata_pos = pos_idx;
                        plot_this = @plot3_param;
                    case 'sum'
                        zdata_pos = this.data.pos_pts.v_sum_pos;
                        plot_this = @plot3_sum;
                    case 'num'
                        zdata_pos = this.data.pos_pts.v_num_pos;
                        plot_this = @plot3_num;
                    otherwise  % assumes z_axis is a parameter name
                        if ismember(this.z_axis, this.data.req_names)
                            pos_idx = this.data.(this.z_axis).pos_pts.idx;
                            zdata_pos  = this.data.(this.z_axis).pos_pts.rob;
                        else
                            zdata_pos  = B.GetParam(this.z_axis, pos_idx);
                        end
                        
                        plot_this = @plot3_param;
                end
                
                switch this.y_axis
                    case 'auto'
                        if numel(this.BrSet.req_monitors)==1||...
                                has_neg&&~has_pos||...
                                has_pos&&~has_neg
                            ydata_pos = this.data.pos_pts.v_sum_pos;
                            if strcmp(this.z_axis, none_z)
                                
                                plot_this = @plot_sum;
                            end
                        else
                            ydata_pos = this.data.neg_pts.v_num_neg;
                            if strcmp(this.z_axis, none_z)
                                plot_this = @plot_num;
                                
                            end
                        end
                        
                    case 'sum'
                        ydata_pos = this.data.pos_pts.v_sum_pos;
                        plot_this = @plot_sum;
                    case 'num'
                        ydata_pos = this.data.pos_pts.v_num_pos;
                        if strcmp(this.z_axis, none_z)
                            plot_this = @plot_num;
                        end
                    otherwise
                        if ismember(this.y_axis, this.data.req_names)
                            pos_idx = this.data.(this.y_axis).pos_pts.idx;
                            ydata_pos  = this.data.(this.y_axis).pos_pts.rob;
                        else
                            ydata_pos  = B.GetParam(this.y_axis, pos_idx);
                        end
                        if strcmp(this.z_axis, none_z)
                            plot_this = @plot_param;
                        end
                end
                
                
                switch this.x_axis
                    case 'idx'
                        xdata_pos = pos_idx;
                    otherwise  % assumes x_axis is a parameter name
                        xdata_pos = B.GetParam(this.x_axis, pos_idx);
                end
                
            end
            
            if has_neg
                neg_idx = this.data.neg_pts.idx;
                
                switch this.z_axis
                    case none_z
                    case 'idx'
                        zdata_neg = neg_idx;
                        plot_this = @plot3_param;
                    case 'sum'
                        zdata_neg = this.data.neg_pts.v_sum_neg;
                        plot_this = @plot3_sum;
                    case 'num'
                        zdata_neg = this.data.neg_pts.v_num_neg;
                        plot_this = @plot3_num;
                    otherwise  % assumes z_axis is a parameter name
                        if ismember(this.z_axis, this.data.req_names)
                            neg_idx = this.data.(this.z_axis).neg_pts.idx;
                            zdata_neg = this.data.(this.z_axis).neg_pts.rob;
                        else
                            zdata_neg = B.GetParam(this.z_axis, neg_idx);
                        end
                        plot_this = @plot3_param;
                        
                end
                
                switch this.y_axis
                    case 'auto'
                        if numel(this.BrSet.req_monitors)==1||...
                                has_neg&&~has_pos||...
                                has_pos&&~has_neg
                            ydata_neg = this.data.neg_pts.v_sum_neg;
                            if strcmp(this.z_axis, none_z)
                                plot_this = @plot_sum;
                            end
                        else
                            ydata_neg = this.data.neg_pts.v_num_neg;
                            if strcmp(this.z_axis, none_z)
                                plot_this=@plot_num;
                            end
                        end
                    case 'sum'
                        ydata_neg = this.data.neg_pts.v_sum_neg;
                        if strcmp(this.z_axis, none_z)
                            plot_this = @plot_sum;
                        end
                    case 'num'
                        ydata_neg = this.data.neg_pts.v_num_neg;
                        if strcmp(this.z_axis, none_z)
                            plot_this = @plot_num;
                        end
                    otherwise
                        if ismember(this.y_axis, this.data.req_names)
                            neg_idx = this.data.(this.y_axis).neg_pts.idx;
                            ydata_neg = this.data.(this.y_axis).neg_pts.rob;
                        else
                            ydata_neg = B.GetParam(this.y_axis, neg_idx);
                        end
                        if strcmp(this.z_axis, none_z)
                            plot_this =  @plot_param;
                        end
                end
                
                
                switch this.x_axis
                    case 'idx'
                        xdata_neg = neg_idx;
                    otherwise  % assumes parameter name
                        xdata_neg = B.GetParam(this.x_axis, neg_idx);
                end
                
                
            end
            
            plot_this();
            
            function plot_param()
                if has_pos&&~isempty(xdata_pos)
                    this.pos_plot = plot(xdata_pos,ydata_pos,'.g', 'MarkerSize', 20);
                end
                hold on;
                if has_neg&&~isempty(xdata_neg)
                    this.neg_plot = plot(xdata_neg,ydata_neg,'.r', 'MarkerSize', 20);
                end
                grid on;
                xlabel(this.x_axis, 'Interpreter', 'None');
                ylabel(this.y_axis, 'Interpreter', 'None');
            end
            
            function plot3_param()
                if has_pos&&~isempty(xdata_pos)
                    this.pos_plot = plot3(xdata_pos,ydata_pos,zdata_pos, '.g', 'MarkerSize', 20);
                end
                hold on;
                if has_neg&&~isempty(xdata_neg)
                    this.neg_plot = plot3(xdata_neg,ydata_neg,zdata_neg,'.r', 'MarkerSize', 20);
                end
                grid on;
                xlabel(this.x_axis, 'Interpreter', 'None');
                ylabel(this.y_axis, 'Interpreter', 'None');
                zlabel(this.z_axis, 'Interpreter', 'None');
            end
            
            
            function plot_sum()
                if has_pos&&~isempty(xdata_pos)
                    this.pos_plot = plot(xdata_pos,ydata_pos,'.g', 'MarkerSize', 20);
                end
                hold on;
                if has_neg&&~isempty(xdata_neg)
                    this.neg_plot = plot(xdata_neg,ydata_neg,'.r', 'MarkerSize', 20);
                end
                grid on;
                xlabel(this.x_axis, 'Interpreter', 'None');
                ylabel('Cumulative satisfactions/violations');
            end
            
            function plot3_sum()
                if has_pos&&~isempty(xdata_pos)
                    this.pos_plot = plot3(xdata_pos,ydata_pos,zdata_pos, '.g', 'MarkerSize', 20);
                end
                hold on;
                if has_neg&&~isempty(xdata_neg)
                    this.neg_plot = plot3(xdata_neg,ydata_neg,zdata_neg,'.r', 'MarkerSize', 20);
                end
                grid on;
                xlabel(this.x_axis, 'Interpreter', 'None');
                ylabel(this.y_axis, 'Interpreter', 'None');
                zlabel('Cumulative satisfactions/violations');
            end
            
            
            function plot_num()
                if has_pos&&~isempty(xdata_pos)
                    ydata_pos = this.data.pos_pts.v_num_pos;
                    this.pos_plot = bar(xdata_pos, ydata_pos ,0.5,'g');
                end
                hold on;
                grid on;
                if has_neg&&~isempty(xdata_neg)
                    ydata_neg = this.data.neg_pts.v_num_neg;
                    this.neg_plot = bar(xdata_neg, ydata_neg ,0.5,'r');
                    %set(gca, 'YLim', [min(ydata_neg)-.1, max(ydata_pos)+.1],  'Ytick', ceil(min(ydata_neg)-.1):1:floor(max(ydata_pos)+.1));
                end
                xlabel(this.x_axis, 'Interpreter', 'None');
                ylabel('Num. requirement falsified/satisfied');
                set(gca, 'XTick', 1:numel(this.data.all_pts.idx) );
            end
            
            
            function plot3_num()
                if has_pos
                    ydata_pos = this.data.pos_pts.v_num_pos;
                    this.pos_plot = plot3(xdata_pos, ydata_pos, zdata_pos,'.g', 'MarkerSize', 20);
                end
                hold on;
                grid on;
                if has_neg
                    ydata_neg = this.data.neg_pts.v_num_neg;
                    this.neg_plot = plot3(xdata_neg, ydata_neg , zdata_neg,'.r', 'MarkerSize', 20);
                    %set(gca, 'YLim', [min(ydata_neg)-.1, max(ydata_pos)+.1],  'Ytick', ceil(min(ydata_neg)-.1):1:floor(max(ydata_pos)+.1));
                end
                xlabel(this.x_axis, 'Interpreter', 'None');
                ylabel(this.y_axis, 'Interpreter', 'None');        
                zlabel('Num. requirement falsified/satisfied');
            end

            h = title('Left click on data to get details, right click to plot signals', 'FontWeight', 'normal', 'FontSize', 10);
            
            
            %% Datacursor mode customization
            cursor_h = datacursormode(this.Fig);
            cursor_h.UpdateFcn = @myupdatefcn;
            cursor_h.SnapToDataVertex = 'on';
            %   datacursormode on
            
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
                    if isfield(this.summary, 'file_names')
                        [~, f ,e ] = fileparts(this.summary.file_names{data_trace_idx_});
                        txt{end+1} = ['File Name: ' f e];
                    end
                    
                    txt{end+1} = '--------------';
                    for irr = 1:numel(this.summary.requirements.names)
                        txt{end+1} = [this.summary.requirements.names{irr} ':' num2str(this.summary.requirements.rob(i_pts_req, irr))];
                    end
                else
                    txt{1} = ['Sample idx:' num2str(i_pts_req)];
                    if isfield(this.summary, 'file_names')
                        [~, f ,e ] = fileparts(this.summary.file_names{data_trace_idx_});
                        txt{end+1} = ['File Name: ' f e];
                    end
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
            for ir = 1:numel(this.data.req_names)
                uimenu(top_y, 'Label', this.data.req_names{ir},'Callback',@(o,e)(this.set_y_axis(this.data.req_names{ir})));
            end
            
            
            top_z = uimenu(cm, 'Label', ['Change z-axis']);
            uimenu(top_z, 'Label', none_z,'Callback',@(o,e)(this.set_z_axis(none_z)));
            uimenu(top_z, 'Label', 'sum','Callback',@(o,e)(this.set_z_axis('sum')));
            for ip = 1:numel(this.data.variables)
                uimenu(top_z, 'Label', this.data.variables{ip},'Callback',@(o,e)(this.set_z_axis(this.data.variables{ip})));
            end
            for ir = 1:numel(this.data.req_names)
                uimenu(top_z, 'Label', this.data.req_names{ir},'Callback',@(o,e)(this.set_z_axis(this.data.req_names{ir})));
            end
            
            set(cursor_h, 'UIContextMenu', cm);
            set(this.ax, 'UIContextMenu', cm);
            set(this.Fig, 'UIContextMenu', cm);
            datacursormode on
            
            
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
        
        function set_z_axis(this, param)
            current_axis = this.z_axis;
            try
                this.z_axis = param;
                cla(this.ax,'reset');
                this.update_plot();
            catch ME
                warning('BreachSamplesPlot:set_axis_fail', 'set_z_axis failed with error %s',   ME.message );
                this.set_z_axis(current_axis)
            end
        end
        
    end
end

