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
        vac_plot
        neg_plot        
        all_plot
    end
    
    methods
        
        function this = BreachSamplesPlot(BrSet, params, varargin)
            % BreachSamplesPlot Initial implementation meant to navigate the summary of a BreachRequirement evaluation
            switch nargin
                case 0
                    return;
            end
            
            opt.Fig = 'create';
            
            this.BrSet = BrSet;
            
            if strcmp(opt.Fig,'create')
                this.Fig = figure;
            end
            
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
            if ~isempty(this.params)
                var = this.params;
                if numel(var)>=1
                    this.x_axis = 'idx';
                    this.y_axis = var{1};
                end
                if numel(var)>=2
                    this.x_axis= var{1}; 
                    this.y_axis = var{2};
                end
                if numel(var)>=3
                    this.z_axis = var{3};
                end
                
            elseif ~this.has_rob()
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
            if this.has_rob()
                rob = this.summary.requirements.rob;
                rob_vac =  this.summary.requirements.rob_vac;
            else
                rob = 0*all_pts'+1;
                rob_vac = rob;
                rob_vac(:) = nan;
            end
            
            this.data.all_pts.idx = all_pts;
            
            %% robustnesses of each requirement            
            if this.has_rob()
                this.data.req_names = {};
                for ir = 1:numel(this.BrSet.req_monitors)
                    name = this.BrSet.req_monitors{ir}.name;
                    this.data.req_names{end+1} = name;
                    vals_req = rob(:,ir);
                    
                    vals_pos = vals_req';
                    vals_pos(vals_req'<=0) = 0;
                    vals_pos(vals_req'<=0) = 0;
                    
                    
                    % vacuous pos
                    vals_vac = rob_vac(:,ir)';
                    vals_vac(vals_pos'<inf) = 0;
                    
                    % falsified requirements
                    vals_neg = vals_req';
                    vals_neg(vals_req'>=0) = 0;
                    
                    % idx pos, neg and vac
                    num_vals_pos = sum(vals_pos>=0&vals_neg==0,1);
                    num_vals_neg = sum(vals_neg<0,1);
                    num_vals_vac = sum(vals_pos==inf,1);
                    
                    
                    idx_pos = num_vals_pos >0;
                    idx_vac = num_vals_vac >0;
                    idx_neg = num_vals_neg >0;
                    
                    this.data.(name).pos_pts.idx= [];
                    this.data.(name).pos_pts.idx_traj = [];
                    this.data.(name).pos_pts.rob = [];
                    
                    this.data.(name).vac_pts.idx= [];
                    this.data.(name).vac_pts.idx_traj = [];
                    this.data.(name).vac_pts.rob = [];
                                        
                    this.data.(name).neg_pts.idx= [];
                    this.data.(name).neg_pts.idx_traj = [];
                    this.data.(name).neg_pts.rob = [];
                    
                    if any(idx_pos)
                        this.data.(name).pos_pts.idx= all_pts(idx_pos);
                        this.data.(name).pos_pts.idx_traj = this.BrSet.P.traj_ref(this.data.(name).pos_pts.idx);
                        this.data.(name).pos_pts.rob = vals_pos(idx_pos);
                    end
                    
                    if any(idx_vac)
                        this.data.(name).vac_pts.idx= all_pts(idx_vac);
                        this.data.(name).vac_pts.idx_traj = this.BrSet.P.traj_ref(this.data.(name).vac_pts.idx);
                        this.data.(name).vac_pts.rob = vals_vac(idx_vac);
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
            vals_pos = rob';
            vals_pos(rob'<=0) = 0;
            
            vals_cum_pos = vals_pos>=0;             
            num_req = size(vals_pos,1);
            for ir = 1:num_req % for each requirequirement
              ifalse = find(vals_pos(ir,:)==0,1); 
              if ifalse
                vals_cum_pos(ir, ifalse:end)=0;                
              end
            end
            vals_cum_neg = 1-vals_cum_pos;
            
            vals_vac = vals_pos; 
            vals_vac(vals_pos<inf)=0;
            
            % falsified requirements
            vals_neg = rob';
            vals_neg(rob'>=0) = 0;
            
            % idx pos and neg
            num_vals_pos = sum(vals_pos>=0&vals_neg==0,1);
            sum_num_vals_pos =sum(vals_cum_pos,1); 
            num_vals_neg = sum(vals_neg<0,1);
            sum_num_vals_neg =sum(vals_cum_neg,1); 
            num_vals_vac = sum(vals_vac>0&vals_neg==0,1);
            idx_pos = num_vals_pos >0;
            idx_neg = num_vals_neg >0;
            idx_vac = num_vals_vac >0;
            
            
            if any(idx_pos)
                this.data.pos_pts.idx= all_pts(idx_pos);
                this.data.pos_pts.idx_traj = this.BrSet.P.traj_ref(this.data.pos_pts.idx);
                this.data.pos_pts.v_sum_pos = sum(vals_pos(:,idx_pos),1);
                this.data.pos_pts.v_num_pos = num_vals_pos(idx_pos);
                this.data.pos_pts.v_sum_num_pos = sum_num_vals_pos(idx_pos);
            end

            if any(idx_vac)
                this.data.vac_pts.idx= all_pts(idx_vac);
                this.data.vac_pts.idx_traj = this.BrSet.P.traj_ref(this.data.vac_pts.idx);
                this.data.vac_pts.v_sum_vac = sum(vals_vac(:,idx_vac),1);
                this.data.vac_pts.v_num_vac = num_vals_vac(idx_vac);
            end

            if any(idx_neg)
                this.data.neg_pts.idx= all_pts(idx_neg);
                this.data.neg_pts.idx_traj = this.BrSet.P.traj_ref(this.data.neg_pts.idx);
                this.data.neg_pts.v_sum_neg = sum(vals_neg(:, idx_neg),1);
                this.data.neg_pts.v_num_neg = -num_vals_neg(idx_neg);
                this.data.neg_pts.v_sum_num_neg = -sum_num_vals_neg(idx_neg);
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
            
            if this.has_rob()
                this.update_req_plot();
            else
                this.update_brset_plot();
            end
            
        end
        
        function b = has_rob(this)
           b = isa(this.BrSet, 'BreachRequirement')&&...           
               isfield(this.summary,'requirements')&&...
               isfield(this.summary.requirements,'rob')&&...
               ~isempty(this.summary.requirements.rob);
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
                        if B.P.DimP>B.P.DimX
                            y_label = B.P.ParamList(B.P.DimX+1);
                        else
                            y_label = '';
                        end
                    else
                        y_label = this.data.variables{1};
                    end
                    ydata = B.GetParam(y_label,idx);
                    if isempty(ydata)
                        ydata = xdata*0;
                    end
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
                txt = regexprep(txt,'_',' ');
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
                if this.BrSet.hasTraj()
                    if isempty(this.idx_tipped)
                        it = 1;
                    else
                        it = this.idx_tipped(1);
                    end
                    sig = this.signature.signals{1};
                    F = BreachSignalsPlot(this.BrSet,sig, it);
                else
                    warning('BreachSamplesPlot:no_signal','No signal(s) to plot.');
                end
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
            none_z = {'none (2D plot)', 'none', '', '2D'};
            
            has_pos = isfield(this.data, 'pos_pts');
            has_vac = isfield(this.data, 'vac_pts');
            has_neg = isfield(this.data, 'neg_pts');
            
            %% positive samples
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
                        plot_this = @plot_num;
                        
                    case 'sum'
                        ydata_pos = this.data.pos_pts.v_sum_pos;
                        plot_this = @plot_sum;
                    case 'num'
                        ydata_pos = this.data.pos_pts.v_num_pos;
                        if any(strcmp(this.z_axis, none_z))
                            plot_this = @plot_num;
                        end
                    case 'sum_num'                        
                        ydata_pos = this.data.pos_pts.v_sum_num_pos;
                        if any(strcmp(this.z_axis, none_z))
                            plot_this = @plot_sum_num;
                        end
                    otherwise
                        if ismember(this.y_axis, this.data.req_names)
                            pos_idx = this.data.(this.y_axis).pos_pts.idx;
                            ydata_pos  = this.data.(this.y_axis).pos_pts.rob;
                        else
                            ydata_pos  = B.GetParam(this.y_axis, pos_idx);
                        end
                        if any(strcmp(this.z_axis, none_z))
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
            
            %% Vacuous samples
            if has_vac
                vac_idx = this.data.vac_pts.idx;
                
                switch this.z_axis
                    case none_z
                    case 'idx'
                        zdata_vac = vac_idx;
                        plot_this = @plot3_param;
                    case 'sum'
                        zdata_vac = this.data.vac_pts.v_sum_vac;
                        plot_this = @plot3_sum;
                    case 'num'
                        zdata_vac = this.data.vac_pts.v_num_vac;
                        plot_this = @plot3_num;
                    otherwise  % assumes z_axis is a parameter name
                        if ismember(this.z_axis, this.data.req_names)
                            vac_idx = this.data.(this.z_axis).vac_pts.idx;
                            zdata_vac  = this.data.(this.z_axis).vac_pts.rob;
                        else
                            zdata_vac  = B.GetParam(this.z_axis, vac_idx);
                        end
                        
                        plot_this = @plot3_param;
                end
                
                switch this.y_axis
                    case 'auto'
                        plot_this = @plot_num; % default
                        
                    case 'sum'
                        ydata_vac = this.data.vac_pts.v_sum_vac;
                        plot_this = @plot_sum;
                    case 'num'
                        ydata_vac = this.data.vac_pts.v_num_vac;
                        if any(strcmp(this.z_axis, none_z))
                            plot_this = @plot_num;
                        end
                    otherwise
                        if ismember(this.y_axis, this.data.req_names)
                            vac_idx = this.data.(this.y_axis).vac_pts.idx;
                            ydata_vac  = this.data.(this.y_axis).vac_pts.rob;
                        else
                            ydata_vac  = B.GetParam(this.y_axis, vac_idx);
                        end
                        if any(strcmp(this.z_axis, none_z))
                            plot_this = @plot_param;
                        end
                end
                
                
                switch this.x_axis
                    case 'idx'
                        xdata_vac = vac_idx;
                    otherwise  % assumes x_axis is a parameter name
                        xdata_vac = B.GetParam(this.x_axis, vac_idx);
                end
                
            end
            
       
            
            %% Negative samples
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
                        plot_this=@plot_num;
                    case 'sum'
                        ydata_neg = this.data.neg_pts.v_sum_neg;
                        if any(strcmp(this.z_axis, none_z))
                            plot_this = @plot_sum;
                        end
                    case 'num'
                        ydata_neg = this.data.neg_pts.v_num_neg;
                        if any(strcmp(this.z_axis, none_z))
                            plot_this = @plot_num;
                        end
                    case 'sum_num'
                        ydata_neg = this.data.neg_pts.v_sum_num_neg;
                        if any(strcmp(this.z_axis, none_z))
                            plot_this = @plot_sum_num;
                        end
                    otherwise
                        if ismember(this.y_axis, this.data.req_names)
                            neg_idx = this.data.(this.y_axis).neg_pts.idx;
                            ydata_neg = this.data.(this.y_axis).neg_pts.rob;
                        else
                            ydata_neg = B.GetParam(this.y_axis, neg_idx);
                        end
                        if any(strcmp(this.z_axis, none_z))
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
            leg = {};                
            plot_this();
            
            function plot_param()
                if has_pos&&~isempty(xdata_pos)
                    this.pos_plot = plot(xdata_pos,ydata_pos,'.g', 'MarkerSize', 20);
                end
                hold on;
                
                if has_vac&&~isempty(xdata_vac)
                    this.vac_plot = plot(xdata_vac,ydata_vac,'.m', 'MarkerSize', 20);
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
                
                if has_vac&&~isempty(xdata_vac)
                    this.vac_plot = plot3(xdata_vac,ydata_vac,zdata_vac, '.m', 'MarkerSize', 20);
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
                
                if has_vac&&~isempty(xdata_vac)
                    this.vac_plot = plot(xdata_vac,ydata_vac,'.m', 'MarkerSize', 20);
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

                if has_vac&&~isempty(xdata_vac)
                    this.vac_plot = plot3(xdata_vac,ydata_vac,zdata_vac, '.m', 'MarkerSize', 20);
                end
                
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
                    bar_pos()
                end
                
                if has_neg&&~isempty(xdata_neg)
                    bar_neg();
                end
                
                grid on;
                xlabel(this.x_axis, 'Interpreter', 'None');
                ylabel('Num. requirement falsified/satisfied');
                set(gca, 'XTick', 1:numel(this.data.all_pts.idx), 'XLim', [0 numel(this.data.all_pts.idx)+1]);
                legend(leg);
            end
            
            function bar_pos()
                if has_vac&&~isempty(xdata_vac)
                    ydata_pos = this.data.pos_pts.v_num_pos;
                    ydata_vac = this.data.vac_pts.v_num_vac;
                    this.pos_plot = bar([xdata_pos xdata_pos(end)+1], [ydata_pos nan] ,0.5,'g');
                    hold on
                    this.vac_plot = bar([xdata_vac xdata_vac(end)+1]+.1, [ydata_vac nan],0.4,'m');
                    leg = {'#Satisfied', '#Vacuous Sat.'};
                else
                    ydata_pos = this.data.pos_pts.v_num_pos;
                    this.pos_plot = bar([xdata_pos xdata_pos(end)+1] , [ydata_pos nan] ,0.5,'g');
                    hold on
                    leg =[leg {'#Satisfied'}];
                end
            end
            
            function bar_neg()
                ydata_neg = this.data.neg_pts.v_num_neg;
                this.neg_plot = bar([xdata_neg xdata_neg(end)+1], [ydata_neg nan] ,0.5,'r');                
                leg = [leg {'#Falsified'}];                
            end            
            
            function plot_sum_num()
                
                if has_pos&&~isempty(xdata_pos)
                    bar_sum_pos()
                end
                
                if has_neg&&~isempty(xdata_neg)
                    bar_sum_neg();
                end
                
                grid on;
                xlabel(this.x_axis, 'Interpreter', 'None');
                ylabel('Num. requirement falsified/satisfied');
                set(gca, 'XTick', 1:numel(this.data.all_pts.idx), 'XLim', [0 numel(this.data.all_pts.idx)+1]);
                legend(leg);
            end

            function bar_sum_pos()
                if has_vac&&~isempty(xdata_vac)
                    ydata_pos = this.data.pos_pts.v_sum_num_pos;
                    ydata_vac = this.data.vac_pts.v_num_vac;
                    this.pos_plot = bar([xdata_pos xdata_pos(end)+1], [ydata_pos nan] ,0.5,'g');
                    hold on
                    this.vac_plot = bar([xdata_vac xdata_vac(end)+1]+.1, [ydata_vac nan],0.4,'m');
                    leg = {'#Satisfied', '#Vacuous Sat.'};
                else
                    ydata_pos = this.data.pos_pts.v_sum_num_pos;
                    this.pos_plot = bar([xdata_pos xdata_pos(end)+1] , [ydata_pos nan] ,0.5,'g');
                    hold on
                    leg =[leg {'#Satisfied'}];
                end
            end
            
            function bar_sum_neg()
                ydata_neg = this.data.neg_pts.v_sum_num_neg;
                this.neg_plot = bar([xdata_neg xdata_neg(end)+1], [ydata_neg nan] ,0.5,'r');                
                leg = [leg {'#Falsified'}];                
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
                end
                xlabel(this.x_axis, 'Interpreter', 'None');
                ylabel(this.y_axis, 'Interpreter', 'None');        
                zlabel('Num. requirement falsified/satisfied');
            end

            stat = this.BrSet.GetStatement();
            st = sprintf('#traces:%g  #req:%g  #false traces:%g  #false req:%g  #vacuous sat:%g',...
                stat.num_traces_evaluated, stat.num_requirements, stat.num_traces_violations, stat.num_total_violations, stat.num_vacuous_sat);
            
            h = title(st, 'FontWeight', 'normal', 'FontSize', 10);
            set(this.Fig,'Name', 'Left click on data to get details, right click to plot signals')
            
            %% Datacursor mode customization
            cursor_h = datacursormode(this.Fig);
            cursor_h.UpdateFcn = @myupdatefcn;
            cursor_h.SnapToDataVertex = 'on';
            
            
            function [txt] = myupdatefcn(obj,event_obj)
                pos = event_obj.Position;
                ipos = find(event_obj.Target.XData==pos(1)&event_obj.Target.YData==pos(2),1);
                is_neg = 0;
                is_pos = 0;
                is_vac = 0;                
                
                if isequal(this.neg_plot, event_obj.Target)
                    is_neg = 1;
                    i_pts_req = neg_idx(ipos);
                elseif isequal(this.pos_plot, event_obj.Target)
                    is_pos = 1;
                    i_pts_req = pos_idx(ipos);                    
                elseif isequal(this.vac_plot, event_obj.Target)
                    is_vac = 1; 
                    i_pts_req = vac_idx(ipos);
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
                        rob = this.summary.requirements.rob(i_pts_req, irr);
                        if (rob==inf&&is_vac)||(rob<0)&&is_neg||(rob>=0)&&(rob<inf)&&is_pos
                            txt{end+1} = [this.summary.requirements.names{irr} ':' num2str(rob)];
                            if isfield(this.summary.requirements, 'rob_vac')
                                rob_vac = this.summary.requirements.rob_vac(i_pts_req, irr);
                                if ~isnan(rob_vac)&&rob==inf
                                    txt{end+1} = [this.summary.requirements.names{irr} ' vacuity:' num2str(rob_vac)];
                                end
                            end
                        end
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
                txt = regexprep(txt,'_',' ');
            end
            
            
            %% Context menu
            cm = uicontextmenu;
            uimenu(cm, 'Label', 'Open signals plot','Callback', @ctxtfn_signals_plot)
            
            if this.has_rob()
                top_diag = uimenu(cm, 'Label', ['Plot diagnostics']);
                for ir = 1:numel(this.summary.requirements.names)
                    uimenu(top_diag,'Label', this.summary.requirements.names{ir},'Callback', @(o,e)ctxtfn_plot_diagnostics(ir, o, e));
                end
            end
            top_x = uimenu(cm, 'Label', ['Change x-axis']);
            uimenu(top_x, 'Label', 'idx','Callback',@(o,e)(this.set_x_axis('idx')));
            
            if numel(this.data.variables)<10
                top_var = uimenu(top_x, 'Label', 'Variable');
                for ip = 1:numel(this.data.variables)
                    uimenu(top_var, 'Label', this.data.variables{ip},'Callback',@(o,e)(this.set_x_axis(this.data.variables{ip})));
                end
            else
                num_var_in_menu = 0;                   
                top_var = [];
                some_var_left = true;
                while some_var_left                    
                    var_idx_start = num_var_in_menu+1;
                    var_idx_end= min(numel(this.data.variables),num_var_in_menu+10);
                    some_var_left = (var_idx_end ~= numel(this.data.variables));
                    top_var(end+1) = uimenu(top_x, 'Label', sprintf(['Variable (%d - %d)'],var_idx_start, var_idx_end));
                    for ip = var_idx_start:var_idx_end
                        uimenu(top_var(end), 'Label', this.data.variables{ip},'Callback',@(o,e)(this.set_x_axis(this.data.variables{ip})));
                    end                    
                    num_var_in_menu = var_idx_end;
                end
                
            end
                                   
            top_y = uimenu(cm, 'Label', ['Change y-axis']);                        
            uimenu(top_y, 'Label', 'Sums of robustness','Callback',@(o,e)(this.set_y_axis('sum')));
            uimenu(top_y, 'Label', 'Number of Sat/False req.','Callback',@(o,e)(this.set_y_axis('num')));
            uimenu(top_y, 'Label', 'Cumulative Numbers of Sat/False req.','Callback',@(o,e)(this.set_y_axis('sum_num')));
      
            if numel(this.data.variables)<10
                top_var = uimenu(top_y, 'Label', 'Variable');
                for ip = 1:numel(this.data.variables)
                    uimenu(top_var, 'Label', this.data.variables{ip},'Callback',@(o,e)(this.set_y_axis(this.data.variables{ip})));
                end
            else
                num_var_in_menu = 0;                   
                top_var = [];
                some_var_left = true;
                while some_var_left                    
                    var_idx_start = num_var_in_menu+1;
                    var_idx_end= min(numel(this.data.variables),num_var_in_menu+10);
                    some_var_left = (var_idx_end ~= numel(this.data.variables));
                    top_var(end+1) = uimenu(top_y, 'Label', sprintf(['Variable (%d - %d)'],var_idx_start, var_idx_end));
                    for ip = var_idx_start:var_idx_end
                        uimenu(top_var(end), 'Label', this.data.variables{ip},'Callback',@(o,e)(this.set_y_axis(this.data.variables{ip})));
                    end                    
                    num_var_in_menu = var_idx_end;
                end
                
            end
            
            
            lim = 6;
            if numel(this.data.req_names)<lim
                top_req_y = uimenu(top_y, 'Label', 'Requirement');
                for ip = 1:numel(this.data.req_names)
                    uimenu(top_req_y, 'Label', this.data.req_names{ip},'Callback',@(o,e)(this.set_y_axis(this.data.req_names{ip})));
                end
            else
                num_req_in_menu = 0;                   
                top_req_y = [];
                some_req_left = true;                
                while some_req_left                    
                    req_idx_start = num_req_in_menu+1;
                    req_idx_end= min(numel(this.data.req_names),num_req_in_menu+lim);
                    some_req_left = (req_idx_end ~= numel(this.data.req_names));
                    top_req_y(end+1) = uimenu(top_y, 'Label', sprintf(['Requirement (%d - %d)'],req_idx_start, req_idx_end));
                    for ip = req_idx_start:req_idx_end
                        uimenu(top_req_y(end), 'Label', this.data.req_names{ip},'Callback',@(o,e)(this.set_y_axis(this.data.req_names{ip})));
                    end                    
                    num_req_in_menu = req_idx_end;
                end                
            end           
            
            top_z = uimenu(cm, 'Label', ['Change z-axis']);
            uimenu(top_z, 'Label', none_z{1},'Callback',@(o,e)(this.set_z_axis(none_z{1})));
            uimenu(top_z, 'Label', 'sum','Callback',@(o,e)(this.set_z_axis('sum')));
            
            if numel(this.data.variables)<10
                top_var_z = uimenu(top_z, 'Label', 'Variable');
                for ip = 1:numel(this.data.variables)
                    uimenu(top_var_z, 'Label', this.data.variables{ip},'Callback',@(o,e)(this.set_z_axis(this.data.variables{ip})));
                end
            else
                num_var_in_menu = 0;                   
                top_var_z = [];
                some_var_left = true;
                while some_var_left                    
                    var_idx_start = num_var_in_menu+1;
                    var_idx_end= min(numel(this.data.variables),num_var_in_menu+10);
                    some_var_left = (var_idx_end ~= numel(this.data.variables));
                    top_var_z(end+1) = uimenu(top_z, 'Label', sprintf(['Variable (%d - %d)'],var_idx_start, var_idx_end));
                    for ip = var_idx_start:var_idx_end
                        uimenu(top_var_z(end), 'Label', this.data.variables{ip},'Callback',@(o,e)(this.set_z_axis(this.data.variables{ip})));
                    end                    
                    num_var_in_menu = var_idx_end;
                end
                
            end
            
            lim = 6; 
            if numel(this.data.req_names)<lim
                top_req_z = uimenu(top_z, 'Label', 'Requirement');
                for ip = 1:numel(this.data.req_names)
                    uimenu(top_req_z, 'Label', this.data.req_names{ip},'Callback',@(o,e)(this.set_z_axis(this.data.req_names{ip})));
                end
            else
                num_req_in_menu = 0;                   
                top_req_z = [];
                some_req_left = true;                
                while some_req_left                    
                    req_idx_start = num_req_in_menu+1;
                    req_idx_end= min(numel(this.data.req_names),num_req_in_menu+lim);
                    some_req_left = (req_idx_end ~= numel(this.data.req_names));
                    top_req_z(end+1) = uimenu(top_z, 'Label', sprintf(['Requirement (%d - %d)'],req_idx_start, req_idx_end));
                    for ip = req_idx_start:req_idx_end
                        uimenu(top_req_z(end), 'Label', this.data.req_names{ip},'Callback',@(o,e)(this.set_z_axis(this.data.req_names{ip})));
                    end                    
                    num_req_in_menu = req_idx_end;
                end                
            end
                                                
            set(cursor_h, 'UIContextMenu', cm);
            set(this.ax, 'UIContextMenu', cm);
            set(this.Fig, 'UIContextMenu', cm);
            datacursormode on
            
            
            function ctxtfn_plot_diagnostics(ir, o,e)
                if isempty(this.idx_tipped)
                    it = 1;
                else
                    it = this.idx_tipped(1);
                end
                F = this.BrSet.PlotDiagnostics(ir, it);
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

