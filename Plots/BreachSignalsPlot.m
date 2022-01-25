classdef BreachSignalsPlot < handle
    
    properties
        BrSet
        Fig
        Axes            
        Summary
    end
        
    properties(SetAccess=protected, GetAccess=public)
       ipts       
       zero_rob_line_name = 'zero robustness line'
    end
    
    methods
        
        function this = BreachSignalsPlot(BrSet, signals, ipts)
            
            switch nargin
                case 0
                    return;
                case 1
                    all_sigs = BrSet.GetSignalList; 
                    signals = all_sigs{1};
                    ipts = 1;
                case 2
                    ipts = 1;
            end
            
            this.BrSet = BrSet;
            this.Summary = BrSet.GetSummary();
            if isa(this.BrSet, 'BreachRequirement')
                this.Summary.num_traces =this.Summary.num_traces_evaluated;
            end
            
            this.Fig  = figure;
            
            set(this.Fig, 'KeyPressFcn', @(o,e)this.key_pressed_callback(o,e));
                        
            this.ipts = ipts;
            
            if ischar(signals)
                signals = {signals};
            end
            
            for is = 1:numel(signals) % plot signals on separate axes
                this.AddSignals(signals{is}, is);
            end                        
            
            this.update_title();
        end
        
        function key_pressed_callback(this, o,e)            
          switch e.Key
              case 'rightarrow'
                if ismember(e.Modifier,'shift')
                    this.next_ipts(ceil(this.Summary.num_traces/10));
                elseif ismember(e.Modifier,'control')
                    this.next_ipts(10);
                else
                    this.next_ipts();                 
                end
              case 'leftarrow' 
                  if ismember(e.Modifier,'shift')
                    this.prev_ipts(ceil(this.Summary.num_traces/10));
                  elseif ismember(e.Modifier,'control')
                    this.prev_ipts(10);
                  else
                    this.prev_ipts();                 
                  end
          end
        end
               
        function ax = AddAxes(this,pos)
            % AddAxes add new axe at specified position
            if nargin==1
                pos = numel(this.Axes)+1;
            end
            num_ax_old = numel(this.Axes);
            if numel(this.Axes)>0
                Xlim = get(this.Axes(1).ax,'XLim');
            end
            for ia = 1:num_ax_old
                if ia < pos
                    subplot(num_ax_old+1, 1, ia, this.Axes(ia).ax)
                else
                    subplot(num_ax_old+1, 1, ia+1, this.Axes(ia).ax)
                end
            end
            
            ax = subplot(num_ax_old+1, 1, pos);
            axes(ax);
            grid on;
            
            if exist('Xlim', 'var')
                set(ax, 'XLim', Xlim);
            end
            
            new_ax_struct.ax = ax;
            new_ax_struct.signals= {};
            new_ax_struct.robs = {};
            this.Axes = [this.Axes(1:pos-1) new_ax_struct this.Axes(pos:end)];
            if isempty(this.Fig)
                this.Fig = figure;
            else
                figure(this.Fig);
            end
            
            if numel(this.Axes)>1
                linkaxes(arrayfun(@(c)(c.ax),this.Axes),'x');
            end
            ax = this.AddAxesContextMenu(ax);
            % default to horizontal zoom mode
            h = zoom;
            set(h,'Motion','horizontal','Enable','off');        
        end
                              
        function DeleteAxes(this, pos)
            % DeleteAxe Remove axe at specified position
            
            num_ax_old = numel(this.Axes);
            this.Axes(pos).ax.delete;
            this.Axes = [this.Axes(1:pos-1) this.Axes(pos+1:end)];
            
            for ia = 1:num_ax_old-1
                subplot(num_ax_old-1, 1, ia, this.Axes(ia).ax)
            end
            
            figure(this.Fig);
            if ~isempty(this.Axes)
                linkaxes(arrayfun(@(c)(c.ax),this.Axes), 'x');
            end
        end
        
        function AddSignals(this,sigs,ax,itraces)
            if ischar(sigs)
                sigs = {sigs};
            end
            
            if ~exist('ax', 'var')||isempty(ax)
                ax = numel(this.Axes)+1;
            end
            
            if ~exist('itraces', 'var')
                itraces = this.ipts;
            elseif strcmp(itraces, 'all')
                itraces = 1:size(this.BrSet.P.pts,2);
            end
            
            if isnumeric(ax)
                if ax==0
                    this.AddAxes(1);                    
                    ax = this.Axes(1);                    
                elseif ax > numel(this.Axes)
                    this.AddAxes();                    
                    ax =this.Axes(end).ax;                    
                else                    
                    ax = this.Axes(ax).ax;                    
                end
            end
            if ~isa(ax, 'matlab.graphics.axis.Axes')
                error('Argument should be an Axes object or an index.')
            end
            
            for is = 1:numel(sigs)
                sig = sigs{is};
                this.plot_signal(sig, ax, itraces);
            end
            this.update_legend(ax);            
        end
                                   
        function AddRobSignals(this,sigs,ax,itraces,ireq)
            if ischar(sigs)
                sigs = {sigs};
            end
            
            if ~exist('ax', 'var')||isempty(ax)
                ax = numel(this.Axes)+1;
            end
            
            if ~exist('itraces', 'var')
                itraces = this.ipts;
            elseif strcmp(itraces, 'all')
                itraces = 1:size(this.BrSet.P.pts,2);
            end
            
            if ~exist('ireq', 'var')
                ireq = 1;
            elseif ischar(ireq)
                [~ , ireq] = this.BrSet.get_req_from_name(ireq);
            end            
            
            if isnumeric(ax)
                if ax==0
                    this.AddAxes(1);                    
                    ax = this.Axes(1);                    
                elseif ax > numel(this.Axes)
                    this.AddAxes();                    
                    ax =this.Axes(end).ax;                    
                else                    
                    ax = this.Axes(ax).ax;                    
                end
            end
            if ~isa(ax, 'matlab.graphics.axis.Axes')
                error('Argument should be an Axes object or an index.')
            end
            
            for is = 1:numel(sigs)
                sig = sigs{is};
                this.plot_rob(sig, ax, itraces,ireq);
            end
            this.update_legend(ax);            
        end                                
        
        function PlotDiagnostics(this, req)
            req.plot_full_diagnostics(this);
            for ia = 1:numel(this.Axes)
                this.update_legend(this.Axes(ia).ax);
            end            
        end
        
        function int_false= HighlightFalse(this, sig, ax,inv)
            if ~exist('ax', 'var')||isempty(ax)
                ax = this.Axes(end).ax;
            end
            if ~exist('inv', 'var')||isempty(inv)
                inv = false;
            end
            
            axes(ax);
            hold on;
            traj = this.BrSet.P.traj{this.ipts};
            idx = FindParam(this.BrSet.P, sig);
            tau = traj.time;
            val = traj.X(idx,:);
            if inv
                val(~isnan(val)) = ~val(~isnan(val));
            end
            int_false = highlight_truth_intervals(tau,val, 'g', 0, 'r', 0.3);
            set(ax,'UserData', int_false);
            this.update_legend(ax);
            
        end
              
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
                        if isequal(c(idx).FaceColor, [1 0 0])
                            status = 0;
                        else
                            status = 1; 
                        end
                    end
                    num_patch = num_patch+1;
                elseif isequal(c(idx).Marker, 'x')
                    lc(end+1) = c(idx);
                    if isequal(c(idx).Color, [1 0 0])
                        st{end+1} = 'Possible worst value';
                    else
                        st{end+1} = 'Possible best value';
                    end
                else
                    lst = get(c(idx), 'DisplayName');
                    if ~ismember(lst, st)&&~isempty(lst)
                        lc(end+1) = c(idx);
                        st{end+1} = lst;
                    end
                end
                    
            end
            
            if num_patch>0
                if status==0
                    set(lc(patch_idx),'DisplayName' ,['Violation Interval(s) (' num2str(num_patch) ')']);
                else
                    set(lc(patch_idx),'DisplayName' ,['Satisfying Interval(s) (' num2str(num_patch) ')']);
                end
                st{patch_idx} = get(lc(patch_idx),'DisplayName'); 
            end
            l= legend(lc,st);
            set(l, 'Interpreter', 'None');
        end
        
        function pos = get_pos_from_ax(this,ax)
            for pos = 1:numel(this.Axes)
                if isequal(this.Axes(pos).ax,ax)
                    return;
                end
            end
            if pos ==  numel(this.Axes)
                error('Axes not found!');
            end
        end
        
        function update_axes(this, ax)
            if nargin==1 
               Axes_to_update = this.Axes; 
            elseif isa(ax, 'matlab.graphics.axis.Axes')
                Axes_to_update = ax;
            elseif isnumeric(ax)
                Axes_to_update = this.Axes(ax).ax;
            end
            
            for A = Axes_to_update
                ax = A.ax;                
                this.plot_signal(A.signals,ax,this.ipts); 
                if isa(this.BrSet, 'BreachRequirement')
                    if ~isempty(A.robs)
                        for ireq=1:numel(this.BrSet.req_monitors)  % TODO: precond_monitors...
                            if ~isempty(A.robs{ireq})
                               this.plot_rob(A.robs{ireq},ax,this.ipts,ireq);
                            end
                        end
                    end
                end
            end                        
        end
        
        function set_ipts(this,ipts)
            this.ipts= ipts;
            this.update_axes();
            this.update_title();
        end

        
        function next_ipts(this, num)
            if nargin<2
                num = 1;
            end
            this.ipts = min(this.ipts+num, this.Summary.num_traces);
            this.update_axes();
            this.update_title();
        end
        
        function prev_ipts(this,num)
            if nargin<2
                num = 1;
            end
            this.ipts = max(1,this.ipts-num);
            this.update_axes();
            this.update_title();
        end        
        
        function update_title(this)
            ttle = [this.BrSet.whoamI ': '];
            if isa(this.BrSet, 'BreachRequirement')
                num_traces= this.Summary.num_traces_evaluated;
                ttle = [ttle 'Trace ' num2str(this.ipts) '/' num2str(num_traces) '. '];
                % checks violations
                
                if  this.Summary.num_violations_per_trace(this.ipts)>0
                    status_vio = ' Status: False';
                else
                    status_vio = ' Status: True';
                end
                                
                vio= find(this.Summary.num_violations_per_trace(this.ipts+1:end),1);
                
                if  ~isempty(vio)                   
                    status_vio = [status_vio '  (Next False idx:' num2str(vio+this.ipts) ').'];
                end
                ttle = [ttle status_vio];
                
            else
                num_traces=  this.Summary.num_traces;
                ttle = [ttle 'Trace ' num2str(this.ipts) '/' num2str(num_traces)];
            end
           set(this.Fig,'Name',ttle, 'NumberTitle', 'off');                 
        end
        
        
    end
    
    methods (Access=protected)
         
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
            
            if isa(this.BrSet,'BreachRequirement')
                trm = uimenu(cm, 'Label', 'Debug requirement','Separator', 'on');
                for ir = 1:numel(this.BrSet.req_monitors)
                    uimenu(trm, 'Label', this.BrSet.req_monitors{ir}.name,'Callback', @(o,e)ctxtfn_plot_full_diag(this.BrSet.req_monitors{ir},o,e));
                end
                rob_mnu = uimenu(cm, 'Label', 'Plot Robustness Signal');
                rob_mnu_map = containers.Map();
                
                for ir = 1:numel(this.BrSet.req_monitors)
                    req = this.BrSet.req_monitors{ir};
                    if isa(req, 'stl_monitor')
                        phi_name = req.name;
                        phi = req.formula;
                        cmu = uimenu(rob_mnu, 'Label',phi_name);
                        this.get_mnu(ax,cmu,phi,ir);
                                                                        
                    end
                end
            end
            
            uimenu(cm, 'Label', 'Add axes above','Separator', 'on', 'Callback', @(o,e)ctxtfn_add_axes_above(ax,o,e));
            uimenu(cm, 'Label', 'Add axes below', 'Callback', @(o,e)ctxtfn_add_axes_below(ax, o,e));
            uimenu(cm, 'Label', 'Reset axes','Separator', 'on', 'Callback', @(o,e)ctxtfn_reset_axes(ax, o,e));
            uimenu(cm, 'Label', 'Delete axes', 'Callback', @(o,e)ctxtfn_delete_axes(ax, o,e));
                                               
            function ctxtfn_add_axes_above(ax, ~,~)
                for ia = 1:numel(this.Axes)
                    if isequal(ax,this.Axes(ia).ax)
                        break;
                    end
                end
                this.AddAxes(ia);
            end
            
            function ctxtfn_add_axes_below(ax, ~,~)
                for ia = 1:numel(this.Axes)
                    if isequal(ax,this.Axes(ia).ax)
                        break;
                    end
                end
                this.AddAxes(ia+1);
            end
            
            function ctxtfn_delete_axes(ax, ~,~)
                for ia = 1:numel(this.Axes)
                    if isequal(ax,this.Axes(ia).ax)
                        break;
                    end
                end
                this.DeleteAxes(ia);
            end
            
            function ctxtfn_reset_axes(ax, ~,~)
                cla(ax);
                pos = this.get_pos_from_ax(ax);
                this.Axes(pos).signals = {};
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
            
            function ctxtfn_plot_full_diag(req, ~,~)                
                this.PlotDiagnostics(req);
            end             
                        
        end        
                        
        function get_mnu(this,ax,mnu_parent, phi,ir,o,e)
            
            if exist('o','var')
                mnu_parent = o;                
                set(mnu_parent, 'Callback', @(o,e)do_nothing());
                i0 = 2;
            else
                i0 = 1;
            end
            
            [st_phis, phi_ids, not_expanded] =tree_disp(phi,0,2);
            for i = i0:numel(st_phis)
                
                phi_str = strtrim(st_phis{i});
                phi_str = strrep(phi_str,'|-','');
                if isempty(not_expanded{i})
                    mm= uimenu(mnu_parent, 'Label',st_phis{i},...
                        'Callback', @(o,e)ctxtfn_plot_rob(phi_str,o,e));
                else
                    mm = uimenu(mnu_parent,'Label',st_phis{i}, ...
                        'Callback', @(o,e)(this.get_mnu(ax,mnu_parent, not_expanded{i},ir,o,e)));
                    uimenu(mm, 'Label',disp(not_expanded{i},2),...
                        'Callback', @(o,e)ctxtfn_plot_rob(disp(not_expanded{i}),o,e));                    
                end
            end
                        
            function ctxtfn_plot_rob(phi_str,~,~)                
                this.plot_rob(phi_str, ax,this.ipts, ir);
                this.update_legend(ax);
            end
            function do_nothing()
            end
        end

        function h = get_line_from_signal_name(this, ax, signame)
            h = [];
            ch = get(ax,'Children');            
            for idx_l =1:numel(ch)
               c = ch(idx_l);
                if strcmp(get(c,'DisplayName'), signame)
                   h = c;
               end
            end            
        end
        
        function l = plot_rob(this, formulas, ax, ipts, ireq)
            % 
            if isempty(formulas)
                l = [];
                return
            end   
            if ~iscell(formulas)
                formulas= {formulas};
            end 
            
            pos = this.get_pos_from_ax(ax);
            if isempty(this.Axes(pos).robs)              
                this.Axes(pos).robs=cell(1, numel(this.BrSet.req_monitors));
            end
                
            this.Axes(pos).robs{ireq} = union(this.Axes(pos).robs{ireq}, formulas);
            
            if ~exist('ipts', 'var')
                ipts= this.ipts;
            end
            
            axes(ax);
            hold on;
            itraj = unique(this.BrSet.P.traj_ref(ipts), 'stable');
            
            time = this.BrSet.P.traj{itraj}.time;
            ch = get(ax,'Children');            
            if isempty(ch)
               set(ax, 'XLimMode', 'auto', 'YLimMode', 'auto');
            end
                                    
            for f = formulas 
                l = this.get_line_from_signal_name(ax, f{1}); % FIXME: here we 
                [t,r]= this.BrSet.GetRobustSat(ireq, itraj, f{1});
                if isempty(l)
                    l = plot(t, r, 'DisplayName', f{1});
                else
                    set(l,'XData',t,'YData',r);
                end
                if isempty(this.get_line_from_signal_name(ax,this.zero_rob_line_name))
                    plot(time, 0*time, 'DisplayName', this.zero_rob_line_name, 'LineStyle', '--', 'Color','r');
                end
            end
            
            
        end
        
        function l = plot_signal(this, sig, ax, ipts)
            if isempty(sig)
                l = [];
                return
            end
                
            pos = this.get_pos_from_ax(ax);
                                    
            this.Axes(pos).signals = union(this.Axes(pos).signals, sig);            
            if ~exist('ipts', 'var')
                ipts= this.ipts;
            end
            
            axes(ax);
            hold on;
            itraj = unique(this.BrSet.P.traj_ref(ipts), 'stable');
            
            if ~iscell(sig)
                sig= {sig};
            end 
            time = this.BrSet.P.traj{itraj}.time;
            ch = get(ax,'Children');            
            if isempty(ch)
               set(ax, 'XLimMode', 'auto', 'YLimMode', 'auto');
            end
                                    
            for s = sig                
                % find out if it's a robustness signal or not
                if isa(this.BrSet,'BreachRequirement')
                    [~, b, ~, bb]= this.BrSet.FindSignalsIdx(s{1});
                else
                    b=true;
                end
                if b||bb                    
                    sig_values = this.BrSet.GetSignalValues(s{1}, itraj);
                    l = this.get_line_from_signal_name(ax, s{1});
                    if isempty(l)
                        l = plot(time , sig_values, 'DisplayName', s{1});
                    else
                        set(l, 'XData',time,'YData',sig_values);
                    end
                    
                else % try robustnes... should be obsolete
                    warning('Not sure what I am doing here, plot_signal but trying robustness signal ?');
                    l = this.get_line_from_signal_name(ax, s{1});
                    [t,r]= this.BrSet.GetRobustSat(1, itraj, s{1});
                    if isempty(l)                   
                        l = plot(t, r, 'DisplayName', s{1});                                                
                    else
                        set(l,'XData',t,'YData',r);
                    end                    
                    if isempty(this.get_line_from_signal_name(ax,this.zero_rob_line_name))
                        plot(time, 0*time, 'DisplayName', this.zero_rob_line_name, 'LineStyle', '--', 'Color','r');
                    end
                    
                end
            end
            
        end
              
    end
end