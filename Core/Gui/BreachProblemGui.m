classdef BreachProblemGui < BreachGuiClass
% A class to run, configure, etc, a BreachProblem
    
    properties
        pb
        pause_requested
        Bsearch
    end
    
    methods
        function this =  BreachProblemGui(pb)
            
            set(this.hdle, 'Name', 'BreachProblem GUI'),
            
      
            %% Linking problem object
            this.pb=pb;
            this.pb.callback_obj = @(o,e)(this.zcallback('pb',o,e));
            
            params = this.pb.params;
            ranges = [this.pb.lb this.pb.ub];
            
            B = BreachSet(params);
            B.SetParamRanges(params, ranges);
            B.SetParam(params,this.pb.x0);
            
            this.Bsearch = B; 

            %% Creates all elements
            
            % objective ax
            this.create_axes('ax_main',2,8);
            this.create_panel('panel_plots', 'Objective', {{'ax_main'}});
            
            % param table
            this.button_reset_domain('create');
            this.button_apply_domain('create');
            this.table_params('create');
            this.checkbox_toggle_edit('create');
            this.create_panel('panel_params', 'Search Domain', ...
                {{'checkbox_toggle_edit','button_reset_domain','button_apply_domain'}; ...
                 {'table_params'}});
        
            % buttons
            this.button_stop('create');
            this.button_resume('create');
            
            this.create_text('text_obj_eval', 'Max number of evaluations',.5,1)
            this.edit_max_obj_eval('create');
            this.create_group('group_obj_eval', {{'text_obj_eval', 'edit_max_obj_eval'}});
            this.create_panel('panel_buttons', 'Controls', {{'group_obj_eval'}; ...
                                                            {'button_resume', 'button_stop'}});
            
            % Layout
            layout = {{'panel_plots'};
            {'panel_params'};
            {'panel_buttons'}};

            this.set_layout(layout);
            this.enable_resizable();
            this.zcallback('pb');
        end
        
        function zcallback(this,id,o,e)
        % callback function, called by the problem at each freq_update    
            switch id
                case 'pb'
                    % update number of obj evaluations
                    panel_plots_title = '  Objective   ';
                    obj_eval = num2str(size(this.pb.obj_log,2));
                    max_eval = num2str(this.pb.max_obj_eval);
                    panel_plots_title = [panel_plots_title obj_eval ' / ' max_eval];  
                    this.set_by_id('panel_plots','Title', panel_plots_title);
                    
                    ax=this.uimap('ax_main').hdle;                    
                    axes(ax);
                    enableDefaultInteractivity(ax);
                    grid on;
                    plot(1:this.pb.nb_obj_eval, this.pb.obj_log); 
                    drawnow
                    if this.pause_requested
                        this.pb.is_paused= 1;
                    end
            end
            
        end
     
        %%  Domain panel         
        function button_reset_domain(this, mode,o, e)
            w= .7;
            h= 1;
            cb = @(o,e)(this.button_reset_domain('callback', o, e));
            switch mode
                case 'create'
                    e = this.create_button('button_reset_domain','Reset Search Domain ', cb, w,h);
                case 'callback'
                    this.table_params('reset'); 
            end
        end

       function button_apply_domain(this, mode,o, e)
            w= .7;
            h= 1;
            cb = @(o,e)(this.button_apply_domain('callback', o, e));
            switch mode
                case 'create'
                    e = this.create_button('button_apply_domain','Apply Changes', cb, w,h);
                case 'callback'
                    this.table_params('apply'); 
            end
        end
        
       function checkbox_toggle_edit(this, mode,o, e)
            w= .6;
            h= 1;
            cb = @(o,e)(this.checkbox_toggle_edit('callback', o, e));
            switch mode
                case 'create'
                    e = this.create_checkbox('checkbox_toggle_edit','Edit Search Domain ', cb, w,h);
                    cb([],[]);
                case 'callback'
                    v = this.get_by_id('checkbox_toggle_edit','Value');
                    if v
                        this.set_by_id('button_apply_domain', 'enable', 'on');
                        this.set_by_id('button_reset_domain', 'enable', 'on');
                        this.set_by_id('table_params','ColumnEditable', true);
                    else
                        this.set_by_id('button_apply_domain', 'enable', 'off');
                        this.set_by_id('button_reset_domain', 'enable', 'off');
                        this.set_by_id('table_params','ColumnEditable', false);
                    end
            end
        end

       function table_params(this, mode,o, e)
            id = 'table_params';
            w = 2;
            h = 4;
            switch mode
                case 'create'
                    cb =@(o,e)(this.table_params('callback',o,e));
                    elem_tble = this.create_table(id, {},cb, w, h);
                    htable = elem_tble.hdle;
                    elem_tble.hdle =  fill_uitable_params(htable,this.Bsearch);
                    this.data_gui.(id).prev_data = get(elem_tble.hdle,'data');
                case 'callback'
                    % sanity test of input
                    idx= e.Indices;
                    data = get(o,'data');
                    try
                       read_uitable_params(h_uitable);
                    catch ME
                       warning('BreachProblemGUI:InvalidEntry','Invalid data: %s', data{idx(1),idx(2)});
                       set(o, 'data',this.data_gui.(id).prev_data)
                    end
                case 'reset'
                    elem_tble = this.uimap('table_params');
                    htable = elem_tble.hdle;
                    elem_tble.hdle =  fill_uitable_params(htable,this.Bsearch);                
                case 'apply'
                    elem_tble = this.uimap('table_params');
                    htable = elem_tble.hdle;
                    [params, p0, domains] = read_uitable_params(htable);
                    this.pb.BrSet.SetParam(params, p0);
                    this.pb.BrSet.SetDomain(params, domains);
                    this.pb.ResetObjective();
            end
            
        end
        
       %% Panel Control
       function button_resume(this, mode,o, e)
            w= 1;
            h= 1;
            cb = @(o,e)(this.button_resume('callback', o, e));
            switch mode
                case 'create'
                    e = this.create_button('button_resume','Start/Resume', cb, w,h);
                case 'callback'
                    this.pause_requested= 0;
                    this.pb.is_paused = 0; 
                    this.pb.solve();
            end
        end

        function button_stop(this, mode,o, e)
            w= 1;
            h= 1;
            cb = @(o,e)(this.button_stop('callback', o, e));
            switch mode
                case 'create'
                    e = this.create_button('button_stop','Stop', cb, w,h);
                case 'callback'
                    this.pause_requested = 1;                                        
            end
        end
        
        function edit_max_obj_eval(this, mode,o, e)
            id = 'edit_max_obj_eval';
            w= .5;
            h= 1;
            cb = @(o,e)(this.edit_max_obj_eval('callback', o, e));
            string = '';
            switch mode
                case 'create'
                    this.create_edit('edit_max_obj_eval',string, cb, w,h);
                    this.edit_max_obj_eval('update');
                case 'callback'
                    max_eval = str2double(this.get_by_id(id, 'String'));
                    this.pb.max_obj_eval = max_eval;
                case 'update'
                    max_eval = num2str(this.pb.max_obj_eval);
                    this.set_by_id(id,'String',max_eval);
            end
        end


       %% Plotting 
     % TODO more options, axes, etc
     %  function reset_obj_plot(this)
     %
     %  end
        function reset_top_fig(this)
            % we enable figure toolbar since we have axes
            this.hdle = figure('Name',this.title, ...                
                'Toolbar', 'figure',....
                'Menubar', 'none',...
                'NumberTitle','off',...
                'SizeChangedFcn', @(o,e)(this.resize_callback(o,e)));
        end
    end

     

        

    
end