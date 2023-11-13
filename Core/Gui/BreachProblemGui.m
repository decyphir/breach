classdef BreachProblemGui < BreachGuiClass
% A class to run, configure, etc, a BreachProblem
    
    properties
        pb
        pause_requested
        Bsearch
    end
    
    methods
        function this =  BreachProblemGui(pb)
            
            set(this.hdle, 'Name', 'BreachProblemGUI'),

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
            id = 'table_params';
            cb =@(o,e)(this.table_params('callback',o,e));      
            this.create_table(id, {},cb, 2, 4);
            this.create_panel('panel_params', 'Search Domain', {{'table_params'}});
            this.table_params('init');
        
            % buttons
            this.create_button('button_stop','Stop');
            this.create_button('button_resume','Start/Resume');
            this.create_panel('panel_buttons', 'Controls', {{'button_resume', 'button_stop'}});

            layout = {{'panel_plots'};
            {'panel_params'};
            {'panel_buttons'}};

            this.set_layout(layout);
            this.enable_resizable();
            
        end
        
        function zcallback(this,id,o,e)
        % callback function, called by the problem at each freq_update    
            switch id
                case 'button_stop'
                    this.pause_requested = 1;                                        
                case 'button_resume'
                    this.pause_requested= 0;
                    this.pb.is_paused = 0; 
                    this.pb.solve();
                case 'pb'
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
     
        %%  Domain table
        function table_params(this, mode,o, e)
            switch mode
                case 'init'
                    e = this.uimap('table_params');
                    htable = e.hdle;
                    e.hdle =  fill_uitable_params(htable,this.Bsearch);
                    set(e.hdle,'ColumnEditable', false);
                case 'callback'
                    % Todo
            end
        end
        
       %% Plotting 
     % TODO more options, axes, etc
     %  function reset_obj_plot(this)
     %
     %  end

    end

     

        

    
end