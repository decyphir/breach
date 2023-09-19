classdef BreachProblemGui < BreachGuiClass
    properties
        pb
        pause_requested
    end
    
    
    methods
        function this =  BreachProblemGui(pb)
            
            this.pb=pb;
            this.pb.callback_obj = @(o,e)(this.zcallback('pb',o,e));
            this.create_axes('ax_main',2,8);
            this.create_button('button_stop','Stop');
            this.create_button('button_resume','Start/Resume');
            layout = {{'ax_main'};
                {'button_resume', 'button_stop'}};
            this.set_layout(layout);
            this.enable_resizable();
            
        end
        
        function zcallback(this,id,o,e)
            
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
                    grid on;
                    plot(1:this.pb.nb_obj_eval, this.pb.obj_log); 
                    drawnow
                    if this.pause_requested
                        this.pb.is_paused= 1;
                    end
            end
            
        end
        
    end
    
end