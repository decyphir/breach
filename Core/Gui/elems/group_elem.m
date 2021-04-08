classdef group_elem < gui_elem
     
    properties
       layout       
       margin_top = .02
       margin_bot = .025
       
    end
    
    methods
        
        function this = group_elem(layout)
            this.layout = layout;
            this.wleft =.0;  % no internal margins for groups
            this.wright =.0;
            this.htop  =.0;
            this.hbot  = .0;
        end
        
                
        function disp_in(this, gui, elem_parent, w, h, wunit, hunit)
                                                
            gui.set_layout(this.layout, elem_parent,w,h, 'absolute');
            
        end
    end
end

