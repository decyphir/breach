classdef gui_elem < handle
    properties
        hdle
        
        h = 1  % exterior size
        w = 1
         
        wleft =.05  % internal margins
        wright =.05
        htop  =.05
        hbot  = .05
    end
    
    methods
                
        function disp_in(this, gui, elem_parent, w, h, wunit, hunit)
            
            he = this.h;
            we = this.w;
            wpos = w;
            hpos = h;
            wsz = we*wunit;
            hsz = he*hunit;
            set(this.hdle,'Parent',elem_parent.hdle, 'Position', [wpos hpos wsz hsz], 'visible', 'on');
                    
        end               
        
    end
    
    
    
end

