classdef popup_elem < gui_elem
    
    methods       
        
        
        function disp_in(this, gui, elem_parent, w, h, wunit, hunit)
        % draw checkbox with margins
            
            we = this.w;
            he = this.h;
                    
            wpos = w+this.wleft*wunit;
            hpos = h+this.hbot*hunit;
            wsz = (we-this.wright-this.wleft)*wunit;
            hsz = (he-this.htop-this.hbot)*hunit;

            set(this.hdle,'Parent',elem_parent.hdle, 'Position', [wpos hpos wsz hsz], 'visible', 'on');
                    
        end
    end
end

