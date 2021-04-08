classdef panel_elem < gui_elem
     
    properties
       layout       
       margin_top = .02
       margin_bot = .025
       
    end
    
    methods
        
        function this = panel_elem(hdle, layout)
            this.hdle = hdle;
            this.layout = layout;
            this.wleft =.025;  % default internal margins smaller for panels
            this.wright =.025;
            this.htop  =.025;
            this.hbot  = .025;
        end
        
                
        function disp_in(this, gui, elem_parent, w, h, wunit, hunit)
            
            we = this.w;
            he = this.h;
                                    
            wpos = w+this.wleft*wunit;
            hpos = h+this.hbot*hunit;
            wsz = (we-this.wright-this.wleft)*wunit;
            hsz = (he-this.htop-this.hbot)*hunit;
            
            set(this.hdle,'Parent',elem_parent.hdle, 'Position', [wpos hpos wsz hsz], 'Visible', 'on', 'FontWeight', 'bold');
            gui.set_layout(this.layout, this,-this.wleft*wunit,-this.hbot*hunit, 'absolute');
            
        end
    end
end

