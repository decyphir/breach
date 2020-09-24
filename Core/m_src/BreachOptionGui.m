classdef BreachOptionGui < handle
    % Draw N buttons + two at the bottom relative to a top figure
    
    properties
        dlg
        dlg_pos
        dlg_sz
        room_increment = 45
        button_pos
        button_sz
        gui_components
        font_name = 'Arial'
        font_size = 10
        max_char = 80
        
        output
        choices
        tips
        title
        button_ok_
        button_cancel_
    end
    
    properties
       error_msg 
    end
    
    methods
        function this = BreachOptionGui(title, options, choices, tips,varargin)
           if ~exist('title', 'var')||isempty(title)
               title = 'Choose Options';
           end
            if nargin>1
                options = varargin2struct_breach(options, varargin{:});
                this.output = options;
                this.choices = choices;
                this.tips = tips;
            else
                this.output = struct;
                this.choices = struct;
                this.tips = struct;
            end
            this.title = title;
            reset_fig(this);
            
        end
        
        function add_options(this,name,value, choice, tip)
            
            if ischar(name)&&~iscell(name)
                name = {name};
                value = {value};
                choice = {choice};
                tip = {tip};
            end
            
            if iscell(name)
                for io = 1:numel(name)
                    this.output.(name{io}) = value{io};
                    this.choices.(name{io}) = choice{io};
                    this.tips.(name{io}) = tip{io};
                end
                
            end
            this.reset_fig();
        end
        
        function merge_options(this, opt_other, choices_other, tips_other)
            fn= fieldnames(opt_other);
            for ifn = 1:numel(fn)
                name = fn{ifn};
                this.output.(name) = opt_other.(name); 
                this.choices.(name) = choices_other.(name); 
                this.tips.(name) = tips_other.(name); 
            end
            this.reset_fig();
        end
                
        function reset_fig(this)
            
            this.gui_components = {};
            close(this.dlg);
            this.dlg = figure('Name',this.title, ...
                'Resize', 'on',...
                'Toolbar', 'none',....
                'Menubar', 'none',...
                'NumberTitle','off');
            pos = get(this.dlg, 'Position');
            this.dlg_pos = pos(1:2);
            this.dlg_sz = [600 0];
            set(this.dlg, 'Position', [this.dlg_pos this.dlg_sz])
            this.add_ok_cancel_buttons();
            
            flds = fieldnames(this.output);
            for ii = numel(flds):-1:1
                opt_name = flds{ii};
                opt_value = this.output.(opt_name);
                opt_choices = this.choices.(opt_name);
                opt_tip = this.tips.(opt_name);
                this.add_uicontrol(opt_name,opt_value, opt_choices, opt_tip);
            end
            
        end
        
        function add_uicontrol(this, name, val, choices, tip)

            if ischar(choices)
                switch choices
                    case 'bool'
                        g = this.add_check_box(name, tip);
                        set(g, 'Value', val);
                        set(g, 'Callback', @check_callback);
                    case 'int'
                        g = this.add_edit_box(num2str(val), tip);
                        txt = sprintf('Enter a value for %s (int)' ,  name);
                        this.add_text_box(txt);
                        set(g, 'String', num2str(val));
                        set(g,'Callback', @edit_num_callback);
                    case 'double'
                        g = this.add_edit_box(num2str(val),tip);
                        txt = sprintf('Enter a value for %s (double)' ,  name);
                        this.add_text_box(txt);
                        set(g, 'String', num2str(val));
                        set(g,'Callback', @edit_num_callback);
                    case 'string'
                        g = this.add_edit_box(val,tip);
                        txt = sprintf('Enter a value for %s (string)' ,  name);
                        this.add_text_box(txt);
                        set(g, 'String', val);
                        set(g,'Callback', @edit_string_callback);
                end
                
            elseif iscell(choices)
                   g = this.add_edit_box(num2str(val),tip);
                   txt = sprintf('Choose a value for %s' ,  name);
                   this.add_text_box(txt);
                   set(g,'Style', 'popupmenu','String', choices);
                   val_idx = find(strcmp(choices, val));
                   set(g, 'Value', val_idx);
                   set(g,'Callback', @menu_callback);
                   
            end        
   
            function menu_callback(o, e)
                val_idx = get(o, 'Value');
                val_new  = choices{val_idx};
                this.output.(name)=val_new;
            end
            
            function edit_num_callback(o, e)
                val_st = get(o, 'String');
                val_new  = str2num(val_st);
                this.output.(name) =val_new;
            end

            function edit_string_callback(o, e)
                val_st = get(o, 'String');
                this.output.(name) =val_st;
            end

            function check_callback(o,e)
                val_new = get(o , 'Value');
                this.output.(name) = val_new; 
            end
            
        end
              
        function g = add_button(this, string, callback,tooltip)
            this.make_room();
            pos = [this.button_pos this.button_sz];
            num_elements = size(this.gui_components,1);
            g = this.create_button(string, callback, pos);
            set(g,'TooltipString', tooltip)
            this.gui_components{num_elements+1,1} = g;
        end
        
        function g=  create_button(this,  string, callback, pos)
  
            g=  uicontrol('Parent',this.dlg,...
                'Style','pushbutton',...
                'String', string,...
                'Units', 'Normalized',...
                'FontName', this.font_name,...
                'FontSize', this.font_size,...
                'Position', pos,...
                'Callback', callback ...
                );
            
        end
        
        function g = add_check_box(this, string, tooltip)
            this.make_room();
            g = this.create_check_box([ 'Enable ' string ' option' ], [this.button_pos this.button_sz]);
            set(g, 'TooltipString', tooltip);
            this.gui_components{end+1,1} = g;
        end
        
        function g  = create_check_box(this, string, pos)
            
            g = uicontrol('Parent',this.dlg,...
                'Style','checkbox',...
                'String', string ,...
                'FontName', this.font_name,...
                'FontSize', this.font_size,...
                'Units', 'Normalized',...
                'Position', pos);
        end

        function g  = create_radio(this, string, pos)
            
            g = uicontrol('Parent',this.dlg,...
                'Style','radiobutton',...
                'String', string ,...
                'FontName', this.font_name,...
                'FontSize', this.font_size,...
                'Units', 'Normalized',...
                'Position', pos);
        end

        
        function g = add_text_box(this, string)
            this.make_room();
            pos = [this.button_pos this.button_sz];
            num_elements = size(this.gui_components,1);
            ax = axes;
            set(ax, 'Position', pos, 'visible', 'off');
            g = text(0.0,0.3, string,'Interpreter','none', 'FontWeight', 'bold');
            set(g,'FontName', this.font_name,...
                'FontSize', this.font_size);
            this.gui_components{num_elements+1,1} = ....
                ax;
        end
        
        function g =add_edit_box(this, string, tooltip)
            this.make_room();
            pos = [this.button_pos this.button_sz];
            num_elements = size(this.gui_components,1);
            g = this.create_edit_box(pos);
            set(g, 'String', string);
            set(g, 'TooltipString', tooltip);
            this.gui_components{num_elements+1,1} = g;
        end
        
        function g = create_edit_box(this,  pos) 
            g=  uicontrol('Parent',this.dlg,...
                'Style','edit',...
                'Units', 'Normalized',...
                'FontName', this.font_name,...
                'FontSize', this.font_size,...
                'Position', pos);
        end 
                
        function add_ok_cancel_buttons(this)
            this.make_room();
            sz = this.get_multiple_sz(2);
            num_elements = size(this.gui_components,1);
            this.button_cancel_ = ....
                uicontrol('Parent',this.dlg,...
                'Style','pushbutton',...
                'Units', 'Normalized',...
                'String', 'Cancel',...
                'FontName', this.font_name,...
                'FontSize', this.font_size,...
                'Callback', @cancel_callback,...
                'Position', sz(1,:));
            this.gui_components{num_elements+1,1} = this.button_cancel_;
            this.button_ok_  = ....
                uicontrol('Parent',this.dlg,...
                'Style','pushbutton',...
                'Units', 'Normalized',...
                'String', 'Ok',...
                'FontName', this.font_name,...
                'FontSize', this.font_size,...
                'Callback', @ok_callback,...
                'Position', sz(2,:));
            this.gui_components{num_elements+1,2} = this.button_ok_;
           
            function  ok_callback(hobj, evt)
                set(hobj, 'String', 'Please wait...')
                drawnow;
                close(this.dlg);
            end
    
            
            function  cancel_callback(hobj, evt)
                this.output=[];
                close(this.dlg);
            end
            
        end
        
        
        %% Auxiliary functions
             
        function txt =truncate_txt(this, txt)
            if numel(txt)>this.max_char
                txt = [txt(1:this.max_char) '...'];
            end
        end
 
               
        function make_room(this)
            d = this.room_increment;
            n = size(this.gui_components, 1)+1;
            this.dlg_sz = this.dlg_sz+ [0 d];
            this.dlg_pos = this.dlg_pos - [0 d/2];
            set(this.dlg, 'Position', [this.dlg_pos this.dlg_sz])
            
            total_ratio = 0.9+0.1*(n-1)/(n+10);
            one_minus_tr = 1-total_ratio;
            
            if n>1
                new_inter = one_minus_tr/(n-1);
                new_height =  (total_ratio-(n-1)*new_inter)/n;
                this.button_sz =   [ 0.9  new_height];
                this.button_pos = [ 0.05  one_minus_tr/2+(n-1)*(new_height+new_inter)];
                
                % Resize and pos other elements
                for ic = 1:size(this.gui_components, 1)
                    for irow = 1:size(this.gui_components, 2)
                        g = this.gui_components{ic,irow};
                        if ~isempty(g)
                            pos = get(g,'Position');
                            pos(2) = one_minus_tr/2+(ic-1)*(new_height+new_inter);
                            pos(4) = this.button_sz(2);
                            set(g,'Position', pos);
                        end
                    end
                end
            else
                this.button_sz =   [ 0.90   0.90];
                this.button_pos = [ 0.05  0.05];
            end
            
        end
        
        function sz = get_multiple_sz(this,nb)
            % get size button to accomodate for nb button in the current row
            sz = [this.button_pos this.button_sz];
            sz(1,3) = (1/nb)*(sz(1, 3)-(nb-1)*sz(1,1)/nb);   % divide x size
            for i =2:nb
                sz(i,:) = sz(i-1,:);
                sz(i, 1) = sz(i-1,1)+sz(1,1)/nb+sz(1,3);
            end
        end

        function pause_gui(this)
            % pause_gui disable all that can be disabled
            for ie = 1:numel(this.gui_components)
                g = this.gui_components{ie};
                try
                    set(g, 'enable', 'off')
                end
            end
        end
        
        function unpause_gui(this)
            % pause_gui enable all that can be enabled 
            for ie = 1:numel(this.gui_components)
                g = this.gui_components{ie};
                try
                    set(g, 'enable', 'on')
                end
            end
            figure(this.dlg);
        end
        
        
    end
end

