classdef BreachGuiClass < handle
    
    %% Functional
    properties
        title
        uimap
        output=1
        error_msg
    end
    
    %% Default data field
    properties
        data_gui
    end
    
    
    %% statics helpers
    methods (Static)
        function [value, cfg] = get_struct_field(cfg, field, value)
            % get the value of a field, set it if it does not exist
            if isfield(cfg,field)
                value = cfg.(field);
            else
                cfg.(field)= value;
            end
        end
    end
    
    methods
        %% Constructors
        function this = BreachGuiClass(title)
            this.uimap = containers.Map();
            
            if ~exist('title', 'var')||isempty(title)
                title = 'GuiClass';
            end
            this.title = title;
            
            %% Creates separators
            hs1= this.create_separator('hsep.05');
            hs1.h = 0.05;
            
            hs2= this.create_separator('hsep.1');
            hs2.h = 0.1;
            
            hs2bis= this.create_separator('hsep.2');
            hs2bis.h = 0.2;
            
            hs25 = this.create_separator('hsep.25');
            hs25.h = 0.25;
            
            hs3 = this.create_separator('hsep.3');
            hs3.h = 0.3;
            
            hs4 = this.create_separator('hsep.4');
            hs4.h = 0.4;
                        
            hs5= this.create_separator('hsep.5');
            hs5.h = 0.5;
            
            ws1= this.create_separator('wsep.05');
            ws1.w = 0.05;
            
            ws2= this.create_separator('wsep.1');
            ws2.w = 0.1;
            
            this.create_separator('wsep.5',  .5,1);
            
            this.create_separator('fullsep', 1, 1);
          %% create scaling invisible button
            button_scale = this.create_button('button_scale');
            this.set_by_id('button_scale', 'enable', 'off', 'Position', [0 0 this.wunit this.hunit]);
            
            
          %% Creates default ok and cancel buttons
            
            button_ok = this.create_button('button_ok','Ok',@(o,e)(this.button_ok_callback(o)));
            button_cancel = this.create_button('button_cancel','Cancel',@(o,e)(this.button_cancel_callback()));
            
            button_ok.w = 0.5;
            button_ok.wleft = 0.025;
            button_cancel.w = 0.5;
            button_cancel.wright = 0.025;
            
            this.create_group('ok_group', {{'button_cancel', 'button_ok'}; {'hsep.25'}});
            this.create_group('group_ok', {{'button_cancel', 'button_ok'}; {'hsep.25'}});
            
            this.reset_top_fig();
            this.disable_resizable();
            %           this.layout = {{'ok_group'}};
            %           this.set_layout();
            
        end
        
        function this = zcallback(this, id, o,e)
            % Default/main callback function - does nothing for now
            disp(['event from elem:' id])
            disp(e)
            disp(o)
        end
        
        function e = create_ui(this,id,g)
            e = gui_elem();
            e.hdle = g;
            this.uimap(id)= e;
        end
        
        function e = create_separator(this, id,w,h)
            if nargin<4
                h=1;
            end
            if nargin<3
                w=1;
            end
            e = gui_elem();
            e.w = w;
            e.h = h;
            this.uimap(id)= e;
            
        end
        
        function e = create_axes(this,id, w,h)
                        
            if nargin<4
                h=4;
            end
            if nargin<3
                w=1;
            end
            
            e = axes_elem();
            e.w = w;
            e.h = h;
            
            e.wleft =.1;  % internal margins
            e.wright =.06;
            e.htop =.5;
            e.hbot = .5;
            
            e.hdle = axes();            
            this.uimap(id)= e;
        end
        
        function e = create_button(this, id,  string, callback, w,h)
            
            if ~exist('callback','var')||isempty(callback)
                callback = @(o,e)(this.zcallback(id,o,e));
            end
            
            if nargin < 3
                string = [id ': Push Button'];
            end
            
            g = uicontrol('Parent',this.hdle,...
                'Style','pushbutton',...
                'String', string,...
                'FontName', this.font_name,...
                'FontSize', this.font_size,...
                'visible', 'off',...
                'Callback', callback ...
                );
            
            e = button_elem();
            if nargin<5
                w =1;
            end
            if nargin<6
                h =1;
            end
            e.w=w;
            e.h=h;
            
            e.hdle = g;
            this.uimap(id)= e;
            
        end
                        
        function e = create_slider(this, id, string, callback, w,h)
            if ~exist('callback','var')||isempty(callback)
                callback = @(o,e)(this.zcallback(id,o,e));
            end
            
            if nargin < 3
                string = [id ': Push Button'];
            end
            
            g = uicontrol('Parent',this.hdle,...
                'Style','slider',...
                'String', string,...
                'FontName', this.font_name,...
                'FontSize', this.font_size,...
                'visible', 'off',...
                'Callback', callback ...
                );
            
            e = slider_elem();
            if nargin<5
                w =1;
            end
            if nargin<6
                h =1;
            end
            e.w=w;
            e.h=h;
            
            e.hdle = g;
            this.uimap(id)= e;
            
        end
        
        function e = create_checkbox(this, id, string, callback, w,h)
                        
            if nargin < 3
                string = [id ': Check Box'];
            end
            
            if ~exist('callback','var')||isempty(callback)
                callback = @(o,e)(this.zcallback(id,o,e));
            end
            
            
            g = uicontrol('Parent',this.hdle,...
                'Style','checkbox',...
                'String', string ,...
                'FontName', this.font_name,...
                'FontSize', this.font_size,...
                'Callback', callback,...
                'visible', 'off');
            
            e = checkbox_elem();
            if nargin<5
                w =1;
            end
            if nargin<6
                h =1;
            end
            e.w=w;
            e.h=h;
            e.hdle = g;
            this.uimap(id)= e;
            
        end
        
        function e = create_radio(this,id, string,callback, w,h)
            
            if nargin < 3
                string = [id ': Check Box'];
            end
            
            if ~exist('callback','var')||isempty(callback)
                callback = @(o,e)(this.zcallback(id,o,e));
            end
            
            
            g = uicontrol('Parent',this.hdle,...
                'Style','radiobutton',...
                'String', string ,...
                'FontName', this.font_name,...
                'FontSize', this.font_size,...
                'Callback', callback,...
                'visible', 'off');
            this.create_ui(id,g);
            
            e = radio_elem();
            if nargin<5
                w =1;
            end
            if nargin<6
                h =1;
            end
            e.w=w;
            e.h=h;
            e.hdle = g;
            this.uimap(id)= e;
            
        end
        
        function e = create_popup(this,id,string,callback, w,h)
            
            if nargin < 3
                string = [id ': Check Box'];
            end
            
            
            if ~exist('callback','var')||isempty(callback)
                callback = @(o,e)(this.zcallback(id,o,e));
            end
            
            
            g = uicontrol('Parent',this.hdle,...
                'Style','popupmenu',...
                'String', string ,...
                'FontName', this.font_name,...
                'FontSize', this.font_size,...
                'Callback', callback,...
                'visible', 'off');
            this.create_ui(id,g);
            
            e = popup_elem();
            if nargin<5
                w =1;
            end
            if nargin<6
                h =1;
            end
            e.w=w;
            e.h=h;
            e.htop = .15;
            
            e.hdle = g;
            this.uimap(id)= e;
            
        end
        
        function e = create_listbox(this,id,string,callback, w,h)
            
            if nargin < 3
                string = {};
            end
            
            if ~exist('callback','var')||isempty(callback)
                callback = @(o,e)(this.zcallback(id,o,e));
            end
                                    
            
            g = uicontrol('Parent',this.hdle,...
                'Style','listbox',...
                'String', string ,...
                'FontName', this.font_name,...
                'FontSize', this.font_size,...
                'Callback', callback,...
                'visible', 'off');
            this.create_ui(id,g);
            
            e = listbox_elem();
            if nargin<5
                w =1;
            end
            if nargin<6
                h = 2;
            end
            e.w=w;
            e.h=h;
            e.hdle = g;
            this.uimap(id)= e;
            
        end
        
        function e = create_table(this,id, data, callback, w,h)
            
            if nargin<3
                data = {'',[]};
            end
            
            if ~exist('callback','var')||isempty(callback)
                callback = @(o,e)(this.zcallback(id,o,e));
            end
            
            
            g =  uitable('Parent',this.hdle,...
                'FontName', this.font_name,...
                'FontSize', this.font_size,...
                'Data',data,...
                'visible', 'off');
            if ~isempty(callback)
                set(g,'CellEditCallback', callback);
            end
            
            e = table_elem();
            if nargin<5
                w =1;
            end
            if nargin<6
                h =1;
            end
            e.w=w;
            e.h=h;
            e.hdle = g;
            this.uimap(id)= e;
            
        end
        
        function e = create_edit(this,id, string, callback, w,h)
            
            if nargin<3
                string = '';
            end
            
            if ~exist('callback','var')||isempty(callback)
                callback = @(o,e)(this.zcallback(id,o,e));
            end
            
            g =  uicontrol('Parent',this.hdle,...
                'Style','edit',...
                'FontName', this.font_name,...
                'FontSize', this.font_size,...
                'String',string,...
                'visible', 'off');
            if ~isempty(callback)
                set(g,'Callback', callback);
            end
            
            e = edit_elem();
            if nargin<5
                w =1;
            end
            if nargin<6
                h =1;
            end
            e.w=w;
            e.h=h;
            e.hdle = g;
            this.uimap(id)= e;
            
        end
        
        function e = create_text(this,id, string, w,h)
            
            if nargin<3
                string = '';
            end
            
            g =  uicontrol('Parent',this.hdle,...
                'Style','text',...
                'FontName', this.font_name,...
                'FontSize', this.font_size,...
                'enable', 'inactive', ...
                'String',string,...
                'HorizontalAlignment', 'left',...
                'visible', 'off'...
                );
            %                'BackgroundColor', [0.831; 0.816;0.784]...
            
            e = text_elem();
            if nargin<4
                w =1;
            end
            if nargin<5
                h =1;
            end
            e.w=w;
            e.h=h;
            e.htop = .15;
            e.hdle = g;
            this.uimap(id)= e;
            
        end
          
        function e = create_panel(this, id, string, layout)
            
            if nargin<3
                string = [id ' panel'];
            end
            
            p = uipanel('Parent',this.hdle,  ...
                'Title', string, ...
                'FontName', this.font_name,...
                'FontSize', this.font_size,...
                'Units', 'pixel',...
                'visible', 'off');
            
            
            
            layout = [{{'hsep.5'}}; layout;  {{'hsep.2'}}];
            [w, h] = get_layout_sz(this, layout);
            
            e = panel_elem(p, layout);
            e.w = w;
            e.h = h;
            this.uimap(id)= e;
            
        end
        
        function e = create_group(this, id, layout)
            [w, h] = get_layout_sz(this, layout);
            e = group_elem(layout);
            e.w = w;
            e.h = h;
            this.uimap(id)= e;
        end
        
        %% Advanced elements
        function e = create_changing_button(this, id,names, callback )
        % button cycling through a list of options
        
        e = create_button(this,id, names{1},@(o,e)(callback_changing));
        
        this.data_gui.(id).Value= 1;
        this.data_gui.(id).String = names{1};
        this.data_gui.(id).StringList =names;
        num_options = numel(names);
        
            function callback_changing()
                idx = this.data_gui.(id).Value;
                if idx==num_options
                    idx = 1;
                else
                    idx=idx+1;
                end           
                this.data_gui.(id).Value=idx;
                this.data_gui.(id).String = this.data_gui.(id).StringList{idx};
                update(); 
            end
            
            function update()                
                set(e.hdle, 'String', this.data_gui.(id).String);
                callback();
            end
        end
                       
        function  e = create_domain_slider(this, id, name, domain, value, callback)
         % Creates a slider from a BreachDomain    
            id_min = [id '_min'];
            id_val = [id '_val'];
            id_max = [id '_max'];
            id_slider = [id '_slider'];
            
            e_min = this.create_edit(id_min, num2str(0), @(o,e)(callback_edit()),  1/4, 1/2);
            e_min.htop = e_min.htop/2;
            e_min.wright=0;
            
            e_val =this.create_edit(id_val, num2str(0.5),@(o,e)(callback_edit()), 1/2, 1/2);
            e_val.htop = e_val.htop/2;
            e_val.wleft = e_val.wleft/2;
            e_val.wright = e_val.wright/2;
            
            e_max = this.create_edit(id_max, num2str(1),@(o,e)(callback_edit()), 1/4, 1/2);
            e_max.wleft = 0;
            e_max.htop = e_max.htop/2;
            
            e_slider= this.create_slider(id_slider, name, @(o,e)(callback_slider()));
            
            this.set_h(id_slider, 1/2);
            this.set_w(id_slider, 1);
            
            %% create panel
            layout_panel =  { ...                
                {id_min, id_val, id_max}; ...
                {id_slider};...                
                };
            
            string = name;
            e = this.create_panel(id,string,  layout_panel);
                       
          %% Setup values of stuff                      
            value = domain.checkin(value);
            this.data_gui.(id).name =name;
            this.data_gui.(id).domain =domain;
            this.data_gui.(id).value = value;           
            this.data_gui.(id).callback = callback;
            
            function callback_edit()                
              
                bug= 0;
                dom = this.data_gui.(id).domain.domain;
                if isempty(dom)
                    old_min = -inf;
                    old_max = inf;                    
                else
                    old_min = dom(1);
                old_max =dom(2);                
                end
                old_value = this.data_gui.(id).value;
                
                min_value  = str2double(get(e_min.hdle, 'String'));                                
                if isnan(min_value)
                    bug= 1;
                    min_value = old_min;
                end                
                max_value  = str2double(get(e_max.hdle, 'String'));
                if isnan(max_value)
                    bug= 1;
                    max_value = old_max;
                end
                val  = str2double(get(e_val.hdle, 'String'));                   
                if isnan(val)
                    bug= 1;
                    val = old_value;
                end

                if max_value-min_value<0
                    bug = 1;
                    if old_max ~= max_value
                        max_value = min_value;
                    else
                        min_value= max_value;
                    end
                end
                
                if bug||~(isequal(old_min,min_value)&&isequal(old_max,max_value)&&isequal(old_value,val))
                    if min_value == -inf && max_value== inf
                        this.data_gui.(id).domain.domain = [];
                    else
                        this.data_gui.(id).domain.domain = [min_value, max_value];
                    end
                    this.data_gui.(id).value = this.data_gui.(id).domain.checkin(val);
                    update(~bug&&~isequal(this.data_gui.(id).value,old_value ));
                end
            end                        
            
            function callback_slider()
                
                dom = this.data_gui.(id).domain;
                val=  get(e_slider.hdle,'Value');
                if strcmp(dom.type,'enum')
                  m =str2double( get(e_min.hdle,'String'));
                  M =str2double( get(e_max.hdle,'String'));
                  enum = dom.enum;                  
                  enum = enum(enum>=m);
                  enum = enum(enum<=M);
                   val = enum(val);
                else
                    val = this.data_gui.(id).domain.checkin(val);
                end
                
                old_value = this.data_gui.(id).value;
                
                if ~isequal(old_value , val)
                    this.data_gui.(id).value = val;
                    update(true);
                end
            end
            
            function update(call)
                      this.update_domain_slider(id, call);
            end
                
        end
               
        function update_domain_slider(this, id, call)
                % if call is true, call user provided callback of the slider. 
                % reads Breach Domain and value. value is in domain already.
                name =   this.data_gui.(id).name;
                dom = this.data_gui.(id).domain.domain;
                val = this.data_gui.(id).value;                
                
                id_min = [id '_min'];
                id_val = [id '_val'];
                id_max = [id '_max'];
                id_slider = [id '_slider'];
                
                e = this.uimap(id);
                e_slider = this.uimap(id_slider);
                e_min = this.uimap(id_min);
                e_max = this.uimap(id_max);
                e_val = this.uimap(id_val);

                if isempty(dom)
                    min = '-inf';
                    max = 'inf';
                    min_slider= val-10;
                    max_slider = val+10;
                elseif (dom(1) ==-inf) && dom(2)<inf
                    min = '-inf';
                    max = num2str(dom(2));
                    min_slider= val-100;
                    max_slider = dom(2);                                          
                elseif (dom(2) ==inf) && dom(1)>-inf
                    min = num2str(dom(1));
                    max = 'inf';                   
                    min_slider= dom(1);
                    max_slider = val+100;                                          
                else
                    min = num2str(dom(1));
                    max = num2str(dom(2));
                    min_slider= dom(1);
                    max_slider = dom(2);                                          
                end
                                                        
                switch this.data_gui.(id).domain.type
                    case 'bool'
                        min = 0;
                        max = 1;
                                                
                        set(e_slider.hdle,'Max',1);
                        set(e_slider.hdle,'Min',0);                        
                        set(e_slider.hdle,'SliderStep',[1 1]);
                        set(e_slider.hdle,'Value', val);                                      
                    
                    
                    case 'int'
                        
                        min = ceil(str2num(min));
                        max = floor(str2num(max));
                        min_slider= ceil(min_slider);                        
                        max_slider = floor(max_slider);                        
                        if min_slider == max_slider
                            max_slider = max_slider+1;
                        end
                        
                        set(e_slider.hdle,'Min', min_slider);
                        set(e_slider.hdle,'Max', max_slider);
                        set(e_slider.hdle,'Value', val);                       
                        set(e_slider.hdle,'SliderStep',[1 10]/(max_slider - min_slider));
                    
                    case 'enum'
                        enum = this.data_gui.(id).domain.enum;
                        enum = enum(enum>=str2double(min));
                        enum = enum(enum<=str2double(max));
                        if isempty(enum)
                            enum = nan;                                                    
                        end
                        min_slider = 1;
                        max_slider = numel(enum);
                        val_slider = find(val==enum,1);
                        if isempty(val_slider)
                            val_slider = 1;
                            val = nan;
                        end
                        if min_slider==max_slider
                            max_slider=min_slider+1;
                        end
                        
                        set(e_slider.hdle,'Min', min_slider);
                        set(e_slider.hdle,'Max', max_slider);
                        set(e_slider.hdle,'Value', val_slider);
                        set(e_slider.hdle,'SliderStep',[1 10]/(max_slider - min_slider));                                           
                    otherwise
                        if min_slider == max_slider
                            max_slider = max_slider+10*eps;
                        end                        
                        
                        set(e_slider.hdle,'Min', min_slider);
                        set(e_slider.hdle,'Max', max_slider);
                        set(e_slider.hdle,'Value', val);                                        
                end
                set(e_min.hdle,'String', num2str(min));
                set(e_max.hdle,'String', num2str(max));
                set(e_val.hdle,'String', num2str(val));
                
                set(e.hdle,'Title', [name   '  (' this.data_gui.(id).domain.short_disp() ')']); 
                if call
                      this.data_gui.(id).callback();
                end
        end
       
        function e = create_radio_panel(this,id,title, list, cb_fun) 
            this.data_gui.(id).list = list;
            num_radios = numel(list);
            
            this.data_gui.(id).idx_selected = zeros(1,num_radios);
            this.data_gui.(id).selected='';
            
           
            layout_panel = cell(num_radios,1); 
            id_radios = cell(num_radios);
            for ir = 1:num_radios
                id_radios{ir} = [id '_radio' num2str(ir)];
                this.create_radio(id_radios{ir}, list{ir}, @(o,e)(callback_radio(ir)));                
                layout_panel{ir,1} = id_radios(ir);
            end
                                    
            e = this.create_panel(id,title,layout_panel);
            callback_radio(1);
            
            function callback_radio(idx_r)
                old_selected = this.data_gui.(id).idx_selected;
                this.data_gui.(id).idx_selected = zeros(1,num_radios);
                this.data_gui.(id).idx_selected(idx_r) = 1;
                this.data_gui.(id).selected= list{idx_r};                            
                if ~isequal(old_selected, this.data_gui.(id).idx_selected)
                     update();         
                     drawnow;
                     cb_fun(idx_r);
                end
            end            
            
            function update()
                for ir = 1:num_radios                    
                    id_radio =[id '_radio' num2str(ir)];
                    if this.data_gui.(id).idx_selected(ir)
                        this.set_by_id(id_radio, 'Value',1);
                    else
                        this.set_by_id(id_radio, 'Value',0);
                    end
                end                
                                
            end
            
        end
       
        
        %% Callbacks
        function button_ok_callback(this, hobj)
            set(hobj, 'String', 'Please wait...')
            drawnow;
            close(this.hdle);
        end
        
        function button_cancel_callback(this)
            this.output=[];
            close(this.hdle);
        end
       
    end
    
    %% Layout
    properties
        layout   % cell of cells of dimension n x 1 where n is number of rows
        
        font_name = 'Arial'
        font_size = 10
        max_char = 80
        
        hdle  % handle to GUI top dialog window
        
        wunit = 400
        hunit = 40
                
    end
    
    methods
        function resize_callback(this,o,e)
            % Empty for now
        end
        
        function reset_top_fig(this)
            this.hdle = figure('Name',this.title, ...                
                'Toolbar', 'none',....
                'Menubar', 'none',...
                'NumberTitle','off',...
                'SizeChangedFcn', @(o,e)(this.resize_callback(o,e)));
        end
        
        function set_layout(this,layout, parent, w0,h0,fit)
            %% Arguments
            if nargin < 6
                fit = 'fit'; % fits parent figure to layout size and sets units to Normalized, no rescaling
            end
            
            if nargin<5
                h0 = 0;
            end
            
            if nargin<4
                w0 = 0;
            end
            
            if nargin<3
                parent = this;
            end
            
            if nargin<2
                layout = this.layout;
            end
            
            
            %% Get scale
            if strcmp(fit, 'fit')
                if isempty(this.layout)
                    scale = [1 1];
                else
                    pos = get(this.hdle, 'Position');
                    [w, h] = this.get_layout_sz(this.layout);
                    scale = [pos(3)/(w*this.wunit) pos(4)/(h*this.hunit)];
                    this.wunit = this.wunit*scale(1);
                    this.hunit = this.hunit*scale(2);
                end
            end
            
            this.set_for_all('Units', 'pixels');
            
            wunit = this.wunit;
            hunit = this.hunit;
            
            w = w0;
            hcur = h0;
            
            wmax = w;
            hmax = hcur;
            
            % absolute placement
            for ih = size(layout,1):-1:1  % along rows
                
                rowcurr = layout{ih};
                
                % get hmax for current row
                hmax = 0;
                for ei = 1:numel(rowcurr)
                    try
                        e = this.uimap(rowcurr{ei});
                    catch
                        error('element %s not found',rowcurr{ei});
                    end
                    hmax = max([hmax, e.h]);
                end
                
                % place elements of current row, aligned at top
                w=w0;
                for ei = 1:numel(rowcurr)
                    e = this.uimap(rowcurr{ei});
                    he = e.h;   % height of current elememt
                    h = hcur+(hmax-he)*hunit;
                    e.disp_in(this, parent, w, h, wunit, hunit);
                    w = w + e.w*wunit;
                end
                hcur = hcur+ hmax*hunit;
                wmax = max([wmax, w]);
            end
            hmax = hcur;
            
            switch fit
                case 'init'
                    % fit parent size, preserve scale and top position
                    old_pos = get(parent.hdle,'Position');
                    pos = old_pos;
                    pos(4) = hmax;
                    pos(3) = wmax;
                    
                    % Make sure position of top left corner does not change
                    hpos_old = old_pos(2)+old_pos(4); % pixel h position
                    pos(2) = hpos_old-hmax;
                    
                    set(parent.hdle, 'Position', pos);
                    this.layout = layout;
                    
                case 'fit'
                    
                    % fit parent size, preserve scale and top position
                    old_pos = get(parent.hdle,'Position');
                    pos = old_pos;
                    pos(4) = hmax;
                    pos(3) = wmax;
                    
                    % Make sure position of top left corner does not change
                    hpos_old = old_pos(2)+old_pos(4); % pixel h position
                    pos(2) = hpos_old-hmax;
                    set(parent.hdle, 'Position', pos);
                    this.layout = layout;
                    
            end
            
        end
        
        function set_for_all(this, varargin)
            values = this.uimap.values();
            for ei = 1:numel(values)
                e = values{ei};
                if ~isempty(e.hdle)
                    set(e.hdle,varargin{:});
                end
            end
        end
        
        function set_for_all_layout(this,layout, varargin)
            for ih = size(layout,1):-1:1  % along rows
                rowcurr = layout{ih};
                for ei = 1:numel(rowcurr)
                    e = this.uimap(rowcurr{ei});
                    if ~isempty(e.hdle)
                        set(e.hdle,varargin{:});
                    end
                end
            end
        end
        
        function enable_resizable(this)            
            set(this.hdle,'Resize', 'on');                    
            this.set_for_all('Units', 'Normalized');            
        end
        
        function disable_resizable(this)
            set(this.hdle,'Resize', 'off');        
        end
        
        function set_by_id(this, id, varargin)
            e = this.uimap(id);
            set(e.hdle, varargin{:});
        end
        
        function out= get_by_id(this, id, varargin)
            e = this.uimap(id);
            out = get(e.hdle, varargin{:});
        end
        
        function [w, h] = get_layout_sz(this, layout)
            % get layout size in layout units
            w=0;
            h=0;
            for ih = size(layout,1):-1:1  % along rows
                
                rowcurr = layout{ih};
                
                % get hmax for current row
                hmax = 0;
                wrow = 0;
                for ei = 1:numel(rowcurr)
                    try
                        e = this.uimap(rowcurr{ei});
                    catch
                        error('Element %s not found',rowcurr{ei});
                    end
                    hmax = max([hmax, e.h]);
                    wrow = wrow+e.w;
                end
                w = max([w, wrow]);
                h = h+hmax;
            end
        end
        
        function set_w(this,id,w)
            e = this.uimap(id);
            e.w = w;
            this.uimap(id)=e;
        end
        
        function set_h(this,id,h)
            e = this.uimap(id);
            e.h = h;
            this.uimap(id)=e;
        end
        
        %% Misc
        function txt =truncate_txt(this, txt)
            if numel(txt)>this.max_char
                txt = [txt(1:this.max_char) '...'];
            end
        end
        
    end
    
    %% Templates
    
    methods
        
        
        
        
        
        %% A template for a button. Replace truc with some label.
        function this = button_template(this, mode,name,string,callback,w,h)
            % Names
            label = ['button_' name];
            if ~exist('string','var')||isempty('string')
                if isfield(this.data_gui, name)
                    string = this.data_gui.(name).string;
                else
                    w = name;
                end
                string  = name;
            end
            
            % sizes
            if ~exist('w', 'var')||isempty(w)
                if isfield(this.data_gui, name)
                    w = this.data_gui.(name).w;
                else
                    w = 1;
                end
                
            end
            if ~exist('h', 'var')||isempty(h)
                if isfield(this.data_gui, name)
                    h = this.data_gui.(name).h;
                else
                    h = 1;
                end
            end
            
            % Resolve call
            switch mode
                case 'create'
                    this.create_button(label, string, callback, w,h);
                case 'callback'
                    callback(o,e); % for manual call of the callback function I guess
                case 'update'
                    if isfield(this.data_gui, name)
                        this.set_by_id(label, 'String', this.data_gui.(name).string);
                    end
            end
            
            
        end
        
        %% Template for a slider
        function this = slider_template(this, mode)
            % Global variables
            label = 'button_truc';
            string  = 'Truc';
            
            % sizes
            if ~exist('w', 'var')||isempty(w)
                w = 1;
            end
            if ~exist('h', 'var')||isempty(h)
                h = 1;
            end
            
            % Resolve call
            switch mode
                case 'create'
                    create_this()
                case 'callback'
                    callback_this();
                case 'update'
                    update_this();
            end
            
            % Sub functions
            function create_this() % create and initializes elements
                this.data_gui.truc.num_clicked = 0;
                this.create_button(label, string, @(hObj,evt) (button_truc(this, 'callback')), w,h);
            end
            
            function callback_this()
                % modifies this.data_gui
                this.data_gui.truc.num_clicked = this.data_gui.truc.num_clicked+1;
                
                % calls relevant update functions
                update_this();
            end
            
            function update_this()  % modifies itself based on data_gui
                this.set_by_id(label, 'String', sprintf('Truc: %g', this.data_gui.truc.num_clicked))
            end
            
        end
        
        
    end
    
    
    
end

