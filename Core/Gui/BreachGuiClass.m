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
            
            hs25 = this.create_separator('hsep.25');
            hs25.h = 0.25;
            
            hs3= this.create_separator('hsep.5');
            hs3.h = 0.5;
            
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
        font_size = 11
        max_char = 80
        
        hdle  % handle to GUI top dialog window
        
        wunit = 500
        hunit = 50
                
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

