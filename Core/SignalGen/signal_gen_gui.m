function varargout = signal_gen_gui(varargin)
% SIGNAL_GEN_GUI MATLAB code for signal_gen_gui.fig
%      SIGNAL_GEN_GUI, by itself, creates a new SIGNAL_GEN_GUI or raises the existing
%      singleton*.
%
%      H = SIGNAL_GEN_GUI returns the handle to a new SIGNAL_GEN_GUI or the handle to
%      the existing singleton*.
%
%      SIGNAL_GEN_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIGNAL_GEN_GUI.M with the given input arguments.
%
%      SIGNAL_GEN_GUI('Property','Value',...) creates a new SIGNAL_GEN_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before signal_gen_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to signal_gen_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help signal_gen_gui

% Last Modified by GUIDE v2.5 11-Jun-2019 17:16:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @signal_gen_gui_OpeningFcn, ...
    'gui_OutputFcn',  @signal_gen_gui_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before signal_gen_gui is made visible.
function signal_gen_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to signal_gen_gui (see VARARGIN)

% Set fonts and size depending on system
if ismac
    FONT=12;
    %POS = [50 10 200 50];
    handles.TBL_SZ = {200 120 120 150 80} ;
else
    FONT=10;
    %POS = [50 10 200 50];
    handles.TBL_SZ = {300 80 80 100 100 100};
end

hfn = fieldnames(handles);
for ifn = 1:numel(hfn)
    try
        set(handles.(hfn{ifn}), 'FontSize', FONT);
    end
end
%set(handles.main, 'Position',POS);

handles.select_cells = [];
handles.constraints_map = containers.Map('KeyType', 'int32','ValueType','any');

set(handles.checkbox_enveloppe, 'Value', 0);
set(handles.text_max_time_enveloppe, 'Visible', 'off');
set(handles.edit_max_time_enveloppe, 'Visible', 'off');
handles.max_time_enveloppe = 3;

% get signal names
if isa(varargin{1}, 'BreachOpenSystem')
    handles.B = varargin{1};
    handles.IG = handles.B.InputGenerator.copy();
    %recover domains
    handles.IG.Domains = handles.B.GetDomain(handles.IG.P.ParamList);    
    signal_names = handles.B.Sys.InputList;
    if numel(varargin)>=2&&(isstruct(varargin{2})||ischar(varargin{2}))
        cfg_in = varargin{2};
        % Add reading contraints    
        if isfield(cfg_in,'constraints_cfg')
            data = {};
            for ic = 1:numel(cfg_in.constraints_cfg)
                f = cfg_in.constraints_cfg{ic};
                data{ic,1} = f.id;
                data{ic,2} = f.expr;
                handles.constraints_map(ic) = struct('id',  f.id, 'expr',f.expr);
            end
            if ~isempty(data)
                set(handles.table_constraints, 'Data', data);
            end
        end
    end
            
elseif isstruct(varargin{1})||ischar(varargin{1})  % configuration struct
    cfg_in = varargin{1};
    handles.B = [];
    handles.IG = ReadInputGenCfg(varargin{1});
    signal_names = handles.IG.GetSignalList();
    if isfield(varargin{1}, 'sim_time')
        handles.time = varargin{1}.sim_time;
    end
    % Add reading contraints
    if isfield(cfg_in,'constraints_cfg')
        data = {};
        for ic = 1:numel(cfg_in.constraints_in)
            f = cfg_in.constraints_in{ic};
            data{ic}{1} = f.id;
            data{ic}{2} = f.expr; 
            handles.constraints_map(ic) = struct('id',  f.id, 'expr',f.expr);
        end
        if ~isempty(data)
            set(handles.table_constraints, 'Data', data);
        end
    end        
end
set(handles.popupmenu_signal_name, 'String', signal_names);

%
signal_types= {
    'constant_signal_gen',...
    'step_signal_gen',...
    'fixed_cp_signal_gen',...
    'var_cp_signal_gen',...
    'pulse_signal_gen',...
    'random_signal_gen'...
    'exponential_signal_gen'...
    'sinusoid_signal_gen'...
    'spike_signal_gen'...
    'from_file_signal_gen',...  % done from main gui now.
    };
set(handles.popupmenu_signal_gen_type, 'String', signal_types);

% Assign default signal gen
handles.signal_gen_map = containers.Map();
for isig= 1:numel(signal_names)
    c = signal_names{isig};
    handles.signal_gen_map(c)=constant_signal_gen({c});

    % try to import from B - works when one sg for one signal (TODO: generalize)
    try
        sg = handles.IG.GetSignalGenFromSignalName(c);
        if numel(sg.signals)==1
            handles.signal_gen_map(c)=sg;
            sg_class = class(sg);
            idx = find(strcmp(signal_types, sg_class));
            if isig == 1
                set(handles.popupmenu_signal_gen_type,'Value', idx);
            end
            sg.p0 = handles.IG.GetParam(sg.params,1);
        else
            sgs = sg.split();
            for isg =1:numel(sgs)               
                if strcmp(sgs{isg}.signals{1}, c)
                    sg = sgs{isg};
                    handles.signal_gen_map(c)= sg;
                    sg_class = class(sg);
                    idx = find(strcmp(signal_types, sg_class));
                    if isig == 1
                        set(handles.popupmenu_signal_gen_type,'Value', idx);
                    end
                    break;
                end
            end
        end
    catch
        handles.signal_gen_map(c)=constant_signal_gen({c});
    end
    
end

% Init time
if ~isfield(handles, 'time')
    if ~isempty(handles.B)
        handles.time = get_time_string(handles.B.GetTime());
    else
        handles.time = 0:.01:1;
    end
end

set( handles.edit_time, 'String', get_time_string(handles.time));
% Choose default command line output for signal_gen_gui
signal_gens= handles.signal_gen_map.values;
handles.output = BreachSignalGen(signal_gens);


% update config and params
update_config(handles);

% Update table (fix domains checkin)
update_uitable(handles);
sg = get_current_sg(handles);
[sg.params, sg.p0, sg.params_domain] = read_uitable_params(handles.uitable_params);

% Init plot
update_plot(handles);

% Update handles structure
guidata(hObject, handles);

function st = dbl2str(x)
st = num2str(x, '%0.5g');

function time_string = get_time_string(time)
if ischar(time)
    time_string = time;
elseif isnumeric(time)
    if isscalar(time)
        time_string = ['[0 ' dbl2str(time) ']'];
    elseif numel(time)==2
        time_string = ['[' dbl2str(time(1)) ' ' dbl2str(time(2)) ']'];
    elseif max(diff(diff(time)))<1000*eps
        time_string = ['0:' dbl2str(time(2)-time(1)) ':' dbl2str(time(end))];
    else
        time_string = ['[' num2str(time) ']'] ;
    end
end

% --- Outputs from this function are returned to the command line.
function varargout = signal_gen_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout{1} = handles;


% --- Executes on selection change in popupmenu_signal_gen_type.
function popupmenu_signal_gen_type_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_signal_gen_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get signal name and current generator
sig_name = get_current_signal(handles);
sg = get_current_sg(handles);

% assign new signal generator
idx = get(hObject,'Value');
classes = get(hObject,'String');
class_name = classes{idx};
 
niou_sg =  eval([class_name '({ sig_name });']);

cfg = sg.getSignalGenArgs();
niou_cfg = niou_sg.getSignalGenArgs();

% restore common generator parameters such as method, number of cp
cfg_val = {};
for ic =1:numel(niou_cfg);
   io = find(strcmp(niou_cfg{ic}, cfg),1);
    if ~isempty(io) % found in previous signal gen
        cfg_val{ic} =  sg.(niou_cfg{ic});
    else
        cfg_val{ic} =  niou_sg.(niou_cfg{ic});
    end   
end

niou_sg =  eval([class_name '({ sig_name }, cfg_val{:});']);

% restore what is possible for p0
for ip = 1:numel(niou_sg.params)
    p = niou_sg.params{ip};
    ipold= find(strcmp(p, sg.params),1);
    if ~isempty(ipold)
        niou_sg.p0(ip) = sg.p0(ipold);
        niou_sg.params_domain(ip) = sg.params_domain(ipold);
    end
end

handles.signal_gen_map(sig_name) = niou_sg;

% update config and params
update_config(handles);

% update plot
update_plot(handles);

% update stuff
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_signal_gen_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_signal_gen_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_signal_name.
function popupmenu_signal_name_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_signal_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String')); %returns popupmenu_signal_gen_type contents as cell array
sig_name = contents{get(hObject,'Value')}; %returns selected item from popupmenu_signal_gen_type
sg = handles.signal_gen_map(sig_name);
sg_class = class(sg);

classes = get(handles.popupmenu_signal_gen_type, 'String');
idx = find(strcmp(classes, sg_class));

if isempty(idx)
    classes = [classes {sg_class}];
    set(handles.popupmenu_signal_gen_type, 'String', classes);
    idx = numel(classes);
end
set(handles.popupmenu_signal_gen_type,'Value', idx);

% update config and params
update_config(handles);

% update plot
update_plot(handles);

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_signal_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_signal_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    time_st=get(hObject, 'String');
    time = evalin('base', time_st);
    handles.time = time_st;
    update_plot(handles);
catch
    set(hObject, 'String',handles.time);
end

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_ok.
function pushbutton_ok_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sgs = handles.signal_gen_map.values;
handles.output.InitSignalGen(sgs);
if ~isempty(handles.B)
    handles.B.SetInputGen(handles.output);
end
close(handles.main);


% --- Executes when entered data in editable cell(s) in uitable_params.
function uitable_params_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable_params (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

try 
    sg = get_current_sg(handles);
    [sg.params, sg.p0, sg.params_domain] = read_uitable_params(hObject);
    update_plot(handles);
    guidata(hObject, handles);
catch
    warndlg('Invalid input', 'Problem');
end


% --- Executes when entered data in editable cell(s) in uitable_config.
function uitable_config_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable_config (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)

try
    
    sg = get_current_sg(handles);
    content = get(hObject, 'Data');
    
    cfg_params = sg.getSignalGenArgs();
    
    % update param struct
    if  ~isempty(cfg_params)
        for ip = 1:numel(cfg_params)
            content{ip, 1} = cfg_params{ip};
            old_val = sg.(cfg_params{ip});
            if iscell(old_val)
                val = cell(1,1);
                val{1} = content{ip,2};
            else
                val= content{ip,2};
            end
            args_val{ip}= val;
        end
    end
    
    % create new signal gen
    sg_name = class(sg);
    sig_name = get_current_signal(handles);
    niou_sg = eval([sg_name '(sig_name, args_val{:});']);
    
    
    % restore parameter values
    for ip = 1:numel(niou_sg.params)
        p = niou_sg.params{ip};
        ipold= find(strcmp(p, sg.params),1);
        if ~isempty(ipold)
            niou_sg.p0(ip) = sg.p0(ipold);
            niou_sg.params_domain(ip) = sg.params_domain(ipold);            
        end        
    end
                
    niou_sg.signals_domain = sg.signals_domain;    
    handles.signal_gen_map(sig_name)= niou_sg;
    
    update_config(handles);
    update_plot(handles);
    guidata(hObject, handles);
    
catch
    update_config(handles);
    update_plot(handles);
    guidata(hObject, handles);
    
end

%% update functions
function update_config(handles)
% update config parameters

sg = get_current_sg(handles);
handles.uitable_config= update_cfg_uitable(sg, handles.uitable_config);
set(handles.uitable_config, 'ColumnWidth', {250 250});

function h_uitable = update_cfg_uitable(sg, h_uitable)
set(h_uitable,'RowName',{});
set(h_uitable,'ColumnName',{'Name','Value'});
set(h_uitable,'ColumnEditable', [false true]);

cfg_params = sg.getSignalGenArgs();
content = {'',''};
if  ~isempty(cfg_params)
    content = cell(1,1);
    for ip = 1:numel(cfg_params)
        par = cfg_params{ip};
        content{ip, 1} = par;
        val = sg.(cfg_params{ip});
        if iscell(val)&&~isempty(val)
            content{ip,2} = val{1} ;
        elseif isempty(val)
            content{ip,2} = '';
        else
            content{ip,2} = val;
        end
    end
end

set(h_uitable, 'Data', content);

function update_plot(handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = handles.axes1;
axes(a);
set(handles.text_computing, 'visible', 'on');
drawnow; 

% get current signal name
popup_sel_index = get(handles.popupmenu_signal_name, 'Value');
sig_names = get(handles.popupmenu_signal_name, 'String');
sig_name = sig_names{popup_sel_index};

% fetch current generator
sg = handles.signal_gen_map(sig_name);
cla;
legend off;
% compute and plot signal
time = evalin('base',handles.time);
sg.plot(sig_name, time);

title(sig_name, 'Interpreter', 'None');
grid on;
set(gca, 'FontSize',8)

% compute enveloppe
if get(handles.checkbox_enveloppe, 'Value')
    if ~isempty(handles.constraints_map)
        constraints = handles.constraints_map.values;
        sg.plot_enveloppe(sig_name,time, 'max_time', handles.max_time_enveloppe, 'constraints', constraints);
    else
        sg.plot_enveloppe(sig_name,time, 'max_time', handles.max_time_enveloppe);
    end
end
update_uitable(handles);
update_constraints_table(handles);
set(handles.text_computing, 'visible', 'off');
drawnow; 




function update_uitable(handles)
sg = get_current_sg(handles);
fill_uitable_params(handles.uitable_params, sg.params, sg.p0, sg.params_domain);
set(handles.uitable_params, 'ColumnWidth', handles.TBL_SZ);

%% helpers
function sig_name = get_current_signal(handles)
contents = cellstr(get(handles.popupmenu_signal_name,'String')); %returns popupmenu_signal_gen_type contents as cell array
sig_name = contents{get(handles.popupmenu_signal_name,'Value')}; %returns selected item from popupmenu_signal_gen_type

function sg = get_current_sg(handles)
sig_name = get_current_signal(handles);
sg = handles.signal_gen_map(sig_name);

% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on key press with focus on uitable_params and none of its controls.
function uitable_params_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to uitable_params (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

if (isa(eventdata, 'matlab.ui.eventdata.UIClientComponentKeyEvent'))
    switch eventdata.Key
        case 'return'
            if strcmp(eventdata.Modifier, 'shift')
                val = inputdlg('Enter value');
                if ~isempty(val)&&~isempty(handles.select_cells)
                    tdata = get(hObject,'Data');
                    for irow = 1:size(handles.select_cells)
                        num_val = str2num(val{1});
                        if numel(num_val)>1
                            tdata{handles.select_cells(irow, 1),handles.select_cells(irow, 2)} = val{1};
                        else
                            tdata{handles.select_cells(irow, 1),handles.select_cells(irow, 2)} = num_val;
                        end
                    end
                    set(hObject,'Data',tdata);
                    
                    sg = get_current_sg(handles);
                    [sg.params, sg.p0, sg.params_domain] = read_uitable_params(hObject);
                    
                    update_plot(handles);
                    guidata(hObject, handles);
                end
            end
            
    end
end


% --- Executes when selected cell(s) is changed in uitable_params.
function uitable_params_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable_params (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

handles.select_cells = eventdata.Indices;
guidata(hObject,handles);


% --- Executes on button press in button_cancel.
function button_cancel_Callback(hObject, eventdata, handles)
handles.output.signalGenerators = {};
close(handles.main);


% --- Executes when entered data in editable cell(s) in table_constraints.
function table_constraints_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_constraints (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

idx = eventdata.Indices;
data = get(hObject, 'Data');
if idx(2) == 1
    id = eventdata.EditData;
    expr = data{idx(1), 2};
else
    id = data{idx(1),1};
    expr = eventdata.EditData;    
end

[handles, data, status] = check_constraint(handles, id, expr, data, idx(1));
if status
    set(hObject, 'Data', data);
    update_constraints_table(handles);
end
guidata(hObject,handles);

function [handles, data,status] = check_constraint(handles, id, expr, data, irow)

status =1;
try
    if ~isempty(id)&&~isempty(expr)
        phi = STL_Formula(id,expr);
        [sigs,params] =  STL_ExtractSignals(phi);
        time = evalin('base',handles.time);
        S = BreachSignalGen(handles.signal_gen_map.values);
        Usigs = S.GetSignalList();
        Uparams = S.GetParamList();
        
        undef_sig = setdiff(sigs, Usigs);
        if ~isempty(undef_sig)
            errordlg(sprintf('The constraint [ %s ] involves signal %s unrelated to the input signals. Please double-check for potential typo.', expr, undef_sig{1}), 'Error');
            error('undef_sig');
        end
        
        undef_param = setdiff(params, Uparams);
        if ~isempty(undef_param)
            errordlg(sprintf('The constraint [ %s ] is syntactically correct but involve parameter %s unrelated to the input parameters. Please double-check for potential typo.', expr, undef_param{1}), 'Error');
            error('undef_param');
        end
        
        S.Sim(time);
        v = S.CheckSpec(phi);
        data{irow, 3} = v;
        handles.constraints_map(irow) = struct('id',  id, 'expr',expr);    
    else
        if handles.constraints_map.isKey(irow)
            handles.constraints_map.remove(irow);
        end
        data{irow, 3} = [];
    end
catch ME
    status = 0;
    if ~strcmp(ME.message, 'undef_sig')&&~strcmp(ME.message, 'undef_param')
        errordlg(sprintf('Error Message: %s', ME2st(ME)), 'Problem with Expression');
    end
    data{irow, 3} = NaN;
    set(handles.table_constraints, 'Data', data);
end


function update_constraints_table(handles)
    old_data = get(handles.table_constraints,'Data');
    data = old_data;
    for irow = 1:size(data, 1)
        id = old_data{irow,1};
        expr = data{irow,2};
        [handles, data] = check_constraint(handles, id, expr, data, irow);
    end
    
    if ~(isempty(data{end,1})&&isempty(data{end,2}))
        data{end+1,1} = '';
    end
    set(handles.table_constraints, 'Data', data);
    if ~isequal(data, old_data)
         update_plot(handles);
    end


% --- Executes on button press in checkbox_enveloppe.
function checkbox_enveloppe_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_enveloppe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_plot(handles);
if get(hObject, 'Value')
    set(handles.edit_max_time_enveloppe, 'Visible', 'on');
    set(handles.text_max_time_enveloppe, 'Visible', 'on');
else
    set(handles.edit_max_time_enveloppe, 'Visible', 'off');
    set(handles.text_max_time_enveloppe, 'Visible', 'off');
end

function edit_max_time_enveloppe_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max_time_enveloppe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_max_time_enveloppe as text
%        str2double(get(hObject,'String')) returns contents of edit_max_time_enveloppe as a double
try
    handles.max_time_enveloppe = str2double(get(hObject,'String'));
    guidata(hObject, handles);
    update_plot(handles);
end





% --- Executes during object creation, after setting all properties.
function edit_max_time_enveloppe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max_time_enveloppe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
