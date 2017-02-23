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

% Last Modified by GUIDE v2.5 21-Feb-2017 11:54:35

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

% get signal names

if isa(varargin{1}, 'BreachOpenSystem')
  handles.B = varargin{1};
  signal_names = handles.B.Sys.InputList; 
else
    handles.B = [];
    signal_names = varargin{1};
end
set(handles.popupmenu_signal_name, 'String', signal_names);

% Assign default signal gen 
handles.signal_gen_map = containers.Map();
for isig= 1:numel(signal_names)
    c = signal_names{isig};
    handles.signal_gen_map(c)=constant_signal_gen({c}); 
end

% 
signal_types= {
 'constant_signal_gen',...
 'step_signal_gen',...
 'fixed_cp_signal_gen',...
 'var_cp_signal_gen',...
 'pulse_signal_gen',...
 };
set(handles.popupmenu_signal_gen_type, 'String', signal_types);

% Init time
handles.time = 0:.1:10;

% Choose default command line output for signal_gen_gui
signal_gens= handles.signal_gen_map.values;
handles.output = BreachSignalGen(signal_gens);

% update config and params
update_config(handles);

% Update handles structure
guidata(hObject, handles);

% Init plot
update_plot(handles);

% UIWAIT makes signal_gen_gui wait for user response (see UIRESUME)
% uiwait(handles.main);


% --- Outputs from this function are returned to the command line.
function varargout = signal_gen_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in popupmenu_signal_gen_type.
function popupmenu_signal_gen_type_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_signal_gen_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% get signal name 
sig_name = get_current_signal(handles);

% assign new signal generator
idx = get(hObject,'Value');
classes = get(hObject,'String');
class_name = classes{idx};
handles.signal_gen_map(sig_name) = eval([class_name '({sig_name});']);

% update config and params
update_config(handles);

% update stuff
guidata(hObject,handles);

% update plot 
update_plot(handles);

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

guidata(hObject,handles);

% update plot
update_plot(handles);

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


% --- Executes on selection change in popupmenu_config_param.
function popupmenu_config_param_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_config_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[~, cfg_val]= get_current_cfg(handles);
if isnumeric(cfg_val)
   cfg_val = num2str(cfg_val); 
end
set(handles.edit_cfg_val, 'String', cfg_val);
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function popupmenu_config_param_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_config_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_param.
function popupmenu_param_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[~, pval]= get_current_param(handles);
set(handles.edit_param_val, 'String', num2str(pval));
guidata(hObject,handles); 


% --- Executes during object creation, after setting all properties.
function popupmenu_param_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_cfg_val_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cfg_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get list of arguments 
args = get(handles.popupmenu_config_param, 'String');
idx_arg = get(handles.popupmenu_config_param,'Value');
if isempty(args{1})
    return 
else
 sg = get_current_sg(handles);   
 args_val = cell(1, numel(args));
 for ia = 1:numel(args)
    args_val{ia} = sg.(args{ia});
 end
 niou_arg_val = get(hObject,'String');   
 if isnumeric(args_val{idx_arg})
    niou_arg_val =  str2num(niou_arg_val);   
 end
 args_val{idx_arg} = niou_arg_val;
 
 sg_name = class(sg);
 sig_name = get_current_signal(handles);
 niou_sg = eval([sg_name '(sig_name, args_val{:});']);
 if isequal(size(sg.p0),size(niou_sg.p0))
    niou_sg.p0 = sg.p0;
 end
 
 handles.signal_gen_map(sig_name)= niou_sg;
 
 % update 
 update_config(handles);
 guidata(hObject,handles);
 update_plot(handles);  
 
end


% Hints: get(hObject,'String') returns contents of edit_cfg_val as text
%        str2double(get(hObject,'String')) returns contents of edit_cfg_val as a double


% --- Executes during object creation, after setting all properties.
function edit_cfg_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cfg_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_param_val_Callback(hObject, eventdata, handles)
% hObject    handle to edit_param_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sg = get_current_sg(handles);
pname =get_current_param(handles);
pval = str2num(get(hObject,'String'));
sg.set_param(pname,pval);

update_config(handles);
guidata(hObject,handles);
update_plot(handles);
uicontrol(handles.popupmenu_param);


% --- Executes during object creation, after setting all properties.
function edit_param_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_param_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.time = str2num(get(hObject,'String'));
update_plot(handles);
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
handles.B.RunGUI;
close(handles.main);



%% update functions
function update_config(handles)

% update config parameters
pop_cfg = handles.popupmenu_config_param;
sg = get_current_sg(handles);
cfg_param = sg.getSignalGenArgs();
if  isempty(cfg_param)
cfg_param =  {''};
end
set(pop_cfg, 'String', cfg_param);
set(pop_cfg, 'Value', 1);

% update parameters
pop_param = handles.popupmenu_param;
params = sg.params;
if  isempty(params)
    params =  {''};
end

% display pname: pval 
for ip = 1:numel(params)
    pname = params{ip};
    pval  = sg.get_param(pname); 
    params{ip}= [pname ': ' num2str(pval)];
end

set(pop_param, 'String', params);
set(pop_param, 'Value', 1);

% update cfg and param values
[~, cfg_val] = get_current_cfg(handles);
set(handles.edit_cfg_val, 'String', cfg_val);

[~, param_val] = get_current_param(handles);
set(handles.edit_param_val, 'String', param_val);

function update_plot(handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);

% get current signal name
popup_sel_index = get(handles.popupmenu_signal_name, 'Value');
sig_names = get(handles.popupmenu_signal_name, 'String');
sig_name = sig_names{popup_sel_index};
    
% fetch current generator
sg = handles.signal_gen_map(sig_name);

% compute and plot signal
X = sg.computeSignals(sg.p0, handles.time);
plot(handles.time, X );        
title(sig_name);    
grid on;
set(gca, 'FontSize',8)

%% helpers
function sig_name = get_current_signal(handles)
contents = cellstr(get(handles.popupmenu_signal_name,'String')); %returns popupmenu_signal_gen_type contents as cell array
sig_name = contents{get(handles.popupmenu_signal_name,'Value')}; %returns selected item from popupmenu_signal_gen_type

function sg = get_current_sg(handles)
sig_name = get_current_signal(handles);
sg = handles.signal_gen_map(sig_name);

function [cfg_name, cfg_val] = get_current_cfg(handles)
sg = get_current_sg(handles);
cfg_content = get(handles.popupmenu_config_param, 'String');
cfg_idx = get(handles.popupmenu_config_param, 'Value');
cfg_name = cfg_content{cfg_idx};
if ~isempty(cfg_name)
    cfg_val = sg.(cfg_name);
else
    cfg_val = [];
end

function [param_name, param_val] = get_current_param(handles)
sg = get_current_sg(handles);
param_content = get(handles.popupmenu_param, 'String');
param_idx = get(handles.popupmenu_param, 'Value');

param_name = param_content{param_idx};
sp = strsplit(param_name,':');
param_name = sp{1};
param_val = sg.get_param(param_name);


% --- Executes on key press with focus on popupmenu_param and none of its controls.
function popupmenu_param_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_param (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

if (isa(eventdata, 'matlab.ui.eventdata.UIClientComponentKeyEvent'))
    switch eventdata.Key
     case 'return'
      uicontrol(handles.edit_param_val);
    end
end
