function varargout = BreachGui(varargin)
% BREACHGUI M-file for BreachGui.fig
%      BREACHGUI, by itself, creates a new BREACHGUI or raises the existing
%      singleton*.
%
%      H = BREACHGUI returns the handle to a new BREACHGUI or the handle to
%      the existing singleton*.
%
%      BREACHGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BREACHGUI.M with the given input arguments.
%
%      BREACHGUI('Property','Value',...) creates a new BREACHGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BreachGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BreachGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BreachGui

% Last Modified by GUIDE v2.5 17-Apr-2017 11:54:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @BreachGui_OpeningFcn, ...
    'gui_OutputFcn',  @BreachGui_OutputFcn, ...
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


% --- Executes just before BreachGui is made visible.
function BreachGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BreachGui (see VARARGIN)

handles = info(handles, 'Starting Breach.... (Memory tip: his first name is Millard)');

% Set fonts and size depending on system
if ismac
    FONT=12;
    POS = [50 10 200 50];
    handles.TBL_SZ = {200 120 120 150 80} ;
else
    FONT=10;
    POS = [50 10 200 50];
    handles.TBL_SZ = {300 120 150 200 120} ;
end

hfn = fieldnames(handles);
for ifn = 1:numel(hfn)
    try 
        set(handles.(hfn{ifn}), 'FontSize', FONT);
    end
end
set(handles.breach, 'Position',POS);


crd = pwd;
set(hObject, 'Name', ['Breach (' crd  ')']);

%% used to memorize previous entries
handles.last_options =[];
handles.check_prop_options =[];

%% Init system panel
BrGUI = varargin{1};
Sys = BrGUI.Sys;

if isfield(BrGUI,  'SimInputsOnly')
    set(handles.checkbox_input_only,'Value',BrGUI.SimInputsOnly);
end

if (isfield(Sys,'tspan'))
    if isscalar(Sys.tspan)
        handles.last_tspan = ['[0 ' dbl2str(Sys.tspan) ']'];
    elseif numel(Sys.tspan)==2
        handles.last_tspan = ['[' dbl2str(Sys.tspan(1)) ' ' dbl2str(Sys.tspan(2)) ']'];
    else
        handles.last_tspan = ['0:' dbl2str(Sys.tspan(2)-Sys.tspan(1)) ':' dbl2str(Sys.tspan(end))];
    end
else
    handles.last_tspan = '';
end

handles.figp=[];

%% Init working sets panel

% Find out who's in the workspace 
ws_var = evalin('base', 'who');
for iv= 1:numel(ws_var)
    % is this a BreachSet
    BB__ = evalin('base', ws_var{iv});
    if isa(BB__, 'BreachSet') % found one, keep it
        handles.working_sets.(ws_var{iv}) = BB__;
        if BrGUI== BB__ % this is the caller
            handles.current_set = ws_var{iv};
        end
    end    
end

fnames = fieldnames(handles.working_sets);
igui = find(strcmp(fnames, handles.current_set)); 
set(handles.working_sets_listbox, 'Value', igui);

handles = update_working_sets_panel(handles);

%% Init properties panel
handles.idx_prop= 1;
handles.current_prop = ''; 
handles = update_properties_panel(handles);

%% Init modif panel
handles.current_pts = 1;
handles.refine = 1;
handles.plot_proj = [];
handles.refine_all = 0;
handles.halton = 0;
handles.refine_args = 0;
handles.show_params = BrGUI.P.ParamList;

handles.selected_param=1;
handles.selected_varying_param = 1;

% Init param pts plot
handles.current_plot{1} =[];
handles.current_plot{2} =[];
handles.current_plot{3} =[];
handles.current_marked = [];

handles = update_selected_domain(handles);
handles = update_modif_panel(handles);
handles = info(handles, 'Ready.');

% Choose default command line output for BreachGui
handles.output = BrGUI;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BreachGui wait for user response (see UIRESUME)
% uiwait(handles.breach);


function handles = update_selected_domain(handles)
Br = handles.working_sets.(handles.current_set);
params = Br.GetBoundedDomains();
handles.selected_params = intersect(handles.show_params, params);



% --- Outputs from this function are returned to the command line.
function varargout = BreachGui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function edit_epsi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_epsi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in button_remove_set.
function button_remove_set_Callback(hObject, eventdata, handles)
try
    old_name = handles.current_set;
    ind = handles.selected_varying_param;
    fn = fieldnames(handles.working_sets);
    st = fn{get(handles.working_sets_listbox,'Value')};
    val = get(handles.working_sets_listbox,'Value');
    
    if (val>1)
        val = val-1;
        set(handles.working_sets_listbox,'Value', val);
        handles.current_set = fn{val};
    else
        handles.current_set = fn{val+1};
    end
    
    handles.working_sets = rmfield(handles.working_sets,st);
    evalin('base', ['clear ' old_name]);
    handles =update_working_sets_panel(handles);
    handles= update_properties_panel(handles);
    handles =update_modif_panel(handles);
    guidata(hObject,handles);
    
end
% hObject    handle to button_remove_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on selection change in working_sets_listbox.
function working_sets_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to working_sets_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = get(hObject,'String');
fn = fieldnames(handles.working_sets);
set_name = fn{get(hObject,'Value')};
handles.current_set = set_name;
handles.current_pts = 1;
Br = handles.working_sets.(handles.current_set);
handles.show_params = Br.P.ParamList; 

handles = update_working_sets_panel(handles);
handles = update_modif_panel(handles);
handles= update_properties_panel(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function working_sets_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to working_sets_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in button_compute_traj.
function button_compute_traj_Callback(hObject, eventdata, handles)
% hObject    handle to button_compute_traj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Br = handles.working_sets.(handles.current_set);
tspan = inputdlg('Enter tspan (Format: [ti tf] or [t0 t1 t2 ... tn] or ti:dt:tf)','Compute trajectories', 1, {handles.last_tspan});
if isempty(tspan)
    return;
end
handles.last_tspan = tspan{1};
tspan = eval(tspan{1});

handles= info(handles,'Computing trajectories...');
Br.Sim(tspan);
handles= info(handles,'Computing trajectories... Done.');

handles = update_working_sets_panel(handles);
handles = update_modif_panel(handles);
guidata(hObject, handles);
h = BreachTrajGui(Br, handles);


% --- Executes during object creation, after setting all properties.
function edit_tspan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tspan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_default_param_Callback(hObject, eventdata, handles)
% hObject    handle to edit_default_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Br = handles.working_sets.(handles.current_set);
val = eval(get(hObject,'String'));
param = get(handles.edit_param_name,'String');
ii = FindParam(Br.P, param);

if ~any(Br.P.dim==ii)
    
    Br.P = ...
        SetParam(Br.P,param, val);
    
    handles.selected_param = FindParam(Br.P, param);
     
    if isfield(Br.P, 'traj')
        if(handles.selected_param <= Br.P.DimP)
            [~,~,~,~,tokens] = regexp(handles.last_tspan,'([0-9eE\+-\.]+):([0-9eE\+-\.]+):([0-9eE\+-\.]+)');
            tspan = str2double(tokens{1}{1}):str2double(tokens{1}{2}):str2double(tokens{1}{3});
            handles = info(handles,'Updating trajectories...');
            Br.P = ComputeTraj(Br.Sys,rmfield(Br.P,'traj'), tspan);
            handles = info(handles,'Updating trajectories... Done.');
        end
    end
    
    if isfield(Br.P, 'props')
        
        props = Br.P.props;
        props_values = Br.P.props_values;
        
        P0 = rmfield(Br.P,'props');
        P0 = rmfield(P0,'props_names');
        P0 = rmfield(P0,'props_values');
        
        for ii = 1:numel(props)
            phi = props(ii);
            tau = props_values(ii).tau;
            if isfield(P0,'traj')
                handles = info(handles, 'Computing satisfaction of formula...');
                P0 = SEvalProp(Br.Sys, P0, phi, tau);
                handles = info(handles, 'Computing satisfaction of formula... Done.');
            else
                PO = SPurge_props(P0);
                handles = info(handles, 'Parameters changed, recompute trajectories to re-evaluate properties.');
            end
        end
        Br.P = P0;
    end
    
    handles = update_modif_panel(handles);
    handles = update_properties_panel(handles);   
    guidata(hObject,handles);
    
else
    handles = info(handles, 'To change this value, use Modif current sample');
end

% --- Executes during object creation, after setting all properties.
function edit_default_param_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_default_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_num_samples_Callback(hObject, eventdata, handles)
% hObject    handle to edit_num_samples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Br = handles.working_sets.(handles.current_set);
st = get(hObject,'String');
if strcmp(st, 'all')
    handles.sample_arg_num_samples = 'all';
else
    try
        handles.sample_arg_num_samples = str2num(st);
    catch
        handles.sample_arg_num_samples = 0;
        return
    end
end

handles = run_sample_domain(handles);
handles =update_modif_panel(handles);
guidata(hObject,handles);

function handles  = run_sample_domain(handles)
Br = handles.working_sets.(handles.current_set);
handles = update_sample_args(handles);

Br.SampleDomain(handles.selected_params, ...
    handles.sample_arg_num_samples,...
    handles.sample_arg_method,...
    handles.sample_arg_multi);

function handles = update_sample_args(handles)
 val = get(handles.popup_sample_option,'Value'); 
     switch(val)
         case 1
             handles.sample_arg_multi = 'replace';
         case 2
             handles.sample_arg_multi = 'append';
         case 3
             handles.sample_arg_multi = 'combine';
     end
     
     val = get(handles.popup_method,'Value'); 
     switch(val)
         case 1
             handles.sample_arg_method = 'grid';
         case 2
             handles.sample_arg_method = 'corners';
         case 3
             handles.sample_arg_method = 'rand';
         case 4
             handles.sample_arg_method = 'quasi';
         case 5
             handles.sample_arg_method = 'all';           
     end
     
     
     





% --- Executes during object creation, after setting all properties.
function edit_num_samples_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_num_samples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in button_new_set.
function button_new_set_Callback(hObject, eventdata, handles)
% hObject    handle to button_new_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

   Br = handles.working_sets.(handles.current_set);
   names = fieldnames(handles.working_sets);
    new_name =  genvarname('Br',names);
    Br_new = Br.copy();
    Br_new.Reset();
    Br_new.ResetSelected();
    handles.working_sets = setfield(handles.working_sets, new_name, Br_new);
    assignin('base', new_name, handles.working_sets.(new_name));

    handles = update_working_sets_panel(handles);
    guidata(hObject,handles);

% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function test_memory_Callback(hObject, eventdata, handles)
% hObject    handle to test_memory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgbox('Who was the 13th president of the United States of America ?', ...
    'Memory Test')


% --------------------------------------------------------------------
function menu_load_working_set_Callback(hObject, eventdata, handles)
% hObject    handle to menu_load_working_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

  Br = handles.working_sets.(handles.current_set);
[FileName,PathName] = uigetfile('*.mat','Load Parameter Set...');
if(FileName==0)
    return;
end

handles.working_sets_file_name = [PathName, FileName];
try
    handles.working_sets = evalin('base', ['load(''' handles.working_sets_file_name ''')']);
    
    % Find out who's in the workspace
    ws_var = evalin('base', 'who');
    for iv= 1:numel(ws_var)
        % is this a BreachSet
        BB__ = evalin('base', ws_var{iv});
        if isa(BB__, 'BreachSet') % found one, keep it
            handles.working_sets.(ws_var{iv}) = BB__;
            if Br== BB__ % this is the caller
                handles.current_set = ws_var{iv};
            end
        end
    end
    
    fnames = fieldnames(handles.working_sets);
    igui = find(strcmp(fnames, handles.current_set));
    set(handles.working_sets_listbox, 'Value', igui);

    handles = update_working_sets_panel(handles);
    guidata(hObject,handles);
catch err
    warndlg(['Problem loading: ' err.message] );
    rethrow(err);
end


% --------------------------------------------------------------------
function menu_save_as_Callback(hObject, eventdata, handles)
% hObject    handle to menu_save_as (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    [FileName,PathName] = uiputfile('*.mat','Save Parameter Set As...');
    if(FileName==0)
        return;
    end
    
    handles.working_sets_file_name = [PathName, FileName];
    ws = handles.working_sets;
    handles = info(handles, 'Saving parameter set...');
    save(handles.working_sets_file_name, '-struct','ws');
    handles = info(handles, 'Saving parameter set... Done.');
    handles = update_working_sets_panel(handles);
    guidata(hObject,handles);
catch err
    warndlg(['Problem saving: ' err.message] );
    rethrown(err);
end


% --------------------------------------------------------------------
function quit_Callback(hObject, eventdata, handles)
% hObject    handle to quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(handles.breach);

% --- Executes on button press in button_copy_set.
function button_copy_set_Callback(hObject, eventdata, handles)
% hObject    handle to button_copy_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = get(handles.working_sets_listbox, 'Value');
names = fieldnames(handles.working_sets);
set_name = names{val};
new_name = genvarname(set_name,names);
Br = handles.working_sets.(handles.current_set);
Brcopy = Br.copy();
handles.working_sets = setfield(handles.working_sets, new_name, Brcopy);
assignin('base', new_name, Brcopy);
handles = update_working_sets_panel(handles);

guidata(hObject,handles);

function edit_rename_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    Br = handles.working_sets.(handles.current_set);
    old_name = handles.current_set;
    new_name= get(hObject,'String');
    if (isfield(handles.working_sets,new_name))
        handles = info(handles, 'Name exists already.');
        return
    end
    handles.working_sets = setfield(handles.working_sets,new_name, ...
        Br);
    handles.working_sets = rmfield(handles.working_sets, handles.current_set);
    handles.current_set = new_name;
    handles = update_working_sets_panel(handles);
    handles = update_modif_panel(handles);
    handles= update_properties_panel(handles);
    evalin('base', ['clear ' old_name]);
    guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_rename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menu_files_Callback(hObject, eventdata, handles)
% hObject    handle to menu_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider_current_pts_Callback(hObject, eventdata, handles)
% hObject    handle to slider_current_pts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = get(hObject,'Value');
handles.current_pts = round(val);
handles = update_modif_panel(handles);

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider_current_pts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_current_pts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in popup_pts3.
function popup_pts3_Callback(hObject, eventdata, handles)
% hObject    handle to popup_pts3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = get(hObject,'String');
strpop = contents{get(hObject,'Value')};
handles.current_plot_pts{3} = {strpop};
handles = plot_pts(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popup_pts3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_pts3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_pts2.
function popup_pts2_Callback(hObject, eventdata, handles)
% hObject    handle to popup_pts2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_pts2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_pts2

contents = get(hObject,'String');
strpop = contents{get(hObject,'Value')};
handles.current_plot_pts{2} = {strpop};
handles = plot_pts(handles);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popup_pts2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_pts2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_pts1.
function popup_pts1_Callback(hObject, eventdata, handles)
% hObject    handle to popup_pts1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = get(hObject,'String');
strpop = contents{get(hObject,'Value')};
handles.current_plot_pts{1} = {strpop};
handles = plot_pts(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popup_pts1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_pts1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function menu_save_Callback(hObject, eventdata, handles)
% hObject    handle to menu_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ws = handles.working_sets; 
    
try
    handles = info(handles, ['Saving parameter set to ' handles.working_sets_file_name '...']);
    save(handles.working_sets_file_name, '-struct', 'ws');
catch
    [FileName,PathName] = uiputfile('*.mat','Save Parameter Set As...');
    if(FileName==0)
        return;
    end
    handles.working_sets_file_name = [PathName  FileName];
    handles = info(handles, ['Saving parameter set to ' handles.working_sets_file_name '...']);
    save(handles.working_sets_file_name, '-struct', 'ws');
    handles = update_working_sets_panel(handles);
    
end
handles = info(handles, ['Saving parameter set to ' handles.working_sets_file_name '... Done']);
guidata(hObject, handles);

% --- Executes on button press in button_make_pts_set.
function button_make_pts_set_Callback(hObject, eventdata, handles)
% hObject    handle to button_make_pts_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);

names = fieldnames(handles.working_sets);
new_name = genvarname('P',names);
ipts = find(Br.P.selected);
if isempty(ipts)
    return
end

Ptmp = Sselect(Br.P,ipts);

handles.working_sets = setfield(handles.working_sets, new_name, Ptmp);
handles = update_working_sets_panel(handles);
guidata(hObject,handles);

% --- Executes on selection change in listbox_prop.
function listbox_prop_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_prop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

idx_prop = get(hObject,'Value');
fnames = fieldnames(handles.properties);
handles.current_prop = fnames{idx_prop};
prop = handles.properties.(handles.current_prop);

info_msg = disp(prop,0);
def_par = get_params(prop);
fnames = fieldnames(def_par);

ndf = numel(fnames);
if ndf>0
    info_msg = [info_msg '   |  default params:'];
end
for idf =1:ndf
    info_msg = [info_msg ' ' fnames{idf} '=' num2str(def_par.(fnames{idf}))];
end

set(handles.text_info, 'String', info_msg);
handles = plot_pts(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function listbox_prop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_prop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in button_check_property.
function button_check_property_Callback(hObject, eventdata, handles)
Br = handles.working_sets.(handles.current_set);
prop = handles.properties.(handles.current_prop);

handles = info(handles, 'Computing satisfaction of formula...');
Br.CheckSpec(prop); 
handles = info(handles, 'Computing satisfaction of formula... Done.');
handles = update_working_sets_panel(handles);
handles = update_modif_panel(handles);
handles = update_properties_panel(handles);
guidata(hObject,handles);


% --- Executes on button press in checkbox_selected.
function checkbox_selected_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);

val = get(hObject,'Value');
Br.P.selected(handles.current_pts)=val;
if isempty(find(Br.P.selected))
    set(handles.button_plot_selected, 'String', 'Plot trajs')
else
    set(handles.button_plot_selected, 'String', 'Plot select trajs')
end

guidata(hObject, handles);


% --- Executes on button press in button_edit_prop.
function button_edit_prop_Callback(hObject, eventdata, handles)
% hObject    handle to button_edit_prop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);

prop = handles.properties.(handles.current_prop);

prompt={'Edit formula expression:'};
name='Edit property';
numlines=1;
defaultanswer={disp(prop,-1)};
opts.Resize='on';
opts.WindowStyle='normal';
str=inputdlg(prompt,name,numlines,defaultanswer, opts);

if isempty(str)
    return;
end

try
    id = get_id(prop);
    phi_st = str{1};
    PHI = eval(['STL_Formula(''' id ''',''' phi_st ''')']);
    eval([get_id(PHI) '=PHI']);
    Br.AddSpec(PHI);
    catch
    s = lasterror;
    warndlg(['Invalid formula: ' s.message] );
    error(s);
    return
end


handles = update_properties_panel(handles);
guidata(hObject, handles);

function handles = update_working_sets_panel(handles)

%% Set title
str_name = sprintf('Workspace');
if(numel(str_name)>35)
    str_name = [str_name(1:26) '...' str_name(end-5:end)];
end
set(handles.working_sets_panel,'Title',str_name );


%% Filter out invalid entries     
fn = fieldnames(handles.working_sets);
for ii=1:numel(fn)
    if ~isa(handles.working_sets.(fn{ii}),'BreachSet')
        handles.working_sets = rmfield(handles.working_sets, fn{ii});
    end
end
fn = fieldnames(handles.working_sets);

for ii=1:numel(fn)
    pref = '';
    
    if isfield(handles.working_sets.(fn{ii}).P,'traj_ref')
        if all(handles.working_sets.(fn{ii}).P.traj_ref~=0)
            pref = [pref,'*']; %#ok<AGROW>
        elseif any(handles.working_sets.(fn{ii}).P.traj_ref~=0)
            pref = [pref,'+']; %#ok<AGROW>
        end
    elseif isfield(handles.working_sets.(fn{ii}).P,'traj');
        pref = [pref, '*']; %#ok<AGROW>
    end
    
    fn{ii} = [pref, fn{ii}];
    
end

set(handles.working_sets_listbox,'String', fn);
set(handles.edit_rename,'String',handles.current_set);


function handles = update_properties_panel(handles)
Br = handles.working_sets.(handles.current_set);

str_name = ['Requirements of ' handles.current_set ];
if (numel(str_name)>35)
    str_name= [str_name(1:30) '...' str_name(end-5:end)];
end
set(handles.panel_properties,'Title',str_name );

handles.properties = struct;
specs = Br.Specs;
list_of_their_names = specs.keys();

for iprop = 1:numel(list_of_their_names)
    this_name= list_of_their_names{iprop};
    handles.properties.(this_name) = specs(this_name);
end

fnames = fieldnames(handles.properties);

content={};
for i = 1:numel(fnames)
    st = disp(handles.properties.(fnames{i}),-1);
    iprop = find_prop(fnames{i},  Br.P);
    if (iprop)
        content = {content{:}, ['*' fnames{i} ': ' st]};
    else
        content = {content{:}, [fnames{i} ': ' st]};
    end
end

set(handles.listbox_prop,'String', content);

if ~isfield(handles, 'current_prop')&&~isempty(fnames)
   handles.current_prop = fnames{1};     
   set(handles.listbox_prop, 'Value',1); 
end

val = get(handles.listbox_prop,'Value');
val = min(val, numel(fnames));
if val>0
set(handles.listbox_prop, 'Value',val);
end

function handles = update_modif_panel(handles)

Br = handles.working_sets.(handles.current_set);

%% Current system panel title
if isfield(Br, 'mdl_name')
    str_name = ['Model: ' Br.mdl_name];
else
  str_name =   Br.Sys.name;
end
if (numel(str_name)>35)
    str_name= [str_name(1:30) '...' str_name(end-5:end)];
end
set(handles.system_panel,'Title',str_name );

%% 
modif_panel_title = Br.disp();
set(handles.modif_param_panel,'Title', modif_panel_title);

%% parameters listbox
nb_pts = size( Br.P.pts,2);
if ~isfield(Br.P,'selected')
    Br.P.selected=zeros(1,nb_pts);
end

nb_pts = size(Br.P.pts,2);
DimX = Br.P.DimX;
DimP = Br.P.DimP;
ParamList = Br.P.ParamList;
dim = Br.P.dim;

%% fix selected field if first time in GUI 
if ~isfield(Br.P, 'selected')
    val = get(handles.working_sets_listbox, 'Value');
    names = fieldnames(handles.working_sets);
    set_name = names{val};
    handles.working_sets.(set_name).ResetSelected();
end

handles.current_pts = min(handles.current_pts, nb_pts);
k = handles.current_pts;

%% fill uitable
h_tb= handles.uitable_params;
current_pts = Br.GetParam(handles.show_params,k);   
domains = Br.GetDomain(handles.show_params);
idx = FindParam(Br.P, handles.show_params);
is_signal = find(idx<=DimX);

fill_uitable_params(h_tb, handles.show_params, current_pts, domains, is_signal); % last argument tells which is a signal 
set(h_tb, 'ColumnEditable', [false, true, true, true, true]);
set(h_tb, 'ColumnWidth', handles.TBL_SZ);

handles = update_selected_domain(handles);

%%  slider

tt=strcat('Sample  ',num2str(handles.current_pts),'/',num2str(nb_pts));
set(handles.text_pts, 'String', tt );
set(handles.slider_current_pts, 'Value',handles.current_pts,'Min', .9999, 'Max',nb_pts, 'SliderStep', [1/nb_pts 10/nb_pts]);
set(handles.checkbox_selected, 'Value', Br.P.selected(handles.current_pts));

%% button plot all/selected

if isempty(find(Br.P.selected, 1))
    set(handles.button_plot_selected, 'String', 'Plot traj(s)')
else
    set(handles.button_plot_selected, 'String', 'Plot selected')
end

%% menu for param pts plot

% plot parameter that are visible in the table and bounded
bnd_params = Br.GetBoundedDomains();
show_params = handles.show_params;
params_to_plot = intersect(bnd_params, show_params);
dim =FindParam(Br.P, params_to_plot);

if isempty(dim)
    dim_show = FindParam(Br.P, show_params);
    if isempty(dim_show)
        dim=1;
    else
        idim = find(dim_show>Br.P.DimX);
        if isempty(idim)
            dim = dim_show(1);
        else
            dim = dim_show(idim);
        end
    end
end

set(handles.popup_pts1,'String',ParamList);
handles.current_plot_pts{1} = ParamList(dim(1));
set(handles.popup_pts1,'Value', dim(1));
set(handles.popup_pts2,'String',{'', ParamList{:}});

if(numel(dim)>=2)
    handles.current_plot_pts{2} = ParamList(dim(2));
    set(handles.popup_pts2,'Value', dim(2)+1);
else
    handles.current_plot_pts{2} = '';
    set(handles.popup_pts2,'Value', 1);
end

set(handles.popup_pts3,'String',{'', ParamList{:}});
if(numel(dim)>=3)
    handles.current_plot_pts{3} = ParamList(dim(3));
    set(handles.popup_pts3,'Value', dim(3)+1);
else
    handles.current_plot_pts{3} = '';
    set(handles.popup_pts3,'Value', 1);
end

st_sample= get_sample_string(handles);
handles = info(handles, st_sample);

handles =  plot_pts(handles);

%% Plot function
function handles= plot_pts(handles)
Br = handles.working_sets.(handles.current_set);

select_style = {'sk', 'MarkerSize', 20};

if ~isfield(Br.P,'selected')
    Br.P.selected = zeros(size(Br.P.pts,2))
end

axes(handles.axes_pts);
hold off;
cla;
legend off;
set(gca,'YLimMode','auto');

param_to_plot = handles.current_plot_pts{1};
ipts = [handles.current_pts find(Br.P.selected)];
ipts = unique(sort(ipts));

% TODO (or not) if first param is a signal, plot traces (why not?)

% if FindParam(Br.P,  param_to_plot)<= Br.P.DimX
%     try
%         Br.PlotSignals(param_to_plot, ipts);
% %        Br.GetSig
%         return;
%     end 
% end

if ~strcmp(param_to_plot,'')
    if ~strcmp(handles.current_plot_pts{2},'')
        param_to_plot = {param_to_plot{:} handles.current_plot_pts{2}};
    end
    
    if ~strcmp(handles.current_plot_pts{3},'')
        param_to_plot = {param_to_plot{:} handles.current_plot_pts{3}};
    end
end

i0=0;

% Plot domain
try 
    Br.PlotDomain(param_to_plot);
end

if ~isfield(Br.P, 'props_values')
    Br.PlotParams(param_to_plot)
    Br.PlotParams(param_to_plot,handles.current_pts,{'ok', 'MarkerSize',14});
    
    %S = DiscrimPropValues(Br.P);
    %SplotPts(S,param_to_plot);
else
    for i=1:numel(Br.P.props_names)
        if strcmp(Br.P.props_names{i}, handles.current_prop)
            i0 = i;
        end
    end
    if i0
        val = cat(1, Br.P.props_values(i0,:).val);
        val = val(:,1);
        iparam = FindParam(Br.P,param_to_plot);
        switch (numel(iparam))
            
            case 1
                x = Br.P.pts(iparam(1),:);
                y = zeros(size(x));
                scatter(x,y, 30, val, 'filled');
                
            case 2
                x = Br.P.pts(iparam(1),:);
                y = Br.P.pts(iparam(2),:);
                scatter(x,y, 30, val, 'filled');
                
            case 3
                x = Br.P.pts(iparam(1),:);
                y = Br.P.pts(iparam(2),:);
                z = Br.P.pts(iparam(3),:);
                scatter3(x,y, z, 30, val, 'filled');
        end
        
        stats = sprintf(['val:' dbl2str(val(handles.current_pts)) ...
            '\n#True:' dbl2str(numel(find(val>0))) '/' dbl2str(numel(val))  ...
%            '\nMean:' dbl2str(mean(val)) ...
%          '\nstd:' dbl2str(std(val)) ...
%        '\nMax:' dbl2str(max(val)) ...
%        '\nMin:' dbl2str(min(val))
          ]);
        
        legend(stats);
        clim = sym_clim(val);
        if diff(clim)>0
            set(gca, 'CLim',clim );
        end
        colormap([ 1 0 0; 0 1 0 ]);
        hold on;
        SplotPts(Br.P, param_to_plot,handles.current_pts,{'ok', 'MarkerSize',14});
    else
        SplotPts(Br.P, param_to_plot);
        SplotPts(Br.P, param_to_plot,handles.current_pts,{'ok', 'MarkerSize',14});
        S = DiscrimPropValues(Br.P);
        SplotPts(S,param_to_plot);
    end
    % Plot selected
    
    if any(Br.P.selected)
        hold on;
        kn = find(Br.P.selected);
        SplotPts(Br.P, param_to_plot,kn,select_style);
    end    
    
    
    
end

% --------------------------------------------------------------------
function menu_traj_and_sensi_Callback(hObject, eventdata, handles)
% hObject    handle to menu_traj_and_sensi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_compute_traj_Callback(hObject, eventdata, handles)
% hObject    handle to menu_compute_traj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);

tspan = inputdlg('Enter tspan (Format: [ti tf] or [t0 t1 t2 ... tn] or ti:dt:tf)','Compute trajectories', 1, {handles.last_tspan});
if isempty(tspan)
    return;
end
handles.last_tspan = tspan{1};
tspan = eval(tspan{1});

if (exist([pwd filesep 'cvm'])~=3)
    CompileSystem(Br.Sys);
end

uipause;
set(handles.text_info, 'Simulation running ...');
guidata(hObject, handles);

Br.P = ComputeTraj(Br.Sys, Br.P, tspan);
handles = update_working_sets_panel(handles);
uiresume;
set(handles.text_info, 'Simulation done.');
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_purge_traj_Callback(hObject, eventdata, handles)
% hObject    handle to menu_purge_traj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);

Br.P = ...
    SPurge(Br.P);
handles = update_working_sets_panel(handles);
handles = update_modif_panel(handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function menu_refine_prop_bound_Callback(hObject, eventdata, handles)
% hObject    handle to menu_refine_prop_bound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);

restore_traj =0;
if isfield( Br.Sys, 'type')
    if strcmp(Br.Sys.type,'traces')
        if isfield(Br.P, 'traj')
            traj(1) = Br.P.traj{1};
            restore_traj =1;
        end
    end
end

%  options dialog box
info = 'Refine boundary between regions satisfying different properties';
answers = { handles.current_prop,'','0','3'};
ins = inputdlg({'Enter array of properties (e.g. [phi1, phi2])', 'Enter tspan for trajectories (empty: tspan of computed trajs)','Enter time for property checking', 'Number of iterations'}, info,1, answers );

S = Br.P;
if isempty(ins)
    return;
end

load(handles.properties_file_name);

if isempty(ins{1})
    prop_name = handles.current_prop;
    prop = handles.properties.(prop_name);
else
    prop = eval(ins{1});
end

if isempty(ins{2})
    tspan = S.traj{1}.time;
else
    tspan = eval(ins{2});
end

tprop = eval(ins{3});
nb_iter = eval(ins{4});

S = rmfield(S,'selected');
Sf = SFindPropBoundary(Br.Sys,S, prop, tspan,tprop, nb_iter);

% TOFIX: hasty patch
Sf = rmfield(Sf,'selected');

Br.P = Sf;
handles = update_modif_panel(handles);
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_plot_property_val_Callback(hObject, eventdata, handles)
% hObject    handle to menu_plot_property_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);

try
    close(handles.figp)
end
handles.figp = figure;
hold on;
num_params = numel(handles.selected_params);
Br.PlotRobustMap(handles.properties.(handles.current_prop), handles.selected_params{1:max(num_params,2)});
guidata(hObject,handles);

% --------------------------------------------------------------------
function menu_falsify_Callback(hObject, eventdata, handles)
% hObject    handle to menu_falsify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);

kn = find(Br.P.selected);
if isempty(kn)
    kn = handles.current_pts;
end
Psel = Sselect( Br.P, kn);
prop_name = handles.current_prop;
prop = handles.properties.(prop_name);

BrSys = Br.copy();
BrSys.SetP(Psel);

if isfield(handles, 'falsif_pb')
    options = handles.falsif_pb.solver_options;
    handles.falsif_pb = FalsificationProblem(BrSys, prop);
else
    handles.falsif_pb = FalsificationProblem(BrSys, prop);
    options = handles.falsif_pb.solver_options;
end

options = structdlg(options);
if isempty(options)
    return;
end
handles.falsif_pb.solver_options = options;
handles.falsif_pb.log_traces=1;

handles = info(handles, 'Optimizing...');

handles.falsif_pb.solve();

BrFalse = handles.falsif_pb.GetBrSet_False();
BrLogged = handles.falsif_pb.GetBrSet_Logged();

if isempty(BrFalse)
    handles = info(handles, 'No falsifying input found');
else
    handles = info(handles, 'Done.');
    new_name = ['BrFalse_' get_id(prop)];
    handles.working_sets = setfield(handles.working_sets, new_name, BrFalse);
    assignin('base', new_name, BrFalse);
end
new_name2 = ['BrLog_' get_id(prop)];
handles.working_sets = setfield(handles.working_sets, new_name2, BrLogged);
assignin('base', new_name, BrFalse);

handles = update_working_sets_panel(handles);
handles = update_modif_panel(handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function menu_select_Callback(hObject, eventdata, handles)
% hObject    handle to menu_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_select_prop_gt_Callback(hObject, eventdata, handles)
% hObject    handle to menu_select_prop_gt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);

try
    
    iprop = find_prop(handles.current_prop,Br.P );
    valst = inputdlg('greater than ?');
    if isempty(valst)
        return;
    end
    
    val_threshold = eval(valst{1});
    
    if iprop
        val = cat(1,Br.P.props_values(iprop,:).val);
        val = val(:,1);
        Br.P.selected = (val>=val_threshold)';
    end
    handles = update_modif_panel(handles);
    guidata(hObject,handles);
    
catch
    s = lasterror;
    warndlg(['Problem selection: ' s.message] );
    error(s);
    return
end

% --------------------------------------------------------------------
function menu_select_prop_abs_st_Callback(hObject, eventdata, handles)
% hObject    handle to menu_select_prop_abs_st (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);

try
    iprop = find_prop(handles.current_prop,Br.P);
    valst = inputdlg('smaller than ?');
    if isempty(valst)
        return;
    end
    val_threshold = eval(valst{1});
    if iprop
        val = cat(1,handles.current_prop,Br.P.props_values(iprop,:).val);
        val = val(:,1);
        Br.P.selected = (abs(val)<=val_threshold)';
    end
    
    handles = update_modif_panel(handles);
    guidata(hObject,handles);
catch
    s = lasterror;
    warndlg(['Problem during selection: ' s.message] );
    error(s);
    return
end


% --------------------------------------------------------------------
function menu_inverse_selection_Callback(hObject, eventdata, handles)
% hObject    handle to menu_inverse_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);

Br.P.selected = ~ Br.P.selected;
handles = update_modif_panel(handles);
guidata(hObject,handles);


% --------------------------------------------------------------------
function plot_zero_contour_Callback(hObject, eventdata, handles)
% hObject    handle to plot_zero_contour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);
try
    axes(handles.axes_pts);
    hold on;
    
    iprop = find_prop(handles.current_prop,  Br.P);
    if iprop
        switch numel(Pf.dim)
            case 2
                val = cat(1, Pf.props_values(iprop,:).val);
                Z = val(:,1);
                QuickContourSf(Pf,Z)
        end
    end
catch
    s = lasterror;
    warndlg(['Problem with plot_zero_contour: ' s.message] );
    error(s);
    return
end


% --------------------------------------------------------------------
function menu_unselect_Callback(hObject, eventdata, handles)
% hObject    handle to menu_unselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);
Br.P.selected = 0* ...
    Br.P.selected;
handles = update_modif_panel(handles);

guidata(hObject,handles);


% --------------------------------------------------------------------
function menu_select_range_st_Callback(hObject, eventdata, handles)
% hObject    handle to menu_select_range_st (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);
try
    P = Br.P;
    deltast = inputdlg('Enter range:');
    if isempty(deltast)
        return;
    end
    delta = eval(deltast{:});
    
    if numel(delta)==1
        delta = repmat(delta, [size(P.epsi,1) 1]);
    end
    
    Dp = repmat(delta,[1 size(P.epsi,2)]);
    kn = find(sum(P.epsi<=Dp,1));
    Br.P.selected(kn)=1;
    
    handles = update_modif_panel(handles);
    guidata(hObject,handles);
    
catch
    s = lasterror;
    warndlg(['Problem in selection: ' s.message] );
    error(s);
    return
end


% --- Executes on button press in button_new_prop.
function button_new_prop_Callback(hObject, eventdata, handles)
% hObject    handle to button_new_prop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);

PHI_= Br.AddSpecGUI();

if isa(PHI_,'STL_Formula')
    eval([get_id(PHI_) '=PHI_']);
    handles = update_properties_panel(handles);
    guidata(hObject, handles);
else
    info(handles,'No new formula was defined.');
end

function i = find_prop(st, P)
i=0;

if ~isfield(P, 'props_names')
    return;
else
    props_names =  P.props_names;
end

for k = 1:numel(props_names)
    if strcmp(st,props_names{k})
        i = k;
        return;
    end
end

% --- Executes on button press in button_del_property.
function button_del_property_Callback(hObject, eventdata, handles)
% hObject    handle to button_del_property (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    
    fn = fieldnames(handles.properties);
    
    if numel(fn)==1
        return
    end
    
    st = fn{get(handles.listbox_prop,'Value')};
    val = get(handles.listbox_prop,'Value');
    
    if(val>1)
        val = val-1;
        set(handles.listbox_prop,'Value', val);
        handles.current_prop = fn{val};
    else
        handles.current_prop = fn{val+1};
    end
    
    handles.properties = rmfield(handles.properties,st);
    
    Br = handles.working_sets.(handles.current_set);
    Br.P = SPurge_Props(Br.P);
    Br.Specs.remove(st);
    
    handles = update_properties_panel(handles);
    guidata(hObject, handles);
    
catch
    s = lasterror;
    warndlg(['Problem deleting property: ' s.message] );
    error(s);
    return
end

% --------------------------------------------------------------------
function menu_load_properties_Callback(hObject, eventdata, handles)
% hObject    handle to menu_load_properties (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Br = handles.working_sets.(handles.current_set);   
    [FileName,PathName,FilterIndex] = uigetfile('*.mat; *.stl','Load Properties Set...');
    
    if (FileName ==0)
        return
    end
    
    [~,~,ext] = fileparts(FileName);
    
    if strcmp(ext,'.mat')
        props = load(FileName);
        fnames = fieldnames(props);
        
        % find properties in the file loaded
        nprops = [];
        
        for j = 1:numel(fnames)
            if isa(props.(fnames{j}), 'STL_Formula')
                nprops = [ nprops props.(fnames{j}) ];
            end
        end
        
        if isempty(nprops)
            warndlg('No properties in this file');
            return
        else
            handles.properties = nprops;
        end
        handles.current_prop = fnames{1};
        handles.idx_prop = 1;
        
    elseif strcmp(ext,'.stl')
        [prop_names, props] = STL_ReadFile([PathName FileName]);
        for iprop = 1:numel(props) 
            Br.AddSpec(props{iprop});
        end
    end
    
    set(handles.listbox_prop, 'Value', 1);
    handles = update_properties_panel(handles);
    guidata(hObject,handles);
    
    

% --------------------------------------------------------------------
function menu_save_properties_Callback(hObject, eventdata, handles)
% hObject    handle to menu_save_properties (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    FileName = uiputfile('*.mat','Save Properties Set As...');
    
    if(FileName==0)
        return
    end
    
    handles.properties_file_name = FileName;
    ws = handles.properties;
    save(FileName, '-struct','ws');
    handles = update_properties_panel(handles);
    guidata(hObject,handles);
    
catch
    s = lasterror;
    warndlg(['Problem saving property: ' s.message] );
    error(s);
    return
end

% --------------------------------------------------------------------
function menu_param_sets_Callback(hObject, eventdata, handles)
% hObject    handle to menu_param_sets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_properties_Callback(hObject, eventdata, handles)
% hObject    handle to menu_properties (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in button_go_refine.
function button_go_refine_Callback(hObject, eventdata, handles)
% hObject    handle to button_go_refine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


Br = handles.working_sets.(handles.current_set);
try
    
    if(handles.refine_args==0)
        return;
    end
    
    restore_traj =0;
    if isfield( Br.Sys, 'type')
        if strcmp(Br.Sys.type,'traces')
            if isfield(Br.P, 'traj')
                traj = Br.P.traj;
                restore_traj =1;
            end
        end
    end
    
        
    ipts = find(Br.P.selected);
    if isempty(ipts)
        ipts = handles.current_pts;
    end
    
    Ptmp = Sselect(Br.P,ipts);
    if(handles.halton)
        Pr =QuasiRefine(Ptmp, handles.refine_args);
    else
        Pr = Refine(Ptmp, handles.refine_args);
    end
    nnb_pts = size(Pr.pts,2);
    Pr.selected = zeros(1, nnb_pts);
    
    nipts = find(~Br.P.selected);
    if(ipts == handles.current_pts)
        nipts = nipts(nipts~=handles.current_pts);
    end
    
    if isempty(nipts)
        Pf = Pr;
    else
        P = Sselect(Br.P, nipts);
        P = SPurge(P);
        P = SPurge_props(P);
        Pf = SConcat(Pr,P);
    end
    
    if(restore_traj)
        Pf.traj = repmat(traj(1), [1 size(Pf.pts,2)]);
        for j = 1:numel(Pf.traj)
            Pf.traj{j}.param = Pf.pts(:,j);
        end
    end
    
    Br.P = Pf;
       
    set(handles.edit_num_samples, 'String','');
    set(hObject,'String','Go (0 new samples)');
    handles.refine_args = 0 ;
    
    handles = update_working_sets_panel(handles);
    handles = update_modif_panel(handles);
    guidata(hObject,handles);
    
catch
    s = lasterror;
    set(handles.edit_num_samples, 'String','');
    set(hObject,'String','Go (0 new samples)');
    handles.refine_args = 0 ;
    handles = update_working_sets_panel(handles);
    handles = update_modif_panel(handles);
    warndlg(['Problem refining: ' s.message] );
end

% --- Executes on button press in button_explore_traj.
function button_explore_traj_Callback(hObject, eventdata, handles)
% hObject    handle to button_explore_traj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);
h = BreachTrajGui(Br, handles);

% --------------------------------------------------------------------
function System_Callback(hObject, eventdata, handles)
% hObject    handle to System (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button_save.
function button_save_Callback(hObject, eventdata, handles)
% hObject    handle to button_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ws = handles.working_sets; %#ok<NASGU>
save(handles.working_sets_file_name, '-struct', 'ws');


% --- Executes on button press in button_del_selected.
function button_del_selected_Callback(hObject, eventdata, handles)
% hObject    handle to button_del_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);

try
    
    % deals with the case where no dynamics is used (only traces)
    restore_traj = 0;
    if isfield(Br.Sys, 'type')
        if strcmp(Br.Sys.type,'traces')
            if isfield(Br.P, 'traj')
                traj = Br.P.traj;
                restore_traj =1;
            end
        end
    end
    
    nb_pts = size(Br.P.pts,2);
    if(nb_pts == 1) % cannot delete a unique sample
        return;
    end
    
    ipts = find(Br.P.selected);
    if isempty(ipts)
        return;
    end
    
    nipts = find(~Br.P.selected);
    if ipts == handles.current_pts
        nipts = nipts(nipts~=handles.current_pts);
    end
    
    P = Sselect(Br.P, nipts);
    
    Br.P = P;
    handles = update_modif_panel(handles);
    guidata(hObject,handles);
    
catch
    s = lasterror;
    warndlg(['Problem deleting selected: ' s.message] );
end


% --------------------------------------------------------------------
function menu_max_robust_satisfaction_Callback(hObject, eventdata, handles)
% hObject    handle to menu_max_robust_satisfaction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);

try
    
    P = Br.P;
    info = 'Find max robust satisfaction ';
    
    dim = P.dim;
    pts = P.pts(:,handles.current_pts);
    epsi = P.epsi(:,handles.current_pts);
    
    lbstr = '';
    ustr  = '';
    
    for i=1:numel(pts)
        lbstr = ['' num2str(pts(dim(i)+epsi(1,i)))];
        upstr = ['' num2str(pts(dim(i)+epsi(2,i)))];
    end
    
    answers = {'',lbstr,upstr};
    ins = inputdlg({'Enter tspan for trajectories (default: computed trajectory)', ...
        'Enter lower bounds', ...
        'Enter upper bounds'}, info,1, answers );
    if isempty(ins)
        return;
    end
    
    if isempty(ins{1})
        tspan = S.traj{1}.time;
    else
        tspan = eval(ins{1});
    end
    
    opt.lbound = eval(ins{2});
    opt.ubound = eval(ins{3});
    
    prop_name = handles.current_prop;
    prop = handles.properties.(prop_name);
    
    Sopt = SOptimProp(Br.Sys, P, prop , tspan, lbound, ubound);
    
    % add optimized parameter set
    
    val = get(handles.working_sets_listbox, 'Value');
    names = fieldnames(handles.working_sets);
    set_name = names{val};
    new_name = genvarname(set_name,names);
    handles.working_sets = setfield(handles.working_sets, new_name, Sopt);
    handles = update_modif_panel(handles);
    
catch
    s = lasterror;
    warndlg(['Problem deleting selected: ' s.message] );
    error(s);
end

guidata(hObject,handles);


% --------------------------------------------------------------------
function menu_remove_prop_Callback(hObject, eventdata, handles)
% hObject    handle to menu_remove_prop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);

Br.P = SPurge_props(Br.P);
handles = update_working_sets_panel(handles);
handles = update_properties_panel(handles);
handles = update_modif_panel(handles);
guidata(hObject,handles);


% --- Executes on button press in button_plot_prop.
function button_plot_prop_Callback(hObject, eventdata, handles)
% hObject    handle to button_plot_prop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


prop = handles.properties.(handles.current_prop);
title = 'Plot property options';
prompt = {'Enter tspan for trajectory computation (empty = existing traj), ',...
    'Enter tspan for property evaluation (empty = same as above)' ...
    'Levels of subformulas to plot'};
opt.tspan_traj = [];
opt.tspan_prop_eval = [];
opt.break_lev = inf;

%[opt, handles.check_prop_options] = structdlg(opt, title, prompt, handles.check_prop_options);

if isempty(opt)
    return;
end

Ptmp = Sselect(Br.P, handles.current_pts);

if ~isempty(opt.tspan_traj)
    Pphi = SPurge(Ptmp);
    Pphi = ComputeTraj(Br.Sys, Pphi, opt.tspan_traj);
else
    Pphi = Ptmp;
end

if isempty(opt.tspan_prop_eval)
    opt.tspan_prop_eval = Pphi.traj{Pphi.traj_ref(1)}.time;
end

handles = info(handles, 'Computing satisfaction of formula...');
Pphi = SEvalProp(Br.Sys, Pphi, prop, opt.tspan_prop_eval, 1, opt.break_lev);
PplotFormula(Br.Sys, Pphi, prop, 1, opt.break_lev);
handles = info(handles, 'Computing satisfaction of formula... Done.');
guidata(hObject,handles);


% --- Executes on button press in button_break_prop.
function button_break_prop_Callback(hObject, eventdata, handles)
% hObject    handle to button_break_prop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Br = handles.working_sets.(handles.current_set);
prop = handles.properties.(handles.current_prop);
props = STL_Break(prop);

for i = 1:numel(props)
    PHI_ = props(i);
    eval([get_id(PHI_) '=PHI_']);
    handles.properties.(get_id(props(i))) = props(i);
    Br.AddSpec(props(i));
end

handles = update_properties_panel(handles);

guidata(hObject, handles);


% --- Executes on button press in button_plot_selected.
function button_plot_selected_Callback(hObject, eventdata, handles)
% hObject    handle to button_plot_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);

try
    
    % getting parameters
    if isfield(Br.Sys, 'time_mult')
        Br.P.time_mult = Br.Sys.time_mult;
    end
    
    ipts = find(Br.P.selected);
    
    if isfield(Br.P, 'plot_args')
        args = Br.P.plot_args;
        args = plots_args_gui(Br.P, args);
    else
        args = plots_args_gui(Br.P, []);
    end
    
    Br.P.plot_args = args;
    
    if isempty(args)
        return;
    end
    
    h = figure;
    assignin('base','h_breach_', h);
    
    if(args.phase_portrait == 1)
        SplotTraj(Br.P, args.iX ,ipts ,[]);
    else
        
        SplotVar(Br.P, args.iX,ipts ,[], args.same_axe);
    end
    
catch
    s = lasterror;
    handles = info(handles,['Problem plotting trajectories: did you run simulation first?'] );
    guidata(hObject, handles);
end

return


% --- Executes on key press with focus on breach and none of its controls.
function breach_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to breach (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


function st = dbl2str(x)
st = num2str(x, '%0.5g');

% --------------------------------------------------------------------
function menu_select_satisfied_Callback(hObject, eventdata, handles)
% hObject    handle to menu_select_satisfied (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    Br = handles.working_sets.(handles.current_set);
    iprop = find_prop(handles.current_prop,Br.P);
    val_threshold = 0;
    if iprop
        val = cat(1,Br.P.props_values(iprop,:).val);
        val = val(:,1);
        Br.P.selected = (val>=val_threshold)';
    end
    
    handles = update_modif_panel(handles);
    guidata(hObject,handles);
    
function edit_param_name_Callback(hObject, eventdata, handles)
% hObject    handle to edit_param_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function edit_param_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_param_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function h = info(h,msg)
% INFO write the message into the information panel.
set(h.text_info, 'String', msg);
drawnow();


% --- Executes during object creation, after setting all properties.
function system_panel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to system_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --------------------------------------------------------------------
function create_input_parameter_set_Callback(hObject, eventdata, handles)
% hObject    handle to create_input_parameter_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Br = handles.working_sets.(handles.current_set);
idx = Br.GetParamsInputIdx();
Pnew = CreateParamSet(Br.Sys, idx);
names = fieldnames(handles.working_sets);
new_name = genvarname('Pin',names);
handles.working_sets = setfield(handles.working_sets, new_name, Pnew);
%set(handles.working_sets_listbox,'Value', numel(names)+1);
handles = update_working_sets_panel(handles);
guidata(hObject,handles);


% --------------------------------------------------------------------
function menu_param_synth_Callback(hObject, eventdata, handles)
% hObject    handle to menu_param_synth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);

kn = find(Br.P.selected);
if isempty(kn)
    kn = handles.current_pts;
end
Psel = Sselect( Br.P, kn);
prop_name = handles.current_prop;
prop = handles.properties.(prop_name);

BrSys = Br.copy();
BrSys.SetP(Psel);

handles.param_synth_pb = ParamSynthProblem(BrSys, prop);

handles = info(handles, 'Optimizing...');

handles.param_synth_pb.solve();

BrSynth = handles.param_synth_pb.GetBrSet_Best();
BrLogged = handles.param_synth_pb.GetBrSet_Logged();

% TODO what if no synth param exist?
names = fieldnames(handles.working_sets);
handles = info(handles, 'Done.');
PSynth = BrSynth.P;
PSynth.epsi = 0*PSynth.epsi;
new_name = genvarname(['PSynth_' get_id(prop)],names);
handles.working_sets = setfield(handles.working_sets, new_name, PSynth);

PLog =  BrLogged.P;
PLog.epsi = 0*PLog.epsi;
new_name2 = genvarname(['PLog_' get_id(prop)],names);
handles.working_sets = setfield(handles.working_sets, new_name2, PLog);

handles = update_working_sets_panel(handles);
handles = update_modif_panel(handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function menu_reset_param_files_Callback(hObject, eventdata, handles)
% hObject    handle to menu_reset_param_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);

filename = [Br.Sys.name '_param_sets.mat'];
backup = ['.' Br.Sys.name '_param_sets.mat.bck'];
title = ['Reset param set file'];
answ = questdlg('Are you sure?', title);

if strcmp(answ,'Yes')
    movefile(filename, backup);
    Br.RunGUI;
end


% --------------------------------------------------------------------
function menu_reset_prop_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_reset_prop_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Br = handles.working_sets.(handles.current_set);
filename = [Br.Sys.Dir filesep Br.Sys.name '_properties.mat'];
backup = [Br.Sys.Dir filesep '.' Br.Sys.name '_properties.mat.bck'];
title = ['Reset prop file'];
answ = questdlg('Are you sure?', title);

if strcmp(answ,'Yes')
    movefile(filename, backup);
    Br.RunGUI;
end


% --- Executes on button press in pushbutton_input_gen.
function pushbutton_input_gen_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_input_gen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Br = handles.working_sets.(handles.current_set);
Br.SetInputGenGUI;


% --- Executes on button press in checkbox_input_only.
function checkbox_input_only_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_input_only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);

Br.ResetSimulations();
Br.P = SPurge(Br.P);
Br.SimInputsOnly = get(hObject,'Value');
guidata(hObject,handles);


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
Br = handles.working_sets.(handles.current_set);

[params, p0, domains] = read_uitable_params(hObject); 
idx_params = find(strcmp(params, handles.selected_params) );
Br.SetDomain( handles.selected_params,domains(idx_params));

idx = FindParam(Br.P, handles.selected_params);
if isempty(domains(idx_params).domain)  
    Br.P.pts(idx, : ) = p0(idx_params);
else
    Br.P.pts(idx, handles.current_pts ) = p0(idx_params);
end
handles = update_modif_panel(handles);
guidata(hObject,handles);




% --- Executes on button press in button_all.
function button_all_Callback(hObject, eventdata, handles)
% hObject    handle to button_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.show_params = handles.working_sets.(handles.current_set).P.ParamList;
handles = update_modif_panel(handles);
guidata(hObject,handles);


% --- Executes on button press in button_sys_params.
function button_sys_params_Callback(hObject, eventdata, handles)
% hObject    handle to button_sys_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);
handles.show_params = Br.GetParamsSysList();

handles = update_modif_panel(handles);
guidata(hObject,handles);



% --- Executes on button press in button_inputs.
function button_inputs_Callback(hObject, eventdata, handles)
% hObject    handle to button_inputs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);
idx = Br.GetParamsInputIdx();
handles.show_params = Br.P.ParamList(idx);

handles = update_modif_panel(handles);
guidata(hObject,handles);


% --- Executes on button press in button_prop_param.
function button_prop_param_Callback(hObject, eventdata, handles)
% hObject    handle to button_prop_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);
handles.show_params = Br.GetPropParamList();

handles = update_modif_panel(handles);
guidata(hObject,handles);



function edit_filter_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Br = handles.working_sets.(handles.current_set);
st = get(hObject,'String'); 
if ~isempty(st)
   handles.show_params = Br.P.ParamList(cellfun(@(c)(~isempty(c)),  regexp( Br.P.ParamList,st)));
end

handles = update_modif_panel(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_filter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_domain.
function button_domain_Callback(hObject, eventdata, handles)
% hObject    handle to button_domain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Br = handles.working_sets.(handles.current_set);
handles.show_params = Br.GetBoundedDomains();

handles = update_modif_panel(handles);
guidata(hObject,handles);




% --- Executes on selection change in popup_sample_option.
function popup_sample_option_Callback(hObject, eventdata, handles)
% hObject    handle to popup_sample_option (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 val = get(hObject,'Value'); 
     switch(val)
         case 1
             handles.sample_arg_multi = 'replace';
         case 2
             handles.sample_arg_multi = 'append';
         case 3
             handles.sample_arg_multi = 'combine';
     end
     
     handles = info(handles, st_info);
     
guidata(hObject,handles);


% --- Executes on selection change in popup_method.
function popup_method_Callback(hObject, eventdata, handles)
% hObject    handle to popup_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

     val = get(hObject,'Value'); 
     switch(val)
         case 1
             handles.sample_arg_method = 'grid';
         case 2
             handles.sample_arg_method = 'corners';
         case 3
             handles.sample_arg_method = 'rand';
         case 4
             handles.sample_arg_method = 'quasi';
     end
     
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popup_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected cell(s) is changed in uitable_params.
function uitable_params_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable_params (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

idx = eventdata.Indices;
if ~isempty(idx)
    handles.selected_params = handles.show_params(idx(:,1));
    st_sample= get_sample_string(handles);
    handles = info(handles, st_sample);
    guidata(hObject,handles);
end

function st_sample = get_sample_string(handles)
st_sample = ['Sample domain {' get_domain_string(handles) '}'];
  

function st_dom = get_domain_string(handles)
st_dom =  cell2mat(cellfun(@(c) ( [c ' ' ] ) , handles.selected_params, 'UniformOutput', false ));


% --- Executes on button press in button_reset.
function button_reset_Callback(hObject, eventdata, handles)
% hObject    handle to button_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = handles.working_sets.(handles.current_set);
Br.ResetParamSet();
handles = update_working_sets_panel(handles);
handles = update_modif_panel(handles);
guidata(hObject,handles);


