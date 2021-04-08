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

% Last Modified by GUIDE v2.5 19-Jan-2018 12:50:58

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

handles = info(handles, 'Starting Breach....');

% Set fonts and size depending on system
if ismac
    FONT=12;
    %POS = [50 10 200 50];
    handles.TBL_SZ = {200 120 120 150 80} ;
else
    FONT=10;
    %POS = [50 10 200 50];
    handles.TBL_SZ = {300 120 150 200 120} ;
end

hfn = fieldnames(handles);
for ifn = 1:numel(hfn)
    try
        set(handles.(hfn{ifn}), 'FontSize', FONT);
    end
end
% set(handles.breach, 'Position',POS);

crd = pwd;
set(hObject, 'Name', ['Breach (' crd  ')']);

%% used to memorize previous entries
handles.last_options =[];
handles.check_prop_options =[];
handles.current_set = [];

%% Init stuff
if numel(varargin)>=1
    BrGUI = varargin{1};
else
    BrGUI = [];
end

%% Init properties panel
handles.idx_prop= 1;
handles.current_prop = '';
handles.properties = struct;

%% Init modif panel
handles.current_pts = 1;
handles.refine = 1;
handles.refine_all = 0;
handles.halton = 0;
handles.refine_args = 0;
handles.select_cells = [];
handles.current_plot_pts = {};
handles.working_sets = struct;

% Init param pts plot
handles.current_plot{1} =[];
handles.current_plot{2} =[];
handles.current_plot{3} =[];
handles.current_marked = [];
handles.figp=[];

% Find out who's in the workspace
[handles, BrGUI] = get_param_sets(handles,BrGUI);
handles = update_all(handles, BrGUI);

handles = info(handles, 'Ready.');

% Choose default command line output for BreachGui
handles.output = handles.breach;

% Update handles structure
guidata(hObject, handles);

function handles = update_all(handles, BrGUI)

    fnames = fieldnames(handles.working_sets);
    if ~isempty(fnames)
        if isempty(BrGUI)
            BrGUI = handles.current_set;
        end
        
        igui = find(strcmp(fnames, handles.current_set));
        set(handles.working_sets_listbox, 'Value', igui);
        
        Sys = BrGUI.Sys;
        if (isfield(Sys,'tspan'))
            set( handles.edit_time, 'String', get_time_string(Sys.tspan));
        else
            set( handles.edit_time, 'String', get_time_string(0:.01:1));
        end
        
        handles.show_params = BrGUI.P.ParamList;
        
        handles = update_working_sets_panel(handles);
        handles = update_modif_panel(handles);
        handles = update_properties_panel(handles);
        handles = set_default_plot(handles);
        handles = plot_pts(handles);
        
   end

function [handles, BrGUI] = get_param_sets(handles, BrGUI)
% find all param sets in workspace
ws_var = evalin('base', 'who');
for iv= 1:numel(ws_var)
    % is this a BreachSet?
    BB__ = evalin('base', ws_var{iv});
    if isa(BB__, 'BreachSet')&&(~isequal(ws_var{iv}, 'ans')) % found one, keep it, exclude ans
        handles.working_sets.(ws_var{iv}) = BB__;
        if isempty(BrGUI)
            BrGUI = BB__;
            handles.current_set = ws_var{iv};
        else
            if BrGUI== BB__ % this is the caller
                handles.current_set = ws_var{iv};
            end
        end
    elseif isa(BB__,'BreachProblem')
        
%         set_name  = [ws_var{iv} '__BrSet'];
%         handles.working_sets.(set_name) = BB__.BrSet;
%         if isempty(BrGUI)
%             BrGUI = BB__.BrSet;
%             handles.current_set = set_name;
%         end
        
        if isprop(BB__, 'BrSet_Best')&& ~isempty(BB__.BrSet_Best)
            set_name  = [ws_var{iv} '__Best'];
            handles.working_sets.(set_name) = BB__.GetBrSet_Best;
            if isempty(BrGUI)
                BrGUI = BB__.GetBrSet_Best;
                handles.current_set = set_name;
            end
        end
        
        if isprop(BB__, 'BrSet_Logged')&& ~isempty(BB__.BrSet_Logged)
            set_name  = [ws_var{iv} '__Logged'];
            handles.working_sets.(set_name) = BB__.BrSet_Logged;
            if isempty(BrGUI)
                BrGUI = BB__.BrSet_Logged;
                handles.current_set = set_name;
            end
        end
        
    end
end

function h_scat = get_scatter_handle()
% might get useless when I implement a new class for updatable plots..
ch = get(gca,'Children');
h_scat = [];
for ic = 1:numel(ch)
    if  isa(ch(ic),'matlab.graphics.chart.primitive.Scatter') || isa(ch(ic),'matlab.graphics.chart.primitive.Line')
        h_scat = ch(ic);
        return
    end
end

function brush_callback(handles,e,d)
set_selected_from_brush(handles);
set_brush_style(e,d);

function set_brush_style(~,~)
% capture brushed point and customize  the markers
h_scat = get_scatter_handle();
if ~isempty(h_scat)
    hB = h_scat.BrushHandles;  % handles for brush
    hM  = hB.Children;  % handle for brush markers
    if isempty(hM) % try waiting to get the marker...
        
        pause(0.01);
        hM  = hB.Children;  % handle for brush markers
        if isempty(hM)
            return % give up - hopefully nothing bad happened
        end
    end
end
hM(1).Style= 'square';
hM(1).FaceColorData(1:4) = [0 0 0 64]';
hM(1).Size = 14;

function handles = update_selected_domain(handles)
Br = get_current_set(handles);
if ~isempty(Br)
    params = Br.GetBoundedDomains();
    handles.selected_params = intersect(handles.show_params, params);
    if isempty(handles.selected_params)
        if numel(handles.show_params) >=1
            handles.selected_params = handles.show_params(1);
        end
    end
end

% --- Outputs from this function are returned to the command line.
function varargout = BreachGui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if ~isempty(handles)
    varargout{1} = handles.output;
end

% --- Executes on button press in button_remove_set.
function button_remove_set_Callback(hObject, eventdata, handles)
    Br = get_current_set(handles);
    if ~isempty(Br)
        
        old_name = handles.current_set;
        fn = fieldnames(handles.working_sets);
        if numel(fn)>1
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
            
            
            handles.show_params = Br.P.ParamList;
            
            handles =update_working_sets_panel(handles);
            handles= update_properties_panel(handles);
            handles =update_modif_panel(handles);
            handles = set_default_plot(handles);
            handles = plot_pts(handles);
            guidata(hObject,handles);
        end
    end

% --- Executes on selection change in working_sets_listbox.
function working_sets_listbox_Callback(hObject, eventdata, handles)
contents = get(hObject,'String');
if ~isempty(contents)
    fn = fieldnames(handles.working_sets);
    set_name = fn{get(hObject,'Value')};
    handles.current_set = set_name;
    handles.current_pts = 1;
    Br = get_current_set(handles);
    if ~isempty(Br)
        handles.show_params = Br.P.ParamList;
        
        handles = update_working_sets_panel(handles);
        handles= update_properties_panel(handles);
        handles = update_modif_panel(handles);
        
        handles = set_default_plot(handles);
        handles = plot_pts(handles);
        guidata(hObject, handles);
    end
end

function handles = set_default_plot(handles)
Br = get_current_set(handles);
if ~isempty(Br)
    handles.current_plot_pts = Br.GetBoundedDomains();
    if isempty(handles.current_plot_pts)&&Br.P.DimP>Br.P.DimX % not like I'm happy with that
        handles.current_plot_pts = Br.P.ParamList{Br.P.DimX+1};
    end
end

% --- Executes during object creation, after setting all properties.
function working_sets_listbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in button_compute_traj.
function button_compute_traj_Callback(hObject, eventdata, handles)

Br = get_current_set(handles);
if ~isempty(Br)
    time_str  = get(handles.edit_time, 'String');
    tspan = eval(time_str);
    if isempty(tspan)
        handles = info(handles, 'Enter a valid simulation time, e.g., 0:.01:10');
        return;
    end
    
    handles= info(handles,'Computing trajectories...');
    Br.Sim(tspan);
    handles= info(handles,'Computing trajectories... Done.');
    
    handles = update_working_sets_panel(handles);
    handles = update_modif_panel(handles);
    guidata(hObject, handles);
    h = BreachTrajGui(Br, handles);
end

function edit_num_samples_Callback(hObject, eventdata, handles)
% hObject    handle to edit_num_samples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Br = get_current_set(handles);
if ~isempty(Br)
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
    handles = update_sample_args(handles);
    guidata(hObject,handles);
end

function handles  = run_sample_domain(handles)
Br = get_current_set(handles);
if ~isempty(Br)
    handles = update_sample_args(handles);
    
    Br.SampleDomain(handles.selected_params, ...
        handles.sample_arg_num_samples,...
        handles.sample_arg_method,...
        handles.sample_arg_multi);
end

function handles = update_sample_args(handles)
Br = get_current_set(handles);
if ~isempty(Br)
    
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
            handles.sample_arg_method = 'quasi-random';
    end
    
    % checks if any selected param is actually a signal
    params = handles.selected_params;
    isSig = Br.isSignal(params);
    if any(Br.isSignal(params))
        isSig = params(isSig>0);
        handles = info(handles,  ['Cannot sample ' isSig{1} ', as this is a signal.']);
        set(handles.button_sample, 'Enable', 'off');
        return;
    end
    
    % checks whether everybody is a variables
    domains = Br.GetDomain(params);
    variables = {};
    for id = 1:numel(domains)
        if  ~isempty(domains(id).domain)
            variables = [variables params(id)];
        end
    end
    
    params = setdiff(params,variables);
    if (~isempty(params))
        set(handles.button_sample, 'Enable', 'off');
        msg = ['Parameter ' params{1} ' is not a variable (empty domain), cannot sample.'];
        handles = info(handles,  msg);
        return
    end
    
    num_args = get(handles.edit_num_samples, 'String');
    if ~isempty(num_args)
        set( handles.button_sample, 'Enable', 'on');
        msg = ['Ready to sample domain of { ' get_domain_string(handles) '}'];
        handles = info(handles,  msg);
    else
        set( handles.button_sample, 'Enable', 'off');
        msg = ['Enter a  number of samples. It can be a scalar: int or the keyword ''all'' or an array or cell of scalar and keyword.' ];
        handles = info(handles,  msg);
    end
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

Br = get_current_set(handles);
if ~isempty(Br)
    names = fieldnames(handles.working_sets);
    new_name =  genvarname('Br',names);
    Br_new = Br.copy();
    Br_new.ResetParamSet();
    Br_new.ResetSelected();
    handles.working_sets = setfield(handles.working_sets, new_name, Br_new);
    assignin('base', new_name, handles.working_sets.(new_name));
    
    handles = update_working_sets_panel(handles);
    handles = update_properties_panel(handles);
    handles = update_modif_panel(handles);
    
    guidata(hObject,handles);
end

% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_load_working_set_Callback(hObject, eventdata, handles)
% hObject    handle to menu_load_working_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

load_breachset(hObject,handles);

function load_breachset(hObject, handles)

    [FileName,PathName] = uigetfile('*.mat','Load Parameter Set...');
    if(FileName==0)
        return;
    end
    
    handles.working_sets_file_name = [PathName, FileName];
    Br = get_current_set(handles);
    %handles.working_sets = evalin('base', ['load(''' handles.working_sets_file_name ''')']);
    evalin('base', ['load(''' handles.working_sets_file_name ''')']);
    handles = get_param_sets(handles, Br);
    handles = update_all(handles,  Br);
    
    guidata(hObject,handles);
   
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
    BreachSave(handles.working_sets_file_name);
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
Br = get_current_set(handles);
if ~isempty(Br)
    val = get(handles.working_sets_listbox, 'Value');
    names = fieldnames(handles.working_sets);
    set_name = names{val};
    new_name = genvarname(set_name,names);
    Brcopy = Br.copy();
    handles.working_sets = setfield(handles.working_sets, new_name, Brcopy);
    assignin('base', new_name, Brcopy);
    handles = update_working_sets_panel(handles);
    guidata(hObject,handles);
end

function edit_rename_Callback(hObject, eventdata, handles)
Br = get_current_set(handles);
if ~isempty(Br)
    
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
    assignin('base', new_name, Br);
    guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function edit_rename_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function menu_files_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function menu_save_Callback(hObject, eventdata, handles)
ws = handles.working_sets;
try
    handles = info(handles, ['Saving parameter set to ' handles.working_sets_file_name '...']);
    BreachSave(handles.working_sets_file_name);
catch
    [FileName,PathName] = uiputfile('*.mat','Save Parameter Set As...');
    if(FileName==0)
        return;
    end
    handles.working_sets_file_name = [PathName  FileName];
    handles = info(handles, ['Saving parameter set to ' handles.working_sets_file_name '...']);
    BreachSave(handles.working_sets_file_name);
    handles = update_working_sets_panel(handles);
    
end
handles = info(handles, ['Saving parameter sets to ' handles.working_sets_file_name '... Done']);
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_copy_selected_Callback(hObject, eventdata, handles)
% hObject    handle to menu_copy_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Br = get_current_set(handles);
if ~isempty(Br)
    
    names = fieldnames(handles.working_sets);
    new_name = genvarname(handles.current_set,names);
    
    ipts = find(Br.P.selected);
    if isempty(ipts)
        return
    end
    Br_new = Br.copy();
    Br_new.P = Sselect(Br.P, ipts);
    
    handles.working_sets.(new_name) =  Br_new;
    assignin('base', new_name, Br_new);
    handles = update_working_sets_panel(handles);
    guidata(hObject,handles);
end

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
Br = get_current_set(handles);
if ~isempty(Br)&&~isempty(handles.current_prop) 
    if ~Br.hasTraj
        handles = info(handles, 'No traces to check requirement on. Run simulations first.');
        guidata(hObject,handles);
        return;
    end 
    prop = handles.properties.(handles.current_prop);  
    handles = info(handles, 'Computing satisfaction of formula...');
    Br.CheckSpec(prop);
    Br.SortbyRob();
    Br.SortbySat();
    handles = info(handles, 'Computing satisfaction of formula... Done.');
    
    handles = update_working_sets_panel(handles);
    %handles = update_modif_panel(handles);
    plot_pts(handles);
    handles = update_properties_panel(handles);
    guidata(hObject,handles);
end

% --- Executes on button press in button_edit_prop.
function button_edit_prop_Callback(hObject, eventdata, handles)
% hObject    handle to button_edit_prop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = get_current_set(handles);

if ~isempty(Br)
    
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
    
    id = get_id(prop);
    phi_st = str{1};
    PHI = eval(['STL_Formula(''' id ''',''' phi_st ''')']);
    eval([get_id(PHI) '=PHI']);
    Br.AddSpec(PHI);
    
    handles = update_properties_panel(handles);
    guidata(hObject, handles);
end

function handles = update_working_sets_panel(handles)

%% Set title
str_name = sprintf('Workspace');
if(numel(str_name)>35)
    str_name = [str_name(1:26) '...' str_name(end-5:end)];
end
set(handles.working_sets_panel,'Title',str_name );

% Filter out invalid entries
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
Br = get_current_set(handles);
if ~isempty(Br)
    
    str_name = ['Requirements of ' handles.current_set ];
    if (numel(str_name)>35)
        str_name= [str_name(1:30) '...' str_name(end-5:end)];
    end
    set(handles.panel_properties,'Title',str_name );
    
    handles.properties = struct;
    specs = Br.Specs;
    specs_names = specs.keys();
    
    for iprop = 1:numel(specs_names)
        this_name= specs_names{iprop};
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
        handles.current_prop = fnames{val};
    end
end

function handles= fill_uitable(handles)
Br = get_current_set(handles);
if ~isempty(Br)
    
    k = handles.current_pts;
    DimX = Br.P.DimX;
    
    h_tb= handles.uitable_params;
    % filter out invalid parameters (occurs when changing BreachSet to
    % another with different input generator for example
    handles.show_params = intersect(handles.show_params, Br.P.ParamList,'stable'); 
    
    current_pts = Br.GetParam(handles.show_params,k);
    domains = Br.GetDomain(handles.show_params);
    idx = FindParam(Br.P, handles.show_params);
    is_signal = find(idx<=DimX);
    
    fill_uitable_params(h_tb, handles.show_params, current_pts, domains, is_signal); % last argument tells which is a signal
    set(h_tb, 'ColumnEditable', [false, true, true, true, true]);
    set(h_tb, 'ColumnWidth', handles.TBL_SZ);
end

function time_string = get_time_string(time)

if isscalar(time)
    time_string = ['[0 ' dbl2str(time) ']'];
elseif numel(time)==2
    time_string = ['[' dbl2str(time(1)) ' ' dbl2str(time(2)) ']'];
elseif max(diff(diff(time)))<100*eps
    time_string = ['0:' dbl2str(time(2)-time(1)) ':' dbl2str(time(end))];
else
    time_string = num2str(time);
end

function handles = update_modif_panel(handles)

Br = get_current_set(handles);
if ~isempty(Br)
    
    %% Parallel checkbox
    set(handles.button_parallel, 'Value', Br.use_parallel~=0);
    
    %% Title
    modif_panel_title = Br.disp();
    set(handles.modif_param_panel,'Title', modif_panel_title);
    
    %% parameters listbox
    nb_pts = size( Br.P.pts,2);
    if ~isfield(Br.P,'selected')
        Br.P.selected=zeros(1,nb_pts);
    end
    nb_pts = size(Br.P.pts,2);
    
    %% fix selected field if first time in GUI
    if ~isfield(Br.P, 'selected')
        val = get(handles.working_sets_listbox, 'Value');
        names = fieldnames(handles.working_sets);
        set_name = names{val};
        handles.working_sets.(set_name).ResetSelected();
    end
    
    handles.current_pts = min(handles.current_pts, nb_pts);
    
    %% fill uitable
    handles = fill_uitable(handles);
    handles = update_selected_domain(handles);
    st_sample= get_sample_string(handles);
    handles = info(handles, st_sample);
    set(handles.breach,'UserData', Br);
end

%% Plot function
function handles= plot_pts(handles)
Br = get_current_set(handles);
if ~isempty(Br)
    if ~isfield(Br.P,'selected')
        Br.P.selected = zeros(size(Br.P.pts,2));
    end
    
    %% Determines projection
    params_to_plot = handles.current_plot_pts;
    if isempty(params_to_plot);
        return;
    end
    %% Reset axes
    
    axes(handles.axes_pts);
    %hold off;
    cla;
    legend off;
    title('');
    hold on;
    set(gca, 'XTickMode','auto', 'YTickMode', 'auto', 'ZTickMode', 'auto',...
        'XLimMode', 'auto', 'YLimMode', 'auto', 'ZLimMode', 'auto');
    xlabel('');   ylabel('');   zlabel(''); 
    if numel(params_to_plot)>=3
        view(3);
    else
        view(2);
    end
    
    %% Plots domain - only if all parameters are variables
    var = Br.GetVariables();
    var_dom_to_plot = intersect(params_to_plot,var,'stable');
    if isequal(var_dom_to_plot, params_to_plot)
        hold on;
         Br.PlotDomain(var_dom_to_plot);
    end

    %end
    %% Determines property
    spec_id = handles.current_prop;
    spec = [];
    if ~isempty(spec_id)
        if isfield(Br.P, 'props_names')&&any(strcmp(Br.P.props_names, handles.current_prop))
            spec = handles.properties.(spec_id);
        end
    end
    
    if isempty(spec)
        Br.PlotParams(params_to_plot);
        title('');
    else
        Br.PlotSatParams(spec, params_to_plot);
    end
    
    set_brush_from_selected(handles);
end

function set_brush_from_selected(handles)
Br = get_current_set(handles);
if ~isempty(Br)
    
    if any(Br.P.selected)
        kn = find(Br.P.selected);
        num_axes = min(numel(handles.current_plot_pts), 3);
        pval_select = Br.GetParam(handles.current_plot_pts(1:num_axes), kn);
        h_scat = get_scatter_handle();
        switch num_axes
            case 1
                data = get(h_scat, 'XData');
            case 2
                data = [get(h_scat, 'XData'); get(h_scat, 'YData')];
            case 3
                data = [get(h_scat, 'XData'); get(h_scat, 'YData');get(h_scat, 'ZData')];
        end
        % Find data in select
        if ~isempty(data)
            h_scat.BrushData  = double(ismember(data', pval_select','rows')');
            set_brush_style();
        end
    end
end

function set_selected_from_brush(handles)
Br = get_current_set(handles);

if ~isempty(Br)
    
    num_axes = min(numel(handles.current_plot_pts), 3);
    pval = Br.GetParam(handles.current_plot_pts(1:num_axes));
    
    h_scat = get_scatter_handle();
    switch num_axes
        case 1
            data = get(h_scat, 'XData');
        case 2
            data = [get(h_scat, 'XData'); get(h_scat, 'YData')];
        case 3
            data = [get(h_scat, 'XData'); get(h_scat, 'YData');get(h_scat, 'ZData')];
    end
    
    % Find pval in brush data in select
    brushed_data = data(:, logical(h_scat.BrushData));
    Br.P.selected =  double(ismember(pval', brushed_data','rows')');
end

% --------------------------------------------------------------------
function menu_plot_property_val_Callback(hObject, eventdata, handles)
% hObject    handle to menu_plot_property_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = get_current_set(handles);
if ~isempty(Br)
    try
        close(handles.figp)
    end
    handles.figp = figure;
    hold on;
    num_params = numel(handles.selected_params);
    Br.PlotRobustMap(handles.properties.(handles.current_prop), handles.selected_params(1:min(num_params,2)));
    guidata(hObject,handles);
end

function req = get_current_req_name(handles)
req = handles.current_prop;


% --------------------------------------------------------------------
    function menu_falsify_Callback(hObject, eventdata, handles)
        Br = get_current_set(handles);
        if ~isempty(Br)
            req = get_current_req_name(handles);
            pb = FalsifWizard('ParamSet', handles.current_set,'Requirement', req);
            if ~isempty(pb)
                pb.solve();
                handles = get_param_sets(handles,Br);
                handles = update_working_sets_panel(handles);
                guidata(hObject, handles);
            end
        end

% --------------------------------------------------------------------
function menu_select_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function menu_select_N_worst_Callback(hObject, eventdata, handles)
Br = get_current_set(handles);
if ~isempty(Br)
    valst = inputdlg('Number of samples to select:','',1,{'1'});
    if isempty(valst)
        return;
    else
        Br.SortbyRob();
        Br.SortbySat();
        val = str2num(valst{1});
        val = min(val, Br.GetNbParamVectors());
        Br.P.selected(:,:) = 0;
        Br.P.selected(1:val) = 1;
    end
    handles =  plot_pts(handles);
    guidata(hObject,handles);
end

% --------------------------------------------------------------------
function menu_select_prop_gt_Callback(hObject, eventdata, handles)
Br = get_current_set(handles);
if ~isempty(Br)
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
    handles =  plot_pts(handles);
    guidata(hObject,handles);
end

% --------------------------------------------------------------------
function menu_select_prop_abs_st_Callback(hObject, eventdata, handles)
Br = get_current_set(handles);

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
handles =  plot_pts(handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function menu_inverse_selection_Callback(hObject, eventdata, handles)
Br = get_current_set(handles);
if ~isempty(Br)
    Br.P.selected = ~ Br.P.selected;
    handles =  plot_pts(handles);
    guidata(hObject,handles);
end

% --------------------------------------------------------------------
function menu_unselect_Callback(hObject, eventdata, handles)
Br = get_current_set(handles);
if ~isempty(Br)
    Br.P.selected = 0* Br.P.selected;    
    handles = plot_pts(handles);
    guidata(hObject,handles);
end

% --- Executes on button press in button_new_prop.
function button_new_prop_Callback(hObject, eventdata, handles)
Br = get_current_set(handles);
if ~isempty(Br)
    PHI_= Br.AddSpecGUI();
    
    if isa(PHI_,'STL_Formula')
        eval([get_id(PHI_) '=PHI_']);
        handles = update_properties_panel(handles);
        guidata(hObject, handles);
    else
        info(handles,'No new formula was defined.');
    end
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

fn = fieldnames(handles.properties);

if numel(fn)==0
    return
end

st = fn{get(handles.listbox_prop,'Value')};
val = get(handles.listbox_prop,'Value');

if(val>1)
    val = val-1;
    set(handles.listbox_prop,'Value', val);
    handles.current_prop = fn{val};
elseif numel(fn)==1
    handles.current_prop = '';
    
else
    handles.current_prop = fn{val+1};
end

handles.properties = rmfield(handles.properties,st);

Br = get_current_set(handles);
Br.P = SPurge_props(Br.P);
Br.Specs.remove(st);

handles = update_properties_panel(handles);
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_load_properties_Callback(hObject, eventdata, handles)
% hObject    handle to menu_load_properties (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load_requirement(hObject, handles);

function load_requirement(hObject, handles)
Br = get_current_set(handles);
if ~isempty(Br)
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
end

% --------------------------------------------------------------------
function menu_save_properties_Callback(hObject, eventdata, handles)
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

% --------------------------------------------------------------------
function menu_properties_Callback(hObject, eventdata, handles)

% --- Executes on button press in button_explore_traj.
function button_explore_traj_Callback(hObject, eventdata, handles)
Br = get_current_set(handles);
if ~isempty(Br)
    h = BreachTrajGui(Br, handles);
end

% --- Executes on button press in button_save.
function button_save_Callback(hObject, eventdata, handles)
ws = handles.working_sets; %#ok<NASGU>
save(handles.working_sets_file_name, '-struct', 'ws');

% --------------------------------------------------------------------
function menu_del_select_Callback(hObject, eventdata, handles)
Br = get_current_set(handles);
if ~isempty(Br)
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
end
% --------------------------------------------------------------------
function menu_remove_prop_Callback(hObject, eventdata, handles)
% hObject    handle to menu_remove_prop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = get_current_set(handles);
if ~isempty(Br)
    Br.P = SPurge_props(Br.P);
    handles = update_working_sets_panel(handles);
    handles = update_properties_panel(handles);
    handles = update_modif_panel(handles);
    guidata(hObject,handles);
end

% --- Executes on button press in button_break_prop.
function button_break_prop_Callback(hObject, eventdata, handles)
% hObject    handle to button_break_prop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Br = get_current_set(handles);
%JOHAN UGLY FIX
prop = handles.properties.(handles.current_prop);
%prop = handles.properties.phi;
% END JOHAN UGLY FIX
props = STL_Break(prop);

for i = 1:numel(props)
    PHI_ = props(i);
    this_id = get_id(PHI_);
    this_id = strrep(this_id, '__', '_');
    PHI_ = set_id(PHI_, this_id);
    eval([this_id '=PHI_;']);
    handles.properties.(this_id) = PHI_;
    Br.AddSpec(PHI_);
end

handles = update_properties_panel(handles);
guidata(hObject, handles);

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
Br = get_current_set(handles);
if ~isempty(Br)
    iprop = find_prop(handles.current_prop,Br.P);
    val_threshold = 0;
    if iprop
        val = cat(1,Br.P.props_values(iprop,:).val);
        val = val(:,1);
        Br.P.selected = (val<val_threshold)';
    end   
    handles = plot_pts(handles);
    guidata(hObject,handles);
end

function h = info(h,msg)
% INFO write the message into the information panel.
set(h.text_info, 'String', msg);
drawnow();

% --------------------------------------------------------------------
function menu_param_synth_Callback(hObject, eventdata, handles)
Br = get_current_set(handles);
req = get_current_req_name(handles);
pb = ParamSynthWizard('ParamSet', handles.current_set,'Requirement', req);
if ~isempty(pb)
    pb.solve();
    handles = get_param_sets(handles,Br);
    handles = update_working_sets_panel(handles);
    guidata(hObject, handles);
end



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
Br = get_current_set(handles);
[params, p0, domains] = read_uitable_params(hObject);

% look for parameter to change
idx_params = find(strcmp(params, handles.selected_params(1)) );

% check whether we changed a value or a domain
col = eventdata.Indices(2);
switch col
    case 2  % change value
        idx = FindParam(Br.P,handles.selected_params(1));
        if isempty(domains(idx_params).domain)  % constant parameter, set everybody
            Br.SetParam(idx,  p0(idx_params), true);
        else % set only current_pts, tricky and slightly inefficient- consider having SetParam handling this in the future
            all_values = Br.GetParam(idx);
            all_values( handles.current_pts) = p0(idx_params);
            Br.SetParam(idx, all_values, true);
        end
    case {3,4,5}  % change domain
        Br.SetDomain(handles.selected_params,domains(idx_params));
end
Br.CheckinDomain();
handles = info(handles, 'Parameters changed - rerun simulations and/or check requirements');
handles = update_modif_panel(handles);
guidata(hObject,handles);

% --- Executes on button press in button_all.
function button_all_Callback(hObject, eventdata, handles)
Br = get_current_set(handles);
if ~isempty(Br)
    handles.show_params = Br.P.ParamList;
    handles = fill_uitable(handles);
    guidata(hObject,handles);
end

% --- Executes on button press in button_sys_params.
function button_sys_params_Callback(hObject, eventdata, handles)
Br = get_current_set(handles);
if ~isempty(Br)
    handles.show_params = Br.GetSysParamList();
    handles = fill_uitable(handles);
    guidata(hObject,handles);
end

% --- Executes on button press in button_inputs.
function button_inputs_Callback(hObject, eventdata, handles)
Br = get_current_set(handles);
if ~isempty(Br)
    idx = Br.GetParamsInputIdx();
    handles.show_params = Br.P.ParamList(idx);
    handles = fill_uitable(handles);
    guidata(hObject,handles);
end

% --- Executes on button press in button_prop_param.
function button_prop_param_Callback(hObject, eventdata, handles)
Br = get_current_set(handles);
if ~isempty(Br)
    handles.show_params = Br.GetPropParamList();
    handles = fill_uitable(handles);
    guidata(hObject,handles);
end

function edit_filter_Callback(hObject, eventdata, handles)
Br = get_current_set(handles);
if ~isempty(Br)
    st = get(hObject,'String');
    if ~isempty(st)
        handles.show_params = Br.P.ParamList(cellfun(@(c)(~isempty(c)),  regexp( Br.P.ParamList,st)));
    end
    
    handles = fill_uitable(handles);
    guidata(hObject,handles);
end

% --- Executes during object creation, after setting all properties.
function edit_filter_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_domain.
function button_domain_Callback(hObject, eventdata, handles)
Br = get_current_set(handles);
if ~isempty(Br)
    handles.show_params = Br.GetBoundedDomains();
    
    handles = fill_uitable(handles);
    guidata(hObject,handles);
end

% --- Executes on selection change in popup_sample_option.
function popup_sample_option_Callback(hObject, eventdata, handles)

val = get(hObject,'Value');
switch(val)
    case 1
        handles.sample_arg_multi = 'replace';
    case 2
        handles.sample_arg_multi = 'append';
    case 3
        handles.sample_arg_multi = 'combine';
end

handles = update_sample_args(handles);
guidata(hObject,handles);


% --- Executes on selection change in popup_method.
function popup_method_Callback(hObject, eventdata, handles)

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
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected cell(s) is changed in uitable_params.
function uitable_params_CellSelectionCallback(hObject, eventdata, handles)
idx = eventdata.Indices;
handles.select_cells = idx;
if ~isempty(idx)
    handles.selected_params = handles.show_params(idx(:,1));
    st_sample= get_sample_string(handles);
    st = [st_sample ' - press ''ctrl-p'' to change plot axis.'];
    handles = info(handles, st);
    guidata(hObject,handles);
end

function st_sample = get_sample_string(handles)
Br = get_current_set(handles);
if ~isempty(Br)
    params = handles.selected_params;
    isSig = Br.isSignal(params);
    signals = params(isSig);
    params = params(~isSig);
    st_sample = ['Selected: ' ];
    
    if ~isempty(signals)
        st_signals  = cell2mat(cellfun(@(c) ( [c ' ' ] ) , signals, 'UniformOutput', false ));
        st_sample = [st_sample 'Signals: { ' st_signals '} '];
        set(handles.button_sample, 'Enable', 'off');
    end
    
    if ~isempty(params)
        
        domains = Br.GetDomain(params);
        variables = {};
        for id = 1:numel(domains)
            if  ~isempty(domains(id).domain)
                variables = [variables params(id)];
            end
        end
        if ~isempty(variables)
            st_variables = cell2mat(cellfun(@(c) ( [c ' ' ] ) ,variables, 'UniformOutput', false ));
            st_sample  = [st_sample ' Variables { ' st_variables '} '];
        end
        
        params = setdiff(params,variables);
        if ~isempty(params)
            st_params  = cell2mat(cellfun(@(c) ( [c ' ' ] ) , params, 'UniformOutput', false ));
            st_sample  = [st_sample 'Parameters { ' st_params '}'];
            set(handles.button_sample, 'Enable', 'off');
        end
    end
end


function st_dom = get_domain_string(handles)
st_dom =  cell2mat(cellfun(@(c) ( [c ' ' ] ) , handles.selected_params, 'UniformOutput', false ));

% --- Executes on button press in button_reset.
function button_reset_Callback(hObject, eventdata, handles)
Br = get_current_set(handles);
if ~isempty(Br)
    Br.ResetParamSet();
    handles = update_working_sets_panel(handles);
    handles = update_properties_panel(handles);
    handles = update_modif_panel(handles);
    handles = set_default_plot(handles);
    handles = plot_pts(handles);
    guidata(hObject,handles);
end

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
                    
                    Br = get_current_set(handles);
                    [params, p0, domains] = read_uitable_params(hObject);
                    
                    for ip = 1:numel(handles.selected_params)
                        idx_param = find(strcmp(params, handles.selected_params(ip)) );
                        Br.SetDomain(handles.selected_params,domains(idx_param));
                        idx_P = FindParam(Br.P, handles.selected_params);
                        if isempty(domains(idx_param).domain)
                            Br.P.pts(idx_P, : ) = p0(ip);
                        else
                            Br.P.pts(idx_P, handles.current_pts ) = p0(ip);
                        end
                    end
                    Br.CheckinDomain();
                    handles = info(handles, 'Parameters changed - rerun simulations and/or check requirements');
                    handles = plot_pts(handles);
                    guidata(hObject,handles);
                    
                end
            end
        case 'p'
            if strcmp(eventdata.Modifier, 'control')
                handles.current_plot_pts = handles.selected_params;
                handles = plot_pts(handles);
                guidata(hObject,handles);
            end
    end
end

function breach_WindowButtonDownFcn(hObject, eventdata, handles)

% --- Executes on mouse motion over figure - except title and menu.
function breach_WindowButtonMotionFcn(hObject, eventdata, handles)

% --------------------------------------------------------------------
function brushtoggle_OnCallback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function brushtoggle_ClickedCallback(hObject, eventdata, handles)
delete(findall(gca,'Type','hggroup','HandleVisibility','off'));

if isequal(hObject.State, 'on')
    br = brush( handles.breach);
    set(br,'Enable','on','Color', [0.8 0.8 0.8],  'ActionPostCallback', @(e,d) brush_callback(handles,e,d));
    %  set(br,'Enable','on','Color', [0.8 0.8 0.8]);
else
    brush( handles.breach, 'off' );
end

function handles = disable_brush(handles)
handles.brushtoggle.State = 'Off';
brush( handles.breach, 'off' );

% --------------------------------------------------------------------
function uitoggletool5_OnCallback(hObject, eventdata, handles)
handles = disable_brush(handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function uitoggletool1_OnCallback(hObject, eventdata, handles)
handles = disable_brush(handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function uitoggletool2_OnCallback(hObject, eventdata, handles)
handles = disable_brush(handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function uitoggletool3_OnCallback(hObject, eventdata, handles)
handles = disable_brush(handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function uitoggletool4_OnCallback(hObject, eventdata, handles)
handles = disable_brush(handles);
guidata(hObject,handles);


% --- Executes on button press in button_load_req.
function button_load_req_Callback(hObject, eventdata, handles)
load_requirement(hObject,handles);

% --------------------------------------------------------------------
function uipushtool_load_ClickedCallback(hObject, eventdata, handles)
load_breachset(hObject,handles);

% --------------------------------------------------------------------
function uipushtool_export_ClickedCallback(hObject, eventdata, handles)
ws = handles.working_sets;
try
    handles = info(handles, ['Saving parameter set to ' handles.working_sets_file_name '...']);
    %   save(handles.working_sets_file_name, '-struct', 'ws');
    BreachSave(handles.working_sets_file_name);
catch
    [FileName,PathName] = uiputfile('*.mat','Save Parameter Set As...');
    if(FileName==0)
        return;
    end
    handles.working_sets_file_name = [PathName  FileName];
    handles = info(handles, ['Saving parameter set to ' handles.working_sets_file_name '...']);
    BreachSave(handles.working_sets_file_name);
    handles = update_working_sets_panel(handles);
end

handles = info(handles, ['Saving parameter sets to ' handles.working_sets_file_name '... Done']);
guidata(hObject, handles);


% --- Executes on button press in button_sample.
function button_sample_Callback(hObject, eventdata, handles)
try
    handles = run_sample_domain(handles);
    Br = get_current_set(handles);
    if ~isempty(Br)
        modif_panel_title = Br.disp();
        set(handles.modif_param_panel,'Title', modif_panel_title);
        plot_pts(handles);
        guidata(hObject,handles);
    end
catch
    [msg, msgid] = lasterr;
    handles = info(handles,['Error: ' msgid '--' msg]);
    guidata(hObject,handles);
end

% --------------------------------------------------------------------
function menu_inputs_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function menu_import_inputs_Callback(hObject, eventdata, handles)

Br = get_current_set(handles);
if ~isempty(Br)
    [filenames, paths] = uigetfile( ...
        {  '*.mat','MAT-files (*.mat)'}, ...
        'Pick one or more files', ...
        'MultiSelect', 'on');
    
    if isequal(filenames,0) % cancel
        return
    end
    if ~iscell(filenames)
        filenames= {filenames};
    end
    
    files = cellfun( @(c)([ paths c  ] ), filenames,'UniformOutput',false);
    signals = Br.GetSignalList();
    input_signals = signals(Br.GetInputSignalsIdx);
    BInputData = BreachImportData(files, input_signals);
    
    if numel(BInputData.GetParamList)>1
        params_all = BInputData.GetParamList;
        params = select_cell_gui(params_all(2:end), params_all(2:end), 'Select system parameters to import from files');    
        if isequal(params,0)
            return
        end
        BInputData = BreachImportData(files, input_signals, params);
    end
    
    Br.SetInputGen(BInputData);
    Br.use_precomputed_inputs = true;
    
    idx = Br.GetParamsInputIdx();
    handles.show_params = Br.P.ParamList(idx);
    handles = update_modif_panel(handles);
    guidata(hObject, handles);
end

% --------------------------------------------------------------------
function menu_export_signals_Callback(hObject, eventdata, handles)
% hObject    handle to menu_export_signals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = get_current_set(handles);
if ~isempty(Br)
    pat= uigetdir;
    if isequal(pat,0)
        return;
    end
    Br.SaveSignals([], pat );
end

% --------------------------------------------------------------------
function open_mdl_Callback(hObject, eventdata, handles)
% hObject    handle to open_mdl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Br = get_current_set(handles);
if isa(Br,'BreachSimulinkSystem')
    Br.OpenMdl();
end


% --------------------------------------------------------------------
function open_breach_mdl_Callback(hObject, eventdata, handles)
Br = get_current_set(handles);
if isa(Br,'BreachSimulinkSystem')
    Br.OpenBreachMdl();
end


% --------------------------------------------------------------------
function reset_mdl_Callback(hObject, eventdata, handles)
Br = get_current_set(handles);
Bnew = BreachSimulinkWizard(Br.mdl.name, Br, 'VariableName',handles.current_set);
if ~isempty(Bnew)
    handles.working_sets.(handles.current_set) = Bnew;
    handles = update_working_sets_panel(handles);
    handles = update_modif_panel(handles);
    guidata(hObject, handles);
end

% --- Executes on button press in button_parallel.
function button_parallel_Callback(hObject, eventdata, handles)
Br = get_current_set(handles);
if ~isempty(Br)
    if Br.use_parallel
        Br.StopParallel();
    else
        Br.SetupParallel();
    end
end

function Br = get_current_set(handles)
Br = [];
try
    Br = handles.working_sets.(handles.current_set);
end

% --------------------------------------------------------------------
function menu_set_input_gen_Callback(hObject, eventdata, handles)

Br = get_current_set(handles);
if ~isempty(Br)
    hsi = Br.SetInputGenGUI;
    waitfor(hsi);
    idx = Br.GetParamsInputIdx();
    handles.show_params = Br.P.ParamList(idx);
    handles = update_modif_panel(handles);
    handles = set_default_plot(handles);
    handles = plot_pts(handles);
    guidata(hObject,handles);
end

function edit_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    Br = get_current_set(handles);
    if ~isempty(Br)
        old_time_string = get_time_string(Br.GetTime());
        st = get(hObject,'String'); % returns contents of edit_time as text
        if ~isequal(old_time_string, st)
            Br.ResetSimulations();
        end
        time = eval(get(hObject,'String'));
        Br.SetTime(time);
        handles = update_working_sets_panel(handles);
        set(handles.button_compute_traj, 'Enable','on');
        guidata(hObject,handles);
    end
catch
    handles = info(handles, 'Invalid simulation time.' );
    set(handles.button_compute_traj, 'Enable','off');
    guidata(hObject,handles);
end

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


% --------------------------------------------------------------------
function menu_clear_traces_Callback(hObject, eventdata, handles)
% hObject    handle to menu_clear_traces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Br = get_current_set(handles);
if ~isempty(Br)
    Br.ResetSimulations();
    handles = update_working_sets_panel(handles);
    guidata(hObject, handles);
end

% --------------------------------------------------------------------
function menu_import_traces_Callback(hObject, eventdata, handles)

new_name = evalin('base','matlab.lang.makeUniqueStrings(''Bimport'', who)');
B = BreachImportData();
if isa(B.signalGenerators{1}, 'constant_signal_gen') % canceled
    return;
end

all_signals = B.GetSignalList();
signals = select_cell_gui(all_signals ,all_signals, 'Select signals from list');
if ~isempty(signals)&&~isequal(signals,0)
    B = BreachImportData(B.fname, signals);
else
    return;
end
assignin('base', new_name,B);
handles = get_param_sets(handles,B);
handles.show_params = B.P.ParamList;
handles = update_working_sets_panel(handles);
handles = update_modif_panel(handles);
handles =set_default_plot(handles);
handles = plot_pts(handles);

guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_reqmining_Callback(hObject, eventdata, handles)
Br = get_current_set(handles);
if ~isempty(Br)
    req = get_current_req_name(handles);
    if isempty(req)
    pb = ReqMiningWizard('ParamSet', handles.current_set);
    else    
    pb = ReqMiningWizard('ParamSet', handles.current_set,'Requirement', req);
    end
    if ~isempty(pb)
        pb.solve();
        % Extract
        pb_name = pb.whoamI;
        sol_name = [pb_name '_result'];
        log_name = [pb_name '_falsif_log'];
        Bres = pb.synth_pb.GetBrSet_Best();
        Blog = pb.falsif_pb.GetBrSet_Logged();
        assignin('base', sol_name, Bres);
        assignin('base', log_name, Blog);
        handles = get_param_sets(handles,Br);
        handles.current_set = sol_name;
        handles.show_params = Bres.GetVariables();
        handles.current_plot_pts = handles.show_params;
        handles = update_working_sets_panel(handles);
        handles = update_modif_panel(handles);
        handles = plot_pts(handles);
        guidata(hObject, handles);
    end
end


% --------------------------------------------------------------------
function menu_export_to_excel_Callback(hObject, eventdata, handles)

    Br = get_current_set(handles);

    if ~isempty(Br)
        opt.FileName = 'Results.xlsx';
        choices.FileName = 'string';
        tips.FileName = 'Choose a name for Excel file.';
        
        gu = BreachOptionGui('Export to Excel sheet', opt, choices, tips);
        uiwait(gu.dlg);
        if ~isempty(gu.output)
            if isa(Br, 'BreachImportData')
               Br.ExportToExcel('FileName', gu.output.FileName);
            elseif isa(Br.InputGenerator, 'BreachImportData')
                % Let's cheat
                Bi = Br.InputGenerator.copy();
                Bi.P = Br.P;  % arrrg.
                Bi.ExportToExcel('FileName', gu.output.FileName);
            elseif isa(Br, 'BreachSimulinkSystem')
               Br.ExportToExcel(gu.output.FileName);
            else
                handles = info(handles, 'Export to Excel not supported for this type of system.');
                guidata(hObject, handles);
            end
        end
    end
    




% --- Executes during object creation, after setting all properties.
function breach_CreateFcn(hObject, eventdata, handles)
% hObject    handle to breach (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
