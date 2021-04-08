function varargout = BreachTrajGui(varargin)
% BREACHTRAJGUI M-file for BreachTrajGui.fig
%      BREACHTRAJGUI, by itself, creates a new BREACHTRAJGUI or raises the existing
%      singleton*.
%
%      H = BREACHTRAJGUI returns the handle to a new BREACHTRAJGUI or the handle to
%      the existing singleton*.
%
%      BREACHTRAJGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BREACHTRAJGUI.M with the given input arguments.
%
%      BREACHTRAJGUI('Property','Value',...) creates a new BREACHTRAJGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BreachTrajGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BreachTrajGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BreachTrajGui

% Last Modified by GUIDE v2.5 18-Apr-2017 15:20:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',     mfilename,  ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @BreachTrajGui_OpeningFcn, ...
    'gui_OutputFcn',  @BreachTrajGui_OutputFcn, ...
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


% --- Executes just before BreachTrajGui is made visible.
function BreachTrajGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BreachTrajGui (see VARARGIN)

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
%set(handles.main, 'Position',POS);


%% load parameter set
handles.BrSys = varargin{1};
handles_main = varargin{2};
handles.protected_names = fieldnames(handles_main.working_sets);
handles.properties= handles_main.properties;

if isfield(handles_main, 'TrajSet')
    handles.TrajSet = handles_main.TrajSet;
else
    handles.TrajSet= handles_main.working_sets.(handles_main.current_set);
end

if (isfield(handles.BrSys.P, 'traj'))
    if ~isfield(handles.BrSys.P, 'traj_ref')
        handles.BrSys.P.traj_ref = 1:numel(handles.BrSys.P.traj);
    end
    handles.traj_ref = handles.BrSys.P.traj_ref;
end

%%  init things
if (isfield(handles.BrSys.Sys,'time_mult'))
    time_mult = handles.BrSys.Sys.time_mult;
else
    time_mult=1;
end
handles.BrSys.P.time_mult = time_mult;

handles.has_dynamics=1;
if(isfield(handles.BrSys.Sys, 'type'))
    if (strcmp(handles.BrSys.Sys.type,'traces'))
        handles.has_dynamics = 0;
    end
end

handles.auto_recompute = 1;
set(handles.recompute_auto,'Value',1);
handles.traj_opt = {'b'};
plist = handles.BrSys.P.ParamList(1:handles.BrSys.P.DimX);

pnames = fieldnames(handles.properties);
plist = {plist{:} pnames{:}};

while (numel(plist)<3)
    plist = {plist{:} plist{end}};
end

handles.var_list = plist;
handles.change_tspan_for_all = 0;

set(handles.edit_min_slider_param, 'String', num2str(0));
set(handles.edit_max_slider_param, 'String', num2str(0));

handles.exported_axes=[];
handles.selected_prop =0;

% menu for axes 1
set(handles.param1,'String',plist);
handles.current_var{1,1} = plist(1);
handles.current_var{1,2} = '';
handles.current_var{1,3} = '';

set(handles.param1,'Value', 1);

plist2 = {'' plist{:}};
handles.param_list=handles.BrSys.P.ParamList(:);

set(handles.param12,'String',plist2);
handles.current_var{2,1} = plist2(1);
set(handles.param12,'Value', 1);

set(handles.param13,'String',plist2);
handles.current_var{3,1} = plist2(1);
set(handles.param13,'Value', 1);

handles.current_plot{1} =[];
handles.current_plot{2} =[];
handles.current_plot{3} =[];
handles.plot_tout = 0;

% menu for axes 2

if numel(plist)>=2
    set(handles.param2,'String',plist);
    handles.current_var{2,1} = plist(2);
    set(handles.param2,'Value', 2);
    
    set(handles.param22,'String',plist2);
    handles.current_var{2,2} = plist2(1);
    set(handles.param22,'Value', 1);
    
    set(handles.param23,'String',plist2);
    handles.current_var{2,3} = plist2(1);
    set(handles.param23,'Value', 1);
    
else
    set(handles.param2,'String',{''})
    handles.current_var{2,1} = '';
    set(handles.param2,'Value',1);
    set(handles.param22,'String',{''})
    set(handles.param22,'Value',1);
    set(handles.param23,'String',{''})
    set(handles.param23,'Value',1);
end

% menu for axes 3
if numel(plist)>=3
    
    set(handles.param3,'String',plist);
    handles.current_var{3,1} = plist(3);
    set(handles.param3,'Value', 3);
    
    set(handles.param32,'String',plist2);
    handles.current_var{3,2} = plist2(1);
    set(handles.param32,'Value', 1);
    
    set(handles.param33,'String',plist2);
    handles.current_var{3,3} = plist2(1);
    set(handles.param33,'Value', 1);
    
else
    
    handles.current_var{3,1} = '';
    set(handles.param3,'String',{''})
    set(handles.param3,'Value',1);
    set(handles.param32,'String',{''})
    set(handles.param32,'Value',1);
    set(handles.param33,'String',{''})
    set(handles.param33,'Value',1);
    
end

% menu for param pts plot

set(handles.popup_pts1,'String',handles.param_list);
handles.current_plot_pts{1} = handles.param_list(handles.BrSys.P.dim(1));
set(handles.popup_pts1,'Value', handles.BrSys.P.dim(1));

set(handles.popup_pts2,'String',{'', handles.param_list{:}});
if numel(handles.BrSys.P.dim)>=2
    handles.current_plot_pts{2} = handles.param_list(handles.BrSys.P.dim(2));
    set(handles.popup_pts2,'Value', 3);
else
    handles.current_plot_pts{2} = '';
    set(handles.popup_pts2,'Value', 1);
end

set(handles.popup_pts3,'String',{'', handles.param_list{:}});
if numel(handles.BrSys.P.dim)>=3
    handles.current_plot_pts{3} = handles.param_list(handles.BrSys.P.dim(3));
    set(handles.popup_pts3,'Value', 4);
else
    handles.current_plot_pts{3} = '';
    set(handles.popup_pts3,'Value', 1);
end

% init pts_axes

h = figure;
axes;
handles.pts_axes = gca;
handles.pts_figure = h;
set(handles.pts_figure,'Visible','off');

% slider

handles.current_pts = 1;
handles.nb_pts = size(handles.BrSys.P.pts,2);
tt=strcat('pts ',num2str(handles.current_pts),'/',num2str(handles.nb_pts));
handles.slider_param_coeff =1;
set(handles.text_pts, 'String', tt );
set(handles.slider2, 'Value',1,'Min', .9999, 'Max',handles.nb_pts, 'SliderStep', [1/handles.nb_pts 10/handles.nb_pts]);

% param list

handles.current_change_param = 0;
handles.tspan= [];

handles = update_listbox_param(handles);
handles = UpdatePlots(handles);

% Choose default command line output for BreachTrajGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BreachTrajGui wait for user response (see UIRESUME)
% uiwait(handles.main);

% --- Outputs from this function are returned to the command line.
function varargout = BreachTrajGui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in param1.
function param1_Callback(hObject, eventdata, handles)
% hObject    handle to param1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.plot_tout)
    handles.plot_tout =0;
    handles.traj_opt = {'b'};
    set(handles.plot_tout_button, 'Value', 0);
end

contents = get(hObject,'String');
strpop = contents{get(hObject,'Value')};
handles.current_var{1,1} = {strpop};
new_plot = plot_param(handles,1);
handles.current_plot{1} = new_plot;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function param1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to param1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    val = get(hObject,'Value');
    handles.current_pts = round(val);
    tt=strcat('pts ',num2str(handles.current_pts),'/',num2str(handles.nb_pts));
    set(handles.text_pts, 'String', tt );
    
    %content = get(handles.listbox, 'String');
    %Pf = handles.BrSys.P;
    
    %for i=1:numel(Pf.dim)
    %  st = Pf.ParamList{Pf.dim(i)};
    %  st = strcat(st, ':',' ',dbl2str(Pf.pts(Pf.dim(i), handles.current_pts)));
    %  content{i+2} = st;
    %  content{Pf.dim(i)+6} = st;
    %  handles.current_varying_param{i} = st;
    %end
    
    %set(handles.listbox, 'String', content);
    
    handles = update_listbox_param(handles,0);
    handles = UpdatePlots(handles);
    guidata(hObject,handles);
end


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in param2.
function param2_Callback(hObject, eventdata, handles)
% hObject    handle to param2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.plot_tout)
    handles.plot_tout =0;
    handles.traj_opt = {'b'};
    set(handles.plot_tout_button, 'Value', 0);
end
contents = get(hObject,'String');
strpop = contents{get(hObject,'Value')};
handles.current_var{2,1} = {strpop};
new_plot = plot_param(handles,2);
handles.current_plot{2} = new_plot;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function param2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to param2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in param3.
function param3_Callback(hObject, eventdata, handles)
% hObject    handle to param3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.plot_tout)
    handles.plot_tout =0;
    handles.traj_opt = {'b'};
    set(handles.plot_tout_button, 'Value', 0);
end

contents = get(hObject,'String');
strpop = contents{get(hObject,'Value')};
handles.current_var{3,1} = {strpop};
new_plot = plot_param(handles,3);
handles.current_plot{3} = new_plot;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function param3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to param3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on button press in plot_tout_button.
function plot_tout_button_Callback(hObject, eventdata, handles)
% hObject    handle to plot_tout_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.plot_tout = get(hObject,'Value');
if (handles.plot_tout == 0)
    handles.traj_opt = {'b'};
else
    handles.traj_opt = {'r', 'LineWidth',4};
end
handles = plot_tout(handles);
handles = UpdatePlots(handles);
guidata(hObject,handles);


% --- Executes on selection change in listbox.
function listbox_Callback(hObject, eventdata, handles)
% hObject    handle to listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    val = get(hObject,'Value');
    nbd = numel(handles.BrSys.P.dim);
    nbp = numel(handles.BrSys.P.ParamList);
    if (val>2&&val<=nbd+2)
        
        val = handles.BrSys.P.dim(val-2);
        valpts = handles.BrSys.P.pts(val,handles.current_pts);
        
    elseif (val>nbd+5)&&(val<=nbd+5+nbp)
        
        val = val - (5+nbd);
        valpts = handles.BrSys.P.pts(val,handles.current_pts);
        
    elseif (handles.indices_selected_prop(val))
        
        pnames = fieldnames(handles.properties);
        prop = handles.properties.(pnames{handles.indices_selected_prop(val)});
        
        if (handles.indices_selected_prop(val-1)~=handles.indices_selected_prop(val))
            handles.selected_prop = handles.indices_selected_prop(val);
            handles.current_change_param = 0;
            handles = update_listbox_param(handles);
            guidata(hObject,handles);
            return
        else
            handles.selected_prop = 0;
        end
        
        par = get_params(prop);
        try
            valpts = min(par.(handles.names_selected_param{val}),1e99);
        catch
            valpts = 0;
        end
        
    else
        
        handles.current_change_param = 0;
        handles = update_listbox_param(handles);
        guidata(hObject,handles);
        
        return
    end
    
    if (valpts>0)
        minv = 0;
        maxv = 2*valpts;
    elseif (valpts<0)
        minv = 2*valpts;
        maxv = 0;
    else
        minv=-1;
        maxv=1;
    end
    
    lbda = maxv-minv;
    handles.slider_param_coeff = lbda;
    set(handles.slider_param, 'Value', valpts/lbda,'Min', minv/lbda, 'Max',maxv/lbda);
    handles.current_change_param = val;
    handles = update_listbox_param(handles);
    
    guidata(hObject,handles);
    
catch
    s = lasterror;
    BreachReport(s);
    return
end



% --- Executes during object creation, after setting all properties.
function listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in param12.
function param12_Callback(hObject, eventdata, handles)
% hObject    handle to param12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = get(hObject,'String');
strpop = contents{get(hObject,'Value')};

if (handles.plot_tout)
    handles.plot_tout =0;
    handles.traj_opt = {'b'};
    set(handles.plot_tout_button, 'Value', 0);
end

handles.current_var{1,2} = {strpop};

new_plot = plot_param(handles,1);
handles.current_plot{1} = new_plot;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function param12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to param12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in param22.
function param22_Callback(hObject, eventdata, handles)
% hObject    handle to param22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.plot_tout)
    handles.plot_tout =0;
    handles.traj_opt = {'b'};
    set(handles.plot_tout_button, 'Value', 0);
end

contents = get(hObject,'String');
strpop = contents{get(hObject,'Value')};
new_plot = plot_param(handles,2);
handles.current_plot{2} = new_plot;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function param22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to param22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in param32.
function param32_Callback(hObject, eventdata, handles)
% hObject    handle to param32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.plot_tout)
    handles.plot_tout =0;
    handles.traj_opt = {'b'};
    set(handles.plot_tout_button, 'Value', 0);
end

contents = get(hObject,'String');
strpop = contents{get(hObject,'Value')};
new_plot = plot_param(handles,3);
handles.current_plot{3} = new_plot;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function param32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to param32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in param13.
function param13_Callback(hObject, eventdata, handles)
% hObject    handle to param13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = get(hObject,'String');
strpop = contents{get(hObject,'Value')};

if (handles.plot_tout)
    handles.plot_tout =0;
    handles.traj_opt = {'b'};
    set(handles.plot_tout_button, 'Value', 0);
end

new_plot = plot_param(handles,1);
handles.current_plot{1} = new_plot;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function param13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to param13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in param23.
function param23_Callback(hObject, eventdata, handles)
% hObject    handle to param23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.plot_tout)
    handles.plot_tout =0;
    handles.traj_opt = {'b'};
    set(handles.plot_tout_button, 'Value', 0);
end

contents = get(hObject,'String');
strpop = contents{get(hObject,'Value')};

new_plot = plot_param(handles,2);
handles.current_plot{2} = new_plot;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function param23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to param23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in param33.
function param33_Callback(hObject, eventdata, handles)
% hObject    handle to param33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.plot_tout)
    handles.plot_tout =0;
    handles.traj_opt = {'b'};
    set(handles.plot_tout_button, 'Value', 0);
end

contents = get(hObject,'String');
strpop = contents{get(hObject,'Value')};

new_plot = plot_param(handles,3);
handles.current_plot{3} = new_plot;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function param33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to param33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in export_button1.
function export_button1_Callback(hObject, eventdata, handles)
% hObject    handle to export_button1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nb_ax=  numel(handles.exported_axes)+3+1;
handles.current_var{nb_ax,1} = handles.current_var{1,1};
handles.current_var{nb_ax,2} = handles.current_var{1,2};
handles.current_var{nb_ax,3} = handles.current_var{1,3};

handles.current_plot{nb_ax} = [];
figure;
axes;
ax = gca;
handles.exported_axes = [handles.exported_axes ax];
if (handles.plot_tout)
    param_to_plot = handles.current_var{1,1};
    if (~strcmp(param_to_plot,''))
        if (~strcmp(handles.current_var{1,2},''))
            param_to_plot = {param_to_plot{:} handles.current_var{1,2}};
        end
        
        if (~strcmp(handles.current_var{1,3},''))
            param_to_plot = {param_to_plot{:} handles.current_var{1,3}};
        end
    end
    SplotTraj(handles.BrSys.P, param_to_plot);
end

new_plot = plot_param(handles,nb_ax);
handles.current_plot{nb_ax} = new_plot;
grid on;
guidata(hObject,handles);

% --- Executes on button press in export_button2.
function export_button2_Callback(hObject, eventdata, handles)
% hObject    handle to export_button2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nb_ax=  numel(handles.exported_axes)+3+1;
handles.current_var{nb_ax,1} = handles.current_var{2,1};
handles.current_var{nb_ax,2} = handles.current_var{2,2};
handles.current_var{nb_ax,3} = handles.current_var{2,3};

%handles.current_sensi{nb_ax,1} = handles.current_sensi{2,1};
%handles.current_sensi{nb_ax,2} = handles.current_sensi{2,2};
%handles.current_sensi{nb_ax,3} = handles.current_sensi{2,3};
handles.current_plot{nb_ax} = [];

%handles.plot_sensi = [handles.plot_sensi handles.plot_sensi(2)];
figure;
axes;
ax = gca;
handles.exported_axes = [handles.exported_axes ax];
if (handles.plot_tout)
    param_to_plot = handles.current_var{2,1};
    if (~strcmp(param_to_plot,''))
        if (~strcmp(handles.current_var{2,2},''))
            param_to_plot = {param_to_plot{:} handles.current_var{2,2}};
        end
        
        if (~strcmp(handles.current_var{2,3},''))
            param_to_plot = {param_to_plot{:} handles.current_var{2,3}};
        end
    end
    SplotTraj(handles.BrSys.P, param_to_plot);
end

new_plot = plot_param(handles,nb_ax);
handles.current_plot{nb_ax} = new_plot;
grid on;
guidata(hObject,handles);


% --- Executes on button press in export_button3.
function export_button3_Callback(hObject, eventdata, handles)
% hObject    handle to export_button3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nb_ax=  numel(handles.exported_axes)+3+1;
handles.current_var{nb_ax,1} = handles.current_var{3,1};
handles.current_var{nb_ax,2} = handles.current_var{3,2};
handles.current_var{nb_ax,3} = handles.current_var{3,3};

handles.current_plot{nb_ax} = [];

figure;
axes;
ax = gca;
handles.exported_axes = [handles.exported_axes ax];
if (handles.plot_tout)
    param_to_plot = handles.current_var{3,1};
    if (~strcmp(param_to_plot,''))
        if (~strcmp(handles.current_var{3,2},''))
            param_to_plot = {param_to_plot{:} handles.current_var{3,2}};
        end
        
        if (~strcmp(handles.current_var{3,3},''))
            param_to_plot = {param_to_plot{:} handles.current_var{3,3}};
        end
    end
    SplotTraj(handles.BrSys.P, param_to_plot);
end

new_plot = plot_param(handles,nb_ax);
handles.current_plot{nb_ax} = new_plot;
grid on;
guidata(hObject,handles);


function edit_change_param_Callback(hObject, eventdata, handles)
% hObject    handle to edit_change_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% try

ind_p = handles.current_change_param;
if (ind_p)
    new_val = str2double(get(hObject,'String'));
    lbda = handles.slider_param_coeff;
    maxv = get(handles.slider_param,'Max')*lbda;
    minv = get(handles.slider_param,'Min')*lbda;
    
    if (new_val>maxv)
        set(handles.slider_param, 'Value', new_val/lbda,'Max',new_val/lbda);
    elseif (new_val<minv)
        set(handles.slider_param, 'Value', new_val/lbda,'Min',new_val/lbda);
    else
        set(handles.slider_param, 'Value', new_val/lbda);
    end
    
    handles = update_listbox_param(handles,1);
    
    guidata(hObject,handles);
    return;
    
elseif (handles.selected_prop)
    
    new_val = get(hObject,'String');
    pnames = fieldnames(handles.properties);
    prop = handles.properties.(pnames{handles.selected_prop});
    prop_params = get_params(prop);
    if (~strcmp( new_val , disp(prop,0)  ))
        
        handles.properties.(pnames{handles.selected_prop}) = STL_Formula(get_id(prop), new_val);
        
        handles.properties.(pnames{handles.selected_prop}) = set_params(handles.properties.(pnames{handles.selected_prop}), prop_params);
        handles = update_listbox_param(handles,1);
        handles = UpdatePlots(handles);
        guidata(hObject,handles);
        
    end
    
end


% --- Executes during object creation, after setting all properties.
function edit_change_param_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_change_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_change_tspan_Callback(hObject, eventdata, handles)
% hObject    handle to edit_change_tspan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    
    handles.tspan= eval(get(hObject,'String'));
    
    if (handles.auto_recompute)
        handles= update_trajectories(handles);
    end
    
    guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_change_tspan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_change_tspan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_pts_checkbox.
function plot_pts_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to plot_pts_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = get(hObject,'Value');
if (val)
    try
        set(handles.pts_figure,'Visible','on');
    catch
        handles= init_pts_figure(handles);
        set(handles.pts_figure,'Visible','on');
    end
    handles = plot_pts(handles);
else
    try
        set(handles.pts_figure,'Visible','off');
    end
end

guidata(hObject,handles);

% --- Executes on button press in for_all_checkbox.
function for_all_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to for_all_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.change_tspan_for_all  = get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in recompute_auto.
function recompute_auto_Callback(hObject, eventdata, handles)
% hObject    handle to recompute_auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.auto_recompute = get(hObject,'Value');
guidata(hObject,handles);


% --- Executes on selection change in popup_pts1.
function popup_pts1_Callback(hObject, eventdata, handles)
% hObject    handle to popup_pts1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = get(hObject,'String');
strpop = contents{get(hObject,'Value')};
handles.current_plot_pts{1} = {strpop};
handles=  plot_pts(handles);

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


% --- Executes on selection change in popup_pts2.
function popup_pts2_Callback(hObject, eventdata, handles)
% hObject    handle to popup_pts2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = get(hObject,'String');
strpop = contents{get(hObject,'Value')};
handles.current_plot_pts{2} = {strpop};
handles=  plot_pts(handles);

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


% --- Executes on selection change in popup_pts3.
function popup_pts3_Callback(hObject, eventdata, handles)
% hObject    handle to popup_pts3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = get(hObject,'String');
strpop = contents{get(hObject,'Value')};
handles.current_plot_pts{3} = {strpop};
handles=  plot_pts(handles);

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


% --- Executes on button press in button_go.
function button_go_Callback(hObject, eventdata, handles)
% hObject    handle to button_go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = update_trajectories(handles);
guidata(hObject,handles);

function handles = UpdatePlots(handles)

handles = plot_pts(handles);
if (~isfield(handles.BrSys.P,'traj'))
    return
end

new_plot = plot_param(handles,1);
handles.current_plot{1} = new_plot;

new_plot = plot_param(handles,2);
handles.current_plot{2} = new_plot;

new_plot = plot_param(handles,3);
handles.current_plot{3} = new_plot;

for i=1:numel(handles.exported_axes)
    new_plot = plot_param(handles,3+i);
    handles.current_plot{3+i} = new_plot;
end

function handles = init_pts_figure(handles)

h = figure;
axes;
handles.pts_axes = gca;
handles.pts_figure = h;

function handles= plot_pts(handles)

if (get(handles.plot_pts_checkbox,'Value'))
    try
        axes(handles.pts_axes);
    catch
        handles = init_pts_figure(handles);
        axes(handles.pts_axes);
    end
    cla;
    param_to_plot = handles.current_plot_pts{1};
    
    if (~strcmp(param_to_plot,''))
        if (~strcmp(handles.current_plot_pts{2},''))
            param_to_plot = {param_to_plot{:} handles.current_plot_pts{2}};
        end
        
        if (~strcmp(handles.current_plot_pts{3},''))
            param_to_plot = {param_to_plot{:} handles.current_plot_pts{3}};
        end
    end
    
    %    SplotPts(handles.BrSys.P,param_to_plot);
    S = DiscrimPropValues(handles.BrSys.P);
    SplotPts(S,param_to_plot);
    SplotBoxPts(handles.BrSys.P, param_to_plot,handles.current_pts,'+k','r',.1);
    
end

function new_plot= plot_param(handles,ax)

if (isfield(handles.BrSys.Sys,'time_mult'))
    time_mult = handles.BrSys.Sys.time_mult;
else
    time_mult=1;
end

if (~isfield(handles.BrSys.P,'traj'))
    
    new_plot = [];
    return
end

switch (ax)
    case 1
        axes(handles.axes1);
    case 2
        axes(handles.axes2);
    case 3
        axes(handles.axes3);
    otherwise
        try
            axes(handles.exported_axes(ax-3))
        catch
            new_plot= [];
            return;
        end
end

param_to_plot = handles.current_var{ax,1};
pnames = fieldnames(handles.properties);
nprop = find_prop(param_to_plot{1}, pnames);

if (nprop)  % plot values for a property
    prop = handles.properties.(pnames{nprop});
    cla;legend('off');
    hold on;
    grid on;
    st = param_to_plot{1};
    ipts = handles.current_pts;
    %      phi_val = handles.BrSys.P.props_values(nprop,ipts).val;
    %      phi_tspan = handles.BrSys.P.props_values(nprop,ipts).tspan;
    
    colors = {'b','r'};
    
    prop = handles.properties.(pnames{nprop});
    cla;legend('off');
    hold on;
    grid on;
    st = param_to_plot{1};
    ipts = handles.current_pts;
    phi_tspan = handles.BrSys.P.traj{handles.traj_ref(ipts)}.time;
    
    Ptmp = Sselect(handles.BrSys.P, ipts);
    
    phi_val = STL_Eval(handles.BrSys.Sys,prop,Ptmp,handles.BrSys.P.traj{handles.traj_ref(ipts)},phi_tspan);
    ylabel(get_id(prop),'Interpreter','none');
    xlabel(['time']);
    
    %new_plot = plot(phi_tspan*time_mult, phi_val);
    hold on;
    
    % JOHAN CHANGE
    new_plot = stairs(phi_tspan*time_mult, phi_val);
    %plot([phi_tspan(1) phi_tspan(end)]*time_mult, [0 0],'-k');
    stairs(phi_tspan*time_mult, (phi_val>=0)*max(abs(phi_val))/2,'-r');
    % lgh = legend(short_disp(prop,100));
    prop_type = get_type(prop);
    prop_id = get_id(prop);
    children = get_children(prop);
    
    child1_id = get_id(children{1});
    child1_id = strrep(child1_id, '__', '_');
    
    if numel(children) > 1
        child2_id = get_id(children{2});
        child2_id = strrep(child2_id, '__', '_');
    end
    
    if strcmp(prop_type,'always')
        legendText1 = ['always(' child1_id ')'];
        legendText2 = ['min_[t+a, t+b](' child1_id ')'];
    elseif strcmp(prop_type,'or')
        legendText1 = [child1_id ' or ' child2_id];
        legendText2 = ['max(' child1_id ', ' child2_id ')'];
    elseif strcmp(prop_type,'and')
        legendText1 = [child1_id ' and ' child2_id];
        legendText2 = ['min(' child1_id ', ' child2_id ')'];
    elseif strcmp(prop_type,'not')
        legendText1 = ['not(' child1_id ')'];
        legendText2 = ['-' child1_id ''];
    elseif strcmp(prop_type,'predicate')
        legendText1 = get_st(prop);
        legendText2 = get_st(prop);
    elseif strcmp(prop_type,'=>')
        legendText1 = [child1_id ' => ' child2_id];
        legendText2 = ['max( - ' child1_id ', ' child2_id ')'];
    else
        legendText1 = 'Legend text missing!';
        legendText2 = 'Robust text missing!';
    end
    lgh = legend(legendText1,legendText2);
    % END JOHAN CHANGE
    
    set(lgh,'Interpreter','none');
    
    new_plot = [];
    
else  % plot values for a variable
    
    if (~strcmp(handles.current_var{ax,2},''))
        param_to_plot = {param_to_plot{:} handles.current_var{ax,2}};
    end
    
    if (~strcmp(handles.current_var{ax,3},''))
        param_to_plot = {param_to_plot{:} handles.current_var{ax,3}};
    end
    
    if (handles.plot_tout)
        if (~isempty(handles.current_plot{ax}))
            for i = 1:numel(handles.current_plot{ax})
                set(handles.current_plot{ax}(i), 'XData',[], 'YData', [], 'ZData', []);
            end
        end
    else
        cla;legend('off');
    end
    
    SplotTraj(handles.BrSys.P, param_to_plot,handles.current_pts, handles.traj_opt);
    children = get(gca, 'Children');
    if (numel(param_to_plot)>1)
        new_plot=  children(1:2);
    else
        new_plot=  children(1);
    end
end % if (nprop)

function handles= plot_tout(handles)

if (~isfield(handles.BrSys.P,'traj'))
    return
end

axes(handles.axes1); cla;legend('off');
handles.current_plot{1} = [];
axes(handles.axes2); cla;legend('off');
handles.current_plot{2} = [];
axes(handles.axes3); cla;legend('off');
handles.current_plot{3} = [];
for i = 1:numel(handles.exported_axes)-3
    try
        axes(handles.exported_axes(i))
    catch
        continue;
    end
    cla; legend('off');
    handles.current_plot{i+3}=[];
end

pnames = fieldnames(handles.properties);

if (handles.plot_tout == 1)
    
    % Axaes 1
    axes(handles.axes1);
    param_to_plot = handles.current_var{1,1};
    
    if find_prop(param_to_plot{1}, pnames)==0
        
        if (~strcmp(handles.current_var{1,2},''))
            param_to_plot = {param_to_plot{:} handles.current_var{1,2}};
        end
        
        if (~strcmp(handles.current_var{1,3},''))
            param_to_plot = {param_to_plot{:} handles.current_var{1,3}};
        end
        
        SplotTraj(handles.BrSys.P, param_to_plot);
    end
    
    % Axes2
    axes(handles.axes2);
    param_to_plot = handles.current_var{2,1};
    if find_prop(param_to_plot{1}, pnames)==0
        if (~strcmp(param_to_plot,''))
            if (~strcmp(handles.current_var{2,2},''))
                param_to_plot = {param_to_plot{:} handles.current_var{2,2}};
            end
            if (~strcmp(handles.current_var{2,3},''))
                param_to_plot = {param_to_plot{:} handles.current_var{2,3}};
            end
            SplotTraj(handles.BrSys.P, param_to_plot);
        end
    end
    % Axes 3
    axes(handles.axes3);
    param_to_plot = handles.current_var{3,1};
    
    if find_prop(param_to_plot{1}, pnames)==0
        if (~strcmp(param_to_plot,''))
            if (~strcmp(handles.current_var{3,2},''))
                param_to_plot = {param_to_plot{:} handles.current_var{3,2}};
            end
            
            if (~strcmp(handles.current_var{3,3},''))
                param_to_plot = {param_to_plot{:} handles.current_var{3,3}};
            end
        end
        SplotTraj(handles.BrSys.P, param_to_plot);
    end
    % Other Axes
    
    for i = 4:numel(handles.exported_axes)+3
        try
            axes(handles.exported_axes(i))
        catch
            continue;
        end
        param_to_plot = handles.current_var{i,1};
        if find_prop(param_to_plot{1}, fnames)==0
            if (~strcmp(param_to_plot,''))
                if (~strcmp(handles.current_var{i,2},''))
                    param_to_plot = {param_to_plot{:} handles.current_var{i,2}};
                end
                
                if (~strcmp(handles.current_var{i,3},''))
                    param_to_plot = {param_to_plot{:} handles.current_var{i,3}};
                end
            end
            SplotTraj(handles.BrSys.P, param_to_plot);
        end
    end
end

% --- Executes on slider movement.
function slider_param_Callback(hObject, eventdata, handles)
% hObject    handle to slider_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles = update_listbox_param(handles,1);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider_param_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function edit_max_slider_param_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max_slider_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    
    maxv = str2double(get(hObject,'String'));
    lbda = handles.slider_param_coeff;
    
    val  = get(handles.slider_param,'Value')/lbda;
    minv = get(handles.slider_param,'Min')/lbda;
    old_maxv = get(handles.slider_param,'Max')/lbda;
    old_diff = old_maxv-minv;
    changed=0;
    
    if (maxv<minv)
        minv = maxv-old_diff;
    end
    
    if (val>maxv)
        val = maxv;
        changed=1;
    end
    
    lbda = maxv-minv;
    handles.slider_param_coeff = lbda;
    
    set(handles.slider_param, 'Value', val/lbda,'Min', minv/lbda, 'Max',maxv/lbda);
    
    handles = update_listbox_param(handles,changed);
    guidata(hObject, handles);
    
catch
    s = lasterror;
    warndlg(['Problem edit_max_slider: ' s.message] );
    error(s);
    return
end

% --- Executes during object creation, after setting all properties.
function edit_max_slider_param_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max_slider_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_min_slider_param_Callback(hObject, eventdata, handles)
% hObject    handle to edit_min_slider_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    lbda = handles.slider_param_coeff;
    minv = str2double(get(hObject,'String'));
    val  = get(handles.slider_param,'Value')*lbda;
    maxv = get(handles.slider_param,'Max')*lbda;
    old_minv = get(handles.slider_param,'Min')*lbda;
    old_diff = maxv-old_minv;
    changed=0;
    
    if (minv>maxv)
        maxv = minv+old_diff;
    end
    
    if (minv>val)
        val = minv;
        changed=1;
    end
    
    lbda = maxv-minv;
    handles.slider_param_coeff = lbda;
    set(handles.slider_param, 'Value', val/lbda,'Min', minv/lbda, 'Max',maxv/lbda);
    
    handles = update_listbox_param(handles,changed);
    guidata(hObject, handles);
    
catch
    s = lasterror;
    warndlg(['Problem edit_min_slider: ' s.message] );
    error(s);
    return
end


% --- Executes during object creation, after setting all properties.
function edit_min_slider_param_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_min_slider_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles = update_trajectories(handles)

if isempty(handles.tspan)
    try
        tspan = handles.BrSys.P.traj{handles.current_pts}.time;
    catch
        tspan= [0 1];
    end
    
else
    tspan = handles.tspan;
end

% check if we need to recompute everything
if (get(handles.for_all_checkbox,'Value')||~isfield(handles.BrSys.P,'traj'))
    
    handles.BrSys.ResetSimulations();
    handles.BrSys.Sim(tspan);
    traj_ref = handles.BrSys.P.traj_ref;
    
    if isfield(handles.BrSys.P,'props')
        for i=1:numel(handles.BrSys.P.props)
            prop = handles.BrSys.P.props(i);
            for j= 1:size(handles.BrSys.P.pts,2)
                traj = handles.BrSys.P.traj{traj_ref(j)};
                P = Sselect(handles.BrSys.P,j);
                handles.BrSys.P.props_values(i,j).val = ...
                    STL_Eval(handles.BrSys.Sys,props, P, traj,handles.BrSys.P.props_values(i,j).tspan);
            end
        end
    end
    new_plot = plot_param(handles,1);
    handles.current_plot{1} = new_plot;
    
    new_plot = plot_param(handles,2);
    handles.current_plot{2} = new_plot;
    
    new_plot = plot_param(handles,3);
    handles.current_plot{3} = new_plot;
    
else % only one traj needs to be computed
    
    Btmp = handles.BrSys.copy();
    Btmp.P = Sselect(handles.BrSys.P,handles.current_pts);
    Btmp.ResetSimulations();
    Btmp.Sim(tspan);
    
    traj_ref = handles.BrSys.P.traj_ref;
    handles.BrSys.P.traj{traj_ref(handles.current_pts)} = Btmp.P.traj{1};
    handles.BrSys.P.Xf(:,traj_ref(handles.current_pts)) = Btmp.P.traj{1}.X(:,end);
    
    % This is needed if, e.g., ComputeTraj called an init_fun which changed
    % some other values in Ptmp.pts
    handles.BrSys.P.pts(:,handles.current_pts) = Btmp.P.pts;
    
    new_plot = plot_param(handles,1);
    handles.current_plot{1} = new_plot;
    
    new_plot = plot_param(handles,2);
    handles.current_plot{2} = new_plot;
    
    new_plot = plot_param(handles,3);
    handles.current_plot{3} = new_plot;
    
end

function handles =  update_listbox_param(handles, changed)

Br = handles.BrSys;

% update values and plots

if nargin==1
    changed =0;
end

lbda = handles.slider_param_coeff;
ind_p = handles.current_change_param;

if (ind_p)
    if (changed)
        new_val = get(handles.slider_param,'Value')*lbda;
        if (ind_p<=numel(Br.P.ParamList))
            if (get(handles.for_all_checkbox,'Value'))
                Br.P.pts(ind_p, :) = new_val;
                if (ind_p<=Br.P.DimP)
                    for i = 1:numel(Br.P.traj)
                        Br.P.traj{i}.param(ind_p) = new_val;
                    end
                end
            else
                Br.P.pts(ind_p, handles.current_pts) = new_val;
                if (ind_p<=Br.P.DimP)
                    Br.P.traj{handles.current_pts}.param(ind_p) = new_val;
                end
                
            end
        else
            iprop = handles.indices_selected_prop(ind_p);
            par = handles.names_selected_param{ind_p};
            if (new_val>1e98)
                new_val = inf;
            end
            pnames = fieldnames(handles.properties);
            handles.properties.(pnames{iprop}) = set_params(handles.properties.(pnames{iprop}),par,new_val);
        end
        
        if (handles.auto_recompute)
            if (ind_p<=Br.P.DimP)
                handles = update_trajectories(handles);
            end
            handles = UpdatePlots(handles);
        end
    end
    minv = get(handles.slider_param,'Min')*lbda;
    maxv = get(handles.slider_param,'Max')*lbda;
    val = get(handles.slider_param,'Value')*lbda;
    
    set(handles.edit_min_slider_param, 'String', dbl2str(minv));
    set(handles.edit_max_slider_param, 'String', dbl2str(maxv));
    set(handles.edit_change_param, 'String', dbl2str(val));
    
elseif (handles.selected_prop)
    pnames = fieldnames(handles.properties);
    set(handles.edit_change_param, 'String', disp(handles.properties.(pnames{handles.selected_prop}),0));
else
    set(handles.edit_min_slider_param, 'String', '');
    set(handles.edit_max_slider_param, 'String', '');
    set(handles.edit_change_param, 'String', '');
    
    
end

% update listbox display

content = {'Varying parameters' '-------------------'};
 handles.current_varying_param = {};
variables = Br.GetVariables();
for i=1:numel(variables)
    st = variables{i};
    val = Br.GetParam(st);
    st = strcat(st, ':',' ',dbl2str(val(handles.current_pts)));
    handles.current_varying_param{i} = st;
end

content = {content{:} handles.current_varying_param{:} '' 'Systems and props parameters' '-------------------'};

for i=1:numel(Br.P.ParamList)
    st = Br.P.ParamList{i};
    st = strcat(st, ':',' ',dbl2str(Br.P.pts(i, handles.current_pts)));
    content = {content{:} st};
end

content = {content{:} ''};

pnames = fieldnames(handles.properties);
for i= 1:numel(pnames)
    
    prop = handles.properties.(pnames{i});
    st_prop = disp(prop,-1);
    content = {content{:}  [get_id(prop) ': '  st_prop]  '-------------------'};
    
    indx = numel(content);
    handles.indices_selected_prop(indx-1:indx) = i;
    
    phi_params = get_params(prop);
    param_names = fieldnames(phi_params);
    param_values = struct2cell(phi_params);
    
    for j=1:numel(param_names)
        if ~strcmp(param_names{j},'fn')
            try
                st = strcat(param_names{j}, ':',' ',dbl2str(param_values{j}));
            catch
                st = strcat(param_names{j}, ':',' ',class(param_values{j}));
            end
            
            content = {content{:} st};
            indx = numel(content);
            handles.indices_selected_prop(indx) = i;
            handles.names_selected_param{indx} = param_names{j};
        end
    end
end

set(handles.listbox,'String', content);

function i = find_prop(st, props_names)
i=0;
for k = 1:numel(props_names)
    if strcmp(st,props_names{k})
        i = k;
        return;
    end
end

function expr = evalin_caller(n,expr)

for i=0:n
    expr = ['evalin(''caller'',' '''' regexprep(expr,'''','''''') ''');'];
end

function st = dbl2str(x)
st = num2str(x, '%0.5g');
