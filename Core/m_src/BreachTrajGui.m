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

% Last Modified by GUIDE v2.5 08-Aug-2012 16:56:10

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

% load parameter set  

handles_main = varargin{2};    
handles.protected_names = fieldnames(handles_main.working_sets);
handles.working_sets_filename = handles_main.working_sets_file_name;
handles.Sys= handles_main.Sys;
handles.properties= handles_main.properties;
handles.TrajSet= handles_main.working_sets.(handles_main.current_set);

if (isfield(handles.TrajSet, 'traj'))
  if ~isfield(handles.TrajSet, 'traj_ref')
    handles.TrajSet.traj_ref = 1:numel(handles.TrajSet.traj);   
  end
  handles.traj_ref = handles.TrajSet.traj_ref;    
end

%  init things

if (isfield(handles.Sys,'time_mult'))
  time_mult = handles.Sys.time_mult;
else
  time_mult=1;
end
handles.TrajSet.time_mult = time_mult;      

handles.has_dynamics=1;
if(isfield(handles.Sys, 'type'))
  if (strcmp(handles.Sys.type,'traces'))
    handles.has_dynamics = 0;
  end
end

handles.auto_recompute = 1;
set(handles.recompute_auto,'Value',1);
handles.traj_opt = {'b'};
plist = handles.TrajSet.ParamList(1:handles.TrajSet.DimX);

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
handles.param_list=handles.TrajSet.ParamList(:);

set(handles.param12,'String',plist2);
handles.current_var{2,1} = plist2(1);
set(handles.param12,'Value', 1);

set(handles.param13,'String',plist2);
handles.current_var{3,1} = plist2(1);
set(handles.param13,'Value', 1);

handles.current_sensi{1,2} = '';
handles.current_sensi{1,3} = '';
handles.current_sensi{2,2} = '';
handles.current_sensi{2,3} = '';
handles.current_sensi{3,2} = '';
handles.current_sensi{3,3} = '';

handles.plot_sensi(1) = 0;
handles.plot_sensi(2) = 0;
handles.plot_sensi(3) = 0;

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
handles.current_plot_pts{1} = handles.param_list(handles.TrajSet.dim(1));
set(handles.popup_pts1,'Value', handles.TrajSet.dim(1));

set(handles.popup_pts2,'String',{'', handles.param_list{:}});
if numel(handles.TrajSet.dim)>=2
  handles.current_plot_pts{2} = handles.param_list(handles.TrajSet.dim(2));
  set(handles.popup_pts2,'Value', 3);  
else
  handles.current_plot_pts{2} = '';
  set(handles.popup_pts2,'Value', 1);  
end

set(handles.popup_pts3,'String',{'', handles.param_list{:}});
if numel(handles.TrajSet.dim)>=3
  handles.current_plot_pts{3} = handles.param_list(handles.TrajSet.dim(3));
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
handles.nb_pts = size(handles.TrajSet.pts,2);
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
% uiwait(handles.figure1);

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
  
% Hints: contents = get(hObject,'String') returns param1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from param1


% --- Executes during object creation, after setting all properties.
function param1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to param1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
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
    %Pf = handles.TrajSet;

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

  catch 
    s = lasterror;
    warndlg(['Problem slider 2: ' s.message] );
    error(s);    
    return
  
  end


% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


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


% Hints: contents = get(hObject,'String') returns param2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from param2


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

% Hints: contents = get(hObject,'String') returns param3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from param3


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

% Hint: place code in OpeningFcn to populate axes1


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3

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

% Hint: get(hObject,'Value') returns toggle state of plot_tout_button


% --- Executes on selection change in listbox.
function listbox_Callback(hObject, eventdata, handles)
% hObject    handle to listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 try
   val = get(hObject,'Value');
   nbd = numel(handles.TrajSet.dim);
   nbp = numel(handles.TrajSet.ParamList);
   if (val>2&&val<=nbd+2)
     
     val = handles.TrajSet.dim(val-2);
     valpts = handles.TrajSet.pts(val,handles.current_pts);
   
   elseif (val>nbd+5)&&(val<=nbd+5+nbp)
   
     val = val - (5+nbd);
     valpts = handles.TrajSet.pts(val,handles.current_pts);
   
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


% Hints: contents = get(hObject,'String') returns listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox


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

if (handles.plot_sensi(1))
  handles.current_sensi{1,2} = {strpop};
else
  handles.current_var{1,2} = {strpop};
end

new_plot = plot_param(handles,1);
handles.current_plot{1} = new_plot;

guidata(hObject, handles);


% Hints: contents = get(hObject,'String') returns param12 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from param12


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
  if (handles.plot_sensi(2))
    handles.current_sensi{2,2} = {strpop};
  else
    handles.current_var{2,2} = {strpop};
  end
  new_plot = plot_param(handles,2);
  handles.current_plot{2} = new_plot;
  guidata(hObject, handles);
  
% Hints: contents = get(hObject,'String') returns param22 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from param22


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
  if (handles.plot_sensi(3))
    handles.current_sensi{3,2} = {strpop};
  else
    handles.current_var{3,2} = {strpop};
  end
  new_plot = plot_param(handles,3);
  handles.current_plot{3} = new_plot;
  
  guidata(hObject, handles);


% Hints: contents = get(hObject,'String') returns param32 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from param32


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

if (handles.plot_sensi(1))
  handles.current_sensi{1,3} = {strpop};
else
  handles.current_var{1,3} = {strpop};
end

new_plot = plot_param(handles,1);
handles.current_plot{1} = new_plot;

guidata(hObject, handles);

% Hints: contents = get(hObject,'String') returns param13 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from param13


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

  if (handles.plot_sensi(2))
    handles.current_sensi{2,3} = {strpop};
  else
    handles.current_var{2,3} = {strpop};
  end
  
  new_plot = plot_param(handles,2);
  handles.current_plot{2} = new_plot;
  guidata(hObject, handles);

% Hints: contents = get(hObject,'String') returns param23 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from param23


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
  
  if (handles.plot_sensi(3))
    handles.current_sensi{3,3} = {strpop};
  else
    handles.current_var{3,3} = {strpop};
  end
  
  new_plot = plot_param(handles,3);
  handles.current_plot{3} = new_plot;
  guidata(hObject, handles);

% Hints: contents = get(hObject,'String') returns param33 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from param33


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

  handles.current_sensi{nb_ax,1} = handles.current_sensi{1,1};
  handles.current_sensi{nb_ax,2} = handles.current_sensi{1,2};
  handles.current_sensi{nb_ax,3} = handles.current_sensi{1,3};  
  handles.current_plot{nb_ax} = [];

  handles.plot_sensi = [handles.plot_sensi handles.plot_sensi(1)];
  figure;
  axes;
  ax = gca; 
  handles.exported_axes = [handles.exported_axes ax];
  if (handles.plot_tout&&~handles.plot_sensi(1))
    param_to_plot = handles.current_var{1,1};
    if (~strcmp(param_to_plot,''))
      if (~strcmp(handles.current_var{1,2},''))
        param_to_plot = {param_to_plot{:} handles.current_var{1,2}};
      end
      
      if (~strcmp(handles.current_var{1,3},''))
        param_to_plot = {param_to_plot{:} handles.current_var{1,3}};
      end
    end
    SplotTraj(handles.TrajSet, param_to_plot);    
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

  handles.current_sensi{nb_ax,1} = handles.current_sensi{2,1};
  handles.current_sensi{nb_ax,2} = handles.current_sensi{2,2};
  handles.current_sensi{nb_ax,3} = handles.current_sensi{2,3};  
  handles.current_plot{nb_ax} = [];

  handles.plot_sensi = [handles.plot_sensi handles.plot_sensi(2)];
  figure;
  axes;
  ax = gca; 
  handles.exported_axes = [handles.exported_axes ax];
  if (handles.plot_tout&&~handles.plot_sensi(2))
    param_to_plot = handles.current_var{2,1};
    if (~strcmp(param_to_plot,''))
      if (~strcmp(handles.current_var{2,2},''))
        param_to_plot = {param_to_plot{:} handles.current_var{2,2}};
      end
      
      if (~strcmp(handles.current_var{2,3},''))
        param_to_plot = {param_to_plot{:} handles.current_var{2,3}};
      end
    end
    SplotTraj(handles.TrajSet, param_to_plot);    
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
  
  handles.current_sensi{nb_ax,1} = handles.current_sensi{3,1};
  handles.current_sensi{nb_ax,2} = handles.current_sensi{3,2};
  handles.current_sensi{nb_ax,3} = handles.current_sensi{3,3};  
  handles.current_plot{nb_ax} = [];

  handles.plot_sensi = [handles.plot_sensi handles.plot_sensi(2)];
  figure;
  axes;
  ax = gca; 
  handles.exported_axes = [handles.exported_axes ax];
  if (handles.plot_tout&&~handles.plot_sensi(2))
    param_to_plot = handles.current_var{3,1};
    if (~strcmp(param_to_plot,''))
      if (~strcmp(handles.current_var{3,2},''))
        param_to_plot = {param_to_plot{:} handles.current_var{3,2}};
      end
      
      if (~strcmp(handles.current_var{3,3},''))
        param_to_plot = {param_to_plot{:} handles.current_var{3,3}};
      end
    end
    SplotTraj(handles.TrajSet, param_to_plot);    
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
      
      if (~strcmp( new_val , disp(prop,0)  ))

        handles.properties.(pnames{handles.selected_prop}) = QMITL_Formula(get_id(prop), new_val);
        handles = update_listbox_param(handles,1);               
        handles = UpdatePlots(handles);               
        guidata(hObject,handles);
        
      end  
      
    end
 % catch
 %   s = lasterror;
 %   warndlg(['Problem edit_change_param: ' s.message] );
 %   return
 % end

  
% Hints: get(hObject,'String') returns contents of edit_change_param as text
%        str2double(get(hObject,'String')) returns contents of edit_change_param as a double


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


% --- Executes on button press in button_save_pts.
function button_save_pts_Callback(hObject, eventdata, handles)
% hObject    handle to button_save_pts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  
  try
    
    Ptmp = select(handles.TrajSet,handles.current_pts);
    Ptmp = SPurge(Ptmp);
    
    name = handles.traj_set_name;
    new_name = genvarname(name, handles.protected_names);
    eval([new_name '= Ptmp']);   
    save(handles.working_sets_filename,'-append', new_name);
    guidata(hObject,handles);
    
  catch 
    s = lasterror;
    warndlg(['Problem saving: ' s.message] );
    error(s);
    return
  end
  
  
% --- Executes on button press in button_save_all.
function button_save_all_Callback(hObject, eventdata, handles)
% hObject    handle to button_save_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

  try
    Ptmp = SPurge(handles.TrajSet);
    name = handles.traj_set_name;
    new_name = genvarname(name, handles.protected_names);
    eval([new_name '= Ptmp']);   
    save(handles.working_sets_filename,'-append', new_name);
    
    guidata(hObject,handles);
  catch 
    s = lasterror;
    warndlg(['Problem saving: ' s.message] );
    error(s);
    return
  end

  
  
function edit_change_tspan_Callback(hObject, eventdata, handles)
% hObject    handle to edit_change_tspan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  
  try
  
    handles.tspan= eval(get(hObject,'String'));

    if (handles.auto_recompute)
      handles= update_trajectories(handles);
    end
  
    guidata(hObject,handles);
  catch 
    s = lasterror;
    warndlg(['Problem change_tspan: ' s.message] );
    error(s);
    return
  end

  
% Hints: get(hObject,'String') returns contents of edit_change_tspan as text
%        str2double(get(hObject,'String')) returns contents of edit_change_tspan as a double


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
  
% Hint: get(hObject,'Value') returns toggle state of plot_pts_checkbox




% --- Executes on button press in for_all_checkbox.
function for_all_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to for_all_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

  handles.change_tspan_for_all  = get(hObject,'Value');
  guidata(hObject, handles);
  
% Hint: get(hObject,'Value') returns toggle state of for_all_checkbox


% --- Executes on button press in sensibutton1.
function sensibutton1_Callback(hObject, eventdata, handles)
% hObject    handle to sensibutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = get(hObject,'Value');
handles.plot_sensi(1) = val;
if (val)
  set(handles.param12,'String',{'' handles.param_list{:}});
  set(handles.param12,'Value',1);
  set(handles.param13,'String',{'' handles.param_list{:}});
  set(handles.param13,'Value',1);
  handles.current_sensi{1,2} = '';
  handles.current_sensi{1,3} = '';
else
  set(handles.param12,'String',{'' handles.var_list{:}});
  set(handles.param12,'Value',1);
  set(handles.param13,'String',{'' handles.var_list{:}});
  set(handles.param13,'Value',1);
  handles.current_var{1,2} = '';
  handles.current_var{1,3} = '';
end

new_plot = plot_param(handles,1);
handles.current_plot{1} = new_plot;

guidata(hObject,handles);

% Hint: get(hObject,'Value') returns toggle state of sensibutton1


% --- Executes on button press in sensi_button2.
function sensi_button2_Callback(hObject, eventdata, handles)
% hObject    handle to sensi_button2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

  val = get(hObject,'Value');
  handles.plot_sensi(2) = val;
  if (val)
    set(handles.param22,'String',{'' handles.param_list{:}});
    set(handles.param22,'Value',1);
    set(handles.param23,'String',{'' handles.param_list{:}});
    set(handles.param23,'Value',1);
    handles.current_sensi{2,2} = '';
    handles.current_sensi{2,3} = '';
  else
    set(handles.param22,'String',{'' handles.var_list{:}});
    set(handles.param22,'Value',1);
    set(handles.param23,'String',{'' handles.var_list{:}});
    set(handles.param23,'Value',1);
    handles.current_var{2,2} = '';
    handles.current_var{2,3} = '';
    
  end

  new_plot = plot_param(handles,2);
  handles.current_plot{2} = new_plot;

  guidata(hObject,handles);

  
% Hint: get(hObject,'Value') returns toggle state of sensi_button2


% --- Executes on button press in sensi_button3.
function sensi_button3_Callback(hObject, eventdata, handles)
% hObject    handle to sensi_button3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  
  val = get(hObject,'Value');
  handles.plot_sensi(3) = val;
  if (val)    
    set(handles.param32,'String',{'' handles.param_list{:}});
    set(handles.param32,'Value',1);
    set(handles.param33,'String',{'' handles.param_list{:}});
    set(handles.param33,'Value',1);
    handles.current_sensi{3,2} = '';
    handles.current_sensi{3,3} = '';
  else
    set(handles.param32,'String',{'' handles.var_list{:}});
    set(handles.param32,'Value',1);
    set(handles.param33,'String',{'' handles.var_list{:}});
    set(handles.param33,'Value',1);
    handles.current_var{3,2} = '';
    handles.current_var{3,3} = '';
  
  end
  
  new_plot = plot_param(handles,3);
  handles.current_plot{3} = new_plot;
  
  guidata(hObject,handles);

  % Hint: get(hObject,'Value') returns toggle state of sensi_button3


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
  
% Hints: contents = get(hObject,'String') returns popup_pts1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_pts1


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

  
% Hints: contents = get(hObject,'String') returns popup_pts2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_pts2


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

% Hints: contents = get(hObject,'String') returns popup_pts3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_pts3


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
  if (~isfield(handles.TrajSet,'traj'))
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
        
    %    SplotPts(handles.TrajSet,param_to_plot);  
    S = DiscrimPropValues(handles.TrajSet);
    SplotPts(S,param_to_plot);  
    SplotBoxPts(handles.TrajSet, param_to_plot,handles.current_pts,'+k','r',.1);
  
  end


function new_plot= plot_param(handles,ax)
  
  if (isfield(handles.Sys,'time_mult'))
    time_mult = handles.Sys.time_mult;
  else
    time_mult=1;
  end

  
  if (~isfield(handles.TrajSet,'traj'))
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
    if (handles.plot_sensi(ax)) % plot sensitivity of a property
      prop = handles.properties.(pnames{nprop});
      cla;legend('off');
      hold on;
      grid on;
      st = param_to_plot{1};
      ipts = handles.current_pts;
      %      phi_val = handles.TrajSet.props_values(nprop,ipts).val;
      %      phi_tspan = handles.TrajSet.props_values(nprop,ipts).tspan;
        
      sensi_param = [];
      if (~strcmp(handles.current_sensi{ax,2},''))
        sensi_param = [sensi_param  FindParam(handles.TrajSet,handles.current_sensi{ax,2})];
      end

      if (~strcmp(handles.current_sensi{ax,3},''))
        sensi_param = [sensi_param  FindParam(handles.TrajSet,handles.current_sensi{ax,3})];
      end
      colors = {'b','r'};
      if (numel(sensi_param)>0)

        % if needed, recompute sensitivity of the variable and the property 
        plot_done =0;
        for is=1:numel(sensi_param)
          if isfield(handles.TrajSet.traj(handles.current_pts), 'XS')
            isf =find(handles.TrajSet.traj(handles.current_pts).sensis==sensi_param(is));
            if ~isempty(isf) % ok, no need to recompute
              phi_tspan = handles.TrajSet.traj(handles.traj_ref(ipts)).time;
              [phi_val phivald] = QMITL_EvalSensi(handles.Sys,prop,handles.TrajSet.traj(handles.traj_ref(ipts)),phi_tspan);
              plot(phi_tspan*time_mult, phi_vald,colors{is});            
              plot_done=1;            
            end
          end
        
          % recompute
          if (~plot_done)
              
            if isfield(handles.TrajSet, 'traj')
              tspan = handles.TrajSet.traj(handles.current_pts).time;
            else
              if ~isempty(handles.tspan)
                tspan = handles.tspan;
              else
                tspan = [0:.1:1];
              end
            end
              
            Ptmp = CreateParamSet(handles.Sys,sensi_param(is));
            Ptmp.pts = handles.TrajSet.pts(:,handles.current_pts);
            Pftmp = ComputeTrajSensi(handles.Sys, Ptmp,tspan);          

%            handles.TrajSet.traj(handles.current_pts) = Pftmp.traj;
            handles.TrajSet.Xf(:,handles.current_pts) = Pftmp.Xf;
            
            phi_tspan = handles.TrajSet.traj(handles.traj_ref(ipts)).time;
            [phi_val phi_vald] = QMITL_EvalSensi(prop,Pftmp.traj,phi_tspan);
            plot(phi_tspan*time_mult, phi_vald,colors{is});
          end
        
        end %for        
        ylabel(['sensi(' short_disp(prop,20) ')'],'Interpreter','none');
        legend(['sensi(' short_disp(prop,100) ')'],'Interpreter','none');
        xlabel(['time']);
             
      end %  if (numel(sensi_param)>0)
    else      
      
      prop = handles.properties.(pnames{nprop});
      cla;legend('off');
      hold on;
      grid on;
      st = param_to_plot{1};
      ipts = handles.current_pts;
      phi_tspan = handles.TrajSet.traj(handles.traj_ref(ipts)).time;      
      
      % checks if sensitivities needs be (re)computed for the formula
      
      %% OBSOLETE, TO REMOVE OR UPDATE...
      is = QMITL_ExtractSensi(prop);
      if (~isempty(is))        
        Ptmp = CreateParamSet(handles.Sys,is);
        Ptmp.pts = handles.TrajSet.pts(:,handles.current_pts);
        Pftmp = ComputeTrajSensi(handles.Sys, Ptmp, handles.TrajSet.traj(handles.traj_ref(ipts)).time);          
        phi_val = QMITL_Eval(handles.Sys,prop,Ptmp, Pftmp.traj,phi_tspan);
      else        
        Ptmp = Sselect(handles.TrajSet, ipts);
        phi_val = QMITL_Eval(handles.Sys,prop,Ptmp,handles.TrajSet.traj(handles.traj_ref(ipts)),phi_tspan);   
      end
                 
      
      ylabel(get_id(prop),'Interpreter','none');
      xlabel(['time']);
            
      new_plot = plot(phi_tspan*time_mult, phi_val); 
      hold on; 
      plot([phi_tspan(1) phi_tspan(end)]*time_mult, [0 0],'-k');
      stairs(phi_tspan*time_mult, (phi_val>0)*max(abs(phi_val))/2,'-r');
      lgh = legend(short_disp(prop,100));           
      set(lgh,'Interpreter','none');
      
    end % if (handles.plot_sensi(ax))
    new_plot = [];
    
  else  % plot values for a variable      
      
    if (handles.plot_sensi(ax)) % plot sensitivity of a variable

      cla;legend('off');
      sensi_param = [];
      if (~strcmp(handles.current_sensi{ax,2},''))
        sensi_param = [sensi_param  FindParam(handles.TrajSet,handles.current_sensi{ax,2})];
      end

      if (~strcmp(handles.current_sensi{ax,3},''))
        sensi_param = [sensi_param  FindParam(handles.TrajSet,handles.current_sensi{ax,3})];
      end
 
      if (numel(sensi_param)>0)
        Ptmp = CreateParamSet(handles.Sys,sensi_param);
        Ptmp.pts = handles.TrajSet.pts(:,handles.current_pts);
        if isfield(handles.TrajSet, 'traj')
          tspan = handles.TrajSet.traj(handles.current_pts).time;
        else
          if ~isempty(handles.tspan)
            tspan = handles.tspan;
          else
            tspan = [0:.1:1];
          end
        end
        Pftmp = ComputeTrajSensi(handles.Sys, Ptmp, tspan);

        handles.TrajSet.Xf(:,handles.current_pts) = Pftmp.Xf;
        Pftmp.time_mult = time_mult;
        SplotSensi(Pftmp, param_to_plot, sensi_param);

      end
      new_plot=[];
    
    else
      
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

      SplotTraj(handles.TrajSet, param_to_plot,handles.current_pts, handles.traj_opt);
      children = get(gca, 'Children');
      if (numel(param_to_plot)>1)
        new_plot=  children(1:2);
      else
        new_plot=  children(1);
      end
    end % if (handles.plot_sensi(ax))
  end % if (nprop) 

  
  
function handles= plot_tout(handles)
   
  if (~isfield(handles.TrajSet,'traj'))
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
    
    % Axes 1
    axes(handles.axes1);
    param_to_plot = handles.current_var{1,1};
    
    if (find_prop(param_to_plot{1}, pnames)==0)&&(~handles.plot_sensi(1))
      
      if (~strcmp(handles.current_var{1,2},''))
        param_to_plot = {param_to_plot{:} handles.current_var{1,2}};
      end

      if (~strcmp(handles.current_var{1,3},''))
        param_to_plot = {param_to_plot{:} handles.current_var{1,3}};
      end

      SplotTraj(handles.TrajSet, param_to_plot);
    end
    
    % Axes2
    axes(handles.axes2);
    param_to_plot = handles.current_var{2,1};
    if (find_prop(param_to_plot{1}, pnames)==0)&&(~handles.plot_sensi(2))
      if (~strcmp(param_to_plot,''))
        if (~strcmp(handles.current_var{2,2},''))
          param_to_plot = {param_to_plot{:} handles.current_var{2,2}};
        end
        if (~strcmp(handles.current_var{2,3},''))
          param_to_plot = {param_to_plot{:} handles.current_var{2,3}};
        end
        SplotTraj(handles.TrajSet, param_to_plot);
      end
    end
    % Axes 3
    axes(handles.axes3);
    param_to_plot = handles.current_var{3,1};
  
    if (find_prop(param_to_plot{1}, pnames)==0)&&(~handles.plot_sensi(3))
      if (~strcmp(param_to_plot,''))
        if (~strcmp(handles.current_var{3,2},''))
          param_to_plot = {param_to_plot{:} handles.current_var{3,2}};
        end

        if (~strcmp(handles.current_var{3,3},''))
          param_to_plot = {param_to_plot{:} handles.current_var{3,3}};
        end
      end
      SplotTraj(handles.TrajSet, param_to_plot);    
    end
    % Other Axes

    for i = 4:numel(handles.exported_axes)+3
      try
        axes(handles.exported_axes(i))
      catch        
        continue;
      end
      param_to_plot = handles.current_var{i,1};
      if (find_prop(param_to_plot{1}, fnames)==0)&&(~handles.plot_sensi(i))
        if (~strcmp(param_to_plot,''))
          if (~strcmp(handles.current_var{i,2},''))
            param_to_plot = {param_to_plot{:} handles.current_var{i,2}};
          end
        
          if (~strcmp(handles.current_var{i,3},''))
            param_to_plot = {param_to_plot{:} handles.current_var{i,3}};
          end
        end
        SplotTraj(handles.TrajSet, param_to_plot);    
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
  
  
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


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
% Hints: get(hObject,'String') returns contents of edit_max_slider_param as text
%        str2double(get(hObject,'String')) returns contents of edit_max_slider_param as a double


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
   
  
% Hints: get(hObject,'String') returns contents of edit_min_slider_param as text
%        str2double(get(hObject,'String')) returns contents of edit_min_slider_param as a double


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
       tspan = handles.TrajSet.traj(handles.current_pts).time;       
    catch
      tspan= [0 1];             
    end
    
  else
    tspan = handles.tspan;
  end
  
  % check if we need to recompute everything
  if (get(handles.for_all_checkbox,'Value')||~isfield(handles.TrajSet,'traj'))
    
    handles.TrajSet = ComputeTraj(handles.Sys, handles.TrajSet, tspan);
    traj_ref = handles.TrajSet.traj_ref;
   
    if isfield(handles.TrajSet,'props')
      for i=1:numel(handles.TrajSet.props)
        prop = handles.TrajSet.props(i);
        for j= 1:size(handles.TrajSet.pts,2)          
          traj = handles.TrajSet.traj(traj_ref(j));
          P = Sselect(handles.TrajSet,j);
          handles.TrajSet.props_values(i,j).val = ...
          QMITL_Eval(handles.Sys,props, P, traj,handles.TrajSet.props_values(i,j).tspan);         
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
    
    Ptmp = Sselect(handles.TrajSet,handles.current_pts);
    
    % check if sensi had been computed for this traj
           
    if (handles.has_dynamics==1)
      if isfield(Ptmp.traj, 'sensis')
        Ptmp = SPurge(Ptmp);
        Pftmp = ComputeTrajSensi(handles.Sys, Ptmp, tspan);
      else
        Ptmp = SPurge(Ptmp);
        Pftmp = ComputeTraj(handles.Sys, Ptmp, tspan);
      end
    else
        Pftmp = Ptmp;
    end
    traj_ref = handles.TrajSet.traj_ref;
    handles.TrajSet.traj(traj_ref(handles.current_pts)) = Pftmp.traj;
    handles.TrajSet.Xf(:,traj_ref(handles.current_pts)) = Pftmp.traj.X(:,end);
    
    new_plot = plot_param(handles,1);
    handles.current_plot{1} = new_plot;
    
    new_plot = plot_param(handles,2);
    handles.current_plot{2} = new_plot;
  
    new_plot = plot_param(handles,3);
    handles.current_plot{3} = new_plot;
 
  end
  
function handles =  update_listbox_param(handles, changed)
  
% update values and plots 
  
  if nargin==1
    changed =0;
  end
  
  lbda = handles.slider_param_coeff;
  ind_p = handles.current_change_param;
  
  if (ind_p)
    if (changed)
      new_val = get(handles.slider_param,'Value')*lbda;
      if (ind_p<=numel(handles.TrajSet.ParamList))
        if (get(handles.for_all_checkbox,'Value'))        
          handles.TrajSet.pts(ind_p, :) = new_val;
          if (ind_p<=handles.TrajSet.DimP)
            for i = 1:numel(handles.TrajSet.traj)
              handles.TrajSet.traj(i).param(ind_p) = new_val;
            end
          end
        else
          handles.TrajSet.pts(ind_p, handles.current_pts) = new_val;      
          if (ind_p<=handles.TrajSet.DimP)
            handles.TrajSet.traj(handles.current_pts).param(ind_p) = new_val;
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
        if (ind_p<=handles.TrajSet.DimP)
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

  for i=1:numel(handles.TrajSet.dim)
    st = handles.TrajSet.ParamList{handles.TrajSet.dim(i)};
    st = strcat(st, ':',' ',dbl2str(handles.TrajSet.pts(handles.TrajSet.dim(i), handles.current_pts)));
    handles.current_varying_param{i} = st;
  end
  
  content = {content{:} handles.current_varying_param{:} '' 'Systems and props parameters' '-------------------'};

  for i=1:numel(handles.TrajSet.ParamList)
    st = handles.TrajSet.ParamList{i};
    st = strcat(st, ':',' ',dbl2str(handles.TrajSet.pts(i, handles.current_pts)));
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
