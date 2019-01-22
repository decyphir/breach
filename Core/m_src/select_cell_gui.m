function varargout = select_cell_gui(varargin)
% SELECT_CELL_GUI M-file for select_cell_gui.fig
%      SELECT_CELL_GUI, by itself, creates a new SELECT_CELL_GUI or raises the existing
%      singleton*.
%
%      H = SELECT_CELL_GUI returns the handle to a new SELECT_CELL_GUI or the handle to
%      the existing singleton*.
%
%      SELECT_CELL_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECT_CELL_GUI.M with the given input arguments.
%
%      SELECT_CELL_GUI('Property','Value',...) creates a new SELECT_CELL_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before select_cell_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to select_cell_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help select_cell_gui

% Last Modified by GUIDE v2.5 28-Sep-2018 14:53:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @select_cell_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @select_cell_gui_OutputFcn, ...
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


% --- Executes just before select_cell_gui is made visible.
function select_cell_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to select_cell_gui (see VARARGIN)

if ismac
    FONT=12;
    %POS = [60 10 120 40];
else
    FONT=10;
    %POS = [60 10 120 40];
end

hfn = fieldnames(handles);
for ifn = 1:numel(hfn)
    try
        set(handles.(hfn{ifn}), 'FontSize', FONT);
    end
end
%set(handles.main, 'Position',POS);
 
set(hObject, 'Name', ['Choose from list']);
  
% gui takes three arguments : a parameter set and default options
  content_all = varargin{1};
  content_select = varargin{2};
  
  % foolproof - make sure content_all contains content_select
  content_all = union(content_all,content_select);
  
  prompt = varargin{3};
  set(hObject, 'Name', prompt);

  handles.content_really_all = content_all;
  handles.content_all = content_all;   

  % fill variable list 
  set(handles.listbox_all_variables,'String',content_all);
 
  % fill selected var  
  set(handles.listbox_selected,'String',content_select);
    
  % Choose default command line output for select_cell_gui
  handles.output = content_all;
  
  % Init configurable edit boxes
  set(handles.edit1, 'UserData', containers.Map());
  set(handles.edit2, 'UserData', containers.Map());
  
  % focus on listbox with all 
  uicontrol(handles.listbox_all_variables);

  % 
  handles= update_selected(handles);
  handles = update_custom_edits(handles);
  
  % wait or not (controlled by caller)
  if numel(varargin)>=4
      handles.wait = varargin{4};
  else
      handles.wait = true;
  end
  
  % Update handles structure
  guidata(hObject, handles);
  
% UIWAIT makes select_cell_gui wait for user response (see UIRESUME)
if handles.wait
    uiwait(handles.main);
end



% --- Outputs from this function are returned to the command line.
function varargout = select_cell_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output';
varargout{2} = handles;

% The figure can be deleted now
if handles.wait
    delete(handles.main);
end

% --- Executes on selection change in listbox_all_variables.
function listbox_all_variables_Callback(hObject, eventdata, handles)

  handles= update_selected(handles);
  handles = update_custom_edits(handles);
  guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function listbox_all_variables_CreateFcn(hObject, eventdata, handles)

  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end
  

% --- Executes on selection change in listbox_selected.
function listbox_selected_Callback(hObject, eventdata, handles)

  handles= update_selected(handles);
  handles = update_custom_edits(handles);
  guidata(hObject, handles);
 

% --- Executes during object creation, after setting all properties.
function listbox_selected_CreateFcn(hObject, eventdata, handles)

  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end


% --- Executes on button press in button_add.
function button_add_Callback(hObject, eventdata, handles)

  content_all = get(handles.listbox_all_variables, 'String');      
  content_select = get(handles.listbox_selected, 'String');
  vall = get(handles.listbox_all_variables,'Value');
  vtp = get(handles.listbox_selected,'Value');
  niou_cont= content_all(vall);
  content_select = union(content_select, niou_cont); 
  set(handles.listbox_selected, 'String', content_select);  
  handles= update_selected(handles);
  handles = update_custom_edits(handles);
  guidata(hObject, handles);
  
% --- Executes on button press in button_rem.
function button_rem_Callback(hObject, eventdata, handles)
handles = rem(handles);
handles= update_selected(handles);
handles = update_custom_edits(handles);
guidata(hObject, handles);

function handles= rem(handles)
content_select = get(handles.listbox_selected, 'String');
if numel(content_select)==0
    return;
end

vtp = get(handles.listbox_selected,'Value');
nb = numel(content_select);

new_idx = setdiff(1:nb, vtp);
content_select = content_select(new_idx);

set(handles.listbox_selected, 'String', content_select);
set(handles.listbox_selected, 'Value', max(vtp(1)-1,1));

% --- Executes on button press in button_ok.
function button_ok_Callback(hObject, eventdata, handles)
 
  content_selected = get(handles.listbox_selected, 'String');
  handles.output= content_selected ;
  guidata(hObject, handles);
  uiresume(handles.main);
 
% --- Executes on button press in button_cancel.
function button_cancel_Callback(hObject, eventdata, handles)
  handles.output=0;
  set(hObject, 'UserData', true);
  guidata(hObject, handles);
  uiresume(handles.main);

% --- Executes when user attempts to close main.
function main_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end

function handles = add(handles)
% code of button_add_Callback
content_all = get(handles.listbox_all_variables, 'String');
content_select = get(handles.listbox_selected, 'String');
vall = get(handles.listbox_all_variables,'Value');
vtp = get(handles.listbox_selected,'Value');
niou_cont= content_all(vall);
content_select = union(content_select, niou_cont);
handles= update_selected(handles);
handles = update_custom_edits(handles);

set(handles.listbox_selected, 'String', content_select);


% --- Executes on key press with focus on main and none of its controls.
function main_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to main (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
  
% --- Executes on key press with focus on listbox_all_variables and none of its controls.
function listbox_all_variables_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to listbox_all_variables (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if (isa(eventdata, 'matlab.ui.eventdata.UIClientComponentKeyEvent'))
    switch eventdata.Key
     case 'rightarrow'
     handles= add(handles);
      guidata(hObject, handles);
      
     case 'leftarrow'
      handles= rem(handles);
      guidata(hObject, handles);
         
     case 'return'
      content_selected = get(handles.listbox_selected, 'String');
      handles.output = content_selected ;
  
      guidata(hObject, handles);
      uiresume(handles.main);

     case 'escape'
      handles.output=[];
      uiresume(handles.main);
    %  otherwise 
    %  eventdata      
    end 
end

% --- Executes on button press in pushbutton_add_all.
function pushbutton_add_all_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_add_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

content_all = get(handles.listbox_all_variables,'String');
set(handles.listbox_selected,'String',content_all);
set(handles.listbox_selected,'Value',1);
handles= update_selected(handles);
handles = update_custom_edits(handles);
guidata(hObject, handles);
    
% --- Executes on button press in pushbutton_rem_all.
function pushbutton_rem_all_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_rem_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.listbox_selected,'String',{});
set(handles.listbox_selected,'Value',1);
handles= update_selected(handles);
handles = update_custom_edits(handles);
guidata(hObject, handles);

% --- Executes on key press with focus on listbox_selected and none of its controls.
function listbox_selected_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to listbox_selected (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if (isa(eventdata, 'matlab.ui.eventdata.UIClientComponentKeyEvent'))
    switch eventdata.Key
        case 'rightarrow'
            handles= add(handles);
            
            handles= update_selected(handles);
            handles = update_custom_edits(handles);
            
            guidata(hObject, handles);
            
        case 'leftarrow'
            handles= rem(handles);
            handles= update_selected(handles);
            handles = update_custom_edits(handles);
            guidata(hObject, handles);
            
        case 'return'
            
            content_selected = get(handles.listbox_selected, 'String');
            handles.output = content_selected ;
            handles= update_selected(handles);
            handles = update_custom_edits(handles);
            
            guidata(hObject, handles);
            uiresume(handles.main);
            
        case 'escape'
            handles.output=[];
            uiresume(handles.main);
            %  otherwise
            %  eventdata
    end
end

function edit_search_Callback(hObject, eventdata, handles)
% hObject    handle to edit_search (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

exp = get(hObject, 'String');
if isempty(exp)
    handles.content_all = handles.content_really_all;
    set(handles.listbox_all_variables, 'String', handles.content_all);
    
else
    content_all = get(handles.listbox_all_variables, 'String');
    if ~isempty(content_all)
        val_all = get(handles.listbox_all_variables, 'Value');
        content_all_select= content_all(val_all);
    else
        content_all_select = {};
    end
    s = cellfun(@(c)isempty(c), regexp(handles.content_really_all, exp));
    handles.content_all = handles.content_really_all(~s);
    set(handles.listbox_all_variables, 'Value',1);
    set(handles.listbox_all_variables, 'String', handles.content_all);
    if ~isempty(content_all_select)
    nval_all = find(strcmp(handles.content_all, content_all_select));
    if isempty(nval_all)
        nval_all = 1;
    end
    else
        nval_all =1;
    end
    set(handles.listbox_all_variables, 'Value', nval_all);
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_search_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_search (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

var = handles.selected_var;
 
if ~isempty(var)
    m = get(hObject, 'UserData');
    st=  get(hObject, 'String');
    if ~isempty(st)
        m(var)= strtrim(st);
    end
else
    set(hObject, 'String', '');
end


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

var = handles.selected_var;
if ~isempty(var)
    m = get(hObject, 'UserData');
    st=  get(hObject, 'String');
    if ~isempty(st)
        m(var)= strtrim(st);
    end
else
    set(hObject, 'String', '');
end

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function handles =update_custom_edits(handles)
v= handles.selected_var;
set(handles.edit1, 'String','');
set(handles.edit2, 'String', '');
if ~isempty(v)
    m1 = get(handles.edit1,'UserData');
    if m1.isKey(v)
        set(handles.edit1, 'String', m1(v));
    end
    
    m2 = get(handles.edit2,'UserData');
    if m2.isKey(v)
        set(handles.edit2, 'String', m2(v));
    end
        
end

function handles = update_selected(handles)

% find focus? 
h = gco(handles.main);
if ~strcmp(get(h, 'Style'), 'listbox')||isempty(get(h, 'String'))
    h = handles.listbox_all_variables;
end

val= get(h,'Value');
content_all = get(h, 'String');
if ~isempty(content_all)
    handles.selected_var = content_all{val};
else
    handles.selected_var = '';
end
