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

% Last Modified by GUIDE v2.5 19-Feb-2017 16:45:11

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

  set(hObject, 'Name', ['Choose from list']);
  
% gui takes two arguments : a parameter set and default options
  content_all = varargin{1};
  prompt = varargin{2};
  set(hObject, 'Name', prompt);
  
  
  handles.content_all = content_all;   

  % fill variable list 
  set(handles.listbox_all_variables,'String',content_all);
 
  % fill selected var  
  set(handles.listbox_selected,'String',content_all);
    
  % Choose default command line output for select_cell_gui
  handles.output = content_all;
  
  % focus on listbox with all 
  uicontrol(handles.listbox_all_variables);
  
  % Update handles structure
  guidata(hObject, handles);
  
% UIWAIT makes select_cell_gui wait for user response (see UIRESUME)
  uiwait(handles.main);


% --- Outputs from this function are returned to the command line.
function varargout = select_cell_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output';

% The figure can be deleted now
delete(handles.main);

% --- Executes on selection change in listbox_all_variables.
function listbox_all_variables_Callback(hObject, eventdata, handles)

  val= get(hObject,'Value');
  handles.selected_var = val;  
  guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function listbox_all_variables_CreateFcn(hObject, eventdata, handles)

  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end
  

% --- Executes on selection change in listbox_selected.
function listbox_selected_Callback(hObject, eventdata, handles)
  val= get(hObject,'Value');
  handles.selected_var_to_plot = val;  
  guidata(hObject, handles);
  
% Hints: contents = get(hObject,'String') returns listbox_selected contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_selected


% --- Executes during object creation, after setting all properties.
function listbox_selected_CreateFcn(hObject, eventdata, handles)

  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end


% --- Executes on button press in button_add.
function button_add_Callback(hObject, eventdata, handles)

  content_all = get(handles.listbox_all_variables, 'String');      
  content_to_plot = get(handles.listbox_selected, 'String');
  
  vall = get(handles.listbox_all_variables,'Value');
  vtp = get(handles.listbox_selected,'Value');
  
  niou_cont= content_all{vall};
  
  for i=1:numel(content_to_plot)
    if strcmp(content_to_plot{i}, niou_cont)
      return
    end
  end
  content_to_plot = {content_to_plot{:}, niou_cont}; 
  set(handles.listbox_selected, 'String', content_to_plot);
  
  guidata(hObject, handles);

% --- Executes on button press in button_rem.
function button_rem_Callback(hObject, eventdata, handles)
  
  content_to_plot = get(handles.listbox_selected, 'String');
  vtp = get(handles.listbox_selected,'Value');
  nb = numel(content_to_plot);
  
  if numel(content_to_plot)==1
    return;
  end
 
  if (vtp==1)
    content_to_plot = {content_to_plot{2:end}};
  elseif (vtp == nb) 
    content_to_plot = { content_to_plot{1:vtp-1}};
    set(handles.listbox_selected,'Value',vtp-1);
  else
    content_to_plot = { content_to_plot{1:vtp-1} content_to_plot{vtp+1:end}};
  end
  set(handles.listbox_selected, 'String', content_to_plot);
  
  guidata(hObject, handles);
  
 
% --- Executes on button press in button_ok.
function button_ok_Callback(hObject, eventdata, handles)
 
  content_selected = get(handles.listbox_selected, 'String');
  handles.output= content_selected ;
  guidata(hObject, handles);
  uiresume(handles.main);
 

% --- Executes on button press in button_cancel.
function button_cancel_Callback(hObject, eventdata, handles)

  handles.output=[];
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
     
      % code of button_add_Callback
      content_all = get(handles.listbox_all_variables, 'String');      
      content_selected = get(handles.listbox_selected, 'String');
      
      vall = get(handles.listbox_all_variables,'Value');
      niou_cont= content_all{vall};
      
      for i=1:numel(content_selected)
        if strcmp(content_selected{i}, niou_cont)
          return
        end
      end
      content_selected = {content_selected{:}, niou_cont}; 
      set(handles.listbox_selected, 'String', content_selected);
     
      %set(handles.listbox_all_variables, 'Value',1);
      
      guidata(hObject, handles);
      
     case 'leftarrow'
       
      content_all = get(handles.listbox_all_variables, 'String');      
      content_selected = get(handles.listbox_selected, 'String');
      
      vall = get(handles.listbox_all_variables,'Value');
      niou_cont= content_all{vall};
      
      vtpl =0;
      for i=1:numel(content_selected)
        if strcmp(content_selected{i}, niou_cont)
          vtpl = i;
          break
        end
      end
      
      if (vtpl == 0)
        return;
      end
      if numel(content_selected)==1
        return;
      end
      
      if (vtpl==1)
        content_selected = {content_selected{2:end}};
      else
        content_selected = { content_selected{1:vtpl-1} content_selected{vtpl+1: ...
                            end}};
      end
      set(handles.listbox_selected, 'String', content_selected);
      set(handles.listbox_selected, 'Value', 1);
      
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




% --- Executes on button press in pushbutton_rem_all.
function pushbutton_rem_all_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_rem_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
