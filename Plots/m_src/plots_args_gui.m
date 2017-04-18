function varargout = plots_args_gui(varargin)
% PLOTS_ARGS_GUI M-file for plots_args_gui.fig
%      PLOTS_ARGS_GUI, by itself, creates a new PLOTS_ARGS_GUI or raises the existing
%      singleton*.
%
%      H = PLOTS_ARGS_GUI returns the handle to a new PLOTS_ARGS_GUI or the handle to
%      the existing singleton*.
%
%      PLOTS_ARGS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOTS_ARGS_GUI.M with the given input arguments.
%
%      PLOTS_ARGS_GUI('Property','Value',...) creates a new PLOTS_ARGS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plots_args_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plots_args_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plots_args_gui

% Last Modified by GUIDE v2.5 18-Apr-2017 16:32:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plots_args_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @plots_args_gui_OutputFcn, ...
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


% --- Executes just before plots_args_gui is made visible.
function plots_args_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plots_args_gui (see VARARGIN)


% Set fonts and size depending on system
if ismac
    FONT=12;
    POS = [50 10 200 50];
else
    FONT=10;
    POS = [50 10 200 50];
end

hfn = fieldnames(handles);
for ifn = 1:numel(hfn)
    try 
        set(handles.(hfn{ifn}), 'FontSize', FONT);
    end
end
set(handles.main, 'Position',POS);

  set(hObject, 'Name', ['Choose variables to plot and options']);
  
% gui takes two arguments : a parameter set and default options
  P = varargin{1};
  args = varargin{2};
  
  % fill variable list
  content = {};  
  for i=1:P.DimX
    st = P.ParamList{i}; 
    content = {content{:} st};
  end
  
  set(handles.listbox_all_variables,'String',content);
 
  % fill var to plot llist
 
  if isempty(args)
    set(handles.listbox_var_to_plot,'String',{content{1}});
  else
    set(handles.listbox_var_to_plot,'String', args.iX);
    set(handles.checkbox_same_axe, 'Value', args.same_axe)
    set(handles.checkbox_phase_portrait, 'Value',args.phase_portrait)
  
  end
    
  % Choose default command line output for plots_args_gui
  handles.output = [];
  
  
  % Update handles structure
  guidata(hObject, handles);

  
% UIWAIT makes plots_args_gui wait for user response (see UIRESUME)
  uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = plots_args_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% The figure can be deleted now
delete(handles.figure1);

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
  

% --- Executes on selection change in listbox_var_to_plot.
function listbox_var_to_plot_Callback(hObject, eventdata, handles)
  val= get(hObject,'Value');
  handles.selected_var_to_plot = val;  
  guidata(hObject, handles);
  
% Hints: contents = get(hObject,'String') returns listbox_var_to_plot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_var_to_plot


% --- Executes during object creation, after setting all properties.
function listbox_var_to_plot_CreateFcn(hObject, eventdata, handles)

  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end


% --- Executes on button press in button_add.
function button_add_Callback(hObject, eventdata, handles)

  content_all = get(handles.listbox_all_variables, 'String');      
  content_to_plot = get(handles.listbox_var_to_plot, 'String');
  
  vall = get(handles.listbox_all_variables,'Value');
  vtp = get(handles.listbox_var_to_plot,'Value');
  
  niou_cont= content_all{vall};
  
  for i=1:numel(content_to_plot)
    if strcmp(content_to_plot{i}, niou_cont)
      return
    end
  end
  content_to_plot = {content_to_plot{:}, niou_cont}; 
  set(handles.listbox_var_to_plot, 'String', content_to_plot);
  
  guidata(hObject, handles);

% --- Executes on button press in button_rem.
function button_rem_Callback(hObject, eventdata, handles)
  
  content_to_plot = get(handles.listbox_var_to_plot, 'String');
  vtp = get(handles.listbox_var_to_plot,'Value');
  nb = numel(content_to_plot);
  
  
  if numel(content_to_plot)==1
    return;
  end
  
  
  if (vtp==1)
    content_to_plot = {content_to_plot{2:end}};
  elseif (vtp == nb) 
    content_to_plot = { content_to_plot{1:vtp-1}};
    set(handles.listbox_var_to_plot,'Value',vtp-1);
  else
    content_to_plot = { content_to_plot{1:vtp-1} content_to_plot{vtp+1:end}};
  end
  set(handles.listbox_var_to_plot, 'String', content_to_plot);
  
  guidata(hObject, handles);
  
  
% --- Executes on button press in checkbox_phase_portrait.
function checkbox_phase_portrait_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of checkbox_phase_portrait


% --- Executes on button press in checkbox_same_axe.
function checkbox_same_axe_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of checkbox_same_axe


% --- Executes on button press in button_ok.
function button_ok_Callback(hObject, eventdata, handles)
 
  content_to_plot = get(handles.listbox_var_to_plot, 'String');
  handles.output.iX = content_to_plot ;
  handles.output.same_axe = get(handles.checkbox_same_axe, 'Value'); 
  handles.output.phase_portrait = get(handles.checkbox_phase_portrait, 'Value'); 
  
  guidata(hObject, handles);
  close(handles.figure1);
 

% --- Executes on button press in button_cancel.
function button_cancel_Callback(hObject, eventdata, handles)

  handles.output=[];
  close(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end

  

% --- Executes during object deletion, before destroying properties.
function main_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
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
      content_to_plot = get(handles.listbox_var_to_plot, 'String');
      
      vall = get(handles.listbox_all_variables,'Value');
      niou_cont= content_all{vall};
      
      for i=1:numel(content_to_plot)
        if strcmp(content_to_plot{i}, niou_cont)
          return
        end
      end
      content_to_plot = {content_to_plot{:}, niou_cont}; 
      set(handles.listbox_var_to_plot, 'String', content_to_plot);
     
      %set(handles.listbox_all_variables, 'Value',1);
      
      guidata(hObject, handles);
      
     case 'leftarrow'
       
      content_all = get(handles.listbox_all_variables, 'String');      
      content_to_plot = get(handles.listbox_var_to_plot, 'String');
      
      vall = get(handles.listbox_all_variables,'Value');
      niou_cont= content_all{vall};
      
      vtpl =0;
      for i=1:numel(content_to_plot)
        if strcmp(content_to_plot{i}, niou_cont)
          vtpl = i;
          break
        end
      end
      
      
      if (vtpl == 0)
        return;
      end
      if numel(content_to_plot)==1
        return;
      end
      
      if (vtpl==1)
        content_to_plot = {content_to_plot{2:end}};
      else
        content_to_plot = { content_to_plot{1:vtpl-1} content_to_plot{vtpl+1: ...
                            end}};
      end
      set(handles.listbox_var_to_plot, 'String', content_to_plot);
      set(handles.listbox_var_to_plot, 'Value', 1);
      
      guidata(hObject, handles);      
   
     case 'return'
      
      content_to_plot = get(handles.listbox_var_to_plot, 'String');
      handles.output.iX = content_to_plot ;
      handles.output.same_axe = get(handles.checkbox_same_axe, 'Value'); 
      handles.output.phase_portrait = get(handles.checkbox_phase_portrait, 'Value'); 
  
      guidata(hObject, handles);
      close(handles.figure1);

     case 'escape'
      handles.output=[];
      close(handles.figure1);
    %  otherwise 
    %  eventdata      
    end 
end
