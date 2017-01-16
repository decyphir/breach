function varargout = STL_TemplateGUI(varargin)
% STL_TEMPLATEGUI MATLAB code for STL_TemplateGUI.fig
%      STL_TEMPLATEGUI, by itself, creates a new STL_TEMPLATEGUI or raises the existing
%      singleton*.
%
%      H = STL_TEMPLATEGUI returns the handle to a new STL_TEMPLATEGUI or the handle to
%      the existing singleton*.
%
%      STL_TEMPLATEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STL_TEMPLATEGUI.M with the given input arguments.
%
%      STL_TEMPLATEGUI('Property','Value',...) creates a new STL_TEMPLATEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before STL_TemplateGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to STL_TemplateGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help STL_TemplateGUI

% Last Modified by GUIDE v2.5 13-Jul-2016 10:57:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @STL_TemplateGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @STL_TemplateGUI_OutputFcn, ...
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


% --- Executes just before STL_TemplateGUI is made visible.
function STL_TemplateGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to STL_TemplateGUI (see VARARGIN)

global BreachGlobOpt

handles.custom_phi=0;

% Initialize listbox_stlib
stlib = STL_ReadFile('stlib.stl');
set(handles.listbox_stlib, 'String', stlib);

% update template
handles = update_tmp(handles);

% Initialize formula
handles.STL_Formula = handles.STL_Template; 
set(handles.display_STL_Formula,'String', disp(handles.STL_Formula));
handles.phi_id = 'phi_new';
handles.auto_naming =1;
handles= update_phi(handles);


% Signal names for system
handles.signame_sys = varargin{2};
set(handles.popup_signame_concrete,'String', [{''} varargin{2}]);
set(handles.popup_signame_concrete,'Value',1);

% Choose default command line output for STL_TemplateGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes STL_TemplateGUI wait for user response (see UIRESUME)
uiwait(handles.main);

function handles = update_tmp(handles)
global BreachGlobOpt
 
content = get(handles.listbox_stlib, 'String'); 
nbphi =numel(content);

value = get(handles.listbox_stlib, 'Value');
if value>nbphi
    value=1;
end

st_phi_tmp = content{value};
phi_tmp = BreachGlobOpt.STLDB(st_phi_tmp);
handles.STL_Template = phi_tmp;

sig_templates = STL_ExtractSignals(phi_tmp);
handles.sigmap = containers.Map();
handles.signame_templates = sig_templates;

for sig = sig_templates
    handles.sigmap(sig{1}) = 1;
end
set(handles.popup_signame_template,'String', sig_templates);
set(handles.display_STL_Template,'String', disp(handles.STL_Template,0));


% --- Outputs from this function are returned to the command line.
function varargout = STL_TemplateGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
close(handles.main);

% --- Executes on selection change in popup_signame_template.
function popup_signame_template_Callback(hObject, eventdata, handles)
% hObject    handle to popup_signame_template (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

old_names = get(handles.popup_signame_template,'String' );
new_names = get(handles.popup_signame_concrete,'String' );

old_name = old_names{get(handles.popup_signame_template, 'Value')};
set(handles.popup_signame_concrete,'Value',handles.sigmap(old_name));

% Update handles structure
guidata(hObject, handles);



% Hints: contents = cellstr(get(hObject,'String')) returns popup_signame_template contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_signame_template


% --- Executes during object creation, after setting all properties.
function popup_signame_template_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_signame_template (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in popup_signame_concrete.
function popup_signame_concrete_Callback(hObject, eventdata, handles)
% hObject    handle to popup_signame_concrete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_signame_concrete contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_signame_concrete

% contents = cellstr(get(hObject,'String')) % returns popup_signame_concrete contents as cell array
%  contents{get(hObject,'Value')} returns selected item from popup_signame_concrete

handles = update_phi(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function popup_signame_concrete_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_signame_concrete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function display_STL_Formula_CreateFcn(hObject, eventdata, handles)
% hObject    handle to display_STL_Formula (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in button_cancel.
function button_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to button_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = [];
guidata(hObject, handles);
uiresume(handles.main);


% --- Executes on button press in button_ok.
function button_ok_Callback(hObject, eventdata, handles)
% hObject    handle to button_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = set_id(handles.STL_Formula, handles.phi_id);
guidata(hObject, handles);
uiresume(handles.main);

% --- Executes during object deletion, before destroying properties.
function main_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on selection change in listbox_stlib.
function listbox_stlib_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_stlib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_stlib contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_stlib


set(handles.popup_signame_template,'Value',1);
set(handles.popup_signame_concrete,'Value',1);

handles = update_tmp(handles);
handles = update_phi(handles);
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function listbox_stlib_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_stlib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_id_phi_Callback(hObject, eventdata, handles)
% hObject    handle to edit_id_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_id_phi as text
%        str2double(get(hObject,'String')) returns contents of edit_id_phi as a double
st_id = get(hObject,'String');
handles.phi_id = st_id;
handles.auto_naming=0;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_id_phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_id_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function display_STL_Template_Callback(hObject, eventdata, handles)
% hObject    handle to display_STL_Template (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of display_STL_Template as text
%        str2double(get(hObject,'String')) returns contents of display_STL_Template as a double

st_phi_tmp = get(hObject,'String');
handles.STL_Template = STL_Formula('phi_tmp', st_phi_tmp);
guidata(hObject, handles);


function display_STL_Formula_Callback(hObject, eventdata, handles)
% hObject    handle to display_STL_Formula (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of display_STL_Formula as text
%        str2double(get(hObject,'String')) returns contents of display_STL_Formula as a double

st_phi = get(hObject,'String');
handles.STL_Formula = STL_Formula(handles.phi_id, st_phi); 
guidata(hObject, handles);

function handles = update_phi(handles)

old_names = get(handles.popup_signame_template,'String' );
new_names = get(handles.popup_signame_concrete,'String' );

old_name = old_names{get(handles.popup_signame_template, 'Value')}';

handles.sigmap(old_name) = get(handles.popup_signame_concrete, 'Value');

phi = handles.STL_Template;
% update formula

nb_old = numel(old_names);
for i_old = 1:nb_old;
    old_name = old_names{i_old};
    i_new = handles.sigmap(old_name);
    if i_new>1
      new_name = new_names{i_new};
      phi = STL_RenameSignals(phi, old_name,new_name); 
    end
end
if handles.auto_naming == 1
    handles.phi_id = get_id(phi);
    set(handles.edit_id_phi, 'String', get_id(phi));
    handles.STL_Formula = phi;
else
    handles.STL_Formula=set_id(phi, handles.phi_id);
end
set(handles.display_STL_Formula,'String', disp(phi,0));
