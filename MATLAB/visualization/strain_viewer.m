function varargout = strain_viewer(varargin)
% STRAIN_VIEWER MATLAB code for strain_viewer.fig
%      STRAIN_VIEWER, by itself, creates a new STRAIN_VIEWER or raises the existing
%      singleton*.
%
%      H = STRAIN_VIEWER returns the handle to a new STRAIN_VIEWER or the handle to
%      the existing singleton*.
%
%      STRAIN_VIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STRAIN_VIEWER.M with the given input arguments.
%
%      STRAIN_VIEWER('Property','Value',...) creates a new STRAIN_VIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before strain_viewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to strain_viewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help strain_viewer

% Last Modified by GUIDE v2.5 25-Oct-2017 08:42:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @strain_viewer_OpeningFcn, ...
                   'gui_OutputFcn',  @strain_viewer_OutputFcn, ...
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


% --- Executes just before strain_viewer is made visible.
function strain_viewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to strain_viewer (see VARARGIN)

% Choose default command line output for strain_viewer
handles.output = hObject;
handles.defaultDir = fullfile('C:','Users','dbanco02','CHESS-Research','MATLAB');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes strain_viewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = strain_viewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_load_data.
function button_load_data_Callback(hObject, eventdata, handles)
% hObject    handle to button_load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[fileName, pathName] = uigetfile('*.mat','Select a fit data file',...
                                 handles.defaultDir);
if isstr(fileName)
	% Load file
    load(fullfile(pathName,fileName))
    handles.loaded = 1;
    
    % Save data
    handles.strain = strain;
    
    handles.edit_rowstart.String = '1';
    handles.edit_rowend.String = num2str(size(handles.strain,4));
    handles.r1 = str2num(handles.edit_rowstart.String);
    handles.rend = str2num(handles.edit_rowend.String);
    
    handles.apply_limits = 0;
    
    % Plot images
    handles = update_indices(handles);
    plot_strain(handles)
    
    % Update handles
    guidata(hObject, handles);
end

function clear_axes(handles)
for i = 1:5
        eval(sprintf('cla(handles.axes%i,''reset'')',i))
end

function handles = update_indices(handles)
if handles.radio_exx.Value
    handles.comp_idx = 1;
elseif handles.radio_eyy.Value
    handles.comp_idx = 2;
elseif handles.radio_exy.Value
    handles.comp_idx = 3;
end

if handles.radio_total.Value
    handles.type_idx = 1;
elseif handles.radio_elastic.Value
    handles.type_idx = 2;
elseif handles.radio_plastic.Value
    handles.type_idx = 3;  
end

function plot_strain(handles)
if handles.radio_subtract.Value
    center = squeeze(handles.strain(handles.type_idx,handles.comp_idx,1,handles.r1:handles.rend,:));
else
    center = 0;
end
if handles.apply_limits
    handles.max_limit = str2num(handles.edit_max.String);
    handles.min_limit = str2num(handles.edit_min.String);
elseif handles.radio_subtract.Value
    handles.max_limit = 0.025;
    handles.min_limit = -0.025;
else
    handles.max_limit = 0.025;
    handles.min_limit = -0.025;
end
for i = 1:5
    eval(sprintf('axes(handles.axes%i)',i)) 
    tmp_img = squeeze(handles.strain(handles.type_idx,handles.comp_idx,i,handles.r1:handles.rend,:));
    imshow( (tmp_img-center),'DisplayRange',...
        [handles.min_limit handles.max_limit],'Colormap',jet);
    title(sprintf('Load %i',i))
end
colorbar() 




% --- Executes on button press in radio_total.
function radio_total_Callback(hObject, eventdata, handles)
% hObject    handle to radio_total (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_total

% Plot images
handles = update_indices(handles);
clear_axes(handles)
plot_strain(handles)

% Update handles
guidata(hObject, handles);


% --- Executes on button press in radio_elastic.
function radio_elastic_Callback(hObject, eventdata, handles)
% hObject    handle to radio_elastic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_elastic
% Update cutoff menu

% Plot images
handles = update_indices(handles);
clear_axes(handles)
plot_strain(handles)

% Update handles
guidata(hObject, handles);


% --- Executes on button press in radio_subtract.
function radio_subtract_Callback(hObject, eventdata, handles)
% hObject    handle to radio_subtract (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_subtract

% Plot images
handles = update_indices(handles);
clear_axes(handles)
plot_strain(handles)

% Update handles
guidata(hObject, handles);


% --- Executes on button press in radio_exx.
function radio_exx_Callback(hObject, eventdata, handles)
% hObject    handle to radio_exx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Plot images
handles = update_indices(handles);
clear_axes(handles)
plot_strain(handles)

% Update handles
guidata(hObject, handles);

% --- Executes on button press in radio_eyy.
function radio_eyy_Callback(hObject, eventdata, handles)
% hObject    handle to radio_eyy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Plot images
handles = update_indices(handles);
clear_axes(handles)
plot_strain(handles)

% Update handles
guidata(hObject, handles);

% --- Executes on button press in radio_exy.
function radio_exy_Callback(hObject, eventdata, handles)
% hObject    handle to radio_exy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Plot images
handles = update_indices(handles);
clear_axes(handles)
plot_strain(handles)

% Update handles
guidata(hObject, handles);



function edit_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_max as text
%        str2double(get(hObject,'String')) returns contents of edit_max as a double


% --- Executes during object creation, after setting all properties.
function edit_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_min as text
%        str2double(get(hObject,'String')) returns contents of edit_min as a double


% --- Executes during object creation, after setting all properties.
function edit_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_apply_limits.
function button_apply_limits_Callback(hObject, eventdata, handles)
% hObject    handle to button_apply_limits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.apply_limits = 1;

% Plot images
handles = update_indices(handles);
clear_axes(handles)
plot_strain(handles)

% Update handles
guidata(hObject, handles);

% --- Executes on button press in button_reset_limits.
function button_reset_limits_Callback(hObject, eventdata, handles)
% hObject    handle to button_reset_limits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.apply_limits = 0;

% Plot images
handles = update_indices(handles);
clear_axes(handles)
plot_strain(handles)

% Update handles
guidata(hObject, handles);


% --- Executes on button press in text_rowstart.
function radio_rowstart_Callback(hObject, eventdata, handles)
% hObject    handle to text_rowstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of text_rowstart


% --- Executes on button press in text_rowend.
function radio_rowend_Callback(hObject, eventdata, handles)
% hObject    handle to text_rowend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of text_rowend



function edit_rowstart_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rowstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rowstart as text
%        str2double(get(hObject,'String')) returns contents of edit_rowstart as a double
handles.r1 = str2num(handles.edit_rowstart.String);
% Update handles
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_rowstart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rowstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
    
function edit_rowend_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rowend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rowend as text
%        str2double(get(hObject,'String')) returns contents of edit_rowend as a double
handles.rend = str2num(handles.edit_rowend.String);
% Update handles
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_rowend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rowend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radio_plastic.
function radio_plastic_Callback(hObject, eventdata, handles)
% hObject    handle to radio_plastic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_plastic

% Plot images
handles = update_indices(handles);
clear_axes(handles)
plot_strain(handles)

% Update handles
guidata(hObject, handles);
