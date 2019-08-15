function varargout = spread_visualizer(varargin)
% SPREAD_VISUALIZER MATLAB code for spread_visualizer.fig
%      SPREAD_VISUALIZER, by itself, creates a new SPREAD_VISUALIZER or raises the existing
%      singleton*.
%
%      H = SPREAD_VISUALIZER returns the handle to a new SPREAD_VISUALIZER or the handle to
%      the existing singleton*.
%
%      SPREAD_VISUALIZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPREAD_VISUALIZER.M with the given input arguments.
%
%      SPREAD_VISUALIZER('Property','Value',...) creates a new SPREAD_VISUALIZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spread_visualizer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to spread_visualizer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help spread_visualizer

% Last Modified by GUIDE v2.5 06-Mar-2018 18:29:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spread_visualizer_OpeningFcn, ...
                   'gui_OutputFcn',  @spread_visualizer_OutputFcn, ...
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


% --- Executes just before spread_visualizer is made visible.
function spread_visualizer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spread_visualizer (see VARARGIN)

% Choose default command line output for spread_visualizer
handles.output = hObject;
handles.defaultDir = fullfile('D:','CHESS_data','spread_results');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes spread_visualizer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = spread_visualizer_OutputFcn(hObject, eventdata, handles) 
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
    handles.var_signal = var_signal;
    handles.rel_error = rel_error;
    handles.sparsity = sparsity;
    handles.P = P;
    
    handles.edit_rowstart.String = '1';
    handles.edit_rowend.String = num2str(size(var_signal,4));
    handles.r1 = str2num(handles.edit_rowstart.String);
    handles.rend = str2num(handles.edit_rowend.String);
    
    % Update cutoff menu
    if handles.radio_az.Value
        handles.menu_cutoff.String = update_cutoff_menu(handles.P.var_theta);
    elseif handles.radio_rad.Value
        handles.menu_cutoff.String = update_cutoff_menu(handles.P.var_rad);
    end
    handles.apply_limits = 0;
    
    % Plot images
    plot_images(handles)
    
    % Update handles
    guidata(hObject, handles);
end

function menu = update_cutoff_menu(variances)
menu = {};
for i = 1:numel(variances)
    menu{i} = sprintf('%4d: %4.4f',i,variances(i));
end

function clear_axes(handles)
for i = 1:5
        eval(sprintf('cla(handles.axes%i,''reset'')',i))
end

function plot_images(handles)
if handles.radio_spread.Value
    if  handles.radio_az.Value
        plot_az_spread(handles)
    elseif handles.radio_rad.Value
        plot_rad_spread(handles)
    end
elseif handles.radio_meanvar.Value
    if  handles.radio_az.Value
        plot_az_mean_variance(handles)
    elseif handles.radio_rad.Value
        plot_rad_mean_variance(handles)
    end
elseif handles.radio_sparsity.Value
    plot_az_awmv(handles)
elseif handles.radio_error.Value
    plot_error(handles)
end

function plot_az_awmv(handles)

az_var = squeeze(sum(handles.var_signal(:,:,1:5,:,:),2));
total = squeeze(sum(az_var(:,1:5,:,:),1));
awmv_var = zeros([size(az_var,2),size(az_var,3),size(az_var,4)]);

for i = 1:size(az_var,2)
    for j = 1:size(az_var,3)
        for k = 1:size(az_var,4)
            awmv_var(i,j,k) = sum(sqrt(handles.P.var_theta')./handles.P.dtheta.*squeeze(az_var(:,i,j,k)))/total(i,j,k);
        end
    end
end

if handles.radio_subtract.Value
    center = awmv_var(1,handles.r1:handles.rend,:);
else
    center = 0;
end
if handles.apply_limits
    handles.max_limit = str2num(handles.edit_max.String);
    handles.min_limit = str2num(handles.edit_min.String);
elseif handles.radio_subtract.Value
    handles.max_limit = 1;
    handles.min_limit = -1;
else
    handles.max_limit = 1;
    handles.min_limit = 0;
end
for i = 1:5
    eval(sprintf('axes(handles.axes%i)',i)) 
    handles.viewData{i} = squeeze(awmv_var(i,handles.r1:handles.rend,:)-center);
    imshow(handles.viewData{i},'DisplayRange',...
        [handles.min_limit handles.max_limit],'Colormap',jet)
    title(sprintf('Load %i',i))
end
colorbar()  



function plot_az_spread(handles)
cutoff = handles.menu_cutoff.Value;
total_var = squeeze(sum(sum(handles.var_signal(:,:,1:5,:,:),1),2));
high_var_theta = squeeze(sum(sum(handles.var_signal(cutoff:end,:,1:5,:,:),1),2))./total_var;
if handles.radio_subtract.Value
    center = high_var_theta(1,handles.r1:handles.rend,:);
else
    center = 0;
end
if handles.apply_limits
    handles.max_limit = str2num(handles.edit_max.String);
    handles.min_limit = str2num(handles.edit_min.String);
elseif handles.radio_subtract.Value
    handles.max_limit = 1;
    handles.min_limit = -1;
else
    handles.max_limit = 1;
    handles.min_limit = 0;
end
for i = 1:5
    eval(sprintf('axes(handles.axes%i)',i)) 
    handles.viewData{i} = squeeze(high_var_theta(i,handles.r1:handles.rend,:)-center);
    imshow(handles.viewData{i},'DisplayRange',...
        [handles.min_limit handles.max_limit],'Colormap',jet)
    title(sprintf('Load %i',i))
end
colorbar()  

function plot_rad_spread(handles)
cutoff = handles.menu_cutoff.Value;
total_var = squeeze(sum(sum(handles.var_signal(:,:,1:5,:,:),1),2));
high_var_rad = squeeze(sum(sum(handles.var_signal(:,cutoff:end,1:5,:,:),1),2))./total_var; 
if handles.radio_subtract.Value
    center = high_var_rad(1,handles.r1:handles.rend,:);
else
    center = 0;
end
if handles.apply_limits
    handles.max_limit = str2num(handles.edit_max.String);
    handles.min_limit = str2num(handles.edit_min.String);
elseif handles.radio_subtract.Value
    handles.max_limit = 1;
    handles.min_limit = -1;
else
    handles.max_limit = 1;
    handles.min_limit = 0;
end
for i = 1:5
    eval(sprintf('axes(handles.axes%i)',i))
    handles.viewData{i} = squeeze(high_var_rad(i,handles.r1:handles.rend,:)-center);
    imshow(handles.viewData{i},'DisplayRange',...
           [handles.min_limit handles.max_limit],'Colormap',jet)
    title(sprintf('Load %i',i))
end
colorbar()

function mvar = compute_mean_az_variance(var_signal,var_theta)
az_var_signal = squeeze(sum(var_signal,2));
[~, s2, s3, s4] = size(az_var_signal);
mvar = zeros(s2,s3,s4);
for j = 1:s2
    for k = 1:s3
        for l = 1:s4
            total = sum(az_var_signal(:,j,k,l));
            for i = 1:numel(var_theta)
                mvar(j,k,l) = mvar(j,k,l) +...
                      var_theta(i)*az_var_signal(i,j,k,l)/total;
            end
        end        
    end
end

function mvar = compute_mean_rad_variance(var_signal,var_theta)
az_var_signal = squeeze(sum(var_signal,1));
[~, s2, s3, s4] = size(az_var_signal);
mvar = zeros(s2,s3,s4);
for j = 1:s2
    for k = 1:s3
        for l = 1:s4
            total = sum(az_var_signal(:,j,k,l));
            for i = 1:numel(var_theta)
                mvar(j,k,l) = mvar(j,k,l) +...
                      var_theta(i)*az_var_signal(i,j,k,l)/total;
            end
        end        
    end
end


function plot_az_mean_variance(handles)
mvar = compute_mean_az_variance(handles.var_signal,handles.P.var_theta);
if handles.radio_subtract.Value
    center = mvar(1,handles.r1:handles.rend,:);
else
    center = 0;
end
if handles.apply_limits
    handles.max_limit = str2num(handles.edit_max.String);
    handles.min_limit = str2num(handles.edit_min.String);
elseif handles.radio_subtract.Value
    handles.max_limit = 1;
    handles.min_limit = -1;
else
    handles.max_limit = handles.P.var_theta(end);
    handles.min_limit = handles.P.var_theta(1);
end
for i = 1:5
    eval(sprintf('axes(handles.axes%i)',i))
    handles.viewData{i} = squeeze(mvar(i,handles.r1:handles.rend,:)-center);
    imshow(handles.viewData{i},'DisplayRange',...
           [handles.min_limit handles.max_limit],'Colormap',jet)
    title(sprintf('Load %i',i))
end
colorbar()

function plot_rad_mean_variance(handles)
mvar = compute_mean_rad_variance(handles.var_signal,handles.P.var_rad);
if handles.radio_subtract.Value
    center = mvar(1,handles.r1:handles.rend,:);
else
    center = 0;
end
if handles.apply_limits
    handles.max_limit = str2num(handles.edit_max.String);
    handles.min_limit = str2num(handles.edit_min.String);
elseif handles.radio_subtract.Value
    handles.max_limit = 1;
    handles.min_limit = -1;
else
    handles.max_limit = handles.P.var_rad(end);
    handles.min_limit = handles.P.var_rad(1);
end
for i = 1:5
    eval(sprintf('axes(handles.axes%i)',i))
    handles.viewData{i} = squeeze(mvar(i,handles.r1:handles.rend,:)-center);
    imshow(handles.viewData{i},'DisplayRange',...
           [handles.min_limit handles.max_limit],'Colormap',jet)
    title(sprintf('Load %i',i))
end
colorbar()

function plot_sparsity(handles)
if handles.apply_limits
    handles.max_limit = str2num(handles.edit_max.String);
    handles.min_limit = str2num(handles.edit_min.String);
else
    handles.max_limit = max(handles.sparsity(:));
    handles.min_limit = min(handles.sparsity(:));
end
for i = 1:5
    eval(sprintf('axes(handles.axes%i)',i))
    imshow(squeeze(handles.sparsity(i,handles.r1:handles.rend,:)),'DisplayRange',...
            [handles.min_limit handles.max_limit],'Colormap',jet)
    title(sprintf('Load %i',i))
end
colorbar()

function plot_error(handles)
if handles.apply_limits
    handles.max_limit = str2num(handles.edit_max.String);
    handles.min_limit = str2num(handles.edit_min.String);
else
    handles.max_limit = 1;
    handles.min_limit = 0;
end

for i = 1:5
    eval(sprintf('axes(handles.axes%i)',i))
    imshow(squeeze(handles.rel_error(i,handles.r1:handles.rend,:)),'DisplayRange',...
            [handles.min_limit handles.max_limit],'Colormap',jet)
    title(sprintf('Load %i',i))
end
colorbar()

% --- Executes on selection change in menu_cutoff.
function menu_cutoff_Callback(hObject, eventdata, handles)
% hObject    handle to menu_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menu_cutoff contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menu_cutoff
clear_axes(handles)
plot_images(handles)

% --- Executes during object creation, after setting all properties.
function menu_cutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menu_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radio_az.
function radio_az_Callback(hObject, eventdata, handles)
% hObject    handle to radio_az (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_az
% Update cutoff menu
handles.menu_cutoff.String = update_cutoff_menu(handles.P.var_theta);
handles.menu_cutoff.Value = min(numel(handles.menu_cutoff.String),...
                                handles.menu_cutoff.Value);
% Plot images
clear_axes(handles)
plot_images(handles)

% Update handles
guidata(hObject, handles);


% --- Executes on button press in radio_rad.
function radio_rad_Callback(hObject, eventdata, handles)
% hObject    handle to radio_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_rad
% Update cutoff menu

handles.menu_cutoff.String = update_cutoff_menu(handles.P.var_rad);
handles.menu_cutoff.Value = min(numel(handles.menu_cutoff.String),...
                                handles.menu_cutoff.Value);
% Plot images
clear_axes(handles)
plot_images(handles)

% Update handles
guidata(hObject, handles);


% --- Executes on button press in radio_subtract.
function radio_subtract_Callback(hObject, eventdata, handles)
% hObject    handle to radio_subtract (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_subtract
% Plot images
clear_axes(handles)
plot_images(handles)

% Update handles
guidata(hObject, handles);


% --- Executes on button press in radio_spread.
function radio_spread_Callback(hObject, eventdata, handles)
% hObject    handle to radio_spread (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_spread
clear_axes(handles)
plot_images(handles)

% Update handles
guidata(hObject, handles);

% --- Executes on button press in radio_sparsity.
function radio_sparsity_Callback(hObject, eventdata, handles)
% hObject    handle to radio_sparsity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_sparsity
clear_axes(handles)
plot_images(handles)

% Update handles
guidata(hObject, handles);

% --- Executes on button press in radio_meanvar.
function radio_meanvar_Callback(hObject, eventdata, handles)
% hObject    handle to radio_meanvar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_meanvar
clear_axes(handles)
plot_images(handles)

% Update handles
guidata(hObject, handles);

% --- Executes on button press in radio_error.
function radio_error_Callback(hObject, eventdata, handles)
% hObject    handle to radio_error (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_error
clear_axes(handles)
plot_images(handles)

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
clear_axes(handles)
plot_images(handles)
% Update handles
guidata(hObject, handles);

% --- Executes on button press in button_reset_limits.
function button_reset_limits_Callback(hObject, eventdata, handles)
% hObject    handle to button_reset_limits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.apply_limits = 0;
clear_axes(handles)
plot_images(handles)

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

% % --- Executes on mouse press over figure background.
% function figure1_ButtonDownFcn(hObject, eventdata, handles)
% % hObject    handle to figure1 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% C = get;
% row = round(C(1,2));
% col = round(C(1,1));
% print(sprintf('Val: %2.2f, Row: %i, Col: %i',handles.viewData{5}(row,col),row,col))



