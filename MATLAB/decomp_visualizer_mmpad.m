function varargout = decomp_visualizer_mmpad(varargin)
% DECOMP_VISUALIZER_MMPAD MATLAB code for decomp_visualizer_mmpad.fig
%      DECOMP_VISUALIZER, by itself, creates a new DECOMP_VISUALIZER or raises the existing
%      singleton*.
%
%      H = DECOMP_VISUALIZER returns the handle to a new DECOMP_VISUALIZER or the handle to
%      the existing singleton*.
%
%      DECOMP_VISUALIZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DECOMP_VISUALIZER.M with the given input arguments.
%
%      DECOMP_VISUALIZER('Property','Value',...) creates a new DECOMP_VISUALIZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before decomp_visualizer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to decomp_visualizer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help decomp_visualizer

% Last Modified by GUIDE v2.5 21-Sep-2017 11:12:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @decomp_visualizer_OpeningFcn, ...
                   'gui_OutputFcn',  @decomp_visualizer_OutputFcn, ...
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


% --- Executes just before decomp_visualizer is made visible.
function decomp_visualizer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to decomp_visualizer (see VARARGIN)

% Choose default command line output for decomp_visualizer
handles.output = hObject;
handles.defaultDir = fullfile('D:','MMPAD_data');
handles.loaded = 0;

eval_menu_state(handles,hObject)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes decomp_visualizer wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function eval_menu_state(handles,hObject)
if handles.fix_var_theta_button.Value
    handles.var_theta_menu.Visible = 'on';
else
    handles.var_theta_menu.Visible = 'off'; 
end
if handles.fix_var_rad_button.Value
    handles.var_rad_menu.Visible = 'on';
else
    handles.var_rad_menu.Visible = 'off';
end
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = decomp_visualizer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load_fit_button.
function load_fit_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_fit_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fileName, pathName] = uigetfile('*.mat','Select a fit data file',...
                                 handles.defaultDir);
if isstr(fileName)
    % Add mat file names to menu
    allFiles = dir([pathName,'*.mat']);
    handles.filePath = pathName;
    handles.menu_file.String = natsortfiles({allFiles.name});
    handles.menu_file.Value = find(strcmp(handles.menu_file.String,fileName));
	
    % Load file
    load(fullfile(pathName,fileName))
    handles.loaded = 1;
    
    % Save data
    handles.P = P;
    
    % Compute basis matrices
    handles.A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(handles.P);
    handles.A0_stack = unshifted_basis_matrix_stack_norm2(handles.P);
    
    handles.err = err;
    handles.x_hat = x_hat;  
    handles.fit_image = Ax_ft_2D(handles.A0ft_stack,handles.x_hat)*norm(polar_image(:));
    handles.fit_image = handles.fit_image;
    handles.polar_image = polar_image;
    
    % Update static text
    handles.text_coefs.String = sprintf('Total Coefficients: %2.2e',sum(handles.x_hat(:)>0));
    handles.text_error.String = sprintf('Relative Error: %3.3f',handles.err(end));
    try
        handles.text_alpha.String = sprintf('Alpha: %3.3e',handles.P.alphap);
        handles.text_beta.String = sprintf('Beta: %3.3e',handles.P.betap);
    catch
        handles.text_alpha.String = 'Alpha:';
        handles.text_beta.String = 'Beta:';
    end
    handles.text_img.String = sprintf('Image: %i',handles.P.img);
    try
        handles.text_load_step.String = sprintf('Load Step: %i',handles.P.load_step);
    catch
        handles.text_load_step.String = sprintf('Set: %i',handles.P.set);
        handles.P.load_step = handles.P.set;
    end
    % Update menus
    handles.var_rad_menu.String = update_var_menu(handles.P.var_rad,handles.P.drad);
    handles.var_theta_menu.String = update_var_menu(handles.P.var_theta,handles.P.dtheta);
    
    % Set default region
    [n,m] = size(handles.polar_image);
    if n > m
        handles.polar_image = handles.polar_image';
        handles.fit_image = handles.fit_image';
        [n,m] = size(handles.polar_image);
    end
    x = floor(m/2);
    y = 1;
    Nx = min([200,x]);
    Ny = n;
    handles.rows = y:y+Ny-1;
    handles.cols = x-floor(Nx/2):x+floor(Nx/2)-1;
    handles.region = [handles.cols(1), handles.rows(1), Nx-1, Ny-1]; 
    handles.box_data = 0;
    handles.box_fit = 0;
    
    % Plot images
    plot_full_images(handles)
    handles = init_boxes(handles);
    
    % Update handles
    guidata(hObject, handles);
end

function menu = update_var_menu(variances,dx)
menu = {};
for i = 1:numel(variances)
    menu{i} = sprintf('%4d: %4.4f',i,variances(i)/dx);
end


function plot_full_images(handles) 
if ~handles.loaded
    return
end
% Plot full data iamges
axes(handles.axis_data)
lim2 = max(handles.polar_image(:));
handles.im_data = imshow(handles.polar_image,'DisplayRange',[0 lim2],'Colormap',jet);
set(handles.im_data,'ButtonDownFcn',@box_ButtonDownFcn);
title(sprintf('Original. Load: %i, Image: %i',handles.P.load_step,handles.P.img))

axes(handles.axis_fit)
handles.im_fit = imshow(handles.fit_image,'DisplayRange',[0 lim2],'Colormap',jet);
set(handles.im_fit,'ButtonDownFcn',@box_ButtonDownFcn);
title('Fit Image')
 
axes(handles.axes_variance_dist)
handles.var_dist = compute_var_dist(handles);
lim1 = min(handles.var_dist(:));
lim2 = max(handles.var_dist(:));
handles.im_var_dist = imshow(handles.var_dist,'DisplayRange',[lim1 lim2],'Colormap',jet);
%xlabel('\sigma_\theta','Fontsize',10)
ylabel('\sigma_r','Fontsize',10)
title('Signal Distribution','Fontsize',10)
colorbar()

axes(handles.axes_coef_dist)
handles.coef_dist = compute_coef_dist(handles);
lim1 = min(handles.coef_dist(:));
lim2 = max(handles.coef_dist(:));
handles.im_coef_dist = imshow(handles.coef_dist,'DisplayRange',[lim1 lim2],'Colormap',jet);
%xlabel('\sigma_\theta','Fontsize',10)
ylabel('\sigma_r','Fontsize',10)
title(' Coefficient Distribution','Fontsize',10)
colorbar()

function var_dist = compute_var_dist(handles)
n = size(handles.x_hat,4);
m = size(handles.x_hat,3);
var_dist = zeros(n,m);
for i = 1:n
    for j = 1:m
        signal_slice = Ax_ft_2D(handles.A0ft_stack(:,:,j,i),handles.x_hat(:,:,j,i));
        var_dist(i,j) = sum(signal_slice(:));
    end
end

function coef_dist = compute_coef_dist(handles)
n = size(handles.x_hat,4);
m = size(handles.x_hat,3);
coef_dist = zeros(n,m);
for i = 1:n
    for j = 1:m
        coef_slice = handles.x_hat(:,:,j,i);
        coef_dist(i,j) = sum(coef_slice(:));
    end
end

function handles = init_boxes(handles)
axes(handles.axis_data)
handles.box_data = rectangle('Position',handles.region,'EdgeColor','r','LineWidth',1);

axes(handles.axis_fit)
handles.box_fit = rectangle('Position',handles.region,'EdgeColor','r','LineWidth',1);

% Plot regions corresponding to boxes
plot_regions(handles)

function handles = update_boxes(handles)
axes(handles.axis_data)
delete(handles.box_data)
handles.box_data = rectangle('Position',handles.region,'EdgeColor','r','LineWidth',1);

axes(handles.axis_fit)
delete(handles.box_fit)
handles.box_fit = rectangle('Position',handles.region,'EdgeColor','r','LineWidth',1);

% Plot regions corresponding to boxes
plot_regions(handles)

function clear_region_axes(handles)
for i = 1:15
        eval(sprintf('cla(handles.axes%i,''reset'')',i))
end

function plot_regions(handles)
% Plot depending on the display mode
if handles.button_signal.Value
    plot_signals(handles)
elseif handles.button_coefficients.Value
    plot_coefs(handles)
elseif handles.button_basis.Value
    plot_basis(handles)
end

function plot_signals(handles)
if ~handles.loaded 
    return
end
lim2 = max(handles.polar_image(:));
scale = norm(handles.polar_image(:));
axes(handles.axis_region_data)
imshow(handles.polar_image(handles.rows,handles.cols),'DisplayRange',[0 lim2],'Colormap',jet)
title(sprintf('Original. Load: %i, Image: %i',handles.P.load_step,handles.P.img))
axes(handles.axis_region_fit) 
imshow(handles.fit_image(handles.rows,handles.cols),'DisplayRange',[0 lim2],'Colormap',jet)
title('Fit')

clear_region_axes(handles)
% Fix radial variance
if handles.fix_var_rad_button.Value
    variances = handles.P.var_theta;
    var_idx = handles.var_rad_menu.Value;
    % Show contributions from different coef values (fixed radial variance)
    for i = 1:numel(variances)
        eval(sprintf('axes(handles.axes%i)',i))
        signal_slice = Ax_ft_2D(handles.A0ft_stack(:,:,i,var_idx),handles.x_hat(:,:,i,var_idx))';
        imshow(signal_slice(handles.rows,handles.cols)*scale,'DisplayRange',[0 lim2],'Colormap',jet)
        signal_sum = sum(signal_slice(:));
        title(sprintf('%4d: %6.4f',i,signal_sum));
    end
% Fix theta variance
elseif handles.fix_var_theta_button.Value
    variances = handles.P.var_rad;
    var_idx = handles.var_theta_menu.Value;
    % Show contributions from different coef values (fixed radial variance)
    for i = 1:numel(variances)
        eval(sprintf('axes(handles.axes%i)',i))
        signal_slice = Ax_ft_2D(handles.A0ft_stack(:,:,var_idx,i),handles.x_hat(:,:,var_idx,i))';
        imshow(signal_slice(handles.rows,handles.cols)*scale,'DisplayRange',[0 lim2],'Colormap',jet)
        signal_sum = sum(signal_slice(:));
        title(sprintf('%4d: %6.4f',i,signal_sum));
    end
end  

function plot_coefs(handles)
if ~handles.loaded 
    return
end
lim2 = max(handles.polar_image(:));
axes(handles.axis_region_data)
imshow(handles.polar_image(handles.rows,handles.cols),'DisplayRange',[0 lim2],'Colormap',jet)
title(sprintf('Original. Load: %i, Image: %i',handles.P.load_step,handles.P.img))
axes(handles.axis_region_fit) 
imshow(handles.fit_image(handles.rows,handles.cols),'DisplayRange',[0 lim2],'Colormap',jet)
title('Fit')

clear_region_axes(handles)
maxLim = (max(real(handles.x_hat(:))));
% Fix radial variance
if handles.fix_var_rad_button.Value
    variances = handles.P.var_theta;
    var_idx = handles.var_rad_menu.Value;   
    % Show contributions from different coef values (fixed radial variance)
    for i = 1:numel(variances)
        eval(sprintf('axes(handles.axes%i)',i))
        coef_slice = handles.x_hat(:,:,i,var_idx)';
        imshow(coef_slice(handles.rows,handles.cols),'DisplayRange',[0 maxLim],'Colormap',jet)
        coef_sum = sum(coef_slice(:));
        title(sprintf('%4d: %6.4f',i,coef_sum));
    end
% Fix theta variance
elseif handles.fix_var_theta_button.Value
    variances = handles.P.var_rad;
    var_idx = handles.var_theta_menu.Value;
    % Show contributions from different coef values (fixed radial variance)
    for i = 1:numel(variances)
        eval(sprintf('axes(handles.axes%i)',i))
        coef_slice = handles.x_hat(:,:,var_idx,i)';
        imshow(coef_slice(handles.rows,handles.cols),'DisplayRange',[0 maxLim],'Colormap',jet)
        coef_sum = sum(coef_slice(:));
        title(sprintf('%4d: %6.4f',i,coef_sum));
    end
end  


function plot_basis(handles)
if ~handles.loaded 
    return
end
P = handles.P;
lim2 = max(handles.polar_image(:));
axes(handles.axis_region_data)
imshow(handles.polar_image(handles.rows,handles.cols),'DisplayRange',[0 lim2],'Colormap',jet)
title(sprintf('Original. Load: %i, Image: %i',handles.P.load_step,handles.P.img))
axes(handles.axis_region_fit) 
imshow(handles.fit_image(handles.rows,handles.cols),'DisplayRange',[0 lim2],'Colormap',jet)
title('Fit')

clear_region_axes(handles)
% Fix radial variance
if handles.fix_var_rad_button.Value
    variances = handles.P.var_theta;
    var_idx = handles.var_rad_menu.Value;
    % Show contributions from different coef values (fixed radial variance)
    for i = 1:numel(variances)
        eval(sprintf('axes(handles.axes%i)',i))
        basis = handles.A0_stack(:,:,i,var_idx)';
        basis = basis./max(basis(:));
        basis_shift = shift2D(basis,numel(handles.rows)/2,numel(handles.cols)/2);
        imshow(basis_shift(1:numel(handles.rows),1:numel(handles.cols)),'DisplayRange',[0 1],'Colormap',jet)
        title( ['\sigma_{\theta}/d\theta = ' num2str(sqrt(P.var_theta(i))/P.dtheta),...
            ',    \sigma_r/dr = ', num2str(sqrt(P.var_rad(var_idx))/P.drad)] );
    end
% Fix theta variance
elseif handles.fix_var_theta_button.Value
    variances = handles.P.var_rad;
    var_idx = handles.var_theta_menu.Value;
    % Show contributions from different coef values (fixed radial variance)
    for i = 1:numel(variances)
        eval(sprintf('axes(handles.axes%i)',i))
        basis = handles.A0_stack(:,:,var_idx,i)';
        basis = basis./max(basis(:));
        basis_shift = shift2D(basis,numel(handles.rows)/2,numel(handles.cols)/2);
        imshow(basis_shift(1:numel(handles.rows),1:numel(handles.cols)),'DisplayRange',[0 1],'Colormap',jet)
        title( ['\sigma_{\theta}/d\theta = ' num2str(sqrt(P.var_theta(var_idx))/P.dtheta),...
            ',    \sigma_r/dr = ', num2str(sqrt(P.var_rad(i))/P.drad)] );
    end
end  

% --- Executes on selection change in var_theta_menu.
function var_theta_menu_Callback(hObject, eventdata, handles)
% hObject    handle to var_theta_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns var_theta_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from var_theta_menu
plot_regions(handles)

% --- Executes during object creation, after setting all properties.
function var_theta_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to var_theta_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in var_rad_menu.
function var_rad_menu_Callback(hObject, eventdata, handles)
% hObject    handle to var_rad_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns var_rad_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from var_rad_menu
plot_regions(handles)

% --- Executes during object creation, after setting all properties.
function var_rad_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to var_rad_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fix_var_theta_button.
function fix_var_theta_button_Callback(hObject, eventdata, handles)
% hObject    handle to fix_var_theta_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fix_var_theta_button
eval_menu_state(handles,hObject)

% Update plots
plot_regions(handles)

% --- Executes on button press in fix_var_rad_button.
function fix_var_rad_button_Callback(hObject, eventdata, handles)
% hObject    handle to fix_var_rad_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fix_var_rad_button
eval_menu_state(handles,hObject)

% Update plots
plot_regions(handles)



% --- Executes on mouse press over axes background.
function box_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axis_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get data
axesHandle  = get(hObject,'Parent');
coordinates = get(axesHandle,'CurrentPoint'); 
coordinates = coordinates(1,1:2);


% Save new box coordinates
uiObject = get(axesHandle,'Parent');
handles = guidata(uiObject);
[n,m] = size(handles.polar_image);
handles.region(1) = max(1,min(m-handles.region(3),round(coordinates(1))-floor(handles.region(3)/2))); %x
handles.rows = handles.region(2):handles.region(2)+handles.region(4);
handles.cols = handles.region(1):handles.region(1)+handles.region(3);
    
handles = update_boxes(handles);
plot_regions(handles)

% Update handles
guidata(uiObject, handles);


% --- Executes on button press in button_coefficients.
function button_coefficients_Callback(hObject, eventdata, handles)
% hObject    handle to button_coefficients (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of button_coefficients
plot_regions(handles)


% --- Executes on button press in button_signal.
function button_signal_Callback(hObject, eventdata, handles)
% hObject    handle to button_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of button_signal
plot_regions(handles)


% --- Executes on button press in button_basis.
function button_basis_Callback(hObject, eventdata, handles)
% hObject    handle to button_basis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of button_basis
plot_regions(handles)


% --- Executes on selection change in menu_file.
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menu_file contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menu_file
     % Add mat file names to menu

% Load file
fileName = handles.menu_file.String{handles.menu_file.Value};
load(fullfile(handles.filePath,fileName))
handles.loaded = 1;

% Save data
handles.P = P;

% Compute basis matrices
handles.A0ft_stack = unshifted_basis_matrix_ft_stack_norm2(handles.P);
handles.A0_stack = unshifted_basis_matrix_stack_norm2(handles.P);

handles.err = err;
handles.x_hat = x_hat;  
handles.fit_image = Ax_ft_2D(handles.A0ft_stack,handles.x_hat)*norm(polar_image(:));
handles.fit_image = handles.fit_image;
handles.polar_image = polar_image;

% Update static text
handles.text_coefs.String = sprintf('Total Coefficients: %2.2e',sum(handles.x_hat(:)>0));
handles.text_error.String = sprintf('Relative Error: %3.3f',handles.err(end));
try
    handles.text_alpha.String = sprintf('Alpha: %3.3e',handles.P.alphap);
    handles.text_beta.String = sprintf('Beta: %3.3e',handles.P.betap);
catch
    handles.text_alpha.String = 'Alpha:';
    handles.text_beta.String = 'Beta:';
end
handles.text_img.String = sprintf('Image: %i',handles.P.img);
try
    handles.text_load_step.String = sprintf('Load Step: %i',handles.P.load_step);
catch
    handles.text_load_step.String = sprintf('Set: %i',handles.P.set);
    handles.P.load_step = handles.P.set;
end
% Update menus
handles.var_rad_menu.String = update_var_menu(handles.P.var_rad,handles.P.drad);
handles.var_theta_menu.String = update_var_menu(handles.P.var_theta,handles.P.dtheta);

% Set default region
[n,m] = size(handles.polar_image);
if n > m
    handles.polar_image = handles.polar_image';
    handles.fit_image = handles.fit_image';
    [n,m] = size(handles.polar_image);
end
x = floor(m/2);
y = 1;
Nx = min([200,x]);
Ny = n;
handles.rows = y:y+Ny-1;
handles.cols = x-floor(Nx/2):x+floor(Nx/2)-1;
handles.region = [handles.cols(1), handles.rows(1), Nx-1, Ny-1]; 
handles.box_data = 0;
handles.box_fit = 0;

% Plot images
plot_full_images(handles)
handles = init_boxes(handles);

% Update handles
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function menu_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
