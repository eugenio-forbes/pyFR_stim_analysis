function varargout = voxTool3DCT(varargin)
% VOXTOOL3DCT MATLAB code for voxTool3DCT.fig - Shows the 3d CT plot and
% controls for viewing
%      VOXTOOL3DCT, by itself, creates a new VOXTOOL3DCT or raises the existing
%      singleton*.
%
%      H = VOXTOOL3DCT returns the handle to a new VOXTOOL3DCT or the handle to
%      the existing singleton*.
%
%      VOXTOOL3DCT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VOXTOOL3DCT.M with the given input arguments.
%
%      VOXTOOL3DCT('Property','Value',...) creates a new VOXTOOL3DCT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before voxTool3DCT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to voxTool3DCT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help voxTool3DCT

% Last Modified by GUIDE v2.5 10-Jun-2013 17:28:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @voxTool3DCT_OpeningFcn, ...
                   'gui_OutputFcn',  @voxTool3DCT_OutputFcn, ...
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


% --- Executes just before voxTool3DCT is made visible.
function voxTool3DCT_OpeningFcn(hObject, eventdata, handles, main_handles)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% main_handles   passed to this class from the control figure

% The handles consist of the handles from the control figure and the
% default handles for this class.
handles = concatStructs(handles, main_handles);
handles.selectedPoints = struct();
handles.selectedPoints.coords = [];

% This figure is named CT3D_fig for reference from other figures
handles.CT3D_fig = hObject;
handles.output = handles;

% Save handles to this figure
guidata(hObject, handles);
% Save handles to the control figure
guidata(handles.controlFigure, handles);

% If the coordinates are already defined upon creation of the figure,
% create the pointcloud
if isfield(handles, 'xyz')
    clickA3DPoint([handles.xyz(:,1),  handles.xyz(:,2), handles.xyz(:,3)]', handles);
end

h = rotate3d;
set(h, 'ActionPostCallback',@postRotate);

%Set the keypress function for this figure
set(gcf,'keyPressFcn', {@keypress});

function postRotate(obj, evd)
msgbox(sprintf('The new view is [%d %d].',newView));

% --- Outputs from this function are returned to the command line.
function varargout = voxTool3DCT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles;


% --- Executes on slider movement.
function xSlider1_Callback(hObject, eventdata, handles)
% hObject    handle to xSlider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Delete and replot the pointcloud when the slider is adjusted
delete(findobj(handles.mainAxes,'tag','fullBrain'));
plotWithLimits(handles);
xlim([0,500]);
ylim([0,500]);
zlim([0,600]);

% --- Executes during object creation, after setting all properties.
function xSlider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xSlider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function xSlider2_Callback(hObject, eventdata, handles)
% hObject    handle to xSlider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Delete and replot the pointcloud when the slider is adjusted
delete(findobj(handles.mainAxes,'tag','fullBrain'));
plotWithLimits(handles);
xlim([0,500]);
ylim([0,500]);
zlim([0,600]);

% --- Executes during object creation, after setting all properties.
function xSlider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xSlider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function ySlider2_Callback(hObject, eventdata, handles)
% hObject    handle to ySlider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Delete and replot the pointcloud when the slider is adjusted
delete(findobj(handles.mainAxes,'tag','fullBrain'));
plotWithLimits(handles);
xlim([0,500]);
ylim([0,500]);
zlim([0,600]);

% --- Executes during object creation, after setting all properties.
function ySlider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ySlider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function ySlider1_Callback(hObject, eventdata, handles)
% hObject    handle to ySlider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Delete and replot the pointcloud when the slider is adjusted
delete(findobj(handles.mainAxes,'tag','fullBrain'));
plotWithLimits(handles);
xlim([0,500]);
ylim([0,500]);
zlim([0,600]);

% --- Executes during object creation, after setting all properties.
function ySlider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ySlider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function zSlider1_Callback(hObject, eventdata, handles)
% hObject    handle to zSlider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Delete and replot the pointcloud when the slider is adjusted
delete(findobj(handles.mainAxes,'tag','fullBrain'));
plotWithLimits(handles);
xlim([0,500]);
ylim([0,500]);
zlim([0,600]);

% --- Executes during object creation, after setting all properties.
function zSlider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zSlider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function zSlider2_Callback(hObject, eventdata, handles)
% hObject    handle to zSlider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Delete and replot the pointcloud when the slider is adjusted
delete(findobj(handles.mainAxes,'tag','fullBrain'));
plotWithLimits(handles);
xlim([0,500]);
ylim([0,500]);
zlim([0,600]);

% --- Executes during object creation, after setting all properties.
function zSlider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zSlider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
