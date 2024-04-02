function varargout = thresholdDlg(varargin)
% THRESHOLDDLG MATLAB code for thresholdDlg.fig
%      THRESHOLDDLG, by itself, creates a new THRESHOLDDLG or raises the existing
%      singleton*.
%
%      H = THRESHOLDDLG returns the handle to a new THRESHOLDDLG or the handle to
%      the existing singleton*.
%
%      THRESHOLDDLG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in THRESHOLDDLG.M with the given input arguments.
%
%      THRESHOLDDLG('Property','Value',...) creates a new THRESHOLDDLG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before thresholdDlg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to thresholdDlg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help thresholdDlg

% Last Modified by GUIDE v2.5 02-Jun-2014 15:23:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @thresholdDlg_OpeningFcn, ...
                   'gui_OutputFcn',  @thresholdDlg_OutputFcn, ...
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


% --- Executes just before thresholdDlg is made visible.
function thresholdDlg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to thresholdDlg (see VARARGIN)

% Choose default command line output for thresholdDlg
handles.output = 1;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes thresholdDlg wait for user response (see UIRESUME)
uiwait(handles.threshFig);


% --- Outputs from this function are returned to the command line.
function varargout = thresholdDlg_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
handles.output = str2double(get(handles.threshEdit,'string'));
varargout{1} = handles.output;
delete(hObject)



function threshEdit_Callback(hObject, eventdata, handles)
% hObject    handle to threshEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshEdit as text
%        str2double(get(hObject,'String')) returns contents of threshEdit as a double
handles.output = str2double(get(hObject,'string'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function threshEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.threshFig);
