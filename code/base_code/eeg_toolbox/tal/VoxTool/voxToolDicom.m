function varargout = voxToolDicom(varargin)
% VOXTOOLDICOM MATLAB code for voxToolDicom.fig - Holds the window that
% shows the dicoms
%
%      VOXTOOLDICOM, by itself, creates a new VOXTOOLDICOM or raises the existing
%      singleton*.
%
%      H = VOXTOOLDICOM returns the handle to a new VOXTOOLDICOM or the handle to
%      the existing singleton*.
%
%      VOXTOOLDICOM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VOXTOOLDICOM.M with the given input arguments.
%
%      VOXTOOLDICOM('Property','Value',...) creates a new VOXTOOLDICOM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before voxToolDicom_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to voxToolDicom_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help voxToolDicom

% Last Modified by GUIDE v2.5 10-Jun-2013 17:29:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @voxToolDicom_OpeningFcn, ...
                   'gui_OutputFcn',  @voxToolDicom_OutputFcn, ...
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


% --- Executes just before voxToolDicom is made visible.
function voxToolDicom_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to voxToolDicom (see VARARGIN)

% Choose default command line output for voxToolDicom

% makes the handle to the dicom circle with no size, so that it can be
% changed rather than created on the next click.
handles.dicomCircleHandles = annotation('ellipse',[0,0,0,0],'color','red');

handles.dicomFig = hObject;
handles.output = handles;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes voxToolDicom wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = voxToolDicom_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
