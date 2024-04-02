function varargout = voxToolSlices(varargin)
% VOXTOOLSLICES MATLAB code for voxToolSlices.fig
%      VOXTOOLSLICES, by itself, creates a new VOXTOOLSLICES or raises the existing
%      singleton*.
%
%      H = VOXTOOLSLICES returns the handle to a new VOXTOOLSLICES or the handle to
%      the existing singleton*.
%
%      VOXTOOLSLICES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VOXTOOLSLICES.M with the given input arguments.
%
%      VOXTOOLSLICES('Property','Value',...) creates a new VOXTOOLSLICES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before voxToolSlices_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to voxToolSlices_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help voxToolSlices

% Last Modified by GUIDE v2.5 27-Jun-2013 15:49:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @voxToolSlices_OpeningFcn, ...
                   'gui_OutputFcn',  @voxToolSlices_OutputFcn, ...
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


% --- Executes just before voxToolSlices is made visible.
function voxToolSlices_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to voxToolSlices (see VARARGIN)
set(0,'currentfigure',hObject);
handles.sliceXCircleHandles = annotation('ellipse',[0,0,0,0],'color','red');
handles.sliceYCircleHandles = annotation('ellipse',[0,0,0,0],'color','red');
handles.sliceZCircleHandles = annotation('ellipse',[0,0,0,0],'color','red');


handles.sliceFig = hObject;
% Choose default command line output for voxToolSlices
handles.output = handles;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes voxToolSlices wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = voxToolSlices_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
