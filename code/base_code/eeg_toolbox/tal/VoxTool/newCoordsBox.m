function varargout = newCoordsBox(varargin)
% NEWCOORDSBOX MATLAB code for newCoordsBox.fig - A small box to enter new
% coordinates
%
%      NEWCOORDSBOX, by itself, creates a new NEWCOORDSBOX or raises the existing
%      singleton*.
%
%      H = NEWCOORDSBOX returns the handle to a new NEWCOORDSBOX or the handle to
%      the existing singleton*.
%
%      NEWCOORDSBOX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEWCOORDSBOX.M with the given input arguments.
%
%      NEWCOORDSBOX('Property','Value',...) creates a new NEWCOORDSBOX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before newCoordsBox_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to newCoordsBox_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help newCoordsBox

% Last Modified by GUIDE v2.5 11-Jun-2013 16:05:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @newCoordsBox_OpeningFcn, ...
                   'gui_OutputFcn',  @newCoordsBox_OutputFcn, ...
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


% --- Executes just before newCoordsBox is made visible.
function newCoordsBox_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to newCoordsBox (see VARARGIN)

% Choose default command line output for newCoordsBox
handles.output = hObject;
handles.isOk = false;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes newCoordsBox wait for user response
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = newCoordsBox_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Assigns the outputs to the X and Y coordinates entered, and the figure
% itself
if handles.isOk
    varargout{1} = [get(handles.x,'string'),'x',get(handles.y,'string')];
    varargout{2} = hObject;
else
    varargout{1} = [];
end
    


function x_Callback(hObject, eventdata, handles)
% hObject    handle to x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x as text
%        str2double(get(hObject,'String')) returns contents of x as a double


% --- Executes during object creation, after setting all properties.
function x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_Callback(hObject, eventdata, handles)
% hObject    handle to y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y as text
%        str2double(get(hObject,'String')) returns contents of y as a double


% --- Executes during object creation, after setting all properties.
function y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cancelButton.
function cancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Assign isOk to false to signify that the addition was not accepted
handles.isOk = false;
guidata(hObject,handles);
uiresume(handles.figure1);


% --- Executes on button press in okButton.
function okButton_Callback(hObject, eventdata, handles)
% hObject    handle to okButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Check to see if the values entered are numbers
if(all(isstrprop([get(handles.x,'string'),get(handles.y,'string')],'digit')))
    handles.isOk = true;
    guidata(hObject,handles);
    uiresume(handles.figure1);
else
    msgbox('Entries must be numeric')
end
