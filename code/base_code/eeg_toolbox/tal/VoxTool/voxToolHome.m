function varargout = voxToolHome(varargin)
% VOXTOOLHOME MATLAB code for voxToolHome.fig
%      VOXTOOLHOME, by itself, creates a new VOXTOOLHOME or raises the existing
%      singleton*.
%
%      H = VOXTOOLHOME returns the handle to a new VOXTOOLHOME or the handle to
%      the existing singleton*.
%
%      VOXTOOLHOME('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VOXTOOLHOME.M with the given input arguments.
%
%      VOXTOOLHOME('Property','Value',...) creates a new VOXTOOLHOME or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before voxToolHome_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to voxToolHome_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help voxToolHome

% Last Modified by GUIDE v2.5 22-Jun-2013 15:32:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @voxToolHome_OpeningFcn, ...
                   'gui_OutputFcn',  @voxToolHome_OutputFcn, ...
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


% --- Executes just before voxToolHome is made visible.
function voxToolHome_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to voxToolHome (see VARARGIN)

% Choose default command line output for voxToolHome
handles.homeFigure = hObject;
handles.output = hObject;
handles.isRegistered = false;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes voxToolHome wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = voxToolHome_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function subjNameEdit_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to subjNameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of subjNameEdit as text
%        str2double(get(hObject,'String')) returns contents of subjNameEdit as a double

if ~isempty(get(handles.subjNameEdit,'string'))
    handles.subjName = get(hObject,'string');
    if isempty(get(handles.subjSourceEdit,'string'))
        set(handles.subjSourceEdit,'string',handles.subjName);
        handles.subjSourceName = handles.subjName;
    end
elseif isfield(handles,'subjName')
    handles = rmfield(handles,'subjName');
end
guidata(hObject,handles)
activateButtons(handles);

% --- Executes during object creation, after setting all properties.
function subjNameEdit_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to subjNameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function subjSourceEdit_Callback(hObject, eventdata, handles)
% hObject    handle to subjSourceEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of subjSourceEdit as text
%        str2double(get(hObject,'String')) returns contents of subjSourceEdit as a double
if ~isempty(get(handles.subjSourceEdit,'string'))
    handles.subjSourceName = get(hObject,'string');
elseif isfield(handles,'subjSourceName')
    handles =rmfield(handles,'subjSourceName');
end
guidata(hObject,handles)
activateButtons(handles);

% --- Executes during object creation, after setting all properties.
function subjSourceEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subjSourceEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ctDirEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ctDirEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ctDirEdit as text
%        str2double(get(hObject,'String')) returns contents of ctDirEdit as a double
if ~isempty(get(hObject,'string'))
    handles.ctDir = get(hObject,'string');
    handles = setSubjName(handles);
elseif isfield(handles,'ctDir')
    handles = rmfield(handles,'ctDir');
end
guidata(hObject,handles)
activateButtons(handles);

% --- Executes during object creation, after setting all properties.
function ctDirEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ctDirEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ctDirButton.
function ctDirButton_Callback(hObject, eventdata, handles)
% hObject    handle to ctDirButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'mriDir') || isfield(handles,'ctDir')
    lastDir = getLastUsedDirectory(handles);
    path = uigetdir(lastDir);
else
    path = uigetdir();
end
if path
    set(handles.ctDirEdit,'string',path)
    handles.ctDir = path;
    handles = setSubjName(handles);
end
guidata(hObject,handles)
activateButtons(handles);

function mriDirEdit_Callback(hObject, eventdata, handles)
% hObject    handle to mriDirEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mriDirEdit as text
%        str2double(get(hObject,'String')) returns contents of mriDirEdit as a double
if ~isempty(get(hObject,'string'))
    handles.mriDir = get(hObject,'string');
    handles = setSubjName(handles);
elseif isfield(handles,'mriDir');
    handles = rmfield(handles,'mriDir');
end
guidata(hObject,handles)
activateButtons(handles);

% --- Executes during object creation, after setting all properties.
function mriDirEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mriDirEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in mriDirButton.
function mriDirButton_Callback(hObject, eventdata, handles)
% hObject    handle to mriDirButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'mriDir') || isfield(handles,'ctDir')
    lastDir = getLastUsedDirectory(handles);
    disp('here');
    path = uigetdir(lastDir);
else
    path = uigetdir();
end
if path
    set(handles.mriDirEdit,'string',path)
    handles.mriDir = path;
    handles = setSubjName(handles);
end
activateButtons(handles);
guidata(hObject,handles)


% --- Executes on button press in localizeButton.
function localizeButton_Callback(hObject, eventdata, handles)
% hObject    handle to localizeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.homeFigure,'visible','off')
if isfield(handles,'jacksheet')
vmLoc = voxToolControl(handles.subjName, handles.subjDir, handles.ctDir,...
    handles.jacksheet);
else
    vmLoc = voxToolControl(handles.subjName, handles.subjDir, handles.ctDir);
end
if vmLoc
    handles.voxMother = vmLoc;
    set(handles.voxMotherEdit,'string',handles.voxMother);
end
set(handles.homeFigure,'visible','on');

% --- Executes on button press in registerButton.
function registerButton_Callback(hObject, eventdata, handles)
% hObject    handle to registerButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
convertMotherToVoxCoords(handles.username, handles.voxMother, handles.jacksheet,...
    handles.subjName, handles.subjDir, handles.subjSourceName, ...
    handles.sourceDir)
set(handles.updateButton,'enable','on');
% --- Executes on button press in updateButton.
function updateButton_Callback(hObject, eventdata, handles)
% hObject    handle to updateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateTalEventsDataBase(handles.subjName, handles.subjDir)


function voxMotherEdit_Callback(hObject, eventdata, handles)
% hObject    handle to voxMotherEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of voxMotherEdit as text
%        str2double(get(hObject,'String')) returns contents of voxMotherEdit as a double
if ~isempty(get(hObject,'string'))
    handles.voxMother = get(hObject,'string');
elseif isfield(handles,'voxMother')
    handles = rmfield(handles,'voxMother');
end
guidata(hObject,handles)
activateButtons(handles);

% --- Executes during object creation, after setting all properties.
function voxMotherEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to voxMotherEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in voxMotherButton.
function voxMotherButton_Callback(hObject, eventdata, handles)
% hObject    handle to voxMotherButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'mriDir') || isfield(handles,'ctDir')
    lastDir = getLastUsedDirectory(handles,2);
    disp('here');
    [filename,path] = uigetfile('*.txt','Select VOX Coords Mother',[lastDir,'/VOX_coords_mother.txt']);
else
    [filename,path] = uigetdir('*.txt','Select VOX Coords Mother');
end
if path
    set(handles.voxMotherEdit,'string',[path,filename])
    handles.voxMother = [path,filename];
end
activateButtons(handles);
guidata(hObject,handles)



function subjDirEdit_Callback(hObject, eventdata, handles)
% hObject    handle to subjDirEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of subjDirEdit as text
%        str2double(get(hObject,'String')) returns contents of subjDirEdit as a double
if ~isempty(get(hObject,'string'))
    handles.subjDir = get(hObject,'string');
    if isempty(get(handles.subjNameEdit,'string'))
        subjDirSplit = regexp(handles.subjDir,'/','split');
        handles.subjName = subjDirSplit{length(subjDirSplit)};
        if isempty(handles.subjName)
            handles.subjName = subjDirSplit{length(subjDirSplit)-1};
        end
        set(handles.subjNameEdit,'string',handles.subjName);
        set(handles.subjSourceEdit,'string',handles.subjName);
        handles.subjSourceName = handles.subjName;
    end
    if isempty(get(handles.sourceDirEdit,'string'))
        set(handles.sourceDirEdit,'string',handles.subjDir);
    end
elseif isfield(handles,'subjDir')
    handles = rmfield(handles,'subjDir');
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function subjDirEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subjDirEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'mriDir') || isfield(handles,'ctDir')
    lastDir = getLastUsedDirectory(handles,2);
    disp('here');
    path = uigetdir(lastDir);
else
    path = uigetdir();
end
if path
    set(handles.subjDirEdit,'string',path)
    handles.subjDir = path;
    handles = setSubjName(handles);
end
activateButtons(handles);
guidata(hObject,handles)



function sourceDirEdit_Callback(hObject, eventdata, handles)
% hObject    handle to sourceDirEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sourceDirEdit as text
%        str2double(get(hObject,'String')) returns contents of sourceDirEdit as a double
if ~isempty(get(hObject,'string'))
    handles.sourceDir = get(hObject,'string');
elseif isfield(handles,'sourceDir')
    handles = rmfield(handles,'sourceDir');
end
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function sourceDirEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sourceDirEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in sourceDirButton.
function sourceDirButton_Callback(hObject, eventdata, handles)
% hObject    handle to sourceDirButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function usernameEdit_Callback(hObject, eventdata, handles)
% hObject    handle to usernameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of usernameEdit as text
%        str2double(get(hObject,'String')) returns contents of usernameEdit as a double
if ~isempty(get(hObject,'string'))
    handles.username = get(hObject,'string');
elseif isfield(handles,'username')
    handles = rmfield(handles,'username');
end
guidata(hObject,handles)
activateButtons(handles);

% --- Executes during object creation, after setting all properties.
function usernameEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to usernameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function jacksheetEdit_Callback(hObject, eventdata, handles)
% hObject    handle to jacksheetEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of jacksheetEdit as text
%        str2double(get(hObject,'String')) returns contents of jacksheetEdit as a double
if ~isempty(get(hObject,'string'))
    handles.jacksheet = get(hObject,'string');
elseif isfield(handles,'jacksheet')
    handles = rmfield(handles,'jacksheet');
end
guidata(hObject,handles)
activateButtons(handles);

% --- Executes during object creation, after setting all properties.
function jacksheetEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to jacksheetEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in jacksheetButton.
function jacksheetButton_Callback(hObject, eventdata, handles)
% hObject    handle to jacksheetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'mriDir') || isfield(handles,'ctDir')
    lastDir = getLastUsedDirectory(handles,2);
    [filename,path] = uigetfile('*.txt','Select Jacksheet',[lastDir,'/jacksheet.txt']);
else
    [filename,path] = uigetdir('*.txt','Select Jacksheet');
end
if path
    set(handles.jacksheetEdit,'string',[path,filename])
    handles.jacksheet = [path,filename];
end
activateButtons(handles);
guidata(hObject,handles)
