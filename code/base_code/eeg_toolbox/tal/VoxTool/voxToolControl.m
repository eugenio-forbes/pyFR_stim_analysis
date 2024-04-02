function varargout = voxToolControl(varargin)
% VOXTOOLCONTROL MATLAB code for voxToolControl.fig - the main VoxTool
% window that contains controls assigning names to the electrodes in the 3d
% CT figure.
%
%      VOXTOOLCONTROL, by itself, creates a new VOXTOOLCONTROL or raises the existing
%      singleton*.
%
%      H = VOXTOOLCONTROL returns the handle to a new VOXTOOLCONTROL or the handle to
%      the existing singleton*.
%
%      VOXTOOLCONTROL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VOXTOOLCONTROL.M with the given input arguments.
%
%      VOXTOOLCONTROL('Property','Value',...) creates a new VOXTOOLCONTROL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before voxToolControl_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to voxToolControl_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help voxToolControl

% Last Modified by GUIDE v2.5 15-Aug-2013 17:04:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @voxToolControl_OpeningFcn, ...
                   'gui_OutputFcn',  @voxToolControl_OutputFcn, ...
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


% --- Executes just before voxToolControl is made visible.
function voxToolControl_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to voxToolControl (see VARARGIN)

% Assign this figure to controlFigure so it can be referenced from other
% figures
handles.controlFigure = hObject;

% If arguments are passed in, it can make the 3d CT view automatically.
% Mostly for debugging purposes, so we don't have to wait for the .img to
% be loaded
if nargin>3
    if strcmp(varargin{1},'debug')
        %correct for 0-indexing
        handles.xyz = varargin{2}-1;
        handles.vol = varargin{3};
        handles.conversion = varargin{4};
        handles.subjName = 'TJ039';
        handles.dicomDir = '/Users/iped/Desktop/TJ039/images/2012-05-04_post_implant_CT_axial';
        handles.dicomFiles = dir([handles.dicomDir,'/*.dcm']);
        handles.subjDir = '.';
    else
        disp('HELLO!')
        handles.subjName = varargin{1};
        handles.subjDir = varargin{2};
        handles.dicomDir = varargin{3};
        handles.dicomFiles = dir([handles.dicomDir,'/*.dcm']);
        if length(varargin)>3
            handles.jacksheetDir = varargin{4};
        end
    end
        
else
    handles.subjName = '';
    handles.subjDir = '';
end

% setPoints contains information about the already-assigned electrodes
handles.setPoints = struct();
handles.setPoints.coords = [];
handles.setPoints.names = {};
handles.setPoints.baseNames = {};
handles.setPoints.groupNum = [];
handles.setPoints.type = [];
handles.setPoints.size = [];

% Grouped so that they can disappear when the select multiple button is
% unchecked
handles.multipleSelectUIs = [handles.pointsUpButton, handles.pointsDownButton, handles.pointNumEdit];
% Don't think I ever use this
handles.jacksheetElectrodes = {};
% Electrode names that have already been assigned. Don't think I use it.
% Redundant due to setPoints
handles.assignedElectrodes = {};
% Quick way to determine if the multiple selection box is checked without
% doing a 'get'
handles.selectMultiple = false;
% Don't think I use this
handles.basename = {};
% For future use, in case I want to delete points
handles.removedPoints = [];

% If nothing in the window is selected, the arrow keys should do the same
% thing as they do on the 3d CT window. Allows you to move around the
% selected point without clicking back to the CT.
set(gcf,'keyPressFcn', {@keypress});
set(handles.pointsUpButton,'keyPressFcn', {@keypress});
set(handles.pointsDownButton,'keyPressFcn', {@keypress});

% Close all windows if you close this window
set(gcf,'CloseRequestFcn',@onClose);

% Colormap determines the colors that are chosen for the already-assigned
% electrodes
%colormap_tmp = jet(30);
%handles.colormap = colormap_tmp(randperm(30),:);
handles.colormap = distinguishable_colors(50, [0 0 0;.5 .5 .5; 1 1 1; 0 0 1]);
handles.colormap_labels = cell(50,1);

% If command line arguments have been passed in, make the other two
% windows. Otherwise, it waits for an image file to be assigned.
if nargin>3 && strcmp(varargin{1},'debug')
    CT3D_handles = voxTool3DCT(handles);
    handles = concatStructs(handles, CT3D_handles);
    slice_handles = voxToolSlices;
    handles = concatStructs(handles, slice_handles);
end



%Flag to mark that handles have not been saved
handles.saved = false;

if isfield(handles,'jacksheetDir')
    loadJacksheet(handles.jacksheetDir, handles);
end

%Update handles structure
guidata(hObject, handles);

uiwait(hObject);

% --- Outputs from this function are returned to the command line.
function varargout = voxToolControl_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
handles = guidata(hObject);
if isfield(handles,'voxLoc')
    varargout{1} = handles.voxLoc;
else
    varargout{1} = false;
end
delete(hObject);

% --- Executes on button press in loadImageButton.
function loadImageButton_Callback(hObject, eventdata, handles)
% Loads an .img file for analysis
%
% hObject    handle to loadImageButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,pathname] = uigetfile('*.nii.gz;*.nii;*.bhdr;*.img;*.mgz;*.mgh','Load image file',handles.subjDir);

% filename returns 0 if cancel is hit
if filename
    % make the two windows, set the title string to 'Loading...'
    CT3D_handles = voxTool3DCT(handles);
    handles = concatStructs(handles, CT3D_handles);
    slice_handles = voxToolSlices;
    handles = concatStructs(handles, slice_handles);
    if ~isfield(handles,'subjName')
        splitPath = regexp(pathname,filesep,'split');
        handles.subjName = splitPath{length(splitPath)-2};
    end
    set(handles.titleText,'String','Loading...');
    
    % Pause for a split second to give the loading message time to appear.
    pause on
    pause(.5)
    pause off
    fullFile = [pathname,filesep,filename];
    
    % Data is actually read in
    data=MRIread(fullFile);
    % These are the only important bits from MRIread. TODO: cut down
    % MRIread to only the parts that are necessary?
    handles.vol=data.vol;
    handles.conversion = data.vox2ras1(1:3,1:3);
    handles.mriData = data;
    
    try
        threshLimit = thresholdDlg;
    catch e
        threshLimit = 1;
    end
    % Gets the coordinates of only the spots that have the maximum value in
    % the CT
    [x,y,z] = (ind2sub(size(data.vol), find(data.vol>=threshLimit.*max(max(max(data.vol))))));
    %Correct for 0-indexing
    handles.xyz=[x,y,z]-1;
    
    guidata(hObject, handles);
    %Brings up the image on the display
    load_image(handles);
    handles = guidata(hObject);
    
    %Sets the title to the filename
    set(handles.titleText,'String',filename);
    
    %Dicoms cannot be used until after the img is loaded, or else it throws
    %an error
    set(handles.loadMRIButton, 'enable','on')
end

%Save handles.
guidata(hObject, handles);

% --- Executes on button press in loadMatButton.
function loadMatButton_Callback(hObject, eventdata, handles)
% loads a .mat file with the handles from a previous execution of the
% program
% TODO: Make sure the assigned points still work.
%
% hObject    handle to loadMatButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,pathname] = uigetfile('*.mat');
if filename
    data = load([pathname,'/',filename]);
    handles2 = data.handles;
    handles.xyz = handles2.xyz;
    handles.conversion = handles2.conversion;
        handles.setPoints = handles2.setPoints;
    handles.selectedPoints = handles2.selectedPoints;
    handles.subjName = handles2.subjName;
    
    %TODO: SHOULD LOAD UP THE SAVED POINTS ONTO THE 3DCT AND LIST
    slice_handles = voxToolSlices(handles);
    handles = concatStructs(handles, slice_handles);
    CT3D_handles = voxTool3DCT(handles);
    handles = concatStructs(handles, CT3D_handles);

    set(handles.titleText,'string',handles.subjName);
    try %#ok<TRYNC> 
        % Needs the try statement in case the old structure didn't have an
        % unassigned electrodes field
    handles.unassignedElectrodes = handles2.unassignedElectrodes;
    end
    %***********************
    % TODO: CHECK IF DICOM DIR STILL EXISTS
    %***********************
    handles.dicomDir = handles2.dicomDir;
    handles.dicomFiles = handles2.dicomFiles;
    guidata(hObject, handles);
    set(handles.loadMRIButton,'enable','on');
end

% --- Executes on button press in loadJacksheetButton.
function loadJacksheetButton_Callback(hObject, eventdata, handles)
% Loads a jacksheet with (almost) all of the electrodes to be assigned
%
% hObject    handle to loadJacksheetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%***************************
%    TODO: CHECK IF LIST IS ALREADY POPULATED
%***************************

[filename,pathname] = uigetfile('*.txt');
if filename
    handles = loadJacksheet([pathname '/' filename]);
    handles.jacksheetDir = [pathname '/' filename];
    guidata(hObject,handles);
end

function handles = loadJacksheet(jacksheet, handles)
fid = fopen(jacksheet);
if fid~=-1
    electrodes = textscan(fid, '%s','Delimiter','\n');
    electrodes = electrodes{1};
    handles.jacksheetElectrodes = cell(length(electrodes),1);
    for i=1:length(electrodes)
        splitElectrodeName = regexp(electrodes{i},' ','split');
        handles.jacksheetxElectrodes{i} = splitElectrodeName{2};
    end
    handles.unassignedElectrodes = handles.jacksheetElectrodes;
    unassignedElectrodeStrings = handles.jacksheetElectrodes;
    set(handles.unassignedList, 'String', unassignedElectrodeStrings);
    handles.assignedElectrodes = {};

    guidata(handles.controlFigure, handles);
end

% --- Executes on button press in saveMatButton.
function saveMatButton_Callback(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to saveMatButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file, path] = uiputfile('*.mat','Save file name');
save([path,'/',file],'handles');



function electrodeNameEdit_Callback(hObject, eventdata, handles)
% hObject    handle to electrodeNameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of electrodeNameEdit as text
%        str2double(get(hObject,'String')) returns contents of electrodeNameEdit as a double
assignButton_Callback(hObject,eventdata,handles);



% --- Executes on button press in assignButton.
function assignButton_Callback(hObject, eventdata, handles)
% When the assign button is hit (or enter is pressed in the assignement 
% field), this deals with the actual assigning
%
% hObject    handle to assignButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Does not allow clicks within a second. In a try statement because toc
% fails if there was no previous tic
try %#ok<TRYNC>
    if toc<1
        return
    end
end
tic

% For some reason handles was not reliably passed into this function 
% (once, at least).
handles = guidata(hObject);

newElecBasename = get(handles.electrodeNameEdit,'string');

if isempty(newElecBasename)
    return
end

% Even though it's less efficient, newElecNames has to start out empty or
% else the ismember function complains. I could probably get around this
% with some more work.
newElecNames = {};
% I think this is useless.
newElecNamesCorrected = cell(size(handles.selectedPoints.coords,1),1);

% Loop over the selected points
for e=1:size(handles.selectedPoints.coords,1)
    if handles.selectMultiple
        % Look over all of the numbers up until the number of electrodes
        % that have already been assigned with that basename. Fills in gaps
        % if they exist.
        num = nnz(ismember(handles.setPoints.baseNames,newElecBasename))+length(newElecNames);
        for i=1:num+1
            if ~any(ismember(handles.setPoints.names, [newElecBasename, num2str(i)])) && ...
                    (isempty(newElecNames) || ~any(ismember(newElecNames, [newElecBasename, num2str(i)])))
                newElecNames{e} = [newElecBasename,num2str(i)]; %#ok<AGROW>
                break
            end
        end
    else
        newElecNames{e} = newElecBasename; %#ok<AGROW>
    end
    newElecName = newElecNames{e};
    % Checks to see if the assigned electrode name is actually on the
    % jacksheet
    % TODO: only check this once per basename, instead of for every electrode
    if isfield(handles, 'unassignedElectrodes')
        unAssElecNames = handles.unassignedElectrodes;
        assElecNames = handles.setPoints.baseNames;
        [i,~] = find(strcmp(newElecName, unAssElecNames));

        if isempty(i)
            [i,~] = find(strcmp(newElecName, assElecNames));
            if isempty(i)
                choice = menu(['Electrode ',newElecName,' not on jacksheet. Continue?'],'Yes','No');
            else
                choice = menu(['Electrode ',newElecName,' already assigned. Continue?'],'Yes','No');
            end
            if choice ~=1
                return
            else
                val = newElecName;
            end 
        else
            val = unAssElecNames(i,:);
        end
    else
        % This will only happen if selectMultiple is false, I think.
        i = find(strcmp(newElecName, handles.setPoints.baseNames),1);
        if ~isempty(i)
            errordlg(['Electrode ',newElecName,' is already assigned.']);
            return
        end
        val=newElecName;
    end
    newElecNamesCorrected{e} = val;
    h = findobj(handles.mainAxes,'tag',['active',num2str(e)]);
    delete(h);
end
handles = assignElectrodes(handles, newElecNamesCorrected);

set(handles.electrodeNameEdit,'string','');
handles = clearPoints(handles);
guidata(hObject,handles);

% --- Executes on selection change in assignedList.
function assignedList_Callback(hObject, eventdata, handles)
% Called when new selections are made on the assigned list.
% hObject    handle to assignedList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns assignedList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from assignedList

% Get the numbers of the selected items
selectedNums = get(hObject,'value');

% Delete the previously selected object
h = findobj(handles.mainAxes,'tag','selected');
if ~isempty(h)
    delete(h)
end

% Build up a list of all of the coords that are selected
coords = [];
for i=1:length(selectedNums)
    if ~isempty(handles.setPoints.coords)
        selected = selectedNums(i);
        coords(i,:)=handles.setPoints.coords(selected,:);
    end
end

% Circle the selected coordinates
if ~isempty(coords)
    locatePoint(handles, coords, true, 'red', 'selected')
    uicontrol(hObject)
end

% --- Executes on button press in saveVoxButton.
function saveVoxButton_Callback(hObject, eventdata, handles)
% Saves the VOX_COORDS_MOTHER.txt file
% hObject    handle to saveVoxButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%TODO: INTEGRATE PATH FROM MAIN WINDOW
baseNames = {};
setSizes = [];
realSizes = [];
for i=1:length(handles.setPoints.baseNames)
    baseName = handles.setPoints.baseNames{i};
    if ~any(strcmp(baseNames,baseName))
        baseNames{length(baseNames)+1} = baseName;
        setSizes(length(setSizes)+1) = ...
            handles.setPoints.size(i,1)*handles.setPoints.size(i,2);
        realSizes(length(realSizes)+1) = 0;
    end
    realSizes(strcmp(baseNames, baseName)) = ...
        realSizes(strcmp(baseNames, baseName))+1;
end
if ~all(realSizes==setSizes)
    badBasename = baseNames{realSizes~=setSizes};
    button = questdlg(sprintf('Electrode group %s does not match sizse. Continue anyway?',...
        badBasename),'Electrode size mismatch','Yes','No','No');
    if strcmp(button,'No')
        return
    end

end
    
[filename, path] = uiputfile('VOX_coords_mother.txt','Save file...');
fid = fopen([path,'/',filename],'w');

for i=1:length(handles.setPoints.names)
    name = handles.setPoints.names{i};
    coord = handles.setPoints.coords(i,:);
    x = num2str(coord(1));
    y = num2str(coord(2));
    z = num2str(coord(3));
    type = handles.setPoints.type(i);
    size = handles.setPoints.size(i,:);
    s1 = size(1);
    s2 = size(2);
    % Ugh. Has to be y,x,z to keep it consistent.
    fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%d %d\n', name,y,x,z,type, s1, s2);
end
fclose(fid);

handles.saved = true;
handles.voxLoc = [path,'/',filename];
guidata(hObject,handles);

% --- Executes on button press in loadMRIButton.
function loadMRIButton_Callback(hObject, eventdata, handles)
% Loads the folder with the dicoms in it.
% hObject    handle to loadMRIButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
foldername = uigetdir;
if ~isempty(foldername)
    % I don't think dicomDir is actually used anywhere else
    handles.dicomDir = foldername;
    handles.dicomFiles = dir([handles.dicomDir,'/*.dcm']);
end
guidata(hObject, handles);
guidata(handles.CT3D_fig, handles);

% --- Executes on button press in selectMutlipleCB.
function selectMutlipleCB_Callback(hObject, eventdata, handles)
% Runs when the select multiple checkbox is changed. 
%
% hObject    handle to selectMutlipleCB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of selectMutlipleCB

% change the value of selectMultiple
handles = guidata(hObject);
handles.selectMultiple = get(hObject,'value');

% Make multiple selections on the unassigned list possible/impossible based
% on the value of the checkbox. Also turn off/on the up/down arrows.
if get(hObject,'value')
    set(handles.unassignedList,'max',10);
    set(handles.unassignedList,'min',0);
    set(handles.multipleSelectUIs,'visible','on');
else
    set(handles.unassignedList,'max',0);
    set(handles.unassignedList,'min',10);
    set(handles.multipleSelectUIs,'visible','off');
    handles = clearPoints(handles);
end

% Save this value to both the control figure and CT figure
guidata(hObject,handles);
guidata(handles.CT3D_fig, handles);



function Z1_Callback(hObject, eventdata, handles)
% Runs when the Z coordinate is changed
%
% hObject    handle to Z1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Z1 as text
%        str2double(get(hObject,'String')) returns contents of Z1 as a double
coords = handles.selectedPoints.coords;
currNum =str2double(get(handles.pointNumEdit,'string'));
coords(currNum, 3) = str2double(get(hObject,'string'));
tag = ['active',num2str(currNum)];
locatePoint(handles, coords(currNum,:), true, 'red',tag);



function Y1_Callback(hObject, eventdata, handles)
% Runs when the Y coordinate is changed
%
% hObject    handle to Y1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Y1 as text
%        str2double(get(hObject,'String')) returns contents of Y1 as a double
coords = handles.selectedPoints.coords;
currNum =str2double(get(handles.pointNumEdit,'string'));
coords(currNum, 2) = str2double(get(hObject,'string'));
tag = ['active',num2str(currNum)];
locatePoint(handles, coords(currNum,:), false, 'red',tag);



function X1_Callback(hObject, eventdata, handles)
% Runs when the X coordinate is changed
%
% hObject    handle to X1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of X1 as text
%        str2double(get(hObject,'String')) returns contents of X1 as a double
coords = handles.selectedPoints.coords;
currNum =str2double(get(handles.pointNumEdit,'string'));
coords(currNum, 1) = str2double(get(hObject,'string'));
tag = ['active',num2str(currNum)];
locatePoint(handles, coords(currNum,:), false, 'red',tag);



function pointNumEdit_Callback(hObject, eventdata, handles)
% hObject    handle to pointNumEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pointNumEdit as text
%        str2double(get(hObject,'String')) returns contents of pointNumEdit as a double

%**************************
% TODO: WRITE THIS WHOLE FUNCTION. SHOULDN'T BE HARD
%**************************

function pointNumEdit_CreateFcn(hObject, eventdata, handles)


% --- Executes on button press in pointsUpButton.
function pointsUpButton_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% Callback for the button to shift to the next selected point
%
% hObject    handle to pointsUpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currNum = str2double(get(handles.pointNumEdit,'string'));
updateCoords(handles,currNum+1);
updateDicoms(handles, handles.selectedPoints.coords(currNum+1,:),true);

% --- Executes on button press in pointsDownButton.
function pointsDownButton_Callback(hObject, eventdata, handles)
% Callback for the button to shift to the previously selected point
%
% hObject    handle to pointsDownButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currNum = str2double(get(handles.pointNumEdit,'string'));
updateCoords(handles, currNum-1);
updateDicoms(handles, handles.selectedPoints.coords(currNum-1,:),true);

% --- Executes on button press in clearPointsButton.
function clearPointsButton_Callback(hObject, eventdata, handles)
% Clears the currently selected points
%
% hObject    handle to clearPointsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clearPoints(handles);

% --- Executes on key press with focus on unassignedList and none of its controls.
function unassignedList_KeyPressFcn(hObject, eventdata, handles)
% Runs when you hit enter on the unassigned list
%
% hObject    handle to unassignedList (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

if strcmp(eventdata.Key,'return')
    selectedNums = get(hObject,'value');
    oldContents = get(hObject,'string');
    selectedValues = oldContents(selectedNums);
    if size(selectedValues,1)~=size(handles.selectedPoints.coords,1)
        errordlg('Number of points selected does not equal number of electrodes selected','ERROR')
    else
        handles = assignElectrodes(handles, selectedValues);
        set(handles.unassignedList,'value',1);
    end
end


% --- Executes on selection change in elecSizePopup.
function elecSizePopup_Callback(hObject, eventdata, handles)
% Runs when the popup box choosing the electrode size is changed
%
% hObject    handle to elecSizePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns elecSizePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from elecSizePopup
contents = cellstr(get(hObject,'string'));
% If the user selects 'New', show the box to get new coordinates
if strcmp(contents{(get(hObject,'Value'))},'New')
    [coords, window] = newCoordsBox();
    delete(window);
    contents{size(contents,1)} = coords;
    contents{size(contents,1)+1} = 'New';
    set(hObject,'string',contents);
end

% --- Executes during object creation, after setting all properties.
function elecSizePopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elecSizePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function onClose(src, event)
% On close, try to delete the CT and Dicom figures, then close this window
handles = guidata(src);
 if ~handles.saved
     if ~isfield(handles,'voxLoc')
         choice = menu('VOX coords has not been saved. Continue?','Yes','No');
         if choice~=1
             return
         else
             handles.voxLoc = false;
         end
    else
         choice = menu('VOX coords has changed since last save. Continnue?','Yes','No');
         if choice~=1
             return
         end
    end
 end
try %#ok<TRYNC>
    delete(handles.CT3D_fig);
end
try %#ok<TRYNC>
    delete(handles.sliceFig);
end
uiresume(handles.controlFigure);



% --- Executes on key press with focus on assignedList and none of its controls.
function assignedList_KeyPressFcn(hObject, eventdata, handles)
% Deletes a point from the assigned list when delete is pressed
%
% hObject    handle to assignedList (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

if strcmp(eventdata.Key,'delete') || strcmp(eventdata.Key,'backspace')
    unassignElectrodes(handles, get(hObject,'value'))
    set(hObject,'value',1);
    delete(findobj(handles.mainAxes,'tag','selected'))
end


% --- Executes on button press in saveAndQuitButton.
function saveAndQuitButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveAndQuitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
saveVoxButton_Callback(hObject, eventdata, handles);
handles = guidata(hObject);
 if ~handles.saved
     if ~isfield(handles,'voxLoc')
         choice = menu('VOX coords has not been saved. Continue?','Yes','No');
         if choice~=1
             return
         else
             handles.voxLoc = false;
         end
    else
         choice = menu('VOX coords has changed since last save. Continnue?','Yes','No');
         if choice~=1
             return
         end
    end
 end
try %#ok<TRYNC>
    delete(handles.CT3D_fig);
end
try %#ok<TRYNC>
    delete(handles.sliceFig);
end
uiresume(handles.controlFigure);


% --- Executes on button press in loadVoxCoordsButton.
function loadVoxCoordsButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadVoxCoordsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,pathname] = uigetfile('*.txt','Load vox_coords_mother',handles.subjDir);
if filename
    fid = fopen([pathname, filesep, filename]);
    lines = textscan(fid, '%s\t%d\t%d\t%d\t%s\t%s %s');
    names = lines{1};
    xs = lines{2};
    ys = lines{3};
    zs = lines{4};
    types = lines{5};
    size1 = lines{6};
    size2 = lines{7};
    possibleTypes = get(get(handles.electrodeTypePanel,'children'),'string');
    possibleTypeRBs = get(handles.electrodeTypePanel,'children');

    for i =1:length(names)
        possibleSizes = get(handles.elecSizePopup,'string');
        handles = clearPoints(handles);
        handles.selectedPoints.coords = [ys(i) xs(i) zs(i)];
        foundType = false;
        for t=1:length(possibleTypes)
            this_possibleType = possibleTypes{t};
            if strcmp(upper(this_possibleType(1)),upper(types{i}))
                set(handles.electrodeTypePanel,'SelectedObject', possibleTypeRBs(t))
                foundType = true;
                break
            end
        end
        if ~foundType
            error(['TYPE ' types{i} ' DOES NOT MATCH ANY PRESENT TYPES'])
        end
        foundSize = false;
        for s=1:length(possibleSizes)-1
            this_possibleSize = regexp(possibleSizes{s},'x','split');
            if all([str2double(this_possibleSize{1}) str2double(this_possibleSize{2})] ==  ...
                    [str2double(size1{i}) str2double(size2{i})])
                set(handles.elecSizePopup,'value',s)
                foundSize = true;
                break
            end
        end
        if ~foundSize
            contents = cellstr(get(handles.elecSizePopup,'string'));
            contents{size(contents,1)} = [num2str(size1{i}) 'x' num2str(size2{i})];
            contents{size(contents,1)+1} = 'New';
            set(handles.elecSizePopup, 'string',contents);
            set(handles.elecSizePopup, 'value' , length(contents)-1);
        end
        newElecName = names(i);
        handles = assignElectrodes(handles, newElecName);
    end
end
    


% --- Executes on button press in niigzButton.
function niigzButton_Callback(hObject, eventdata, handles)
% hObject    handle to niigzButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newVol = zeros(size(handles.vol));
maxVal = max(max(max(handles.vol)));
for i=1:size(handles.setPoints.coords,1)
    thesecoords = handles.setPoints.coords(i,:);
    [~,cluster] = getPointCluster(handles, thesecoords');
    for j=1:size(cluster,1)
        newVol(cluster(j,1),cluster(j,2),cluster(j,3))=maxVal;
    end
end
mri = handles.mriData;
mri.vol = newVol;
[file, path] = uiputfile('*.nii.gz','Save file name');
file = [file(1:length(file)-3),'.nii.gz'];
MRIwrite(mri, [path,file])
