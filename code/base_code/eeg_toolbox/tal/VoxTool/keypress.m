function keypress(src, eventData)
%KEYPRESS the actions taken upon a keypress while focused on the CT figure
%  mostly for moving the points around in space with the arrow keys

handles = guidata(gcf);
handles = concatStructs(guidata(handles.controlFigure),handles);
coords = handles.selectedPoints.coords;
if handles.selectMultiple
    currNum = str2double(get(handles.pointNumEdit,'string'));
else
    currNum = 1;
end
% Try to get the current coordinates. If you can't, just leave.
try
    X = coords(currNum,1);
    Y = coords(currNum,2);
    Z = coords(currNum,3);
catch %#ok<CTCH>
    return
end

% Move the coordinates upon clicks of the arrow keys, changing the dicom
% image only if the Z value changes.
if strcmp(eventData.Key,'pagedown')
    Z = Z-1;
    changeDicoms = true;
elseif strcmp(eventData.Key,'pageup')
    Z = Z+1;
    changeDicoms = true;
elseif strcmp(eventData.Key,'uparrow')
    X = X-1;
    changeDicoms = false;
elseif strcmp(eventData.Key,'downarrow')
    X = X+1;
     changeDicoms = false;
elseif strcmp(eventData.Key,'rightarrow')
    Y = Y+1;
     changeDicoms = false;
elseif strcmp(eventData.Key,'leftarrow')
    Y = Y-1;
     changeDicoms = false;
elseif strcmp(eventData.Key,'return')
    uicontrol(handles.electrodeNameEdit);
     changeDicoms = false;
else
    return
end
coords(currNum,1) = X;
coords(currNum,2) = Y;
coords(currNum,3) = Z;
handles.selectedPoints.coords = coords;
updateCoords(handles, currNum);
guidata(src, handles);
tag=['active',num2str(currNum)];
locatePoint(handles, [X Y Z], changeDicoms,'red',tag);

end

