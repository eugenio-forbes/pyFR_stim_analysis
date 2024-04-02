function updateCoords( handles, number )
%UPDATECOORDS Update the coordinates box with the selected coordinate
%   number - the number of the selected coordinate (generally the number
%   specified in the box next to the coordinates

% Get the value that the box used to be set to
oldVal = get(handles.pointNumEdit,'string');

% Sets the box by the coordinates to reflect the right number
set(handles.pointNumEdit,'string',num2str(number));
coords = handles.selectedPoints.coords;
boxHandles = [handles.X1, handles.Y1, handles.Z1];

%Error checking to make sure the coordinate is not set to something weird.
%Also, enables and disables the up and down buttons.
if size(coords, 1)==0
   set(boxHandles, 'string','')
   return
end
if number==size(coords,1)
    set(handles.pointsUpButton,'enable','off');
else
    set(handles.pointsUpButton,'enable','on');
end
if number==1
    set(handles.pointsDownButton,'enable','off');
elseif number>size(coords,1) || number<=0
    return
else
    set(handles.pointsDownButton,...
        'enable','on');
end

%Set each of the coordinate boxes with the right coordinate.
for i=1:3
    set(boxHandles(i),'string',num2str(coords(number, i)));
end

%Delete the old active circle, show the new one.
set(findobj(handles.mainAxes,'Tag',['active',oldVal]),'color',[1 .5 0],...
    'linewidth',1);

set(findobj(handles.mainAxes,'Tag',['active',num2str(number)]),'color','red',...
    'linewidth',2);

%Save data.
guidata(handles.controlFigure, handles);
guidata(handles.CT3D_fig, handles);

