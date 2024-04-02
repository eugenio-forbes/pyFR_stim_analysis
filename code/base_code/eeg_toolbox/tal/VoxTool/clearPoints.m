function handles = clearPoints( handles )
%CLEARPOINTS clears all of the points that are currently selected

%find all of the active coordinates.
coords=  handles.selectedPoints.coords;

%remove the active circles from the CT
for i=1:size(coords,1)
    tag = ['active',num2str(i)];
    h = findobj(handles.mainAxes,'tag',tag);
    delete(h);
end

%remove the active area (selected via clustering) from the CT
h = findobj(handles.mainAxes,'tag','activeArea');
delete(h);

%clear the selected coordinates, update the viewer on the control figure
handles.selectedPoints.coords = [];
updateCoords(handles, 1);
%save to control and CT figures.
guidata(handles.controlFigure,handles);
guidata(handles.CT3D_fig,handles);