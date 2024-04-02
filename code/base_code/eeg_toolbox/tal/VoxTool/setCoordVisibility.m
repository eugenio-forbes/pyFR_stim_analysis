function setCoordVisibility( handles, coords, visible )
%SETCOORDVISIBILITY is supposed to turn on or off points within the point
%cloud. Not yet implemented

%TODO: Implement.

points = handles.pointCloud';
badRows = ismember(points, coords, 'rows');
handles.removedPoints=[handles.removedPoints;points(badRows,:)];
delete(findobj(handles.mainAxes,'tag','fullBrain'));
handles.pointCloud = points(~badRows,:)';
plotWithLimits(handles);
end

