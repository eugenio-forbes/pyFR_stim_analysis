function callbackClickA3DPoint(src, ~, handles)
% CALLBACKCLICK3DPOINT mouse click callback function for CLICKA3DPOINT
%
%   The transformation between the viewing frame and the point cloud frame
%   is calculated using the camera viewing direction and the 'up' vector.
%   Then, the point cloud is transformed into the viewing frame. Finally,
%   the z coordinate in this frame is ignored and the x and y coordinates
%   of all the points are compared with the mouse click location and the 
%   closest point is selected.
%
%   Babak Taati - May 4, 2005
%   revised Oct 31, 2007
%   revised Jun 3, 2008
%   revised May 19, 2009
%   reviesd for VoxTool June 2013

if gca~=handles.mainAxes
    return
end
controlHandles = guidata(handles.controlFigure);
handles = guidata(handles.CT3D_fig);
pointCloud = handles.pointCloud;
point = get(gca, 'CurrentPoint'); % mouse click position
camPos = get(gca, 'CameraPosition'); % camera position
camTgt = get(gca, 'CameraTarget'); % where the camera is pointing to

camDir = camPos - camTgt; % camera direction
camUpVect = get(gca, 'CameraUpVector'); % camera 'up' vector

% build an orthonormal frame based on the viewing direction and the 
% up vector (the "view frame")
zAxis = camDir/norm(camDir);    
upAxis = camUpVect/norm(camUpVect); 
xAxis = cross(upAxis, zAxis);
yAxis = cross(zAxis, xAxis);

rot = [xAxis; yAxis; zAxis]; % view rotation 

% the point cloud represented in the view frame
rotatedPointCloud = rot * pointCloud; 

% the clicked point represented in the view frame
rotatedPointFront = rot * point' ;

% find the nearest neighbour to the clicked point 
pointCloudIndex = dsearchn(rotatedPointCloud(1:2,:)', ... 
    rotatedPointFront(1:2));

% get the selected point
selectedPoint = pointCloud(:, pointCloudIndex); 

% if other points have been selected before, store them to coords
if controlHandles.selectMultiple
    coords = handles.selectedPoints.coords;
else
    coords = [];
end

%get the point that is in the center of the selected cluster
centeredPoint = getPointCluster(handles, selectedPoint);

%check to make sure that points is not already selected
if ~isempty(coords) && any(ismember(coords, centeredPoint,'rows'))
    choice=menu('Two selected electrodes will have same coordinates. Add anyway?','Yes','No');
    if choice~=1
        return
    end
end

%if selecting multiple points, change the color of the previously active
%point from red to orange
if controlHandles.selectMultiple
    old_tag = ['active',num2str(size(coords,1))];
    set(findobj(handles.mainAxes,'Tag',old_tag), 'color',[1,.5,0],'linewidth',1)
end

%append the new point to coords
coords = [coords;centeredPoint(1) centeredPoint(2) centeredPoint(3); ];

%tag it, save it to both the control and CT figures, and circle the point
%on the CT and dicom.
tag = ['active',num2str(size(coords,1))];
handles.selectedPoints.coords = coords;
updateCoords(handles, size(coords,1));
guidata(src, handles);
guidata(handles.controlFigure, handles);
locatePoint(controlHandles, [centeredPoint(1), centeredPoint(2), centeredPoint(3)], true,'red', tag);
