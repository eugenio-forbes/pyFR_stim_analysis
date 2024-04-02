function h = clickA3DPoint(pointCloud, handles, limits)
%CLICKA3DPOINT
%   H = CLICKA3DPOINT(POINTCLOUD) shows a 3D point cloud and lets the user
%   select points by clicking on them. The selected point is highlighted 
%   and its index in the point cloud will is printed on the screen. 
%   POINTCLOUD should be a 3*N matrix, represending N 3D points. 
%   Handle to the figure is returned.
%
%   other functions required:
%       CALLBACKCLICK3DPOINT  mouse click callback function
%       ROWNORM returns norms of each row of a matrix
%       
%   To test this function ... 
%       pointCloud = rand(3,100)*100;
%       h = clickA3DPoint(pointCloud);
% 
%       now rotate or move the point cloud and try it again.
%       (on the figure View menu, turn the Camera Toolbar on, ...)
%
%   To turn off the callback ...
%       set(h, 'WindowButtonDownFcn',''); 
%
%   by Babak Taati
%   http://rcvlab.ece.queensu.ca/~taatib
%   Robotics and Computer Vision Laboratory (RCVLab)
%   Queen's University
%   May 4, 2005 
%   revised Oct 30, 2007
%   revised May 19, 2009
%
%   modified for VoxTool by Isaac Pedisich: June, 2013

axes(handles.mainAxes);
if size(pointCloud, 1)~=3
    error('Input point cloud must be a 3*N matrix.');
end

% show the point cloud
if exist('limits','var')
    for i=1:2
        for j=1:3
            if i==1
                pointCloud(:,pointCloud(j,:)>limits(j,i))=[];
            else
                pointCloud(:,pointCloud(j,:)<limits(j,i))=[];
            end
        end
    end
end
h = gcf;

%Shade all of the points based on their X value. This requires splitting
%the points into as many different series as there are unique X values. 
uniqueXs = unique(pointCloud(1,:));
%shading = bone(length(uniqueXs)*2);
%shading = shading(size(shading,1)/2:size(shading,1),:);
shading = bone(double(int32((length(uniqueXs)*1.5))));
shading = shading(end-length(uniqueXs):end,:);


for i=1:length(uniqueXs)
    thisX = uniqueXs(i);
    thisRow = pointCloud(:,pointCloud(1,:)==thisX);
    %Plot each row, tagging all with 'fullBrain'
    if ~isempty(thisRow)
        j = plot3(thisRow(1,:), thisRow(2,:), thisRow(3,:), '.','color',shading(i,:)); 
        set(j,'tag','fullBrain','color',shading(i,:));
    end
    hold on;
end
cameratoolbar(handles.CT3D_fig,'Show'); % show the camera toolbar

% set the callback, pass handles to the callback function, and save.
handles.pointCloud = pointCloud;
guidata(handles.CT3D_fig, handles);
guidata(handles.controlFigure, handles);
set(gcf, 'WindowButtonDownFcn', {@callbackClickA3DPoint, handles}); 
set(handles.mainAxes,'color','black');
