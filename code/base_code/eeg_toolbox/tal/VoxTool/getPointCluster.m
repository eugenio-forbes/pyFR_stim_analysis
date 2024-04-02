function [new_point, relevant_points_slices] = getPointCluster(handles, point, iterations)
%GETPOINTCLUSTER gets the point in the center of the cluster of points that
%have been clicked on
% inputs:
%    handles - must contain handles.xyz (the pointCloud), 
%              handles.surfaceRB (whether it is a surface electrode),
%              handles.conversion (allows for conversion between slices 
%              and mm)
%
%    point  -  the point at the center of the clustering algorithm
%
%    iterations - The number of iterations the algorithm goes through to
%                 find the cluster.
%
% outputs:
%    newer_point - the point that winds up at the center of mass of the
%                  cluster
%    relevant_points_slices - the points that are in the cluster

% Convert the point and the pointCloud to mms
xyz = handles.xyz*handles.conversion;
point = double(point)'*handles.conversion;

% weight_radius defines the radius in which to look for a cluster
if get(handles.surfaceRB,'value') || get(handles.gridRB,'value')
    weight_radius = 3.25;
else
    weight_radius = 2;
end

if ~exist('iterations','var')
    iterations = 3;
end

% replicate the point for the length of the cluster, then get the distances
% between all points and point. relevant_points_mms gives the points within
% the cluster in mms. Loop as many times as needed.
for i=1:iterations
    rep_point = repmat(point,size(xyz,1),1);
    diff = sqrt(sum((rep_point-xyz).^2,2));
    relevant_points_mms = xyz(diff<weight_radius,:);

    % get the mean position of all of the points in the cluster
    point = round(mean(relevant_points_mms,1));

end

%delete the previously active area
h = findobj(handles.mainAxes, 'tag','activeArea');
if ~isempty(h)
    delete(h)
end

%convert back to CT coords from mms, plot these points, label as the
%activeArea
relevant_points_slices = relevant_points_mms/handles.conversion;
h = plot3(relevant_points_slices(:,1), relevant_points_slices(:,2),...
    relevant_points_slices(:,3), 'b.');
set(h,'tag','activeArea')
uistack(h,'top');

new_point = round(mean(relevant_points_mms,1));
new_point = new_point/handles.conversion;
new_point = round(new_point);
end

