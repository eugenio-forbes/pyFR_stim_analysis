function [relevant_points] = getCluster(xyz, radius, point, iterations)
%[relevant_poitns] = GETPOINTCLUSTER(xyz, radius, point, iterations)
%    gets the point in the center of the cluster of points that
%    have been clicked on
%
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
%    relevant_points - the points that are in the cluster
%
% created 10/13 by Isaac Pedisich (iped@sas.upenn.edu)

point = double(point);


if ~exist('iterations','var')
    iterations = 3;
end

% replicate the point for the length of the cluster, then get the distances
% between all points and point. relevant_points_mms gives the points within
% the cluster in mms. Loop as many times as needed.
for i=1:iterations
    rep_point = repmat(point,size(xyz,1),1);
    diff = sqrt(sum((rep_point-xyz).^2,2));
    relevant_points = xyz(diff<radius,:);

    % get the mean position of all of the points in the cluster
    point = round(mean(relevant_points,1));

end



end

