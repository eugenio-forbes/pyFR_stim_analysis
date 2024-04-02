function [ centroids ] = get_centroids( points, radius )
%[centroids] = GET_CENTROIDS(points, radius) 
%       returns the centers of groups of points
%
%   radius - the number of elements away from a center point that a point
%            may be to consider it a part of the same group
%            (defaults to 5)
%
% created by Isaac Pedisich (iped@sas.upenn.edu)

if isempty(points)
    centroids = [];
    return
end
if ~exist('radius','var')
    radius = 5;
end
% Gets the distances from the main point
center_point = points(1,:);
distances = sqrt(sum((points-repmat(center_point,size(points,1),1)).^2,2));
% Recursively finds the rest of the points
% NOTE: will fail with >500 points, as that is the recursion limit of 
% matlab. Shouldn't be a problem though.
centroids = [mean(points(distances<=radius,:),1);...
    get_centroids(points(distances>radius,:), radius)];

end

