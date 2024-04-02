function [ p ] = plot3_wrapper( x , size, color, marker)
% [p] = PLOT3_WRAPPER(x, size, color, marker): a wrapper for the plot3 function to 
% make it easy to pass in a 3 x N matrix of values 
% and plot them all.
%
% Accepts:
%   x: Nx3 matrix of the points to plot
%   size: defaults to 10
%   color: defaults to black
%   marker: defaults to '.'
%
% Returns 
%   p: the handle to the plot

if ~exist('size','var')
    size = 10;
end
if ~exist('color','var')
    color = 'k';
end
if ~exist('marker','var')
    marker = '.';
end
p = plot3(x(:,1), x(:,2),x(:,3),'.','markersize',size,'color',color,'marker',marker);

end

