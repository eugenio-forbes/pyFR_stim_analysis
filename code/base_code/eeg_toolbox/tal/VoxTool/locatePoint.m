function locatePoint( handles, xyz, dicomChanged ,color, tag)
%LOCATEPOINT Plots a circle on the CT and Dicoms, deleting any points that
%have previously been plotted with the same tag

%Delete the previously plotted point with the same tag
h = findobj(handles.mainAxes,'Tag',tag);
if ~isempty(h)
    delete(h)
end

% Plot the cirle on the CT
h = plot3(xyz(:,1),xyz(:,2), xyz(:,3), 'ro', 'MarkerSize', 20,'linewidth',2,'parent',handles.mainAxes);
set(h,'color',color);
set(h,'Tag',tag); % set its Tag property for later use   

%Plot the circle on the dicoms
updateDicoms(handles, xyz(1,:), dicomChanged);

end

