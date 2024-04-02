function load_image( handles )
%LOAD_IMAGE a very short script to plot the pointcloud on the main axes.
axes(handles.mainAxes)
figure(handles.CT3D_fig);
clickA3DPoint([handles.xyz(:,1),  handles.xyz(:,2), handles.xyz(:,3)]', handles);
%cameratoolbar('Show');
end

