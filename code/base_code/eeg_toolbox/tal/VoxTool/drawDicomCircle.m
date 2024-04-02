function drawDicomCircle( handles,xyz, circleNum, color)
%DRAWDICOMCIRCLE draws the circle on the dicom figure
%   handles - the handles to objects in all figures
%   xy - the x & y coordinates to be circled (in dicom coordinates)
%   circlxeNum - not used, always 1 (Potentially useful for circling more
%                than one item.
%   color - the color of the circle (always red, i believe)

%get the position of the point on the dicom in figure coordinates
% Xylim=ylim(handles.sliceXAxes);
% Yylim=ylim(handles.sliceYAxes);
% Zylim=ylim(handles.sliceZAxes);
set(0,'currentfigure',handles.sliceFig);
figposX = dsxy2figxy(handles.sliceXAxes,[xyz(2)-6, xyz(3)-7, 12, 12]);
figposY = dsxy2figxy(handles.sliceYAxes,[xyz(1)-7,xyz(3)-6, 12, 12]);
figposZ = dsxy2figxy(handles.sliceZAxes,[xyz(2)-6, xyz(1)-7, 12, 12]);
%If it is a new circle, plot it - otherwise just set its position to the
%new location
if circleNum>length(handles.sliceXCircleHandles)
    handles.sliceXCircleHandles(circleNum) = annotation('ellipse',figposX, 'color',color,'linewidth',2);
    handles.sliceYCircleHandles(circleNum) = annotation('ellipse',figposY, 'color',color,'linewidth',2);
    handles.sliceZCircleHandles(circleNum) = annotation('ellipse',figposZ, 'color',color,'linewidth',2);

else
    set(handles.sliceXCircleHandles(circleNum), 'position', figposX, 'color',color,'linewidth',1);
    set(handles.sliceYCircleHandles(circleNum), 'position', figposY, 'color',color,'linewidth',1);
    set(handles.sliceZCircleHandles(circleNum), 'position', figposZ, 'color',color,'linewidth',1);
end
%set dicom image on bottom so the circle can be seen
uistack(handles.sliceXAxes,'bottom');
uistack(handles.sliceYAxes,'bottom');
uistack(handles.sliceZAxes,'bottom');
%Switch the active figure back to the CT
set(0,'currentfigure', handles.CT3D_fig);


