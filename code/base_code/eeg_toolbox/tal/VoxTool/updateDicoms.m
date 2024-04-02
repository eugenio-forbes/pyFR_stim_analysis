function updateDicoms(handles,  xyz , changed)
%UPDATEDICOMS updates the dicoms image with the correct image, and moves
%the circle to the right location
%   xyz - the x,y,z coordinates of the point to be circled
%   changed - whether or not the image needs to be changed (whether the Z
%       component of xyz is different than previously. TODO: could 
%       probably be checked dynamically within this function

% If the new dicom needs to be loaded, load it.
set(0,'currentfigure',handles.sliceFig);
try
set(handles.sliceFig,'currentaxes',handles.sliceZAxes)
imagesc(flipud(squeeze(handles.vol(:,:,xyz(3)+1))));
colormap(bone);
catch e
    disp(e.message)
end
try
set(handles.sliceFig,'currentaxes',handles.sliceYAxes)
imagesc(rot90(squeeze(handles.vol(:,xyz(2),:))));
colormap(bone);
catch e
    disp(e.message)
end
try
set(handles.sliceFig,'currentaxes',handles.sliceXAxes)
imagesc(rot90(squeeze(handles.vol(xyz(1),:,:))));
colormap(bone)
catch e
    disp(e.message)
end

% Draw the circle around the correct x,y point on the dicom.
drawDicomCircle(handles,[xyz], 1, 'red');