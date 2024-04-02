function [ TalXFM, Norig, Torig ] = read_transforms( subject )
%READ_TALXFM Summary of this function goes here
%   Detailed explanation goes here
fid = fopen(fullfile('/data/eeg/freesurfer/subjects',subject,'mri/transforms/talairach.xfm'));
fgetl(fid);fgetl(fid);fgetl(fid);fgetl(fid);fgetl(fid);
TalXFM = reshape(fscanf(fid,' %f %f %f %f'),4,3)';


[~,results] = system(sprintf('mri_info --vox2ras %s',...
    fullfile('/data/eeg/freesurfer/subjects',subject,'/mri/orig.mgz')));
Norig = sscanf(results, '%f %f %f %f',[4,4])';

[~,results] = system(sprintf('mri_info --vox2ras-tkr %s',...
    fullfile('/data/eeg/freesurfer/subjects',subject,'/mri/orig.mgz')));
Torig = sscanf(results, '%f %f %f %f',[4,4])';
end

