% This is a wrapper to be run after VoxTool has been run and after the
% surface reconstruction is completed
function elecs_to_surf_wrapper(subj)
%This function wraps around elecs_to_surfaceRAS.m and converts vox coords
%to surface RAS space
% 

%subj = 'TJ061';
path_vox_coords_mother = ['/data/eeg/' subj '/tal/VOX_coords_mother.txt'];
path_vox_surfRAS = ['/data/eeg/' subj '/tal/VOX_coords_mother_surfRAS.txt'];
imageDir = ['/data/eeg/' subj '/tal/images/combined'];
surfDir = ['/data/eeg/freesurfer/subjects/' subj];
path_ct_combined = fullfile(imageDir,[subj '_CT_combined.img']);
path_reg = fullfile(imageDir,['reg.lta']);
path_surface1 = fullfile(surfDir,'surf','lh.pial');
path_surface2 = fullfile(surfDir,'surf','rh.pial');

% run function
[elec_ras, v1, f1, v2, f2] = elecs_to_surfaceRAS(path_vox_coords_mother, ...
    path_ct_combined, path_reg, path_vox_surfRAS, path_surface1, path_surface2);
