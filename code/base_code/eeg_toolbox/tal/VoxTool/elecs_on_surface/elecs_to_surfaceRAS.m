function [elec_ras, v1, f1, v2, f2] =...
    elecs_to_surfaceRAS(path_vox_coords_mother, ...
    path_ct_combined, path_reg, path_vox_surfRAS, path_surface1, path_surface2)
%[elec_ras, v1, f1, v2, f2] = ELECS_TO_SURFACERAS(path_vox_coords_mother,
%       path_ct_combined, path_reg, out_file, path_surface1, path_surface2) 
%
%   converts a set of coordinates, given by a vox_coords_mother file, 
%   to a set of points in surfaceRAS space
%
%   path_vox_coords_mother - The path to the vox_coords_mother file
%                            containing the points
%   path_ct_combined       - The path to the CT_combined file from which
%                            vox_coords_mother was made
%   path_ct2mri            - The path to the *REWRITTEN* CT2MRI file, without the
%                            extension (it references both .m and .img
%   path_reg               - The path to the reg.lta file created with
%                            tkregister2
%   path_vox_surfRAS       - Optional - The path to the file where the
%                            results should be written. 
%                            False or no input -> no file
%   path_surface(1/2)      - Optional - the path to the surface files for
%                            plotting
%
%                       ----------- OR -----------
%[elec_ras, v1, f1, v2, f2] = ELECS_TO_SURFACERAS(subject)
%
%   performs same operations as above, but fills in filenames with default
%   values.
%
% --------------------------------------
%
% RETURNS: [elec_ras_centroids, vertices1, faces1, vertices2, faces2]
%   elec_ras_centroids    -   The centroids of all of the electrodes in
%                   surfaceRAS space
%   v1/2 -   Only returned if surfaces are passed in:
%                   the location of the vertices for each point on the
%                   surface, in surfaceRAS space
%   f1/2    -   Only returned if surfaces are passed in:
%                   Which vertices connect to one another on the surface
%
%   NOTE: plot elec_ras with plot3_wrapper(elec_ras, size, 'color')
%         plot surfaces with plotsurf_wrapper(vertices1, faces1)
% Written by Isaac Pedisich <iped@sas.upenn.edu> 2/2014
if nargin==1
    subject = path_vox_coords_mother;
    path_vox_coords_mother = fullfile('/data/eeg',subject,'tal/VOX_coords_mother.txt');
    path_ct_combined = fullfile('/data/eeg',subject,...
        ['tal/images/combined/',subject,'_CT_combined.img']);
    path_reg = fullfile('/data/eeg',subject,'tal/images/combined/reg.lta');
    surf_path = fullfile('/data/eeg/freesurfer/subjects',subject,'surf/');
    path_vox_surfRAS = fullfile('/data/eeg',subject,'tal/VOX_coords_mother_surfRAS.txt');
    path_surface1 = fullfile(surf_path,'lh.pial');
    path_surface2 = fullfile(surf_path,'rh.pial');
end

%Open vox_coords_mother
fid = fopen(path_vox_coords_mother);

% Get the coordinates
vox_coords = cell2mat(textscan(fid, '%*s\t%f\t%f\t%f\t%*s\t%*s %*s\n'));
fclose(fid);

% Append 1s to coordinates so it can be affine-transformed
elec_xyz = [vox_coords, ones(length(vox_coords),1)];

% Next lines get the vox2ras_tkr (surfaceRAS?) transformation by executing
% mri_info from the system, then parsing the output
[~,sysOutput] = system(sprintf('mri_info --vox2ras-tkr %s',path_ct_combined));
sysOutput = regexp(sysOutput, '(\n)|(\r)','split');
vox2ras_tkr_matrix = sysOutput(2:5);
vox2ras_tkr_matrix = cellfun(@(x)str2double(regexp(x,'\s+','split')), ...
    vox2ras_tkr_matrix, 'uniformOutput',false);
vox2ras_tkr_matrix = cell2mat(vox2ras_tkr_matrix');
vox2ras_tkr_matrix = vox2ras_tkr_matrix(:,2:5);

% Next lines get the transformation matrix from reg.lta
fid = fopen(path_reg);
for i=1:8
    fgetl(fid);
end
reg_mat = fscanf(fid,'%f',[4,4])';
fclose(fid);

% actually apply the transform! Works! Yay!
elec_ras = reg_mat\vox2ras_tkr_matrix*elec_xyz';
elec_ras = elec_ras(1:3,:)';
disp('plotting')
figure
plot3_wrapper(elec_ras,30,'r');
hold all

if exist('path_surface1','var')
    [v1, f1] = read_surf(path_surface1);
    plotsurf_wrapper(v1, f1+1);
end
if exist('path_surface2','var')
    [v2, f2] = read_surf(path_surface2);
    plotsurf_wrapper(v2, f2+1);
end

% Write out VOX_coords_mother again, but with the new coords, if necessary
if exist('path_vox_surfRAS','var') && ~any(path_vox_surfRAS==false)
    disp('writing file')
    fid = fopen(path_vox_coords_mother);
    vox_coords = textscan(fid, '%s\t%f\t%f\t%f\t%s\t%s %s\n');
    fclose(fid);
    fid = fopen(path_vox_surfRAS,'w');
    for i=1:length(vox_coords{1})
        fprintf(fid,'%s\t%f\t%f\t%f\t%s\t%s %s\n',...
            vox_coords{1}{i},...
            elec_ras(i,1), elec_ras(i,2), elec_ras(i,3), ...
            vox_coords{5}{i},...
            vox_coords{6}{i},...
            vox_coords{7}{i});
    end
    fclose(fid);
end 