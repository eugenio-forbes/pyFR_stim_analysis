function [centroids, labels, vertices1, faces1, vertices2, faces2 ] = ...
    elecs_with_labels_to_surfaceRAS(path_vox_coords_mother,...
    path_ct_combined, path_ct2mri, path_surface1, path_surface2, ...
    outfile)
%[centroids, lables, vertices1, faces1, vertices2, faces2] = 
% ELECS_WITH_LABELS_TO_SURFACERAS(path_vox_coords_mother, path_ct_combined,
%       path_ct2mri, path_surface1, path_surface2, outfile)
%
%   converts a set of coordinates, given by a vox_coords_mother file,
%   to a set of points in surfaceRAS space. Creates separate output files 
%   for each electrode so that they can be identified later. Puts all 
%   electrodes in 'elec' folder, contained in same folder as ct_combined.
%
% INPUTS:
%   path_vox_coords_mother - The path to the vox_coords_mother file
%                            containing the points
%   path_ct_combined       - The path to the CT_combined file from which
%                            vox_coords_mother was made
%   path_ct2mri            - The path to the *REWRITTEN* CT2MRI file, without the
%                            extension (it references both .m and .img
%   path_surface(1/2)      - Optional. The path to the surface files for
%                            plotting
%   outfile                - Optional. If given, writes all electrode and
%                            label information to this location.
%
% RETURNS:
%   centroids   -   The centroids corresponding to each of the electrodes 
%                   in surfaceRAS space
%   labels      -   A cell array of the labels corresponding to each of the
%                   electrodes returned in centroids
%   vertices1/2 -   Only returned if surfaces are passed in:
%                   the location of the vertices for each point on the
%                   surface, in surfaceRAS space
%   faces1/2    -   Only returned if surfaces are passed in:
%                   Which vertices connect to one another on the surface
%
%   NOTE: plot centroids with plot3_wrapper(elec_ras, size, 'color')
%         plot surfaces with plotsurf_wrapper(vertices1, faces1)
%
% Created 10/13 by Isaac Pedisich (iped@sas.upenn.edu)
% Updates:
% AGR 11/19 -- if no outfile, automatically writes text output file to
% 'elecs/[subj]_elecLocs_surfRAS.txt' 
            
%Open vox_coords_mother
fid = fopen(path_vox_coords_mother);

% Get the coordinates
vox_coords = cell2mat(textscan(fid, '%*s\t%d\t%d\t%d\t%*s\t%*s %*s\n'));
frewind(fid)
labels = textscan(fid, '%s\t%*d\t%*d\t%*d\t%*s\t%*s %*s\n');
labels = labels{1};
disp('reading CT combined')
ct_combined = MRIread(path_ct_combined);
disp('thresholding')
[x,y,z] = ind2sub(size(ct_combined.vol), ...
    find(ct_combined.vol==max(max(max(ct_combined.vol)))));

% Subtrace 1 from vertices b/c zero-indexed
ct_combined_xyz = [x,y,z]-1;

[folder,filename,~] = fileparts(path_ct_combined);
subj_name = filename(1:regexp(filename,'_')-1);

elec_folder = [folder, filesep, 'elecs'];
mkdir(elec_folder);
mkdir(fullfile(elec_folder,'CT'));
mkdir(fullfile(elec_folder,'MRI'));
centroids = nan(size(vox_coords,1),3);
disp('done thresholding');
parfor i=1:size(vox_coords,1)
    fprintf(['%d/%d ',labels{i},':'],i,size(vox_coords,1))
    coord=vox_coords(i,:);
    % Make the electrodes_out file for each vox_coord individually.
    % Get the centroid for that vox_coord
    centroids(i,:) = elec_to_surfaceRAS(coord, ct_combined_xyz, ct_combined, ...
        elec_folder, subj_name, path_ct2mri, labels{i});
end

disp('plotting')
figure

% Plot the centroids
plot3_wrapper(centroids,30,'r');
hold all

% plot the labels next to the centroids
for i=1:length(labels)
    text(centroids(i,1), centroids(i,2), centroids(i,3), labels{i});
end

% Plot the surface if it is passed in
if exist('path_surface1','var') && any(path_surface1~=false)
    [vertices1, faces1] = read_surf(path_surface1);
    plotsurf_wrapper(vertices1, faces1+1);
end
if exist('path_surface2','var') && any(path_surface1~=false)
    [vertices2, faces2] = read_surf(path_surface2);
    plotsurf_wrapper(vertices2, faces2+1);
end
%parse outfile
if ~exist('outfile','var') || isempty(outfile)
    outfile = fullfile(folder,'elecs',[subj_name '_elecLocs_surfRAS.txt']);
end
% Create the output file
if exist('outfile','var')
    fid = fopen(outfile,'w');
    for i=1:length(labels)
        fprintf(fid, '%s\t%f\t%f\t%f\n',labels{i},...
            centroids(i,1),centroids(i,2),centroids(i,3));
    end
end

function [centroid] = elec_to_surfaceRAS(elec_coord, ...
    ct_combined_xyz, ct_combined, folder, subj_name, path_ct2mri, label)
%ELECS_TO_SURFACE_RAS converts a set of coordinates, given by a
%   vox_coords_mother file, to a set of points in surfaceRAS space
%
%   path_vox_coords_mother - The path to the vox_coords_mother file
%                            containing the points
%   path_ct_combined       - The path to the CT_combined file from which
%                            vox_coords_mother was made
%   path_ct2mri            - The path to the CT2MRI file, without the
%                            extension (it references both .m and .img
%   path_surface(1/2)      - Optional - the path to the surface files for
%                            plotting
%
% RETURNS:
%   elec_ras    -   The locations of all of the electrodes in surfaceRAS
%                   space
%   vertices1/2 -   Only returned if surfaces are passed in:
%                   the location of the vertices for each point on the
%                   surface, in surfaceRAS space
%   faces1/2    -   Only returned if surfaces are passed in:
%                   Which vertices connect to one another on the surface
%
%   NOTE: plot elec_ras with plot3_wrapper(elec_ras, size, 'color')
%         plot surfaces with plotsurf_wrapper(vertices1, faces1)


new_mri_struct = ct_combined;

% Prepare an empty volume to fill in 
max_val = max(max(max(new_mri_struct.vol)));
new_mri_struct.vol = zeros(size(new_mri_struct.vol));

fprintf('clustering... ');
% Get all of the points in the cluster
points = getCluster(ct_combined_xyz, 5, elec_coord([2,1,3]), 3);


% Get mark all of those points on the volume file
for point=points'
    new_mri_struct.vol(point(1),point(2),point(3))=max_val;
end

elec_CT_loc = [folder, filesep, 'CT',filesep,subj_name, '_',label,'_CT.nii.gz'];
elec_MRI_loc = [folder, filesep, 'MRI',filesep,subj_name, '_',label,'_MRI.nii.gz'];

fprintf('writing... ')
MRIwrite(new_mri_struct, elec_CT_loc);

fprintf('warping... ')
%system(['bash run_applywarp.sh ' elec_CT_loc ' ' elec_MRI_loc ' ' path_ct2mri]);
system(sprintf('applywarp -i %s -o %s -r %s.nii.gz --interp=nn --premat=%s.mat', elec_CT_loc, elec_MRI_loc, path_ct2mri, path_ct2mri));

elec_mri = MRIread(elec_MRI_loc);

fprintf('thresholding... ')
[x,y,z] = ind2sub(size(elec_mri.vol),...
    find(elec_mri.vol==max(max(max(elec_mri.vol)))));
elec_xyz = [y x z, ones(length(x),1)];
elec_ras = elec_mri.vox2ras*elec_xyz';
%elec_ras(3,:) = elec_ras(3,:)*-1; % This is to flip it in the x-coord
elec_ras = elec_ras(1:3,:)' - repmat([-3.01 13.01 7.5], length(elec_ras),1);
centroid = get_centroids(elec_ras);
fprintf('done. \n\r');
