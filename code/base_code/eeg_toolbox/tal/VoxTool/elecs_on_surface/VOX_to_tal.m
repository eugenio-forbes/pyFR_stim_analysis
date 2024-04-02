function tal_coords = VOX_to_tal(subj, CT_coords)
%tal_coords = VOX_TO_TAL(subj, [CT_coords])
% converts coordinates from CT (Voxel) space to tkrRAS surface space on
% the average.
% CT_coords is optional, will be read from VOX_coords_mother if not passed
% in.
if ~exist('CT_coords','var') || isempty(CT_coords)
    fid = fopen(fullfile('/data/eeg/',subj,'tal/VOX_coords_mother.txt'));
    results = fscanf(fid,'%*s\t%f\t%f\t%f\t%*s\t%*d %*d',inf);
    CT_coords = reshape(results, 3, length(results)/3)';
end

CT_coords_1 = [CT_coords, ones(size(CT_coords, 1),1)];

reg_avg_file = fullfile('/data/eeg',subj,'tal/images/combined/reg_avg.lta');
if ~exist(reg_avg_file,'file')
    error('VOX_to_tal:reg_avg_missing','reg_avg.lta is missing for %s',subj);
end

fid = fopen(reg_avg_file);
for i=1:8
    fgetl(fid);
end
lines = fscanf(fid,'%f %f %f %f');
reg_avg = reshape(lines,4,4)';


[~,results] = system(sprintf('mri_info --vox2ras-tkr %s',...
    fullfile('/data/eeg',subj,'/tal/images/combined',[subj '_CT_combined.img'])));
Tmov = sscanf(results(28:end), '%f %f %f %f',[4,4])';

tal_coords = (inv(reg_avg) * Tmov * CT_coords_1')';
tal_coords = tal_coords(:,1:3);

