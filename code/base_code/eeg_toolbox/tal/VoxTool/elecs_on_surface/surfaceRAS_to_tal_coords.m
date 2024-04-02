function [ tal_coords ] = surfaceRAS_to_tal_coords( subj, RAS_coords )
%SURFACERAS_TO_TAL_COORDS(subj) converts a subject's surfaceRAS
% coordinates, as read from VOX_coords_mother_surfRAS.txt, 
% into talairach coordinates

[talXFM, Norig, Torig] = read_transforms(subj);

if ~exist('RAS_coords','var') || isempty(RAS_coords)
    fid = fopen(fullfile('/data/eeg/',subj,'tal/VOX_coords_mother_surfRAS.txt'));
    results = fscanf(fid,'%*s\t%f\t%f\t%f\t%*s\t%*d %*d',inf);
    RAS_coords = reshape(results, 3, length(results)/3)';
end

RAS_coords_1 = [RAS_coords, ones(size(RAS_coords, 1),1)];
tal_coords = (talXFM*Norig*inv(Torig)*RAS_coords_1')';
