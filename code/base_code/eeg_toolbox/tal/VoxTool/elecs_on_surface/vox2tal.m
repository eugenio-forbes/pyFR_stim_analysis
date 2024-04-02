function tal_coords = vox2tal(subj, CT_coords, outfile, reg_file)
%tal_coords = VOX2TAL(subj, [CT_coords, outfile, reg_file])
% converts coordinates from CT (Voxel) space to tkrRAS surface space on
% the average.
%
% subj ............... The subject being analyzed. Used to find location of
%                      CT_combined.img
% CT_coords .......... (Optional) The coordinates to be transformed, or the
%                      path to VOX_coords.txt. Default is subj's VOX_coords
% outfile ............ (Optional) If  exists, the coordinates will be
%                      written out to the file. Otherwise, nothing is
%                      written.
% reg_file ........... The file containing the registration matrix used to
%                      convert from tkr-ras to subject/avg RAS.
% -IMP (MAY 2014)
% iped@sas.upenn.edu

% Get the CT coordinates from the file if it is specified, otherwise get it
% from the default location
if ~exist('CT_coords','var') || isempty(CT_coords) || isa(CT_coords,'char')
    if exist('CT_coords','var') && isa(CT_coords,'char')
        fid = fopen(CT_coords);
    else
        fid = fopen(fullfile('/data/eeg/',subj,'tal/VOX_coords.txt'));
    end
    results = textscan(fid,'%d %f %f %f');
    CT_coords = double([results{2},results{3},results{4}]);
    numbers = results{1};
end
CT_coords_1 = double([CT_coords, ones(size(CT_coords, 1),1)]);

% Error checking...
if ~exist(reg_file,'file')
    error('vox2tal:reg_avg_missing','reg_avg.lta is missing for %s',subj);
end

% Hopefully the transformation matrix always starts on the 9th line.
% If this breaks, check with regular expressions (like below, with
% mri_info)
fid = fopen(reg_file);
for i=1:8
    fgetl(fid);
end
lines = fscanf(fid,'%f %f %f %f');
reg_avg = reshape(lines,4,4)';

% Freesurfer must be on the system path. Add it, if it's not.
currPath = regexp(getenv('PATH'),':','split');
if ~any(strcmp(currPath,'/usr/global/freesurfer/bin'))
    setenv('PATH',[getenv('PATH') ':/usr/global/freesurfer/bin']);
end

% run mri_info on the CT_combined. Get the location at which the
% transformation matrix starts
[~,results] = system(sprintf('mri_info --vox2ras-tkr %s',...
    fullfile('/data/eeg',subj,'/tal/images/combined',[subj '_CT_combined.img'])));
relevant_results_start = regexp(results,'\n(\s+-?[0-9]+\.[0-9]+)+','start');

% Sometimes it has to be run twice, for some strange reason...
if isempty(relevant_results_start)
    [~,results] = system(sprintf('mri_info --vox2ras-tkr %s',...
        fullfile('/data/eeg',subj,'/tal/images/combined',[subj '_CT_combined.img'])));
    relevant_results_start = regexp(results,'\n(\s+-?[0-9]+\.[0-9]+)+','start');
end

% Transformation matrix is always the last thing to be output. Read from the
% start to the end
Tmov = sscanf(results(relevant_results_start:end), '%f %f %f %f',[4,4])';

% error catching (AR)
if isempty(Tmov)
    Tmov = sscanf(results(450:end), '%f %f %f %f',[4,4])';
elseif strcmp(subj,'UP040')
    Tmov = sscanf(results(1:end), '%f %f %f %f',[4,4])';
end

% get the full transform
transform = reg_avg \ Tmov;

% The last column should be all 1s (or at least very close). If it's not,
% something went wrong
tal_coords = (transform * CT_coords_1')';
if ~all(tal_coords(:,4)-1<.001)
    warning('vox2tal:TransformFailed','SOMETHING SEEMS WRONG HERE.')
end

% Drop the last column...
tal_coords = tal_coords(:,1:3);

% If the outfile was specified, write to it
if exist('outfile','var')
    fid = fopen(outfile,'w');
    for elec_i=1:size(tal_coords,1)
        fprintf(fid, '%d %.2f %.2f %.2f\n', numbers(elec_i), ...
            tal_coords(elec_i,1), tal_coords(elec_i,2), tal_coords(elec_i,3));
    end
end
