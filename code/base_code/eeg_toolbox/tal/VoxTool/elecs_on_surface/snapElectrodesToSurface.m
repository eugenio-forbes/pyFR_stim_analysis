function [ snapped_electrodeLocs, orig_electrodeLocs ] =...
    snapElectrodesToSurface(subj, electrodeLocs, isDepth, ...
       avgFlag, outfile, labels, doPlot, correction, auto_correct,...
       minDepth,maxSnap)
%[snapped_electrodeLocs, orig_electrodeLocs] = ...
%   SNAPELECTRODESTOSURFACE(subj, electrodeLocs, isDepth, avgFlag, outfile,
%       labels, doPlot, correction, auto_correct, minDepth, maxSnap)
%   snaps electrodes to the gyri of a surface
%
% INPUTS:
%   subj -------------- the subject, used to get surfaces (optional if avg)
%   electrodeLocs ----- EITHER:
%                       The path to a file with electrode locations in 
%                           format: label \t x \t y \t z \n
%                       OR:
%                       The electrode locations in a matrix, with no label
%   isDepth ----------- Boolean vector, true if that electrode is a depth
%   avgFlag ----------- True if the average surface info should be used
%   outfile ----------- (optional) The file to which to write the 
%                       information. If  no filename is given it will not
%                       write
%   labels ------------ (Only necessary if outfile is present and 
%                       electrodeLocs is not a path)
%                       A cell array (1xN) of the labels corresponding to
%                       each electrode
%   doPlot ------------ (Optional) Show diagnostic plots?
%   correction -------- (Optional) Shift the electrodes by this amount 
%                       before fitting
%   auto_correct ------ (Optional) After snapping, shift the electrodes by 
%                       the mean difference and snap again.
%   minDepth ---------- The minimum depth that the script should consider
%                       a gyrus. -1 = min sulcus, 1 = max gyrus.
%                       Defualts to 0.
%   maxSnap ----------- The maximimum distance an electrode will be snapped
%
% OUTPUTS:
%   snapped_electrodeLocs - the location of each of the electrodes once
%                           they have been snapped to the surface
%
%   orig_electrodeLocs    - the locations of each of the electrodes before
%                           being snapped to the surface
% created 10/13 by Isaac Pedisich (iped@sas.upenn.edu)
%
allSurf = [];
allCurv = [];

hold all;
if ~exist('doPlot','var') || isempty(doPlot)
    doPlot = true;
end
% read in one surface and sulc values
if avgFlag
    fsSubjDir = '/data/eeg/freesurfer/subjects/average';
else
    fsSubjDir = fullfile('/data/eeg/freesurfer/subjects',subj);
end

% Surface files
surfL = fullfile(fsSubjDir,'surf/lh.pial');
surfR = fullfile(fsSubjDir,'surf/rh.pial');
% Gyrus/sulcus files
sulcL = fullfile(fsSubjDir,'surf/lh.sulc');
sulcR = fullfile(fsSubjDir,'surf/rh.sulc');

% Read left surf and sulc files
[v_L, f_L] = read_surf(surfL);
[curv_L] = read_curv(sulcL);
allSurf = v_L;
allCurv = curv_L;
if doPlot
    f=figure;
    plotsurf_wrapper(v_L, f_L+1);
end

% Read right surf and sulc files
[v_R, f_R] = read_surf(surfR);
[curv_R] = read_curv(sulcR);
allSurf = [allSurf; v_R];
allCurv = [allCurv; curv_R];
if doPlot
    plotsurf_wrapper(v_R, f_R+1);
end

% read electrode locs and labels. Labels passed in are ignored/rewritten if
% electrodeLocs and lables are both passed in.
if isstr(electrodeLocs)
    fid = fopen(electrodeLocs);
    electrodeLocs = cell2mat(textscan(fid, '%*s\t%f\t%f\t%f\n'));
    frewind(fid);
    labels = textscan(fid, '%s\t%*f\t%*f\t%*f\n');
    labels = labels{1};
end

% If correction exists, remove that amount from electrodeLocs
if exist('correction','var') && ~isempty(correction)
    electrodeLocs = electrodeLocs...
        -repmat(correction,size(electrodeLocs,1),1);
end

% Get only the outermost part of gyri
if ~exist('minDepth','var') || isempty(minDepth)
    minDepth = 0;
end
gyri = allSurf(allCurv<-minDepth,:);

% Snap electrodes to surface
snapped_electrodeLocs = electrodeLocs;
snapped_electrodeLocs_inds = dsearchn(double(gyri), double(electrodeLocs(~isDepth,:)));
snapped_electrodeLocs(~isDepth,:) = gyri(snapped_electrodeLocs_inds,:);

% If auto_correct
if exist('auto_correct','var') && ~isempty(auto_correct) && auto_correct
    % subtract the mean difference from the electrodelocs to the snapped
    difference = mean(electrodeLocs-snapped_electrodeLocs,1);
    new_electrodeLocs = electrodeLocs-...
        repmat(difference,size(electrodeLocs,1),1);
    % resnap
    snapped_electrodeLocs = dsearchn(double(gyri),double(new_electrodeLocs));
    snapped_electrodeLocs = gyri(snapped_electrodeLocs,:);
end

% plot
if doPlot
    plot3_wrapper(snapped_electrodeLocs, 20, 'red');
%    % plot the electrode names
%     if exist('labels','var') && iscell(labels)
%         for i=1:length(labels)
%                    text(snapped_electrodeLocs(i,1),...
%                        snapped_electrodeLocs(i,2),...
%                        snapped_electrodeLocs(i,3),...
%                        labels{i});
%         end
%     end
end

orig_electrodeLocs = electrodeLocs;

% calculate distances
D= pdist2(snapped_electrodeLocs,orig_electrodeLocs);
y = eye(size(D));y = logical(y);
distances = D(y);

% revert the over-snapped electrodes to their original coords
if exist('maxSnap','var') && ~isempty(maxSnap)
    snapped_electrodeLocs(distances>=maxSnap,:) = orig_electrodeLocs(distances>=maxSnap,:);
end

% If there is an outfile, write the electrode locations to it.
if exist('outfile','var') && any(outfile~=false)
    fid = fopen(outfile,'w');
    for i=1:length(labels)
        fprintf(fid, '%s\t%f\t%f\t%f\n',labels{i},...
            snapped_electrodeLocs(i,1),...
            snapped_electrodeLocs(i,2),...
            snapped_electrodeLocs(i,3));
    end
    fclose(fid);
end


