function [xyz_mov, xyz_orig, xyz_snap, dOrig, dSnap, dSpringPre, dSpringPost] = ...
    energyMinimization( subj, source_subj, electrodeLocs,...
    labels,  isDepth, avgFlag, doPlot, outfile)
% [xyz_mov, xyz_orig, xyz_snap, dOrig, dSnap, dSpringPre, dSpringPost] = 
% ENERGYMINIMIZATION( subj, electrodeLocs, labels, isDepth,...
%          avgFlag, doPlot, outfile)
%   subj           -    
%   electrodeLocs  - Either the path to electrodeLocs or the 3xN matrix
%   labels - I think, always 1:length(electrodeLocs). Used in conjunction
%            with bpTalStruct to determine strips. Not necessary is
%            electrodeLocs is a string
%   isDepth=true are excluded from snapping
%   avgFlag=true uses average surface, otherwise uses indiv.
%   doPlot=true will show animation of snapping
%   outfile existing will write output to file
iterations = 101;

initSnap = true;
snapMod = false;
finalSnap = true;

if ~exist('source_subj','var') || isempty(source_subj)
    source_subj = subj;
end
if ~exist('origWeight','var') || isempty(origWeight )
    origWeight=.1;
end
if ~exist('neighWeight','var') || isempty(neighWeight)
    neighWeight=.2;
end
if ~exist('surfWeight','var') || isempty(surfWeight)
    surfWeight = .2;
end
if ~exist('doPlot','var') || isempty(doPlot)
    doPlot = true;
end
    
if isstr(electrodeLocs)
    fid = fopen(electrodeLocs);
    electrodeLocs = cell2mat(textscan(fid, '%*s\t%f\t%f\t%f\n'));
    %electrodeLocs = electrodeLocs(:,[2,1,3]);
    frewind(fid);
    labels = textscan(fid, '%s\t%*f\t%*f\t%*f\n');
    labels = cellfun(@(x)str2double(x),labels{1});
else
    electrodeLocs = double(electrodeLocs);
end

if ~exist('isDepth') || isempty(isDepth)
    isDepth = false(size(electrodeLocs,1),1);
end
if exist(fullfile('/data/eeg',subj,'tal/bpPairs.mat'),'file')
    load(fullfile('/data/eeg/',subj,'tal/bpPairs.mat'));
else
    allpairs = getBipolarPairs(subj);
end

if avgFlag
    fsSubjDir = '/data/eeg/freesurfer/subjects/average';
else
    fsSubjDir = fullfile('/data/eeg/freesurfer/subjects',source_subj);
end

[vL, fL] = read_surf(fullfile(fsSubjDir,'surf/lh.pial'));
sulcL = read_curv(fullfile(fsSubjDir,'surf/lh.sulc'));
[vR, fR] = read_surf(fullfile(fsSubjDir,'surf/rh.pial'));
sulcR = read_curv(fullfile(fsSubjDir,'surf/rh.sulc'));
vAll = [vL;vR];

good_v_mask = false(length(vAll),1);
for i=1:size(electrodeLocs)
    this_xyz = electrodeLocs(i,:);
    distances = sqrt(sum((vAll-repmat(this_xyz,length(vAll),1)).^2,2));
    good_v_mask = good_v_mask | distances<30;
end

if nnz(good_v_mask)==0
    error('energyMinimization:elecsOffBrain','Electrodes appear to be far off brain');
end

curvAll = [sulcL; sulcR];
gyri = vAll(good_v_mask & curvAll<0,:);

xyz_all = electrodeLocs;

xyz_mov = xyz_all(~isDepth,:);
labels_mov = labels(~isDepth);

if doPlot
    plot3_wrapper(xyz_mov,20,'blue');
end
dOrig = getDistances(allpairs, labels_mov, xyz_mov);

hold off

plot3_wrapper(xyz_mov,20,'blue');
dSnap = getDistances(allpairs, labels_mov, gyri(dsearchn(gyri, xyz_mov),:));
xyz_orig = xyz_mov;
springOrig = @(xyz_current)(xyz_mov-xyz_current)*origWeight;
xyz_snap = gyri(dsearchn(gyri, xyz_mov),:);
closest_indices = dsearchn(gyri, xyz_mov);

plot3_wrapper(xyz_mov,20,'blue');
hold off;
drawnow;

xyz_last_snap = xyz_mov;
if initSnap
    xyz_mov = gyri(dsearchn(gyri, xyz_mov),:);
    closest_indices = dsearchn(gyri, xyz_mov);
end
% BULK OF THE WORK...
for i=1:iterations
        xyz_old = xyz_mov;
    if snapMod~=false && mod(i, snapMod)==0;
        xyz_mov = gyri(dsearchn(gyri, xyz_mov),:);
        closest_indices = dsearchn(gyri, xyz_mov);
        if all(all(xyz_mov==xyz_last_snap))
            break
        end
        xyz_last_snap = xyz_mov;
        
    end

     xyz_new = xyz_mov+springNeighbor(allpairs, labels_mov, xyz_mov, neighWeight) + ...
         springSurface(gyri, xyz_mov, closest_indices, surfWeight) ...
         + springOrig(xyz_mov) ...
         ;

    if all(xyz_new==xyz_old)
       break
    end
    xyz_mov = xyz_new;
    if mod(i,2)==0 && doPlot
        hold off;
        xyz_all(~isDepth,:) = xyz_mov;
        if doPlot
            h = plot3_wrapper(xyz_all,20,'blue');
            drawnow;
        end
    end
    

end
dSpringPre = getDistances(allpairs, labels_mov, xyz_mov);
if finalSnap
    xyz_mov = gyri(dsearchn(gyri, xyz_mov),:);
end
    hold off
    if doPlot
        plot3_wrapper(xyz_mov, 20,'blue');
        hold all
        %plotsurf_wrapper(vR, fR+1);
        plot3_wrapper(xyz_orig, 10, 'red');
    end

dSpringPost = getDistances(allpairs, labels_mov, xyz_mov);


if exist('outfile','var')
    xyz_all(~isDepth,:) = xyz_mov;
    fid = fopen(outfile,'w');
    
    for i=1:size(xyz_all,1)
        fprintf(fid, '%d %.3f %.3f %.3f\n',labels(i), xyz_all(i,1), xyz_all(i,2), xyz_all(i,3));
    end
    fclose(fid);
end


function spring = springNeighbor(bp_manual, elecNums, xyz, neighWeight, eNums)

if ~exist('eNums','var')
    eNums = 1:size(xyz, 1);
end
spring = nan(length(eNums),3);
bp_manual = bp_manual';
for i=1:length(eNums)
    eNum = eNums(i);
    [I,~] = find(bp_manual==elecNums(i));
    index = find(bp_manual==elecNums(i));
    index(I==1) = index(I==1)+1;
    index(I==2) = index(I==2)-1;
    relevant_xyz = xyz(ismember(elecNums, bp_manual(index)),:);
    rep_xyz = repmat(xyz(eNum,:),size(relevant_xyz,1),1);
    xyz_diff = relevant_xyz-rep_xyz;
    distances = sqrt(sum(xyz_diff.^2 ,2));
    forces =  (distances-10)*neighWeight;
    spring(i,:) = sum(repmat(forces,1,3).*(xyz_diff./repmat(distances,1,3)),1);
    spring(isnan(spring)) = 0;
end

function spring = springSurface(vAll, xyz, closest_indices, surfWeight)
xyz_diff = vAll(closest_indices,:)-xyz;
distances = sqrt(sum(xyz_diff.^2,2));
forces = ((distances))*surfWeight;
spring = repmat(forces,1,3).*(xyz_diff./repmat(distances,1,3));
spring(isnan(spring)) = 0;

function distances = getDistances(bp_manual, elecNums, xyz)
distances = [];
bp_manual = bp_manual';
for i=1:size(elecNums,1)
    [I,~] = find(bp_manual==elecNums(i));
    index = find(bp_manual==elecNums(i));
    index(I==1) = index(I==1)+1;
    index(I==2) = index(I==2)-1;
    relevant_xyz = xyz(ismember(elecNums, bp_manual(index)),:);
    rep_xyz = repmat(xyz(i,:),size(relevant_xyz,1),1);
    xyz_diff = relevant_xyz-rep_xyz;
    this_distances = sqrt(sum(xyz_diff.^2,2));
    distances(end+1:end+size(xyz_diff,1),:) = this_distances;
    if any(abs(this_distances-10)>10)
%          for j=find(abs(this_distances-10)>10)
%              hold all
%              plot3([relevant_xyz(j,1), xyz(i,1)],...
%                  [relevant_xyz(j,2), xyz(i,2)],...
%                  [relevant_xyz(j,3), xyz(i,3)], 'ro-')
%          end
%          drawnow;
    end
    if any(distances>20)
 %       disp('bad');
    end
end
