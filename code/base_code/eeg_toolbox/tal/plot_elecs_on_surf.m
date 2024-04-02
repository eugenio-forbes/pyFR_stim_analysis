function [talStruct] = plot_elecs_on_surf(subj)
%This function reads the XYZ coords associated with each contact (no
%bipolar). Each electrode group is associated with a different color. 

%Inputs:
%subj ... 'TJ075'

%Outputs:
%talStruct ... updated tal structure with surfRAS and surfRAS_snapped
                %fields
% To Do:

%Written by Ashwin Ramayya (2/19/2014)
                
%set dirs and paths
docDir  = fullfile('/data/eeg/', subj,'/docs');
talDir = ['/data/eeg/' subj '/tal/'];
surfVOX = fullfile(talDir,'VOX_coords_mother_surfRAS.txt');
freesurferDir = fullfile('/data/eeg/freesurfer/subjects',subj);
surfDir = fullfile(freesurferDir,'surf');
surfL = fullfile(surfDir,'lh.pial');
sulcL = fullfile(surfDir,'lh.sulc');
surfR = fullfile(surfDir,'rh.pial');
sulcR = fullfile(surfDir,'rh.sulc');
jackFileName      = 'jacksheet.txt';
jackFileName_2 = 'jack_sheet.txt';
%error check
if ~exist(surfDir,'dir') || ~exist(surfL,'file') || ~exist(surfR,'file')
    error('Surface Reconstruction (recon all) not done')
end
if ~exist(surfVOX,'file') 
    error('elecs_to_surfRAS.m not done')
end

% get surfRAS electrode locations
fid = fopen(surfVOX);
C = textscan(fid,'%s%d%d%d%s%d%d');
centroids = [C{:,2} C{:,3} C{:,4}]; centroids = double(centroids);
labels = C{:,1};

% get elec groups
for k=1:length(labels)
    thisElecNam     = labels{k};
    name_stripped{k}  = thisElecNam(regexp(thisElecNam,'\D'));
    left_ind(k) = isempty(regexp(thisElecNam,'R'));
end

% snap params 
correction = [];
auto_correct = 1;
minDepth = 0;
radius = 5;

% process left hemisphere: 
% Snap left elecs
if sum(left_ind)>0
[l_snap,~] = snapElectrodesToSurface(centroids(left_ind,:),...
    surfL, sulcL, [], [], [], labels(left_ind), correction, auto_correct,...
    minDepth); view([-90 0]);close
end
if sum(~left_ind)>0
[r_snap,~] = snapElectrodesToSurface(centroids(~left_ind,:),...
    surfR, sulcR, [], [], [], labels(~left_ind), correction, auto_correct,...
    minDepth); close
end
%num elec groups
cmap = colormap(jet(10000));
l_elecGroups = unique(name_stripped(left_ind));
r_elecGroups = unique(name_stripped(~left_ind));
l_idx = randi([1 10000],[1 length(l_elecGroups)]);
r_idx = randi([1 10000],[1 length(r_elecGroups)]);


% plot left elecs
h(1) = swagFig; hold all 
[v,f] =read_surf_wrapper(surfL);plotsurf_wrapper(v,f);view([-90 0])
for i = 1:length(l_elecGroups)
    elecInd = strcmp(l_elecGroups{i},name_stripped(left_ind));
    plot3_wrapper(l_snap(elecInd,:),25,cmap(l_idx(i),:));hold all
end


% plot right elecs
h(2) = swagFig;hold all
[v,f] =read_surf_wrapper(surfR);plotsurf_wrapper(v,f);view([90 0])
for i = 1:length(r_elecGroups)
    elecInd = strcmp(r_elecGroups{i},name_stripped(~left_ind));
    plot3_wrapper(r_snap(elecInd,:),40,cmap(r_idx(i),:));hold all
end
axis('off')
keyboard
camlight

%loop through elecs
