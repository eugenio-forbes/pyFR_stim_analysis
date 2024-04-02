function [talStruct] = updateTalStructWithSurfLocs(subj,bipolFlag)
%This function updates the subject's taliarach structure with the location
%of electrodes in surface RAS space and electrode locations after they have been snapped to the surface
%This should be saved in eeg_toolbox/tal

%Inputs:
%subj ... 'TJ059'
%bipolFlag ... 1 (does bipolar), 0 (no bipolar)

%Outputs:
%talStruct ... updated tal structure with surfRAS and surfRAS_snapped
                %fields
% To Do:
%(x) Incorporate average surface stuff.

%Written by Ashwin Ramayya (11/19/2013)
%Update by AGR (03/28/2014); incorporated snap to avg. surf
% get bipolar structure
try
[talStruct] = getBipolarSubjElecs(subj,bipolFlag);
catch
    %error('Tal structure not done')
    disp('ERRROR!!!! tal struct not done; not updating')
    talStruct = [];
    return
end

%check if it already contains surface fields
% if isfield(talStruct,'x_surfRAS')
%     disp('Surface fields have already been added; returning')
%     return
% end

%set dirs and paths
talDir = ['/data/eeg/' subj '/tal/'];
imageDir = fullfile(talDir,'images');
surfVOX = fullfile(talDir,'VOX_coords_mother_surfRAS.txt');
freesurferDir = fullfile('/data/eeg/freesurfer/subjects',subj);
surfDir = fullfile(freesurferDir,'surf');
lblDir = fullfile(freesurferDir,'label');
surfL = fullfile(surfDir,'lh.pial');
sulcL = fullfile(surfDir,'lh.sulc');
surfR = fullfile(surfDir,'rh.pial');
sulcR = fullfile(surfDir,'rh.sulc');
annotL = fullfile(lblDir,'lh.aparc.annot');
annotR = fullfile(lblDir,'rh.aparc.annot');
surfL_w_elecs = fullfile(surfDir,[subj '_lh_pial_w_elecs.mat']);
surfR_w_elecs = fullfile(surfDir,[subj '_rh_pial_w_elecs.mat']);

%left vs. right 
left_ind = [talStruct.x]<0;
eNames = {talStruct.eNames};

% snap params 
correction = [];
auto_correct = 1;
minDepth = 0;
radius = 5;
maxSnap = 10; %the largest distance an electrode can be snapped

%parse indiv recon
if ~exist(surfDir,'dir') || ~exist(surfL,'file') || ~exist(surfR,'file')
     reconDone = 0;
     disp('Surface Reconstruction (recon all) not done')
elseif ~exist(surfVOX,'file') 
     disp('elecs_to_surfRAS.m not done')
     reconDone = 0;
else
    reconDone = 1;
end

%% PART 1: Process Individual Surface Recon (only if recon is done)
if reconDone == 1
   

    % get surfRAS electrode locations
    fid = fopen(surfVOX);
    C = textscan(fid,'%s%d%d%d%s%d%d');
    centroids = [C{:,2} C{:,3} C{:,4}];
    labels = C{:,1};

    % convert centroids to the bipolar montage
    if bipolFlag
        [XYZ_surfRAS,XYZ_surfRAS_lbl,XYZ] = convToBp_local(talStruct,centroids,labels);
        XYZ_surfRAS = double(XYZ_surfRAS);
    end


    % % check plot
    % figure;plot3_wrapper(XYZ,25,'r')
    figure;
    subplot(1,2,1);
    plot3_wrapper(XYZ_surfRAS(left_ind,:),25,'r');hold all
    [v,f] =read_surf_wrapper(surfL);
    plotsurf_wrapper(v,f);
    subplot(1,2,2);
    plot3_wrapper(XYZ_surfRAS(~left_ind,:),25,'r');hold all
    [v,f] =read_surf_wrapper(surfR);
    plotsurf_wrapper(v,f);
    
%     s = input('CONTINUE?? Y or N','s');
%     if strcmp(s,'n')
%         talStruct = [];
%         return;
%     end


    % process left hemisphere: 
    % Snap left elecs
    [l_elecLocs_surfRAS_snapped,~] = snapElectrodesToSurface(XYZ_surfRAS(left_ind,:),...
        surfL, sulcL, [], [], [], XYZ_surfRAS_lbl(left_ind), correction, auto_correct,...
        minDepth,maxSnap); view([-90 0])
    %get anatomical labes
    [l_anatLbls] = getSurfAnatLabels(XYZ_surfRAS(left_ind,:),surfL,annotL);
    [l_anatLbls_snap] = getSurfAnatLabels(l_elecLocs_surfRAS_snapped,surfL,annotL);


    % write a new surface file labeling appropriate faces to the electrode names (see
    % read_surf_wrapper
    writesurf_w_elecs(surfL,l_elecLocs_surfRAS_snapped,eNames(left_ind),radius,surfL_w_elecs);

    % process right hemisphere:
    [r_elecLocs_surfRAS_snapped,~ ] = snapElectrodesToSurface(XYZ_surfRAS(~left_ind,:),...
       surfR, sulcR,[],[],[], XYZ_surfRAS_lbl(~left_ind), correction, auto_correct,...
        minDepth,maxSnap);view([90 0])
    [r_anatLbls] = getSurfAnatLabels(XYZ_surfRAS(~left_ind,:),surfR,annotR);
    [r_anatLbls_snap] = getSurfAnatLabels(r_elecLocs_surfRAS_snapped,surfR,annotR);
    writesurf_w_elecs(surfR,r_elecLocs_surfRAS_snapped,eNames(~left_ind),radius,surfR_w_elecs);

    % convert to XYZ_surfRAS_snap
    XYZ_surfRAS_snap = nan(size(XYZ));
    XYZ_surfRAS_snap(left_ind,:) = l_elecLocs_surfRAS_snapped;
    XYZ_surfRAS_snap(~left_ind,:) = r_elecLocs_surfRAS_snapped;

    %covert anatomical labels to a big cell array
    surfAnatLabels = cell(1,length(left_ind)); 
    surfAnatLabels_snap = surfAnatLabels;

    surfAnatLabels(left_ind) = l_anatLbls;
    surfAnatLabels(~left_ind) = r_anatLbls;
    surfAnatLabels_snap(left_ind) = l_anatLbls_snap;
    surfAnatLabels_snap(~left_ind) = r_anatLbls_snap;

    % update structure
    [talStruct] = updateTalStruct_local(talStruct,XYZ_surfRAS,XYZ_surfRAS_snap,surfAnatLabels,surfAnatLabels_snap,surfL,surfR,surfL_w_elecs,surfR_w_elecs);
end

%% PART 2: Process Average Surface Recon (regardless of if indiv recon is done)
freesurferDir = fullfile('/data/eeg/freesurfer/subjects','average');
surfDir = fullfile(freesurferDir,'surf');
lblDir = fullfile(freesurferDir,'label');
surfL = fullfile(surfDir,'lh.pial');
sulcL = fullfile(surfDir,'lh.sulc');
surfR = fullfile(surfDir,'rh.pial');
sulcR = fullfile(surfDir,'rh.sulc');
annotL = fullfile(lblDir,'lh.aparc.annot');
annotR = fullfile(lblDir,'rh.aparc.annot');

% get tal coords (always gets x_avgSurf)
xyz_tal = [[talStruct.x_avgSurf]' [talStruct.y_avgSurf]' [talStruct.z_avgSurf]'];
left_ind = xyz_tal(:,1)<0;
    %if reconDone
    %    xyz_tal = surfaceRAS_to_tal_coords(subj,[[talStruct.x_surfRAS]' [talStruct.y_surfRAS]' [talStruct.z_surfRAS]']);
    %elseif ~isnan(talStruct.x_avgSurf(1))
    %    xyz_tal = [[talStruct.x_avgSurf]' [talStruct.y_avgSurf]' [talStruct.z_avgSurf]'];
    %else
    %    xyz_tal = [[talStruct.x]' [talStruct.y]' [talStruct.z]' ];
    %end

% process left hemisphere: 
% Snap left elecs
[l_tal_snap,~] = snapElectrodesToSurface(xyz_tal(left_ind,:),...
    surfL, sulcL, [], [], [], {talStruct(left_ind).eNames}, correction, auto_correct,...
    minDepth,maxSnap); view([-90 0])
%get anatomical labes
[l_anatLbls] = getSurfAnatLabels(l_tal_snap,surfL,annotL);


% process right hemisphere: 
% Snap right elecs
[r_tal_snap,~] = snapElectrodesToSurface(xyz_tal(~left_ind,:),...
    surfR, sulcR, [], [], [], {talStruct(~left_ind).eNames}, correction, auto_correct,...
    minDepth,maxSnap); view([90 0])
%get anatomical labes
[r_anatLbls] = getSurfAnatLabels(r_tal_snap,surfR,annotR);


% concatonate stuff
% convert to XYZ_surfRAS_snap
tal_snap = nan(size(xyz_tal));
tal_snap(left_ind,:) = l_tal_snap;
tal_snap(~left_ind,:) = r_tal_snap;

% calculate tal_snap distances
D = pdist2(xyz_tal,tal_snap);
Y = logical(eye(size(D)));
tal_snap_distances = D(Y);

%covert anatomical labels to a big cell array
tal_anatLbls = cell(1,length(left_ind)); 
tal_anatLbls(left_ind) = l_anatLbls;
tal_anatLbls(~left_ind) = r_anatLbls;

% update tal struct with avg data
[talStruct] = updateTalStructWAvg_local(talStruct,tal_snap,tal_anatLbls,tal_snap_distances,surfL,surfR);

%% save out new tal structure
%keyboard
cd(talDir)
if bipolFlag == 1
    subjTalEvents = talStruct;
    fname = [subj '_talLocs_database_bipol.mat'];
    save(fname,'subjTalEvents')
end
close all



%% %%%%%%%%%%% Functions Below %%%%%%%%%%%%%%%%%%%%

%--------
function [talStruct] = updateTalStructWAvg_local(talStruct,tal_snap,tal_anatLbls,tal_snap_distances,surfL,surfR)
for i = 1:length(talStruct)
    %_surfRAS
    talStruct(i).('x_tal_snap') = tal_snap(i,1);
    talStruct(i).('y_tal_snap') = tal_snap(i,2);
    talStruct(i).('z_tal_snap') = tal_snap(i,3);
    % anatLabels_snapped
    talStruct(i).('anatRegion_tal_snap') = tal_anatLbls{i};
  
    %snap_distance
    talStruct(i).('tal_snap_distance') = tal_snap_distances(i);
    
    %paths to surface files (for convenience)
    talStruct(i).path_surfL = surfL;
    talStruct(i).path_surfR = surfR;
end

%------------------------------------------------------------------
function [talStruct] = updateTalStruct_local(talStruct,XYZ_surfRAS,XYZ_surfRAS_snap,surfAnatLabels,surfAnatLabels_snap,surfL,surfR,surfL_w_elecs,surfR_w_elecs);
for i = 1:length(talStruct)
    %_surfRAS
    talStruct(i).('x_surfRAS') = XYZ_surfRAS(i,1);
    talStruct(i).('y_surfRAS') = XYZ_surfRAS(i,2);
    talStruct(i).('z_surfRAS') = XYZ_surfRAS(i,3);
    % anatLabels_snapped
    talStruct(i).('anatRegion_surfRAS') = surfAnatLabels{i};
  
    %_surfRAS_snapped
    talStruct(i).('x_surfRAS_snap2pial') = XYZ_surfRAS_snap(i,1);
    talStruct(i).('y_surfRAS_snap2pial') = XYZ_surfRAS_snap(i,2);
    talStruct(i).('z_surfRAS_snap2pial') = XYZ_surfRAS_snap(i,3); 
    
    % anatLabels_snapped
    talStruct(i).('anatRegion_surfRAS_snap2pial') = surfAnatLabels_snap{i};
    %talStruct(i).('lobe_surfRAS_snap2pial') = reg2lobe_local(surfAnatLabels{i});
    
    %snap_distance
    talStruct(i).('snap_distance') = pdist([XYZ_surfRAS(i,:);XYZ_surfRAS_snap(i,:)]);
    
    %paths to surface files (for convenience)
    talStruct(i).path_surfL = surfL;
    talStruct(i).path_surfR = surfR;
    talStruct(i).path_surfL_w_elecs = surfL_w_elecs;
    talStruct(i).path_surfR_w_elecs = surfR_w_elecs;
end




%------------------------------------------------------------------
% function[lobe] = reg2lobe_local(reg)
% %http://ftp.nmr.mgh.harvard.edu/fswiki/CorticalParcellation
% switch reg
%     case {'superiorfrontal','rostralmiddlefrontal','caudalmiddlefrontal',...
%             'parsopercularis','parstriangularis','parsorbitalis','lateralorbitofrontal',...
%             'medialorbitofrontal','precentral','paracentral','frontalpole'}
%         lobe = 'frontal';
%     case {'superiorparietal','inferiorparietal','supramarginal','postcentral','precuneus'}
%         lobe = 'parietal';
%     case {'superiortemporal','middletemporal','inferiortemporal',}
%         lobe = 'temporal';
%     case {}
%         lobe = 'occipital';
%     case {}
%         lobe = 'cingulate';
% end

function [centXYZ,centXYZ_lbls,XYZ] = convToBp_local(talStruct,centroids,labels)
%output is a list of centXYZs corresponding with talStruct
for i = 1:length(talStruct)
    % identify electrodes
    c = regexp(talStruct(i).tagName,'-','split');% num2str(talStruct(i).channel(1))];
    elec1 = c{1};
    elec2 = c{2};
    % get the matching centroids
    centXYZ1 = centroids(strcmp(elec1,labels),:);
    centXYZ2 = centroids(strcmp(elec2,labels),:);
    
%     % check whether elec names mate
%     if ~regexp(labels{elec1},talStruct(i).tagName) 
%         error('mismatch in tagName (from struct) and centroid label')
%     end
%     if ~regexp(labels{elec2},talStruct(i).tagName) 
%         error('mismatch in tagName (from struct) and centroid label')
%     end
    
    centXYZ(i,:) = bipolarHalfDistance_local(centXYZ1,centXYZ2);
    centXYZ_lbls{i} = talStruct(i).tagName;
    XYZ(i,:) = [talStruct(i).x talStruct(i).y talStruct(i).z];
end

%------------------------------------------------------------------
function xyz_out = bipolarHalfDistance_local(XYZ1,XYZ2)
  tmp = -XYZ1+XYZ2;
  halftmp = tmp./2;
  xyz_out = XYZ1 + halftmp;