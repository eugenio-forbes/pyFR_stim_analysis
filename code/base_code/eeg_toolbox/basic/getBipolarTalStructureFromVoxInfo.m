function [bpTalStruct, hasMom] = getBipolarTalStructureFromVoxInfo(subj,doTalDaemon,source_subj)
%
% FUNCTION:
%   mkGridStructForBipolar_NEW
% 
% DECRIPTION:
%   For one subject, gets electrodes for bipolar and regular
%   montages. It uses the new (june, 2013) vox coords mother file
%   that comes out of voxTool.
%
% INPUTS:
%   subj............. 'TJ052'
%   optional:
%   doTalDeamon ..... whether or not to compute tal stuff
%   source_subj...... in case you are trying to get e.g. TJ069_1 from TJ069
% OUTPUTS:
%
% NOTES:
%   (1) written by jfburke on 06/2013 (john.fred.burke@gmail.com)
%   (2) Updated by AGR on 06/16/2013:
%                   -fixed:
%                   -changed to loading RAW_coords instead of VOX_coords
%                   -talXYZ(1,2,3) --> talXYZ(2,3,4)
%   (3) updated by AGR on 06/24/2013:
                    %-if it cant find an electrode in vox mother,
                    %it checks if it has been removed from leads.txt
                    %if so it isnt an error because, in the past
                    %only channels in leads.txt were taliarached
%   (4) updated by AGR on 08/20/2013:
                    % now it works for a big grid that has been cut up and
                    % placed in various brain regions
%   (5) updated by AGR on 08/23/2013:
                    % updated to appropriately deal with missing electrodes
                    % in the middle of the grid when they have not been included 
                    %in the recording montage                    
% TO DO:
%   ( ) account for sub-grids of varying sizes
%   (X) save bipolar structure for easy loading
%   (x) inocorporate regular talStructure
%   (X) write separate fcn which makes vox_mother from voxCoords and
%       gridStructForBipolar: see vox2voxMother
%   (x) write a separate function that makes jacksheet from electrodes.m:
%   see electrodes2jacksheet
%   ( ) if missing RAW_Coords, call vox2raw

% am I mounting RHINO?
if ismac
  MOUNT_DIR = '/Volumes/RHINO_root';
else
  MOUNT_DIR = '';
end
if ~exist('badTagNames','var')||isempty(badTagNames)
  badTagNames ={};
end
if ~exist('doTalDaemon','var')||isempty(doTalDaemon)
  doTalDaemon = true;
end
hasMom=true;

% directories
subjDir = fullfile(MOUNT_DIR,'/data/eeg',subj);
docDir  = fullfile(subjDir,'docs');
talDir  = fullfile(subjDir,'tal');
source_talDir = fullfile('/data/eeg',source_subj,'tal');
fsDir   = fullfile(MOUNT_DIR,'/data/eeg/freesurfer/subjects');
fsAvgDir= fullfile(fsDir,'average');
fsSubjDir = fullfile(fsDir,source_subj);

% annotation extension
annotExt  = 'aparc.annot';

% file names
voxMotherFileName = 'VOX_coords_mother.txt';
rawFileName       = 'RAW_coords.txt';
avgSurfFileName   = 'RAW_coords_avgSurf.txt';
avgSurfSnapFileName = 'RAW_coords_avgSurf_snap.txt';
avgSurfEsnapFileName = 'RAW_coords_avgSurf_eSnap.txt';
indivSurfFileName = 'RAW_coords_indivSurf.txt';
indivSurfSnapFileName = 'RAW_coords_indivSurf_snap.txt';
indivSurfEsnapFileName = 'RAW_coords_indivSurf_eSnap.txt';
jackFileName      = 'jacksheet.txt';
jackFileName_2    = 'jack_sheet.txt';
elecFileName      = 'electrodes.m';
leadsFileName     = 'leads.txt';
talStructFileName = [subj '_talLocs_database.mat'];
bpTalStructFileName = [subj '_talLocs_database_bipol.mat'];


% get the vox file
rawFile = fullfile(talDir,rawFileName);
if ~exist(rawFile,'file')
  warning(' no RAW_coords.txt file');
  hasRaw = false;
else
    hasRaw = true;
    fid = fopen(rawFile);
    X = textscan(fid,'%d%d%d%d');
    fclose(fid);
    talElecs = X{1};
    talXYZ   = [X{2} X{3} X{4}]; %changed this to 2 3 4, instead of 1 2 3; %AGR 6/16/2013
end
% get the avg coords
avgSurfFile = fullfile(talDir, avgSurfFileName);
avgSurfSnapFile = fullfile(talDir, avgSurfSnapFileName);
avgSurfEsnapFile = fullfile(talDir, avgSurfEsnapFileName);
 if ~exist(avgSurfFile,'file')
    warning('getBipolarTalStructureFromVoxInfo:NoAvgSurf',...
        '%s has no avgSurf file',subj)
    hasAvg = false;
else
    hasAvg = true;
    
    if ~hasRaw && ~hasAvg
        error('getBipolarTalStructureFromVoxInfo:NoGoodSurfs',...
            '%s has no avgSurf or RAW_coords files',subj);
    end
    fid = fopen(avgSurfFile);
    X = textscan(fid,'%d %d %d %d');
    fclose(fid);
    avgSurfElecs = X{1};
    avgSurfXYZ = [X{2} X{3} X{4}];
    
    if exist(avgSurfSnapFile,'file')
        fid = fopen(avgSurfSnapFile);
        X = textscan(fid,'%d %d %d %d');
        fclose(fid);
        avgSurfSnapElecs = X{1};
        avgSurfSnapXYZ = [X{2} X{3} X{4}];
        hasAvgSnap = true;
    else
        hasAvgSnap = false;
    end
    
    if exist(avgSurfEsnapFile,'file')
        fid = fopen(avgSurfEsnapFile);
        X = textscan(fid,'%d %d %d %d');
        fclose(fid);
        avgSurfEsnapElecs = X{1};
        avgSurfEsnapXYZ = [X{2} X{3} X{4}];
        hasAvgEsnap = true;
    else
        hasAvgEsnap = false;
    end
 end

 % get the indiv coords
indivSurfFile = fullfile(talDir, indivSurfFileName);
indivSurfSnapFile = fullfile(talDir, indivSurfSnapFileName);
indivSurfEsnapFile = fullfile(talDir, indivSurfEsnapFileName);
 if ~exist(indivSurfFile,'file')
    warning('getBipolarTalStructureFromVoxInfo:NoAvgSurf',...
        '%s has no avgSurf file',subj)
    hasIndiv = false;
else
    hasIndiv = true;
    
    fid = fopen(indivSurfFile);
    X = textscan(fid,'%d %d %d %d');
    fclose(fid);
    indivSurfElecs = X{1};
    indivSurfXYZ = [X{2} X{3} X{4}];
    
    if exist(indivSurfSnapFile,'file')
        fid = fopen(indivSurfSnapFile);
        X = textscan(fid,'%d %d %d %d');
        fclose(fid);
        indivSurfSnapElecs = X{1};
        indivSurfSnapXYZ = [X{2} X{3} X{4}];
        hasIndivSnap = true;
    else
        hasIndivSnap = false;
    end
    
    if exist(indivSurfEsnapFile,'file')
        fid = fopen(indivSurfEsnapFile);
        X = textscan(fid,'%d %d %d %d');
        fclose(fid);
        indivSurfEsnapElecs = X{1};
        indivSurfEsnapXYZ = [X{2} X{3} X{4}];
        hasIndivEsnap = true;
    else
        hasIndivEsnap = false;
    end
 end

 
if exist(fullfile(talDir,'bpPairs.mat'),'file')
    load(fullfile(talDir,'bpPairs.mat'))
else
    [allpairs, allTagNames, allGrpNames, allElecType] = getBipolarPairs(subj,source_subj);
end

% reorganize the electrodes to be in order
[~,sortInd] = sort(allpairs(:,1));
allpairs    = allpairs(sortInd,:);
allTagNames = allTagNames(sortInd,:);
allGrpNames = allGrpNames(sortInd,:);
allElecType = allElecType(sortInd);

% now make the events structure
bpTalStruct.subject   = [];
bpTalStruct.channel   = [];
bpTalStruct.tagName   = [];
bpTalStruct.grpName   = [];
bpTalStruct.x         = [];
bpTalStruct.y         = [];
bpTalStruct.z         = [];
bpTalStruct.Loc1      = [];
bpTalStruct.Loc2      = [];
bpTalStruct.Loc3      = [];
bpTalStruct.Loc4      = [];
bpTalStruct.Loc5      = [];
bpTalStruct.Loc6      = [];
bpTalStruct.Montage   = [];
bpTalStruct.eNames    = [];
bpTalStruct.eType     = [];
bpTalStruct.bpDistance = [];

bpTalStruct.avgSurf = struct();
bpTalStruct.avgSurf.x = [];
bpTalStruct.avgSurf.y = [];
bpTalStruct.avgSurf.z = [];
bpTalStruct.avgSurf.bpDistance = [];
bpTalStruct.avgSurf.anatRegion = [];
bpTalStruct.avgSurf.x_snap = [];
bpTalStruct.avgSurf.y_snap = [];
bpTalStruct.avgSurf.z_snap = [];
bpTalStruct.avgSurf.anatRegion_snap = [];
bpTalStruct.avgSurf.distance_snap = [];
bpTalStruct.avgSurf.bpDistance_snap = [];
bpTalStruct.avgSurf.x_eSnap = [];
bpTalStruct.avgSurf.y_eSnap = [];
bpTalStruct.avgSurf.z_eSnap = [];
bpTalStruct.avgSurf.anatRegion_eSnap =[];
bpTalStruct.avgSurf.distance_eSnap = [];
bpTalStruct.avgSurf.bpDistance_eSnap = [];
bpTalStruct.avgSurf.path2surfL = [];
bpTalStruct.avgSurf.path2surfR = [];

bpTalStruct.indivSurf = struct();
bpTalStruct.indivSurf.x = [];
bpTalStruct.indivSurf.y = [];
bpTalStruct.indivSurf.z = [];
bpTalStruct.indivSurf.bpDistance = [];
bpTalStruct.indivSurf.anatRegion =[];
bpTalStruct.indivSurf.x_snap = [];
bpTalStruct.indivSurf.y_snap = [];
bpTalStruct.indivSurf.z_snap = [];
bpTalStruct.indivSurf.anatRegion_snap = [];
bpTalStruct.indivSurf.distance_snap = [];
bpTalStruct.indivSurf.bpDistance_snap =[] ;
bpTalStruct.indivSurf.x_eSnap = [];
bpTalStruct.indivSurf.y_eSnap = [];
bpTalStruct.indivSurf.z_eSnap = [];
bpTalStruct.indivSurf.anatRegion_eSnap =[];
bpTalStruct.indivSurf.distance_eSnap = [];
bpTalStruct.indivSurf.bpDistance_eSnap = [];
bpTalStruct.indivSurf.path2surfL = [];
bpTalStruct.indivSurf.path2surfR = [];

path2indivSurfL = fullfile(fsSubjDir,'surf/lh.pial');
path2indivSurfR = fullfile(fsSubjDir,'surf/rh.pial');

path2avgSurfL = fullfile(fsAvgDir,'surf/lh.pial');
path2avgSurfR = fullfile(fsAvgDir,'surf/rh.pial');

fprintf('adding the talairach labels to the electrodes:\n')
ticker   = 0;
tick_inc = 5;
fprintf('   Wait please: ')
allPairs_mask = true(size(allpairs,1),1);
for i=1:size(allpairs,1)
   if ~any(allpairs(i,1)==talElecs) || ~any(allpairs(i,2)==talElecs)
       allPairs_mask(i) = false;
   end
end
allpairs = allpairs(allPairs_mask, :);
for k=1:size(allpairs,1)
  if k./size(allpairs,1)*100>=ticker
    fprintf(' %d%%', ticker)
    ticker=ticker+tick_inc;
  end
  
  % The set of electrodes (e.g. [1,2])
  thisElec  = allpairs(k,:);
  
  % If RAW_coords exists (it always should), get the bipolar pair
  if hasRaw
      thisXYZ_1 = double(talXYZ(thisElec(1)==talElecs,:));
      thisXYZ_2 = double(talXYZ(thisElec(2)==talElecs,:));
      thisXYZ   = double(bipolarHalfDistance_local(thisXYZ_1,thisXYZ_2));
      rawDist = sqrt(sum((thisXYZ_1-thisXYZ_2).^2));
  else
      thisXYZ = nan(1,3);
      rawDist = nan;
  end
     
  % If the average surface exists
  if hasAvg
      % Get the distance for unsnapped
      thisAvgSurfXYZ_1 = double(avgSurfXYZ(thisElec(1)==avgSurfElecs,:));
      thisAvgSurfXYZ_2 = double(avgSurfXYZ(thisElec(2)==avgSurfElecs,:));
      thisAvgSurfXYZ   = double(bipolarHalfDistance_local(...
          thisAvgSurfXYZ_1, thisAvgSurfXYZ_2));
      avgBpDistance = sqrt(sum((thisAvgSurfXYZ_1-thisAvgSurfXYZ_2).^2));
      
      % Get the bp pair for snapped, if it exists
      if hasAvgSnap
          thisAvgSurfSnapXYZ_1 = double(avgSurfSnapXYZ(thisElec(1)==avgSurfSnapElecs,:));
          thisAvgSurfSnapXYZ_2 = double(avgSurfSnapXYZ(thisElec(2)==avgSurfSnapElecs,:));
          thisAvgSurfSnapXYZ   = double(bipolarHalfDistance_local(...
              thisAvgSurfSnapXYZ_1, thisAvgSurfSnapXYZ_2));
          avgSnapBpDistance = sqrt(sum((thisAvgSurfSnapXYZ_1-thisAvgSurfSnapXYZ_2).^2));
          avgDistance_snap = sqrt(sum((thisAvgSurfXYZ-thisAvgSurfSnapXYZ).^2));
      else
          thisAvgSurfSnapXYZ = nan(1,3);
          avgSnapBpDistance = nan;
          avgDistance_snap = nan;
      end
      
      % Get the bp pair for eSnapped, if it exists
      if hasAvgEsnap
          thisAvgSurfEsnapXYZ_1 = double(avgSurfEsnapXYZ(thisElec(1)==avgSurfEsnapElecs,:));
          thisAvgSurfEsnapXYZ_2 = double(avgSurfEsnapXYZ(thisElec(2)==avgSurfEsnapElecs,:));
          thisAvgSurfEsnapXYZ   = double(bipolarHalfDistance_local(...
              thisAvgSurfEsnapXYZ_1, thisAvgSurfEsnapXYZ_2));
          avgEsnapBpDistance = sqrt(sum((thisAvgSurfEsnapXYZ_1-thisAvgSurfEsnapXYZ_2).^2));
          avgDistance_eSnap = sqrt(sum((thisAvgSurfXYZ-thisAvgSurfEsnapXYZ).^2));
      else
          thisAvgSurfEsnapXYZ = nan(1,3);
          avgEsnapBpDistance = nan;
          avgDistance_eSnap = nan;
      end
  else
      % Otherwise, everything's nans
      thisAvgSurfXYZ = nan(1,3);
      thisAvgSurfSnapXYZ = nan(1,3);
      thisAvgSurfEsnapXYZ = nan(1,3);
      avgBpDistance = nan;
      avgSnapBpDistance = nan;
      avgEsnapBpDistance = nan;
      avgDistance_snap = nan;
      avgDistance_eSnap = nan;
  end

  % If the individual surface exists
  if hasIndiv
      % Get the bipolar pair
      thisIndivSurfXYZ_1 = double(indivSurfXYZ(thisElec(1)==indivSurfElecs,:));
      thisIndivSurfXYZ_2 = double(indivSurfXYZ(thisElec(2)==indivSurfElecs,:));
      thisIndivSurfXYZ   = double(bipolarHalfDistance_local(...
                                  thisIndivSurfXYZ_1, thisIndivSurfXYZ_2));
      indivBpDistance = sqrt(sum((thisIndivSurfXYZ_1-thisIndivSurfXYZ_2).^2));
                         
      % If snapped exists, get that bipolar pair
      if hasIndivSnap
          thisIndivSurfSnapXYZ_1 = double(indivSurfSnapXYZ(thisElec(1)==indivSurfSnapElecs,:));
          thisIndivSurfSnapXYZ_2 = double(indivSurfSnapXYZ(thisElec(2)==indivSurfSnapElecs,:));
          thisIndivSurfSnapXYZ   = double(bipolarHalfDistance_local(...
              thisIndivSurfSnapXYZ_1, thisIndivSurfSnapXYZ_2));
          indivSnapBpDistance = sqrt(sum((thisIndivSurfSnapXYZ_1-thisIndivSurfSnapXYZ_2).^2));
          indivDistance_snap = sqrt(sum((thisIndivSurfXYZ-thisIndivSurfSnapXYZ).^2));
      else
          thisIndivSurfSnapXYZ = nan(1,3);
          indivSnapBpDistance = nan;
          indivDistance_snap = nan;
      end
      
      % If eSnapped exists, get that bipolar pair
      if hasIndivEsnap
          thisIndivSurfEsnapXYZ_1 = double(indivSurfEsnapXYZ(thisElec(1)==indivSurfEsnapElecs,:));
          thisIndivSurfEsnapXYZ_2 = double(indivSurfEsnapXYZ(thisElec(2)==indivSurfEsnapElecs,:));
          thisIndivSurfEsnapXYZ   = double(bipolarHalfDistance_local(...
              thisIndivSurfEsnapXYZ_1, thisIndivSurfEsnapXYZ_2));
          indivEsnapBpDistance = sqrt(sum((thisIndivSurfEsnapXYZ_1-thisIndivSurfEsnapXYZ_2).^2));
          indivDistance_eSnap = sqrt(sum((thisIndivSurfXYZ-thisIndivSurfEsnapXYZ).^2));
      else
          thisIndivSurfEsnapXYZ = nan(1,3);
          indivEsnapBpDistance = nan;
          indivDistance_eSnap = nan;
      end
  else
      % Otherwise, nans.
     thisIndivSurfXYZ = nan(1,3);
     thisIndivSurfSnapXYZ = nan(1,3);
     thisIndivSurfEsnapXYZ = nan(1,3);
     indivBpDistance = nan;
     indivSnapBpDistance = nan;
     indivEsnapBpDistance = nan;
     indivDistance_snap = nan;
     indivDistance_eSnap = nan;
  end

  % Fill in the bipolar structure
  % (anatRegion must be left empty)
  bpTalStruct(k).subject  = subj;
  bpTalStruct(k).channel  = thisElec;
  bpTalStruct(k).tagName  = [allTagNames{k,1} '-' allTagNames{k,2}];
  bpTalStruct(k).grpName = allGrpNames{k};
  bpTalStruct(k).eNames   =  [num2str(thisElec(1,1)) '-' num2str(thisElec(1,2))];
  
  bpTalStruct(k).x        = thisXYZ(1);
  bpTalStruct(k).y        = thisXYZ(2);
  bpTalStruct(k).z        = thisXYZ(3);
  bpTalStruct(k).bpDistance = rawDist;
  
  bpTalStruct(k).avgSurf.x = thisAvgSurfXYZ(1);
  bpTalStruct(k).avgSurf.y = thisAvgSurfXYZ(2);
  bpTalStruct(k).avgSurf.z = thisAvgSurfXYZ(3);
  bpTalStruct(k).avgSurf.bpDistance = avgBpDistance;
  bpTalStruct(k).avgSurf.anatRegion = [];
  
  bpTalStruct(k).avgSurf.x_snap = thisAvgSurfSnapXYZ(1);
  bpTalStruct(k).avgSurf.y_snap = thisAvgSurfSnapXYZ(2);
  bpTalStruct(k).avgSurf.z_snap = thisAvgSurfSnapXYZ(3);
  bpTalStruct(k).avgSurf.bpDistance_snap = avgSnapBpDistance;
  bpTalStruct(k).avgSurf.anatRegion_snap = [];
  bpTalStruct(k).avgSurf.distance_snap = avgDistance_snap;
  
  bpTalStruct(k).avgSurf.x_eSnap = thisAvgSurfEsnapXYZ(1);
  bpTalStruct(k).avgSurf.y_eSnap = thisAvgSurfEsnapXYZ(2);
  bpTalStruct(k).avgSurf.z_eSnap = thisAvgSurfEsnapXYZ(3);
  bpTalStruct(k).avgSurf.bpDistance_eSnap = avgEsnapBpDistance;
  bpTalStruct(k).avgSurf.anatRegion_eSnap = [];
  bpTalStruct(k).avgSurf.distance_eSnap = avgDistance_eSnap;
  
  bpTalStruct(k).avgSurf.path2surfL = path2avgSurfL;
  bpTalStruct(k).avgSurf.path2surfR = path2avgSurfR;
  
  bpTalStruct(k).indivSurf.x = thisIndivSurfXYZ(1);
  bpTalStruct(k).indivSurf.y = thisIndivSurfXYZ(2);
  bpTalStruct(k).indivSurf.z = thisIndivSurfXYZ(3);
  bpTalStruct(k).indivSurf.bpDistance = indivBpDistance;
  bpTalStruct(k).indivSurf.anatRegion = [];
  
  bpTalStruct(k).indivSurf.x_snap = thisIndivSurfSnapXYZ(1);
  bpTalStruct(k).indivSurf.y_snap = thisIndivSurfSnapXYZ(2);
  bpTalStruct(k).indivSurf.z_snap = thisIndivSurfSnapXYZ(3);
  bpTalStruct(k).indivSurf.bpDistance_snap = indivSnapBpDistance;
  bpTalStruct(k).indivSurf.anatRegion_snap =[];
  bpTalStruct(k).indivSurf.distance_snap = indivDistance_snap;
    
  bpTalStruct(k).indivSurf.x_eSnap = thisIndivSurfEsnapXYZ(1);
  bpTalStruct(k).indivSurf.y_eSnap = thisIndivSurfEsnapXYZ(2);
  bpTalStruct(k).indivSurf.z_eSnap = thisIndivSurfEsnapXYZ(3);
  bpTalStruct(k).indivSurf.bpDistance_eSnap = indivEsnapBpDistance;
  bpTalStruct(k).indivSurf.anatRegion_eSnap = [];
  bpTalStruct(k).indivSurf.distance_eSnap = indivDistance_eSnap;
  
  bpTalStruct(k).indivSurf.path2surfL = path2indivSurfL;
  bpTalStruct(k).indivSurf.path2surfR = path2indivSurfR;
    
  if doTalDaemon
    bpTalStruct(k)          = tal2Region(bpTalStruct(k),false);
  end
  
  bpTalStruct(k).eType =  allElecType{k};
  
  if strcmp(allElecType{k},'D')
    bpTalStruct(k).Montage = 'hipp';
  else
    if thisXYZ(1)<0
      bpTalStruct(k).Montage = 'lsag';
    else
      bpTalStruct(k).Montage = 'rsag';
    end
  end
end

% Now we have to snap the avgSurf and indivSurf elecs of all types to
% surface to get the anatLabels associated
for i=1:2
    if (i==1 && hasAvg) || (i==2 && hasIndiv)
        if i==1
            structType = 'avgSurf';
            typeStruct = [bpTalStruct.avgSurf];
            fsTypeDir = fsAvgDir;
        else
            structType = 'indivSurf';
            typeStruct = [bpTalStruct.indivSurf];
            fsTypeDir = fsSubjDir;
        end
        
        XYZ = [[typeStruct.x]',[typeStruct.y]', [typeStruct.z]'];
        XYZ_snap = [[typeStruct.x_snap]', [typeStruct.y]', [typeStruct.z]'];
        XYZ_eSnap = [[typeStruct.x_eSnap]', [typeStruct.y_eSnap]', [typeStruct.z_eSnap]'];
        
        vL_avg = read_surf(fullfile(fsTypeDir,'surf/lh.pial'));
        vR_avg = read_surf(fullfile(fsTypeDir,'surf/rh.pial'));
        v_all = [vL_avg; vR_avg];
        
        [~,labelL, colortableL] = read_annotation(fullfile(fsTypeDir,['label/lh.' annotExt]));
        [~,labelR, colortableR] = read_annotation(fullfile(fsTypeDir,['label/rh.' annotExt]));
        label_all = [labelL; labelR];
        colortable_all = merge_colortables(colortableL, colortableR);
        anatLabels = getSurfAnatLabels(XYZ, v_all, label_all, colortable_all);
        anatLabels_snap = getSurfAnatLabels(XYZ_snap, v_all, label_all, colortable_all);
        anatLabels_eSnap = getSurfAnatLabels(XYZ_eSnap, v_all, label_all, colortable_all);
        
        for i=1:length(anatLabels)
            bpTalStruct(i).(structType).anatRegion = anatLabels{i};
            bpTalStruct(i).(structType).anatRegion_snap = anatLabels_snap{i};
            bpTalStruct(i).(structType).anatRegion_eSnap = anatLabels_eSnap{i};
        end
    end
end

    
fprintf('\n\n')
%------------------------------------------------------------------
function xyz_out = bipolarHalfDistance_local(XYZ1,XYZ2)
  tmp = -XYZ1+XYZ2;
  halftmp = tmp./2;
  xyz_out = XYZ1 + halftmp;

function c1 = merge_colortables(c1, c2)
c1.struct_names(end+1:end+length(c2.struct_names)) = c2.struct_names;
c1.table(end+1:end+length(c2.struct_names),:) = c2.table;
