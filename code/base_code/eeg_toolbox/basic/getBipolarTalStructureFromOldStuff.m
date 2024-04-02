function bpTalStruct = getBipolarTalStructureFromOldStuff(subj);         
%
%
%
%
%
%
%
%
%
%

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

% directories
subjDir = fullfile(MOUNT_DIR,'/data/eeg',subj);
docDir  = fullfile(subjDir,'docs');
talDir  = fullfile(subjDir,'tal');
fsDir   = fullfile(MOUNT_DIR,'/data/eeg/freesurfer/subjects');
fsAvgDir= fullfile(fsDir,'average');
fsSubjDir = fullfile(fsDir,subj);


% annotation extension
annotExt  = 'aparc.annot';

% file names
rawFileName         = 'RAW_coords.txt';
avgSurfFileName     = 'RAW_coords_avgSurf.txt';
jackFileName        = 'jacksheet.txt';
jackFileName_2      = 'jack_sheet.txt';
bpMan_fileName      = 'bp_manual.txt';
gridsFile           = 'grids_struct_for_bipolar.m';
elecFileName        = 'electrodes.m';
leadsFileName       = 'leads.txt';
talStructFileName   = [subj '_talLocs_database.mat'];
bpTalStructFileName = [subj '_talLocs_database_bipol.mat'];
avgSurfFileName   = 'RAW_coords_avgSurf.txt';
avgSurfSnapFileName = 'RAW_coords_avgSurf_snap.txt';
avgSurfEsnapFileName = 'RAW_coords_avgSurf_eSnap.txt';
indivSurfFileName = 'RAW_coords_indivSurf.txt';
indivSurfSnapFileName = 'RAW_coords_indivSurf_snap.txt';
indivSurfEsnapFileName = 'RAW_coords_indivSurf_eSnap.txt';

% get the vox file
rawFile = fullfile(talDir,rawFileName);
if ~exist(rawFile,'file')
  warning('getBipolarTalStructureFromOldStuff:NoAvgSurf', ...
      '%s has no RAW_coords.txt file', subj)  
  hasRaw = false;
else
    hasRaw = true;
    fid = fopen(rawFile);
    X = textscan(fid,'%d%d%d%d');
    fclose(fid);
    talElecs = X{1};
    talXYZ   = [X{2} X{3} X{4}];
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


% get jacksheet
[elecNum elecNam elecNam_stripped] = ...
    getTheContentsOFTheJackFile_local(docDir,jackFileName,jackFileName_2);
if isempty(elecNum);error(' no jacksheet');end

% get the bp manual file
bpMan_file = fullfile(docDir, bpMan_fileName);
if ~exist(bpMan_file,'file')
  error(' no bp manual file') 
end
fid = fopen(bpMan_file);
X   = textscan(fid,'%d%d%s');
fclose(fid);
bpElecs = [X{1} X{2}];
bpMont  = X{3};

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
count=0;
for k=1:size(bpElecs,1)
  if k./size(bpElecs,1)*100>=ticker
    fprintf(' %d%%', ticker)
    ticker=ticker+tick_inc;
  end
  thisElec  = bpElecs(k,:);
  if hasRaw
      thisXYZ_1 = double(talXYZ(thisElec(1)==talElecs,:));
      thisXYZ_2 = double(talXYZ(thisElec(2)==talElecs,:));
      %diffXYZ  = thisXYZ_1-thisXYZ_2
      thisXYZ   = double(bipolarHalfDistance_local(thisXYZ_1,thisXYZ_2));
      %sqrt(sum((diffXYZ).^2))
      rawDist = sqrt(sum((thisXYZ_1-thisXYZ_2).^2));
  else
      thisXYZ = nan(1,3);
      rawDist = nan;
  end
     
  if hasAvg
      thisAvgSurfXYZ_1 = double(avgSurfXYZ(thisElec(1)==avgSurfElecs,:));
      thisAvgSurfXYZ_2 = double(avgSurfXYZ(thisElec(2)==avgSurfElecs,:));
      thisAvgSurfXYZ   = double(bipolarHalfDistance_local(...
          thisAvgSurfXYZ_1, thisAvgSurfXYZ_2));
      avgBpDistance = sqrt(sum((thisAvgSurfXYZ_1-thisAvgSurfXYZ_2).^2));
      
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
      thisAvgSurfXYZ = nan(1,3);
      thisAvgSurfSnapXYZ = nan(1,3);
      thisAvgSurfEsnapXYZ = nan(1,3);
      avgBpDistance = nan;
      avgSnapBpDistance = nan;
      avgEsnapBpDistance = nan;
      avgDistance_snap = nan;
      avgDistance_eSnap = nan;
  end

  if hasIndiv
      thisIndivSurfXYZ_1 = double(indivSurfXYZ(thisElec(1)==indivSurfElecs,:));
      thisIndivSurfXYZ_2 = double(indivSurfXYZ(thisElec(2)==indivSurfElecs,:));
      thisIndivSurfXYZ   = double(bipolarHalfDistance_local(...
                                  thisIndivSurfXYZ_1, thisIndivSurfXYZ_2));
      indivBpDistance = sqrt(sum((thisIndivSurfXYZ_1-thisIndivSurfXYZ_2).^2));
                         
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
     thisIndivSurfXYZ = nan(1,3);
     thisIndivSurfSnapXYZ = nan(1,3);
     thisIndivSurfEsnapXYZ = nan(1,3);
     indivBpDistance = nan;
     indivSnapBpDistance = nan;
     indivEsnapBpDistance = nan;
     indivDistance_snap = nan;
     indivDistance_eSnap = nan;
  end

   thisTag_1_idx = elecNum==thisElec(1);
   thisTag_2_idx = elecNum==thisElec(2);  
   if ~isempty(thisTag_1_idx) && ~isempty(thisTag_2_idx)
     thisTag   = [elecNam{thisTag_1_idx} '-' elecNam{thisTag_1_idx}];
   else
     thisTag   = 'XXX';
   end

  bpTalStruct(k).subject  = subj;
  bpTalStruct(k).channel  = thisElec;
  bpTalStruct(k).tagName  = thisTag;
  bpTalStruct(k).grpName = elecNam_stripped{thisTag_1_idx};
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
  
  bpTalStruct(k).eType =  'S';%allElecType{k};
  
  if strcmp(bpTalStruct(k).eType,'D')
    bpTalStruct(k).Montage = 'hipp';
  else
    if thisXYZ(1)<0
      bpTalStruct(k).Montage = 'lsag';
    else
      bpTalStruct(k).Montage = 'rsag';
    end
  end
end

% Now we have to get the anatLabels associated with all of the electrodes
% of each type. Since it has to dsearchn, it makes more sense to do this
% all at once.
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
% for k=1:size(bpElecs,1)
%   if k./size(bpElecs,1)*100>=ticker
%     fprintf(' %d%%', ticker)
%     ticker=ticker+tick_inc;
%   end
%   thisElec      = bpElecs(k,:);
%   
%   % Deal with RAW_coords
%   if hasRaw
%       thisXYZ_1     = double(talXYZ(thisElec(1)==talElecs,:));
%       if isempty(thisXYZ_1)
%           fprintf('%d not found in VOX.\n',thisElec(1))
%           continue
%       end
%       thisXYZ_2     = double(talXYZ(thisElec(2)==talElecs,:));
%       if isempty(thisXYZ_2)
%           fprintf('%d not found in VOX.\n',thisElec(2))
%           continue
%       end
%       thisXYZ       = double(bipolarHalfDistance_local(thisXYZ_1,thisXYZ_2));
%   else
%       thisXYZ = nan(1,3);
%   end
%   
%   % Deal with RAW_coords_avgSurf
%   if hasAvg
%       thisAvgSurfXYZ_1     = double(avgSurfXYZ(thisElec(1)==avgSurfElecs,:));
%       if isempty(thisAvgSurfXYZ_1)
%           fprintf('%d not found in VOX.\n',thisElec(1))
%           continue
%       end
%       thisAvgSurfXYZ_2     = double(avgSurfXYZ(thisElec(2)==avgSurfElecs,:));
%       if isempty(thisAvgSurfXYZ_2)
%           fprintf('%d not found in VOX.\n',thisElec(2))
%           continue
%       end
%       thisAvgSurfXYZ ...
%           = double(bipolarHalfDistance_local(thisAvgSurfXYZ_1,thisAvgSurfXYZ_2));
%   else
%       thisAvgSurfXYZ = nan(1,3);
%   end
%   
%   
%   thisTag_1_idx = elecNum==thisElec(1);
%   thisTag_2_idx = elecNum==thisElec(2);  
%   if ~isempty(thisTag_1_idx) && ~isempty(thisTag_2_idx)
%     thisTag   = [elecNam{thisTag_1_idx} '-' elecNam{thisTag_1_idx}];
%   else
%     thisTag   = 'XXX';
%   end
%   count = count + 1;
%   bpTalStruct(count).subject  = subj;
%   bpTalStruct(count).channel  = thisElec;
%   bpTalStruct(count).tagName  = thisTag;  
%   bpTalStruct(count).grpName = elecNam_stripped{thisTag_1_idx};
%   bpTalStruct(count).x        = thisXYZ(1);
%   bpTalStruct(count).y        = thisXYZ(2);
%   bpTalStruct(count).z        = thisXYZ(3);  
%   bpTalStruct(count).x_avgSurf= thisAvgSurfXYZ(1);
%   bpTalStruct(count).y_avgSurf= thisAvgSurfXYZ(2);
%   bpTalStruct(count).z_avgSurf= thisAvgSurfXYZ(3);
%   
%   if doTalDaemon
%     bpTalStruct(count)          = tal2Region(bpTalStruct(count),false);
%   end
%   if hasRaw
%     bpTalStruct(count).distance = sqrt(sum((thisXYZ_1-thisXYZ_2).^2));
%   else
%       bpTalStruct(count).distance = sqrt(sum((thisAvgSurfXYZ_1-thisAvgSurfXYZ_2).^2));
%   end
%   bpTalStruct(count).eNames   =  [num2str(thisElec(1,1)) '-' num2str(thisElec(1,2))];
%  
%   if thisXYZ(1)<0
%     bpTalStruct(count).Montage = 'lsag';
%   else
%     bpTalStruct(count).Montage = 'rsag';
%   end
% end
fprintf('\n\n')

%---------------------------------------------------------------------
function  allpairs = make_a_bipolar_motage(startval,numElecs_dim1,numElecs_dim2)
    
  % this is all Aswin G. Ramayya.  This occurred on June 14, 2013.
  % Ashwin is completelty responsible for this code and miustakes
  % therein.  John F. Burke is not reponsible for any errors,
  % mistakes, or general bad things that come out of this code. But
  % John Burke checked thouroughly and has placed his personal
  % stamp of approval, thus accepting all responsibility (but no
  % credit). Ashwin would like to add, spontaneously, that he is a
  % tool. John agrees.
  allpairs = [];
  for j = 1:numElecs_dim2
    for i = 1:numElecs_dim1
      if i < numElecs_dim1
	allpairs(end+1,:) = [i-1 i] + startval + numElecs_dim1*((j-1));
      end
      if j < numElecs_dim2
        allpairs(end+1,:) = [i-1 i+numElecs_dim1-1] + startval + numElecs_dim1*((j-1));
      end
    end
  end
  
        
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % for j=1:numElecs_dim1
%     if j<numElecs_dim1
%       % along first dimension    
%       allpairs(end+1,:) = [j-1 j] + startval;
%     end
%   end
  
%   for j=1:numElecs_dim1
%     for k=1:numElecs_dim2-1
%       % go parallel along multiple rows
%       if j<numElecs_dim1
% 	allpairs(end+1,:) = [k*numElecs_dim1-1 k*numElecs_dim1] + startval + j;
%       end
%     end
%   end
  
%   for j=1:numElecs_dim1
%     for k=1:numElecs_dim2-1
%       % go perpendicularly
%       allpairs(end+1,:) = [(k-1)*numElecs_dim1-1 k*numElecs_dim1-1] + startval + j;
%     end
%   end
  

%------------------------------------------------------------------
function xyz_out = bipolarHalfDistance_local(XYZ1,XYZ2)
  tmp = -XYZ1+XYZ2;
  halftmp = tmp./2;
  xyz_out = XYZ1 + halftmp;

%---------------------------------------------------------------------
function out = getTheContentsOFVoxMom_local(d,f);
  fil = fullfile(d,f);
  fid = fopen(fil,'r');
  if fid==-1
    out = [];
    fprintf('%s does not exist...EXITING\n\n',f)
    return
  end  
  X = textscan(fid,'%s%d%d%d%s%s','delimiter','\t');
  fclose(fid);
  gridNames    = X{1};
  elecXYZ      = [X{2} X{3} X{4}];
  stripOrDepth = X{5};
  geometry_tmp = X{6};
  out          = [];
  for k=1:size(geometry_tmp,1)
    out(k).name = gridNames{k};
    out(k).voxP = elecXYZ(k,:);
    out(k).type = stripOrDepth{k};
    out(k).geom = sscanf(geometry_tmp{k},'%d%d')';
  end
  
%--------------------------------------------------
function out = getTheContentsOFTheElecs_local(d,f);
  elecFile = fullfile(d,f);
  if ~exist(elecFile,'file')  
    out = [];
    fprintf('%s does not exist...EXITING\n\n',f)
    return
  end  
  run(elecFile);
  if exist('r','var')
    out = r;
  elseif exist('grids','var')
    out = grids;
  else
    out = [];
    fprintf('The contects of electrodes.m not right...EXITING\n\n')
    return
  end

%------------------------------------------------------------  
function [num nam nam_stripped] = getTheContentsOFTheJackFile_local(dDir,j1,j2);
  jFil_1 = fullfile(dDir,j1);
  jFil_2 = fullfile(dDir,j2);  
  fid_1  = fopen(jFil_1,'r');
  fid_2  = fopen(jFil_2,'r');
    
  if fid_1~=-1 & fid_2==-1
    %fprintf('Found jack sheet in %s\n\n',jFil_1)
    fid = fid_1;
  elseif fid_1==-1 & fid_2~=-1
    %fprintf('Found jack sheet in %s\n\n',jFil_2)
    fid = fid_2;
  elseif fid_1==-1 & fid_2==-1
    fprintf('Found no jack sheet\n\n')
    num = [];
    nam = [];
    nam_stripped=[];
    return
  else
    error('CHECK THIS OUT....I HAVE NO IDEA WHAT THIS CONDITION MEANS')
  end
    
  X   = textscan(fid,'%f%s');
  num = X{1};
  nam = X{2};
  clear X
    
  Nelec       = size(nam,1);
  elecJackTag = cell(Nelec,1);
  elecJackNum = nan(Nelec,1);
  for k=1:Nelec
    thisElecNam     = nam{k};
    elecJackTag{k}  = thisElecNam(regexp(thisElecNam,'\D'));
    elecJackNum_tmp = str2double(thisElecNam(regexp(thisElecNam,'\d')));
    elecJackNum(k)=elecJackNum_tmp;
  end
  nam_stripped = elecJackTag; 

  
  
function c1 = merge_colortables(c1, c2)
c1.struct_names(end+1:end+length(c2.struct_names)) = c2.struct_names;
c1.table(end+1:end+length(c2.struct_names),:) = c2.table;