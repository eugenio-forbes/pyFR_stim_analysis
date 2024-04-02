function [talStruct, hasMom] = getMonopolarTalStructureFromVoxInfo(subj,doTalDaemon,source_subj)
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

if ~exist('source_subj','var')||isempty(source_subj)
  source_subj = subj;
end

hasMom=true;

% directories
subjDir = fullfile(MOUNT_DIR,'/data/eeg',subj);
subjDir_source =  fullfile(MOUNT_DIR,'/data/eeg',source_subj);
docDir  = fullfile(subjDir,'docs');
talDir  = fullfile(subjDir,'tal');
talDir_source = fullfile(subjDir_source,'tal');
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

fid = fopen(fullfile(talDir_source,voxMotherFileName));
X = textscan(fid,'%s%d%d%d%s%d%d');
elecNames = X{1};
elecTypes = X{5};
fclose(fid);

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

% get jacksheet
[elecNum elecNam elecNam_stripped elecNam_stripped_num] = ...
    getTheContentsOFTheJackFile_local(docDir,jackFileName,jackFileName_2);


% now make the events structure
talStruct.subject   = [];
talStruct.channel   = [];
talStruct.tagName   = [];
talStruct.grpName   = [];
talStruct.x         = [];
talStruct.y         = [];
talStruct.z         = [];
talStruct.Loc1      = [];
talStruct.Loc2      = [];
talStruct.Loc3      = [];
talStruct.Loc4      = [];
talStruct.Loc5      = [];
talStruct.Loc6      = [];
talStruct.Montage   = [];
talStruct.eNames    = [];
talStruct.eType     = [];


talStruct.avgSurf = struct();
talStruct.avgSurf.x = [];
talStruct.avgSurf.y = [];
talStruct.avgSurf.z = [];
talStruct.avgSurf.anatRegion = [];
talStruct.avgSurf.x_snap = [];
talStruct.avgSurf.y_snap = [];
talStruct.avgSurf.z_snap = [];
talStruct.avgSurf.anatRegion_snap = [];
talStruct.avgSurf.distance_snap = [];
talStruct.avgSurf.x_eSnap = [];
talStruct.avgSurf.y_eSnap = [];
talStruct.avgSurf.z_eSnap = [];
talStruct.avgSurf.anatRegion_eSnap =[];
talStruct.avgSurf.distance_eSnap = [];
talStruct.avgSurf.path2surfL = [];
talStruct.avgSurf.path2surfR = [];

talStruct.indivSurf = struct();
talStruct.indivSurf.x = [];
talStruct.indivSurf.y = [];
talStruct.indivSurf.z = [];
talStruct.indivSurf.anatRegion =[];
talStruct.indivSurf.x_snap = [];
talStruct.indivSurf.y_snap = [];
talStruct.indivSurf.z_snap = [];
talStruct.indivSurf.anatRegion_snap = [];
talStruct.indivSurf.distance_snap = [];
talStruct.indivSurf.x_eSnap = [];
talStruct.indivSurf.y_eSnap = [];
talStruct.indivSurf.z_eSnap = [];
talStruct.indivSurf.anatRegion_eSnap =[];
talStruct.indivSurf.distance_eSnap = [];
talStruct.indivSurf.path2surfL = [];
talStruct.indivSurf.path2surfR = [];

path2indivSurfL = fullfile(fsSubjDir,'surf/lh.pial');
path2indivSurfR = fullfile(fsSubjDir,'surf/rh.pial');

path2avgSurfL = fullfile(fsAvgDir,'surf/lh.pial');
path2avgSurfR = fullfile(fsAvgDir,'surf/rh.pial');

fprintf('adding the talairach labels to the electrodes:\n')
ticker   = 0;
tick_inc = 5;
fprintf('   Wait please: ')

if hasRaw
    inMom_elecs = elecNum(ismember(elecNum, talElecs));
elseif hasAvg
    inMom_elecs = elecNum(ismember(elecNum, avgSurfElecs));
elseif hasIndiv
    inMom_elecs = elecNum(ismember(elecNum, indivSurfElecs));
end

leads = load(leadsFileName,'-ascii');
inMom_elecs = inMom_elecs(ismember(inMom_elecs, leads));

for k=1:length(inMom_elecs)
  if k./length(inMom_elecs)*100>=ticker
    fprintf(' %d%%', ticker)
    ticker=ticker+tick_inc;
  end
  thisElec  = inMom_elecs(k);

  if hasRaw
      thisXYZ   = double(talXYZ(thisElec==talElecs,:));

  else
      thisXYZ = nan(1,3);
  end
     
  if hasAvg
      thisAvgSurfXYZ   = double(avgSurfXYZ(thisElec==avgSurfElecs,:));
      
      if hasAvgSnap
          thisAvgSurfSnapXYZ   = double(avgSurfSnapXYZ(thisElec==avgSurfSnapElecs,:));
          avgDistance_snap = sqrt(sum((thisAvgSurfXYZ-thisAvgSurfSnapXYZ).^2));
      else
          thisAvgSurfSnapXYZ = nan(1,3);
          avgDistance_snap = nan;
      end
      
      if hasAvgEsnap
          thisAvgSurfEsnapXYZ   = double(avgSurfEsnapXYZ(thisElec==avgSurfEsnapElecs,:));
          avgDistance_eSnap = sqrt(sum((thisAvgSurfXYZ-thisAvgSurfEsnapXYZ).^2));
      else
          thisAvgSurfEsnapXYZ = nan(1,3);
          avgDistance_eSnap = nan;
      end
  else
      thisAvgSurfXYZ = nan(1,3);
      thisAvgSurfSnapXYZ = nan(1,3);
      thisAvgSurfEsnapXYZ = nan(1,3);
      avgDistance_snap = nan;
      avgDistance_eSnap = nan;
  end

  if hasIndiv
      thisIndivSurfXYZ   = double(indivSurfXYZ(thisElec==indivSurfElecs,:));

      if hasIndivSnap
          thisIndivSurfSnapXYZ   = double(indivSurfSnapXYZ(thisElec==indivSurfSnapElecs,:));
          indivDistance_snap = sqrt(sum((thisIndivSurfXYZ-thisIndivSurfSnapXYZ).^2));
      else
          thisIndivSurfSnapXYZ = nan(1,3);
          indivDistance_snap = nan;
      end
      
      if hasIndivEsnap
          thisIndivSurfEsnapXYZ   = double(indivSurfEsnapXYZ(thisElec==indivSurfEsnapElecs,:));
          indivDistance_eSnap = sqrt(sum((thisIndivSurfXYZ-thisIndivSurfEsnapXYZ).^2));
      else
          thisIndivSurfEsnapXYZ = nan(1,3);
          indivDistance_eSnap = nan;
      end
  else
     thisIndivSurfXYZ = nan(1,3);
     thisIndivSurfSnapXYZ = nan(1,3);
     thisIndivSurfEsnapXYZ = nan(1,3);
     indivDistance_snap = nan;
     indivDistance_eSnap = nan;
  end


  talStruct(k).subject  = subj;
  talStruct(k).channel  = thisElec;
  talStruct(k).tagName  = elecNam{thisElec==elecNum,1};
  talStruct(k).grpName = elecNam_stripped{thisElec==elecNum};
  talStruct(k).eName   =  num2str(thisElec);
  
  talStruct(k).x        = thisXYZ(1);
  talStruct(k).y        = thisXYZ(2);
  talStruct(k).z        = thisXYZ(3);

  talStruct(k).avgSurf.x = thisAvgSurfXYZ(1);
  talStruct(k).avgSurf.y = thisAvgSurfXYZ(2);
  talStruct(k).avgSurf.z = thisAvgSurfXYZ(3);
  talStruct(k).avgSurf.anatRegion = [];
  
  talStruct(k).avgSurf.x_snap = thisAvgSurfSnapXYZ(1);
  talStruct(k).avgSurf.y_snap = thisAvgSurfSnapXYZ(2);
  talStruct(k).avgSurf.z_snap = thisAvgSurfSnapXYZ(3);
  talStruct(k).avgSurf.anatRegion_snap = [];
  talStruct(k).avgSurf.distance_snap = avgDistance_snap;
  
  talStruct(k).avgSurf.x_eSnap = thisAvgSurfEsnapXYZ(1);
  talStruct(k).avgSurf.y_eSnap = thisAvgSurfEsnapXYZ(2);
  talStruct(k).avgSurf.z_eSnap = thisAvgSurfEsnapXYZ(3);
  talStruct(k).avgSurf.anatRegion_eSnap = [];
  talStruct(k).avgSurf.distance_eSnap = avgDistance_eSnap;
  
  talStruct(k).avgSurf.path2surfL = path2avgSurfL;
  talStruct(k).avgSurf.path2surfR = path2avgSurfR;
  
  talStruct(k).indivSurf.x = thisIndivSurfXYZ(1);
  talStruct(k).indivSurf.y = thisIndivSurfXYZ(2);
  talStruct(k).indivSurf.z = thisIndivSurfXYZ(3);
  talStruct(k).indivSurf.anatRegion = [];
  
  talStruct(k).indivSurf.x_snap = thisIndivSurfSnapXYZ(1);
  talStruct(k).indivSurf.y_snap = thisIndivSurfSnapXYZ(2);
  talStruct(k).indivSurf.z_snap = thisIndivSurfSnapXYZ(3);
  talStruct(k).indivSurf.anatRegion_snap =[];
  talStruct(k).indivSurf.distance_snap = indivDistance_snap;
    
  talStruct(k).indivSurf.x_eSnap = thisIndivSurfEsnapXYZ(1);
  talStruct(k).indivSurf.y_eSnap = thisIndivSurfEsnapXYZ(2);
  talStruct(k).indivSurf.z_eSnap = thisIndivSurfEsnapXYZ(3);
  talStruct(k).indivSurf.anatRegion_eSnap = [];
  talStruct(k).indivSurf.distance_eSnap = indivDistance_eSnap;
  
  talStruct(k).indivSurf.path2surfL = path2indivSurfL;
  talStruct(k).indivSurf.path2surfR = path2indivSurfR;
    
  if doTalDaemon
    talStruct(k)          = tal2Region(talStruct(k),false);
  end
  
  talStruct(k).eType =  elecTypes{strcmp(elecNames,elecNam{k})};
  
  if strcmp(talStruct(k).eType,'D')
    talStruct(k).Montage = 'hipp';
  else
    if thisXYZ(1)<0
      talStruct(k).Montage = 'lsag';
    else
      talStruct(k).Montage = 'rsag';
    end
  end
end

% Now we have to snap the avgSurf and indivSurf elecs of all types to
% surface to get the anatLabels associated
for i=1:2
    if (i==1 && hasAvg) || (i==2 && hasIndiv)
        if i==1
            structType = 'avgSurf';
            typeStruct = [talStruct.avgSurf];
            fsTypeDir = fsAvgDir;
        else
            structType = 'indivSurf';
            typeStruct = [talStruct.indivSurf];
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
            talStruct(i).(structType).anatRegion = anatLabels{i};
            talStruct(i).(structType).anatRegion_snap = anatLabels_snap{i};
            talStruct(i).(structType).anatRegion_eSnap = anatLabels_eSnap{i};
        end
    end
end

    
fprintf('\n\n')


% %--------------------------------------------------------------------
% function [allnames_tmp] = pairs2names_local(allpairs_tmp,thisGrpEls,thisGrpNam);
% allnames_tmp = {};
% for i = 1:size(allpairs_tmp,1)
%     for j = 1:size(allpairs_tmp,2)
%         allnames_tmp(i,j) = thisGrpNam(allpairs_tmp(i,j) == thisGrpEls);
%     end
% end
% 
% 
% 
% %---------------------------------------------------------------------
% function  allpairs = make_a_bipolar_motage(startval,numElecs_dim1,numElecs_dim2)
%     
%   % this is all Aswin G. Ramayya.  This occurred on June 14, 2013.
%   % Ashwin is completelty responsible for this code and miustakes
%   % therein.  John F. Burke is not reponsible for any errors,
%   % mistakes, or general bad things that come out of this code. But
%   % John Burke checked thouroughly and has placed his personal
%   % stamp of approval, thus accepting all responsibility (but no
%   % credit). Ashwin would like to add, spontaneously, that he is a
%   % tool. John agrees.
%   allpairs = [];
%   for j = 1:numElecs_dim2
%     for i = 1:numElecs_dim1
%       if i < numElecs_dim1
% 	allpairs(end+1,:) = [i-1 i] + startval + numElecs_dim1*((j-1));
%       end
%       if j < numElecs_dim2
%         allpairs(end+1,:) = [i-1 i+numElecs_dim1-1] + startval + numElecs_dim1*((j-1));
%       end
%     end
%   end
%   
        
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

% %---------------------------------------------------------------------
% function out = getTheContentsOFVoxMom_local(d,f);
%   fil = fullfile(d,f);
%   fid = fopen(fil,'r');
%   if fid==-1
%     out = [];
%     fprintf('%s does not exist...EXITING\n\n',f)
%     return
%   end  
%   X = textscan(fid,'%s%d%d%d%s%s','delimiter','\t');
%   fclose(fid);
%   gridNames    = X{1};
%   elecXYZ      = [X{2} X{3} X{4}];
%   stripOrDepth = X{5};
%   geometry_tmp = X{6};
%   out          = [];
%   for k=1:size(geometry_tmp,1)
%     out(k).name = gridNames{k};
%     out(k).voxP = elecXYZ(k,:);
%     out(k).type = stripOrDepth{k};
%     out(k).geom = sscanf(geometry_tmp{k},'%d%d')';
%   end
%   
% %--------------------------------------------------
% function out = getTheContentsOFTheElecs_local(d,f);
%   elecFile = fullfile(d,f);
%   if ~exist(elecFile,'file')  
%     out = [];
%     fprintf('%s does not exist...EXITING\n\n',f)
%     return
%   end  
%   run(elecFile);
%   if exist('r','var')
%     out = r;
%   elseif exist('grids','var')
%     out = grids;
%   else
%     out = [];
%     fprintf('The contects of electrodes.m not right...EXITING\n\n')
%     return
%   end
% 
% %------------------------------------------------------------  
% function [num nam nam_stripped nam_strimmed_num] = getTheContentsOFTheJackFile_local(dDir,j1,j2);
%   jFil_1 = fullfile(dDir,j1);
%   jFil_2 = fullfile(dDir,j2);  
%   fid_1  = fopen(jFil_1,'r');
%   fid_2  = fopen(jFil_2,'r');
%     
%   if fid_1~=-1 & fid_2==-1
%     %fprintf('Found jack sheet in %s\n\n',jFil_1)
%     fid = fid_1;
%   elseif fid_1==-1 & fid_2~=-1
%     %fprintf('Found jack sheet in %s\n\n',jFil_2)
%     fid = fid_2;
%   elseif fid_1==-1 & fid_2==-1
%     fprintf('Found no jack sheet\n\n')
%     num = [];
%     nam = [];
%     nam_stripped=[];
%     nam_strimmed_num=[];
%     return
%   else
%     error('CHECK THIS OUT....I HAVE NO IDEA WHAT THIS CONDITION MEANS')
%   end
%     
%   X   = textscan(fid,'%f%s');
%   num = X{1};
%   nam = X{2};
%   clear X
%     
%   Nelec       = size(nam,1);
%   elecJackTag = cell(Nelec,1);
%   elecJackNum = nan(Nelec,1);
%   for k=1:Nelec
%     thisElecNam     = nam{k};
%     elecJackTag{k}  = thisElecNam(regexp(thisElecNam,'\D'));
%     elecJackNum_tmp = str2double(thisElecNam(regexp(thisElecNam,'\d')));
%     elecJackNum(k)=elecJackNum_tmp;
%   end
%   nam_stripped = elecJackTag; 
%   nam_strimmed_num = elecJackNum;
  


function c1 = merge_colortables(c1, c2)
c1.struct_names(end+1:end+length(c2.struct_names)) = c2.struct_names;
c1.table(end+1:end+length(c2.struct_names),:) = c2.table;


function [num, nam, nam_stripped, nam_strimmed_num] = getTheContentsOFTheJackFile_local(dDir,j1,j2)
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
    nam_strimmed_num=[];
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
  nam_strimmed_num = elecJackNum;