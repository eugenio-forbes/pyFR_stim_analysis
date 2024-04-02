function convertMotherfromRAW(subj,source_subj, doPlot, force_redo)
% CONVERTMOTHERTOVOXCOORDS(subj, source_subj, [doPlot, force_redo])
%
% FUNCTION:
%   convertMotherToVoxCoords.m
%
% DESCRIPTION:
%   name is self explanatory
%
% INPUTS:
%   subjDir............ this is the locatioin of the subjest's directory
%   vox_mother_loc..... this is the location of the vox mother file
%   jackLoc............ this is the location of the jacksheet
%   subj............... this is the name of the subject
%   source_subj........
%   source_subjDir.....
%   username........... 'jfburke' [your rhino account] use if you
%                       are on a mac
%

if ismac
  fprintf('/n/nTHIS ONLY RUNS ON RHINO/n/n')
  return
else
  mountDir='';
end

if ~exist('doPlot','var') || isempty(doPlot)
  doPlot = 0;
end
if ~exist('source_subj','var')||isempty(source_subj)
  source_subj=subj;
end
if ~exist('force_redo','var') || isempty(force_redo)
    force_redo = false;
end

dataDir = fullfile(mountDir,'/data/eeg');      
subjDir = fullfile(dataDir,subj);
fsDir = fullfile(mountDir,'/data/eeg/freesurfer/subjects');
fsSubjDir = fullfile(fsDir,source_subj);

% dircetories
talDir         = fullfile(subjDir,'tal');
docDir         = fullfile(subjDir,'docs');

if ~exist('source_subjDir','var')||isempty(source_subjDir)
  source_subjDir = subjDir;
end
source_subjDir = fullfile(dataDir,source_subj);

source_talDir  = fullfile(source_subjDir,'tal');
source_docDir  = fullfile(source_subjDir,'doc');

% this is the output file name
voxActualName  = 'VOX_coords.txt';
voxActualFile  =  fullfile(talDir,voxActualName);
voxMotherName  = 'VOX_coords_mother.txt';

% prepare the source and the target files
mni_file       = 'RAW_coords.txt.mni';
tal_file       = 'RAW_coords.txt';
avgSurf_file   = 'RAW_coords_avgSurf.txt';
avgSurfSnap_file = 'RAW_coords_avgSurf_snap.txt';
avgSurfEsnap_file = 'RAW_coords_avgSurf_eSnap.txt';
indivSurf_file = 'RAW_coords_indivSurf.txt';
indivSurfSnap_file = 'RAW_coords_indivSurf_snap.txt';
indivSurfEsnap_file = 'RAW_coords_indivSurf_eSnap.txt';
vox_coords  = fullfile(talDir,voxActualName);
imageroot   = fullfile(source_talDir,'images/combined');
ct_combined = sprintf('%s_CT_combined',source_subj);
ct2mr       = sprintf('%s_CT2MRI',source_subj);
mr_brain    = sprintf('%s_MRI_brain',source_subj);
mr2standard = sprintf('%s_MRI2standard',source_subj);

% use the data in VOX_coords to create a RAW_coords_avgSurf.txt file, which
% contains the coordinates from the electrodes named in the jacksheet in
% the space of the average surface
%reg_avg_file = fullfile(source_talDir, '/images/combined/bbreg_test.dat');
if ~exist(voxActualFile,'file') || force_redo 
    % get the contents of the jacksheet
    jackLoc           = fullfile(docDir,'jacksheet.txt');
    fid_jac           = fopen(jackLoc,'r');
    if fid_jac==-1 
        jackLoc    = fullfile(docDir,'jack_sheet.txt');
        fid_jac    = fopen(jackLoc,'r');
        if fid_jac==-1
            error('convertMotherToVoxCoords:NoJackSheet',...
                'jacksheet does not exist for subject %s',subj);
        end
    end
    jackSheetContents = readTheJackSheet_local(fid_jac);
    fclose(fid_jac);
    
    % get the contents of the of the mother VOX_coords.txt file
    vox_mother_loc  = fullfile(source_talDir,voxMotherName);
    fid_mom         = fopen(vox_mother_loc,'r');
    if fid_mom==-1
        vox2voxMother(subj)
        fid_mom = fopen(vox_mother_loc,'r');
        if fid_mom==-1
            error('convertMotherToVoxCoords:NoVoxCoordsMother',...
                'VOX_coords_mother and VOX_coords do not exist for subject %s',subj);
        end
    end
    voxFileContents = readTheMotherVox_local(fid_mom);
    momTagFullNam   = voxFileContents{1};
    momTagNam       = regexprep(voxFileContents{1},'\d','');
    momTagNum       = regexprep(voxFileContents{1},'\D','');
    momTagXYZ       = [voxFileContents{2} voxFileContents{3} voxFileContents{4}];
    fclose(fid_mom);
    
    % now write the vox file.
    fprintf('\n')
    fprintf('Making the VOX_coords.txt channel: ')
    if ~exist(voxActualFile,'file') || force_redo
        fid_vox = fopen(voxActualFile,'w+');
        if fid_vox==-1;error('fkjndkvjbdkhjfbv');end
        makeTheActualMom(fid_vox,jackSheetContents,voxFileContents);
        fprintf('done\n')
    else
        fprintf('already made\n')
    end
end
reg_avg_file = fullfile(dataDir, subj, 'tal/images/combined/reg_avg.lta');
cd(talDir)
if (~exist(avgSurf_file,'file') || force_redo) && ... 
        (exist(reg_avg_file,'file'))
    
    tal_coords = vox2tal(source_subj,vox_coords, avgSurf_file, reg_avg_file);
    fprintf('    I just made %s\n', avgSurf_file);
    

    isDepth = getIsDepth_local(fullfile(talDir,voxMotherName), ...
        fullfile(docDir, 'jacksheet.txt'),...
        fullfile(vox_coords));
    
    
    % Now snap this to surface using the snap method
    snapElectrodesToSurface(subj, avgSurf_file, isDepth, true, avgSurfSnap_file, [], false, [], [], [], 20);
    fprintf('    I just made %s\n', avgSurfSnap_file);
    % and then again using the spring method
    
    
    energyMinimization(subj, avgSurf_file, [], isDepth, true, false, avgSurfEsnap_file);
    fprintf('    I just made %s\n',avgSurfEsnap_file);
else
    if exist(avgSurf_file,'file') && ~force_redo
        fprintf('    %s already exists\n', avgSurf_file);
    else
        fprintf('    %s does not exist\n',reg_avg_file);
    end
end


% Do the same for RAW_coords_indivSurf.txt
reg_file = fullfile(dataDir, subj, 'tal/images/combined/reg.lta');
cd(talDir)
if (~exist(indivSurf_file,'file') || force_redo)  && exist(reg_file,'file')
    
    tal_coords = vox2tal(source_subj, vox_coords, indivSurf_file, reg_file);
    fprintf('     I just made %s\n', indivSurf_file);
        
    isDepth = getIsDepth_local(fullfile(talDir,voxMotherName),...
        fullfile(docDir, 'jacksheet.txt'),...
        fullfile(vox_coords));

    % Now snap this to surface using the snap method
    snapElectrodesToSurface(subj, indivSurf_file, isDepth, false, indivSurfSnap_file, [], false, [], [], [], 20);
    fprintf('     I just made %s\n',indivSurfSnap_file);
    % and then again using the spring method


    energyMinimization(subj, indivSurf_file, [], isDepth, false, false, indivSurfEsnap_file);
    fprintf('    I just made %s\n',indivSurfEsnap_file);
else
    if exist(indivSurf_file,'file') && ~force_redo
        fprintf('    %s already exists\n',indivSurf_file);
    else
        fprintf('    %s does not exist\n',reg_file);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make the bipolar file from 
bpTalFileName = sprintf('%s_talLocs_database_bipol.mat',subj);
talFileName = sprintf('%s_talLocs_database_monopol.mat',subj);
bpTalFile     = fullfile(talDir,bpTalFileName);
talFile = fullfile(talDir, talFileName);
if exist(bpTalFile,'file')
    system(sprintf('mv %s %s.old',bpTalFile, bpTalFile));
end
[bpTalStruct, hasMom] = getBipolarTalStructureFromVoxInfo(subj,[],source_subj);
if ~hasMom %this should only happen if its an old subject, then we will use the old bp_mat method
    error('NO MOM')
end
cd(talDir)
save(bpTalFile,'bpTalStruct');
%keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check the registration
fprintf('Checking the registration using the bipolar montage:\n')

if  doPlot
    %% plot the histogram
    figure(1);clf
    hist(cat(1,bpTalStruct.bpDistance),length(bpTalStruct))
    xlabel('distance between bp pairs')
    ylabel('number of bp pairs')
    
    
    %% plot electrodes on the brain
    cd(talDir)
    if exist('RAW_coords.txt','file')
    raw_xyz = get_tal_coords('RAW_coords.txt');
    figure;
    tal3d(get_tal_coords('RAW_coords.txt'),[-90 0],[1 0 0],'num',[6], 0,0,1)
    
    %%
    figure
    plot3(raw_xyz(:,2),raw_xyz(:,3),raw_xyz(:,4),'.')
    end

    raw_xyz_surf = get_tal_coords('RAW_coords_avgSurf.txt');
    figure; 
    tal3d(get_tal_coords('RAW_coords_avgSurf.txt'),[-90 0],[0 0 1],'num',10, 0,0,1)
    %% 
    figure; hold all 
    plot3_wrapper(raw_xyz(:,2:4),10,'r')
    plot3_wrapper(raw_xyz_surf(:,2:4),10,'b')


    raw_xyz_surf = get_tal_coords('RAW_coords_avgSurf.txt');
    figure;
    tal3d(get_tal_coords('RAW_coords_avgSurf.txt'),[-90 0],[0 0 1],'num',10, 0,0,1)
    %%
    figure; hold all
    plot3_wrapper(raw_xyz(:,2:4),10,'r')
    plot3_wrapper(raw_xyz_surf(:,2:4),10,'b')
end

%% plot the histogram
% figure(2)
% clf;
% ax=axes;
% for k=1:length(bpTalStruct)
%   line('Parent',ax,'XData',bpTalStruct(k).distance,'YData',0,...
%        'Marker','o');
%   t = text(bpTalStruct(k).distance,0.1,sprintf('%d-%d',bpTalStruct(k).channel));
%   %set(t,'HorizontalAlignment','center');
% end
% xlabel('distance between bp pairs')
% grid on

% xyz = xyz_TAL;
% figure(3);clf;
% res = tal3d([], [-90 0], [1 0 0], 'num', [6], 0, 1); 
% hold on
% set(res.hBrain, 'facealpha', .5);
% fprintf('   done\n\n')

function isDepth = getIsDepth_local(vox_mother_loc, jackLoc, voxLoc)
    fid_mom         = fopen(vox_mother_loc,'r');
    if fid_mom==-1
        error('convertMotherToVoxCoords:NoVoxCoordsMother',...
            'VOX_coords_mother and VOX_coords do not exist for subject %s',subj);
    end
    voxFileContents = readTheMotherVox_local(fid_mom);
    isDepth_mom = strcmp(voxFileContents{5},'D');
    tagNames_mom = voxFileContents{1};
    fclose(fid_mom);
    
    % get the contents of the jacksheet
    fid_jac           = fopen(jackLoc,'r');
    if fid_jac==-1
        jackLoc    = fullfile(docDir,'jack_sheet.txt');
        fid_jac    = fopen(jackLoc,'r');
        if fid_jac==-1
            error('convertMotherToVoxCoords:NoJackSheet',...
                'jacksheet and VOX_coords do not exist for subject %s',subj);
        end
    end
    jackSheetContents = readTheJackSheet_local(fid_jac);
    fclose(fid_jac);
    
    fid_vox = fopen(voxLoc);
    voxContents =  readVox_local(fid_vox);
    tagNums_vox = voxContents{1};
    isDepth = false(length(voxContents(:,1)),1);
    tagNames_jack = cellfun(@(x,y)[x,y],jackSheetContents(:,2), jackSheetContents(:,3),...
        'uniformoutput',false);
    tagNums_jack = cell2mat(jackSheetContents(:,1));
    for i=1:length(tagNums_vox)
        thisName = tagNames_jack{tagNums_jack==tagNums_vox(i)};
        isDepth(i) = isDepth_mom(strcmp(thisName,tagNames_mom));
    end


%------------------------------------------------
function makeTheActualMom(fid_vox,jackSheetContents,voxFileContents);
  momTagFullNam   = voxFileContents{1};
  momTagNam       = regexprep(voxFileContents{1},'\d','');
  momTagNum       = regexprep(voxFileContents{1},'\D','');
  momTagXYZ       = [voxFileContents{2} voxFileContents{3} voxFileContents{4}];
  
  % loop the jacksheet and pull the VOX coord for that electrode
  allTags = {};
  for k = 1:size(jackSheetContents,1);
    thisJackElecNum = jackSheetContents{k,1};
    thisJackTagNam  = jackSheetContents{k,2}; 
    thisJackTagNum  = jackSheetContents{k,3};
    
    tagNamMatches = strcmp(momTagNam,thisJackTagNam);
    tagNumMatches = strcmp(momTagNum,thisJackTagNum);  
    tagMatchIdx   = find(tagNamMatches & tagNumMatches);
    if isempty(tagMatchIdx);
      fprintf('   WARNING!!! %s%s NOT FOUND...NOT WRITING IN VOX COORDS\n',thisJackTagNam,thisJackTagNum)
      continue
    end
    if length(tagMatchIdx)>1;
      error(sprintf('%s%s FOUND TOO MANY TAGS IN VOX MOM FOR THIS ELECTRODE',...
		    thisJackTagNam,thisJackTagNum))
    end
    
    % get the electrode number
    thisElecXYZ = momTagXYZ(tagMatchIdx,:);
    
    % now just print the VOX file
    fprintf(fid_vox,'%d\t%d\t%d\t%d\n',thisJackElecNum, ...
	    thisElecXYZ);
    allTags = cat(1,allTags,[thisJackTagNam thisJackTagNum]);
  end
  pause(.1)

function voxOUT = readVox_local(fid)
    contentStr = '%d%d%d%d';
    voxOUT = textscan(fid, contentStr);
%-------------------------------------------------
function voxOUT = readTheMotherVox_local(fid)
 contentStr = '%s%d%d%d%s%s'; 
 numStuff   = length(regexp(contentStr,'\%')); 
 voxOUT     = textscan(fid,contentStr,'delimiter','\t');
 
%------------------------------------------------
function data = readTheJackSheet_local(fid);
  
  contentStr = '%d%s';
  data       = cell(1,3);
  count      = 0; 
  while true
    thisStr = fgetl(fid);
    if ~ischar(thisStr) || isempty(thisStr)
        break;
    end
    
    % scan the line for the electrode information
    thisLine    = sscanf(thisStr,contentStr);    
    
    % get the electrode number
    thisElecNum = thisLine(1);
    
    % now get the tag name
    thisElecNam = char(thisLine(2:end))';
    numberInds  = regexp(thisElecNam,'\d');
    if length(numberInds)<1;      
      disp(sprintf('\n\n\tWARNING: no numbers fpr %s\n\n',thisElecNam));
      thisTagNam = thisElecNam;
      thisTagNum = 1;
    else
      thisTagNam = thisElecNam(1:numberInds(1)-1);
      thisTagNum = thisElecNam(numberInds(1):end);
    end
    count         = count + 1;
    data{count,1} = thisElecNum;
    data{count,2} = thisTagNam;
    data{count,3} = thisTagNum;    
  end  
  
%------------------------------------------------  
function fid = open_a_file_local(f,d,openStr)

  fil = fullfile(d,f);
  fid = fopen(fil,openStr);
  if fid==-1
    fprintf('  ERROR: % 30.30s .... CANNOT OPEN in .... %s\n',f,d)
    fid=[];
    return
  end

return

%------------------------------------------------
function tal_bp_check_local(subj, subjDir, docsDir)
%
% FUNCTION:
%   tal_bp_check.m
% 
% DESCRIPTION:
%   A way to check that a newly talairached patient has bipolar
%   pairs that are near 10 mm from each other.
%
% INPUT:
%   subj...... 'TJ061'
%
% OUTPUT:
%
% NOTES:
%   (1) written by jfburke 06/2013 (john.fred.burke@gmail.com)
%

% am I mounting RHINO?
% if ismac
%   MOUNT_DIR = '/Volumes/RHINO_root';
% else
%   MOUNT_DIR = '';
% end

% subjDir = fullfile(MOUNT_DIR,'/data/eeg/',subj);
%docsDir = fullfile(subjDir,'docs');
talDir = fullfile(subjDir,'tal');

bp_els  = get_this_file_local(docsDir,'bp_manual.txt','%d%d%s');
bp_els  = [bp_els{1} bp_els{2}];
    
mni_xyz = get_this_file_local(talDir, 'RAW_coords.txt.mni','%d%f%f%f');
allEls  = mni_xyz{1};
xyz     = [mni_xyz{2} mni_xyz{3} mni_xyz{4}];

nBP              = size(bp_els,1);
distBetweenBPels = nan(nBP,1);
for k=1:nBP
  thisBP_1     = bp_els(k,1);
  thisBP_1_ind = find(ismember(allEls,thisBP_1));
  if length(thisBP_1_ind)~=1;error('bad in 1');end
  thisBP_1_xyz = xyz(thisBP_1_ind,:);
  
  thisBP_2 = bp_els(k,2);
  thisBP_2_ind = find(ismember(allEls,thisBP_2));
  if length(thisBP_2_ind)~=1;error('bad in 2');end
  thisBP_2_xyz = xyz(thisBP_2_ind,:);
  
  thisBP_xyz  = bipolarHalfDistance_local(thisBP_1_xyz,thisBP_2_xyz);
  thisBP_dist = sqrt(sum((thisBP_1_xyz-thisBP_2_xyz).^2));

  
  distBetweenBPels(k) = thisBP_dist;
  
end

% plot the histogram
figure(1)
clf
hist(distBetweenBPels,nBP)
% formatPic(gca)
xlabel('distance between bp pairs')
ylabel('number of bp pairs')


% plot the histogram
figure(2)
clf;
ax=axes;
for k=1:nBP
  line('Parent',ax,'XData',distBetweenBPels(k),'YData',0,...
       'Marker','o');
  t = text(distBetweenBPels(k),0.1,sprintf('%d-%d',bp_els(k,:)));
  %set(t,'HorizontalAlignment','center');
end
% formatPic(gca)
xlabel('distance between bp pairs')
grid on

function out = get_this_file_local(d,f,thing)
  fil = fullfile(d,f);
  fid = fopen(fil);  
  if fid==-1
   fprintf('\n\nERROR: %s does not exist in %s\n\n',f,d)
   error('done')
  end  
  out = textscan(fid,thing);
  fclose(fid);
  
  
%------------------------------------------------------------------
function xyz_out = bipolarHalfDistance_local(XYZ1,XYZ2)
  tmp = -XYZ1+XYZ2;
  halftmp = tmp./2;
  xyz_out = XYZ1 + halftmp;


% run the vox coord through the registration

% make the Talairach file, now
%[c,x,y,z] = textread(coordsfile,'%d%n%n%n');
%mni_coords = [x y z];
% convert the coords
%tal_coords = mni2tal(mni_coords);

% backup the old file
%if strcmp(coordsfile(end-3:end),'.mni')
%  backupfile = coordsfile;
%  fprintf('MNI backup file %s exists\n',backupfile);
%  coordsfile = coordsfile(1:end-4);
%else
%  backupfile = [coordsfile '.mni'];
%  fprintf('Backing up MNI file to %s\n',backupfile);
%  movefile(coordsfile,backupfile);
%5end

% write out new file
%outdata = [c tal_coords];
%fid = fopen(coordsfile,'w');
%fprintf(fid,'%d\t%g\t%g\t%g\n',outdata');
%fclose(fid);
%
%fprintf('Wrote new Tal file to %s\n',coordsfile);
