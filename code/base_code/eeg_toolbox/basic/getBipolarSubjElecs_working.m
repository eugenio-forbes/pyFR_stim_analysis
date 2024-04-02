function [subjTalEvents] = getBipolarSubjElecs_working(subj,BIPOLAR_BIT,excludeNonNeuralFlag,excludeEpilepsyFlag);  
%function [subjElec subjXYZ subjMont eNames tagNames lobes bas] = ...
%    getBipolarSubjElecs(subj,BIPOLAR_BIT,excludeNonNeuralFlag,excludeEpilepsyFlag);  
%
% FUNCTION:
%   getBipolarSubjElecs.m
% 
% DECRIPTION:
%   For one subject, gets electrodes for bipolar and regular montages.
%
% INPUTS:
%   subj............. 'UP011'
%   BIPOLAR_BIT...... true OR false
%   allTalEvents..... [optional] talairach events
%
% OUTPUTS:
%
% Notes:
% (1) Written by JFB and AGR 06/2013 based on original by JFB
% (2) Updated 06/24/2013 AGR
%       -outputs bpTalStruct
%       -saves 'bpTalStruct' not 'events'
%       -if tagNames is not a field it remakes it
% To-Do
%   ( ) make sure that it works for bipolflag == 0 even for new subjects
%   ( ) check to make sure the tag names are correct and can handle
%       missing electrodes
%   ( ) Test the flags to see of they work 
%   ( ) Make sure that getbpTalStruct_local has the same output as
%   getBipolarStructFromVoxInfo



% am I mounting RHINO?
if ismac
  MOUNT_DIR = '/Volumes/RHINO_root';
else
  MOUNT_DIR = '';
end

% defauts 
if ~exist('excludeNonNeuralFlag','var')||isempty(excludeNonNeuralFlag)
  excludeNonNeuralFlag=false;
end
if ~exist('excludeEpilepsyFlag','var')||isempty(excludeEpilepsyFlag)
  excludeEpilepsyFlag=false;
end

% HACK IT HERE FOR
% This is pretty bad.  The problem is that CH008b and UP003a
% are unlike all other subjects and this ocde requires a
% directory format like this: '/data/eeg/[subj]'.  So we have
% to do this hack
if strcmp(subj,'CH008b')      
  subj = 'CH008/CH008b';
elseif strcmp(subj,'UP003a')      
  subj = 'UP003/UP003a';
end

if ~exist('BIPOLAR_BIT','var')
    BIPOLAR_BIT = true;
end

% directories
dataDir = fullfile(MOUNT_DIR,'/data/eeg');
subjDir = fullfile(dataDir,subj);
docDir  = fullfile(subjDir,'docs'); 
talDir  = fullfile(subjDir,'tal/talFiles_old');  

% get this subjects info f% run the config and get the tal events
if ~BIPOLAR_BIT
  talFileName = sprintf('%s_talLocs_database_monopol.mat',subj);
  talFile     = fullfile(talDir,talFileName);
else
  
  talFileName = sprintf('%s_talLocs_database_bipol.mat',subj);
  talFile     = fullfile(talDir,talFileName);
end

% load it
subjTalEvents = load(talFile);
names = fieldnames(subjTalEvents);
if length(names)>1
 error('Too many files saved'); 
end
subjTalEvents = subjTalEvents.(names{1});



% remove any electrodes with no neural data. i.e. they are not in:
%   /data/eeg/[subj]/tal/leads.txt
if excludeNonNeuralFlag
  leadsDotTextFile = fullfile(talDir,'leads.txt');
  fid=fopen(leadsDotTextFile,'r');
  if fid==-1
    fprintf('\n\n\tWARNING!!!: There is no leads.txt file here:\n')    
    fprintf('                   cannot exclude non-neural leads\n\n.')    
else
    X         = textscan(fid,'%d'); 
    leads     = X{1};
    fclose(fid);
    els       = cat(1,subjTalEvents.channel);
    tossThese = sum(ismember(els,leads),2)~=size(els,2);
    subjTalEvents(tossThese)=[];
  end
end

% remove any electrodes over epilepsy. i.e. remove any electrodes
% in /data/eeg/[subj]/tal/bad_leads.txt
if excludeEpilepsyFlag
  leadsDotTextFile = fullfile(talDir,'bad_leads.txt');
  fid=fopen(leadsDotTextFile,'r');
  if fid==-1
    fprintf([['  \nWARNING!!!: There is no bad_leads.txt file here:'] ...
	     [ ' trying good_leads\n']])    
    leadsDotTextFile = fullfile(talDir,'good_leads.txt');
    fid=fopen(leadsDotTextFile,'r');
    if fid==-1
      fprintf(['  \nWARNING!!!: There is no good_leads.txt file here:'])    
    else
      X          = textscan(fid,'%d'); 
      gdleads    = X{1};
      fclose(fid);
      els        = cat(1,subjTalEvents.channel);
      tossThese  = sum(ismember(els,gdleads),2)~=size(els,2);
      subjTalEvents(tossThese)=[];
    end
  else
    X          = textscan(fid,'%d'); 
    badleads   = X{1};
    fclose(fid);
    els        = cat(1,subjTalEvents.channel);
    tossThese  = sum(ismember(els,badleads),2)>0;
    subjTalEvents(tossThese)=[];
  end
end

%------------------------------------------------------------
function subjXYZ = getManualBipolElecsXYZ_local(subjElec,subj,subjTalEv)

  allCHANNELS = [subjTalEv.channel]';
  x = [subjTalEv.x]';
  y = [subjTalEv.y]';
  z = [subjTalEv.z]';
  
  subjXYZ  = [];
  for k=1:size(subjElec,1)
    e1 =  subjElec(k,1);
    e2 =  subjElec(k,2);
    e1Ind = find(allCHANNELS==e1);
    e2Ind = find(allCHANNELS==e2);
    thisP1 = [x(e1Ind) y(e1Ind) z(e1Ind)];
    thisP2 = [x(e2Ind) y(e2Ind) z(e2Ind)];
    subjXYZ = cat(1,subjXYZ,bipolarHalfDistance_local(thisP1,thisP2));
  end

%------------------------------------------------------------
function [E M] = rmElecsWithNoNeuralData_local(subj,MOUNT_DIR,E_in,M_in)
  resDataRoot   = fullfile(MOUNT_DIR,'/data/eeg/',subj,'tal'); 
  elFile        = 'leads.txt';
  subj_el_File  = fullfile(resDataRoot,elFile);

  fid = fopen(subj_el_File,'r');
  if fid==-1
    error('leads.txt exist')
  end
  [foo]       = textscan(fid,'%f');
  gdEl        = double(foo{1});
  keepInd_all = ismember(E_in,gdEl);
  keepInd     = keepInd_all(:,1) & keepInd_all(:,2);
  E = E_in(keepInd,:);
  M = M_in(keepInd,:);

%------------------------------------------------------------
function [tagNames eNames] = getECoG_subjTagNames_local(subj,els,MOUNT_DIR);

  % what are bad tag names?
  if ~exist('badNames','var')||isempty(badNames)
    badNames = {'PENN','DRXL','EKG','REF'};
  end
  
  % this is a hack error catch by jfburke.  But, don't blame me,
  % blame the hack way the old jacksheets were recorded (04-2013).
  badSubj = {'BW022','BW023'};
  if ismember(subj,badSubj)
    tagNames=[];    
    return
  end
  
  subjDataDir = fullfile(MOUNT_DIR,'/data/eeg/',subj);
  jacksheet   = fullfile(subjDataDir,'docs','jacksheet.txt');
  jack_sheet  = fullfile(subjDataDir,'docs','jack_sheet.txt');
  
  % open the jacksheet
  fid = fopen(jacksheet);
  if fid==-1;
    %fprintf('No jacksheet in %s',subjDataDir)
    fid = fopen(jack_sheet);
    if fid==-1;
      tagNames=[];    
      return
    end
  end
  [X]          = textscan(fid,'%d%s%s');
  elecNums     = X{1};
  elecTags     = X{2};
  elecTags_HUP = X{3};
  fclose(fid);
  
  % sometimes there is a space that screws us up (HUP jacksheets)
  if length(unique(elecTags))<=1
    elecTags = regexprep(elecTags_HUP,'-REF','');
    elecTags = regexprep(elecTags,'_','');
  end
  
  % is it a bipolar?
  if size(els,2)==1 
    BIPOLAR_BIT=false;
  elseif size(els,2)==2
    BIPOLAR_BIT=true;
  else
    error('band second dimensio size of input electrodes')
  end
  
  % do I want bipolar?-->effects how I get elecs
  if ~BIPOLAR_BIT
    tagNames  = elecTags;  
    goodElecs = true(size(els));
    for k=1:size(els,1)
      eNames{k}         = sprintf('%d',els(k)); 
      thisTagName       = upper(tagNames{k});
      thisTagName_noNum = thisTagName(regexp(thisTagName,'\D'));
      goodElecs(k)      = ~ismember(thisTagName_noNum,badNames);
    end
    
    % get rid of the bad electrode names
    tagNames = tagNames(goodElecs);
  
  else
    
    % now get teh tag and file names    
    tagNames  = cell(size(els,1),1);
    eNames{k} = sprintf('%d-%d',els(k)); 
    for k=1:size(els,1)
      tagNames{k} = sprintf('%s-%s',...
			    elecTags{els(k,1)},...
			    elecTags{els(k,2)});
    end   
  end
  
  % clean the tg names
  if ~isempty(tagNames)
    for k=1:length(tagNames)
      tagNames{k} = regexprep(tagNames{k},'-REF','');
    end   
  end
  
  
function make_AUTO_bipolar_pairs_local(subj)
%
% FUNCTION:
%   make_AUTO_bipolar_pairs.m
% 
% DESCRIPTION:
%   derives bipolar montage automatically
%
% INPUT:
%   fil............. the file to write the bipolar montage to
%
% OUTPUT:
%   allpairs......... the bipolar montage
%
% NOTES:
%   (1) grids_struct: a semicolon-separated list of pairs whose first
%       element is a 1x3 matrix with elements: grid starting #, 
%       grid width, grid height, and whose second element is one of:
%       'left', 'right', 'hipp' or 'inf'
%   (2) written 04/2012 by Aaron Geller (aaron.s.geller@gmail.com)
% 

if ismac
  rootDir='/Volumes/RHINO_root';
else
  rootDir='';
end

grids_struct_FILENAME = 'grids_struct_for_bipolar.m';
fil_OUT_NAME          = 'bp_manual.txt';
docDir                = fullfile(rootDir,'/data/eeg/',subj,'docs');
grids_struct_FILE     = fullfile(docDir,grids_struct_FILENAME);
fil_OUT               = fullfile(docDir,fil_OUT_NAME);
if exist(fil_OUT,'file')
  fprintf('bipolar is made...EXITING\n')
  return
end
if ~exist(grids_struct_FILE,'file')
  fprintf('run mkGridStructForBipolarFromJacksheet.m manually.... EXITING\n')
  return
end
run(grids_struct_FILE)
fid = fopen(fil_OUT, 'w');
if fid==-1
  error('you don''t have permission')
end

for i=1:size(grids_struct,1)
  allpairs = [];
  grid_vec_cell = grids_struct(i,1);
  grid_vec = grid_vec_cell{1};
  startval = grid_vec(1);
  for j=1:grid_vec(2)
    if j<grid_vec(2)
      % along first dimension    
      allpairs(end+1,:) = [j-1 j] + startval;
    end
  end
  
  for j=1:grid_vec(2)
    for k=1:grid_vec(3)-1
      % go parallel along multiple rows
      if j<grid_vec(2)
        allpairs(end+1,:) = [k*grid_vec(2)-1 k*grid_vec(2)] ...
            + startval + j;
      end
    end
  end
  
  for j=1:grid_vec(2)
    for k=1:grid_vec(3)-1
      % go perpendicularly
      allpairs(end+1,:) = [(k-1)*grid_vec(2)-1 k*grid_vec(2)-1] ...
          + startval + j;
    end
  end
  output_mat_to_file(allpairs, fid, char(grids_struct(i,2)));
end
fclose(fid);

function output_mat_to_file(mat, fid, str)
for i=1:size(mat,1)
  thisline = sprintf('%d %d %s\n', mat(i,1), mat(i,2), str);
  fprintf(fid, thisline);
end



%% Graveyard
%   if ~exist(talFile,'file')
%     fprintf('\n')
%     fprintf('%s does not exist.\n',talFile)
%     fprintf('    Loading tal database and saving new....\n',talFile)
%     TALEVENTS=load(fullfile(dataDir,'tal', 'allTalLocs_GM.mat'));
%     allTalEvents=TALEVENTS.events;
%     clear TALEVENTS
%     events = allTalEvents(strcmp({allTalEvents.subject},subj));
%     fprintf('    Found %d electrodes for this patient.\n',length(events));
%     subjTalEvents = events;
%     clear events
%     save(talFile,'subjTalEvents')
%     fprintf('    done\n\n')
%   else
%     subjTalEvents = load(talFile);
%     subjTalEvents = subjTalEvents.talStruct;
%   end
%   
  % make the bp_struct file, if it doesn''t exist
%   if ~exist(talFile,'file')
%     [subjTalEvents hasMom] = getBipolarTalStructureFromVoxInfo(subj); 
%     if ~hasMom || isempty(subjTalEvents)
%       fprintf('Making struct from VOX_coords, jacksheet, and bp_manual.txt...')
%       subjTalEvents = getBipolarTalStructureFromOldStuff(subj);
%       if isempty(subjTalEvents)
% 	error('cannot make the bp structure from the old stuff eithers')
%       end	    
%       pause(.1)
%       fprintf('done\n')
%     end
    
    % this should never happen
    %if isempty(subjTalEvents);error('makes no sense');end
        
    % save the bipolar structure 
%    save(talFile,'subjTalEvents')
%  else
%    clear events; 
%    load(talFile)
%    if exist('events','var') 
%      fprintf('writing over ''events'' with ''subjTalEvents''.\n')
%      subjTalEvents =  events;
%      save(talFile,'subjTalEvents');
%    end    
%     if exist('bpTalStruct','var') 
%       fprintf('writing over ''bpTalStruct'' with ''subjTalEvents''.\n')
%       subjTalEvents =  bpTalStruct;
%       save(talFile,'subjTalEvents') 
%     end
    
% %remake with tagNames as a field
% if ~isfield(subjTalEvents,'tagName')
%   fprintf('\n ADDING TAGNAMES...')
%   els        = cat(1,subjTalEvents.channel);
%   [tagNames] = getECoG_subjTagNames_local(subj,els,MOUNT_DIR);
%   if length(tagNames)~=size(els,1);error('impossible');end
%   for k=1:length(subjTalEvents)
%     subjTalEvents(k).tagName = tagNames{k};
%   end  
%   save(talFile,'subjTalEvents')
%   fprintf('done\n\n')
% end
% 
% %remake with tagNames as a field
% if ~isfield(subjTalEvents,'eNames')
%   fprintf('\n ADDING ENAMES...')
%   for k=1:length(subjTalEvents)
%     thisChan = subjTalEvents(k).channel;    
%     if length(thisChan)==2
%       subjTalEvents(k).eNames = sprintf('%d-%d',thisChan);
%     elseif length(thisChan)==1
%       subjTalEvents(k).eNames = sprintf('%d',thisChan);
%     else
%       error('should never happen')
%     end
%   end
%   save(talFile,'subjTalEvents')
%   fprintf('done\n\n')
% end
%     
% % remove that redundnat isGood field (set the flag).
% if isfield(subjTalEvents,'isGood')
%   subjTalEvents = rmfield(subjTalEvents,'isGood');
% end

