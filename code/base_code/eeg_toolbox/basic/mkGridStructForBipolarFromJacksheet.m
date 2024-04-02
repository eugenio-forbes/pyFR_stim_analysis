function mkGridStructForBipolarFromJacksheet(subj,badTagNames)
%
% FUNCTION:
%   mkGridStructForBipolarFromJacksheet
% 
% DECRIPTION:
%   For one subject, gets electrodes for bipolar and regular
%   montages.  It uses the electrodes.m file to get all the
%   electrodes that need to be included in the jacksheet.
%
% INPUTS:
%   subj............. 'TJ052'
%
% OUTPUTS:
%
% NOTES:
%   (1) written by jfburke on 10/2012 (john.fred.burke@gmail.com)
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
fprintf('\n\n');

% get jacksheet
subjDir = fullfile(MOUNT_DIR,'/data/eeg',subj);
docDir  = fullfile(subjDir,'docs');

% see if you have an electrodes.m file
elecFile = fullfile(docDir,'electrodes.m');
if ~exist(elecFile,'file')  
  fprintf('electrodes.m does not exist...EXITING\n\n')
  return
end

% get the contents of the jackfile
[elecNum elecNam elecNam_stripped] = ...
    getTheContentsOFTheJackFile_local(docDir,'jacksheet.txt','jack_sheet.txt');

% get the electrodes.m contents
run(elecFile);
try
  gridsAndStrips = r; clear r;
catch
  gridsAndStrips = grids; clear grids;
end
  
% define bad elec names
badTagNames_FINAL = {'PENN','REF','EKG','DXL','DRXL','GND','GRD','EOG','ECG','EMG'};
badTagNames_FINAL = cat(2,badTagNames_FINAL,badTagNames);

% get the grid file
gridFile = fullfile(docDir,'grids_struct_for_bipolar.m');
if exist(gridFile)
  fprintf('\n\ngrids_struct_for_bipolar.m exists... EXITING\n\n')
  return
end

% write the grid file
fid = fopen(gridFile,'w'); % <--- not discards contents of current file
if fid==-1
  fprintf('dont have permissions to write...EXITING\n\n')
  return
end

% write the grid file preample
write_preamble_lcoal(fid);

% write the meat of the grid file
for k=1:size(gridsAndStrips,1)
  els                     = [gridsAndStrips(k,1) gridsAndStrips(k,2)];
  theseElecs              = els(1):els(2);
  contacts_ThisTag_totNum = length(theseElecs);
  
  % get the tag names associated with thes electrodes
  if ~isempty(elecNam)
    theseElecsInd = ismember(elecNum,theseElecs);
    theseTags     = elecNam_stripped(theseElecsInd);
    unTags        = unique(theseTags);    
    if length(unTags)>1
      thisTag = 'MULTIPLE_TAGS_FOUND';
    else
      thisTag = unTags{1};        
    end
  else
    thisTag = 'NO_JACKSHEET_FOUND';
  end

  if sum(ismember(upper(badTagNames_FINAL),upper(thisTag)))>0
    fprintf(' %20.20s (%-3.0d to %3.0d)  **NOT** Including in bipolar montage',thisTag,els)
    continue
  else
    fprintf(' %20.20s (%-3.0d to %3.0d)   Including in bipolar montage',thisTag,els)
  end  
  
  % get the hemisphere
  switch upper(thisTag(1))
   case 'R'
    hemStr = 'right';
   case 'L'
    hemStr = 'left';
   otherwise
    hemStr = 'unknown';
  end
    
  if contacts_ThisTag_totNum<=8
    fprintf(' %3d CONTACTS...strip?',contacts_ThisTag_totNum)
    
    % print the line
    fprintf(fid,'    [%-2.3d %-2.2d %-2.2d] ''%s''; ... %% %s: %d-%d: %s\n',...
	    els(1),contacts_ThisTag_totNum,1,...
	    hemStr,thisTag,els(1),els(2),...
	    ['STRIP GUESS MADE AUTOMATICALLY ON ' date]);
  else
    fprintf(' %3d CONTACTS...**GRID**?',contacts_ThisTag_totNum);
    
    % print the line
    fprintf(fid,'    [%-2.3d ???? ????] ''%s''; ... %% %s: %d-%d: %s\n',...
	    els(1),hemStr,thisTag,els(1),els(2),...
	    ['GRID GUESS MADE AUTOMATICALLY ON ' date]);
  end
    
  fprintf('\n')
end
fprintf(fid,'	       };\n');
fclose(fid);
fprintf('\n\n')

function write_preamble_lcoal(fid);
  fprintf(fid,'%%\n');
  fprintf(fid,'%%THIS is used by make_AUTO_bipolar_pairs.m\n');
  fprintf(fid,'%%\n');
  fprintf(fid,'%%grids_struct\n');
  fprintf(fid,'%%  element 1: 1x3 vectore corresponding to:\n');
  fprintf(fid,'%%     1,1 = electrode startig number\n');
  fprintf(fid,'%%     1,2 = # elecs before first break in grid\n');
  fprintf(fid,'%%     1,3 = # repeating elements (above) in grid\n');
  fprintf(fid,'%%  element 2: view\n');
  fprintf(fid,'%%     = ''left'' or ''right'' or ''hipp'' or ''inf''\n');
  fprintf(fid,'%%\n');
  fprintf(fid,'%%examples:\n');
  fprintf(fid,'%%      8x1 strip starting at 64\n');
  fprintf(fid,'%%        = [64 8 1], ''left''\n');
  fprintf(fid,'%%      8x8 grid starting at 1\n');
  fprintf(fid,'%%        = [1 8 8], ''left''\n');
  fprintf(fid,'%%      8x1 strip starting at 36 half on\n');
  fprintf(fid,'%%      left view, half on inferior view\n');
  fprintf(fid,'%%        = [36 4 1], ''left'';...\n');
  fprintf(fid,'%%          [40 4 1], ''inf'';...\n');
  fprintf(fid,'%%\n');
  fprintf(fid,'\n');
  fprintf(fid,'grids_struct = { ...\n');

  
function [num nam nam_stripped] = getTheContentsOFTheJackFile_local(dDir,j1,j2);
  jFil_1 = fullfile(dDir,j1);
  jFil_2 = fullfile(dDir,j2);  
  fid_1  = fopen(jFil_1,'r');
  fid_2  = fopen(jFil_2,'r');
  
  if fid_1~=-1 & fid_2==-1
    fprintf('Found jack sheet in %s\n\n',jFil_1)
    fid = fid_1;
  elseif fid_1==-1 & fid_2~=-1
    fprintf('Found jack sheet in %s\n\n',jFil_2)
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