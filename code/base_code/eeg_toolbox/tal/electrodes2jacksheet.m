function electrodes2jacksheet(subj)
%this function creates a jacksheet.txt from an electrodes.m file for
%subjects who don't have one. It prompts the user to enter electrode
%lbels for each group of electrodes
%%%inputs
%subj = 'UP035';

%written by AGR 06/17/2013 using code from JFB



% am I mounting RHINO?
if ismac
  MOUNT_DIR = '/Volumes/RHINO_root';
else
  MOUNT_DIR = '';
end
if ~exist('badTagNames','var')||isempty(badTagNames)
  badTagNames ={};
end

% directories
subjDir = fullfile(MOUNT_DIR,'/data/eeg',subj);
docDir  = fullfile(subjDir,'docs');
talDir  = fullfile(subjDir,'tal');

% file names
jackFileName      = 'jacksheet.txt';
jackFileName_2    = 'jack_sheet.txt';
elecFileName      = 'electrodes.m';

%check if jacksheet exists (if exists, move it to jacksheet_orig and display warning)
cd(docDir);
if exist(jackFileName,'file') || exist(jackFileName_2,'file')
   %!mv jacksheet.txt jacksheet_orig.txt 
   %disp ('WARNING: Original jacksheet file is copied to jacksheet_orig.txt')
   error('Jacksheet already exists') 
end
% get the elec matrix and elec labels from the electrodes.m file
[gridsAndStrips,lbls] = getTheContentsOFTheElecs_local(docDir,elecFileName);

%write the jacksheet
fid = fopen(jackFileName,'w');
allNeuralElecs = [];
for k=1:size(gridsAndStrips,1)
  numElecs =  1+gridsAndStrips(k,2)-gridsAndStrips(k,1);
  for e = 1:numElecs
     fprintf(fid,'%d %s%d\n',gridsAndStrips(k,1)-1+e,lbls{k},e)
  end
end

%keyboard

%Functions below:
%--------------------------------------------------
function [out,lbls] = getTheContentsOFTheElecs_local(d,f);
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
  
 % get lbls
 lbls = cell(size(out,1),1);
 cd(d)
 fid = fopen(f);
 X = textscan(fid,'%s%s','delimiter','%');
 count = 0;
 for i = 1:size(X{2},1)
    if isempty(X{2}{i});continue; end
    count = count+1;
    Y = textscan(X{2}{i},'%s%s','delimiter','_');
    lbls{count} = Y{1}{1};
 end
 