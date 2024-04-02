function vox2voxMother(subj)
%This function creates a vox mother file for subjects who have a voxCoords
%and a grids_struct_for_bipolar. Assumes that the electrode numbering in
%VOX_coords is same as jacksheet (which it should be).
%Inputs
%subj = UP035
%written by AGR 06/17/2013 using code from JFB

% am I mounting RHINO?
if ismac
  MOUNT_DIR = '/Volumes/RHINO_root';
else
  MOUNT_DIR = '';
end

% directories
subjDir = fullfile(MOUNT_DIR,'/data/eeg',subj);
docDir  = fullfile(subjDir,'docs');
talDir  = fullfile(subjDir,'tal');

% file names
jackFileName      = 'jacksheet.txt';
jackFileName_2    = 'jack_sheet.txt';
voxFileName      = 'VOX_coords.txt';
gridStructFileName = 'grids_struct_for_bipolar.m';
voxMomFileName      = 'VOX_coords_mother.txt';

%check if voxMom already written
cd(talDir)
if exist(voxMomFileName,'file'); error('vox mother already exists'); end

% get the vox file
voxFileName = fullfile(talDir,voxFileName);
if ~exist(voxFileName,'file')
  error(' no VOX_coords.txt file') 
end
fid = fopen(voxFileName);
X = textscan(fid,'%d%d%d%d');
fclose(fid);
voxElecs = X{1};
voxXYZ   = [X{2} X{3} X{4}]; 

% get jacksheet
[elecNum elecNam elecNam_stripped] = ...
    getTheContentsOFTheJackFile_local(docDir,jackFileName,jackFileName_2);
if isempty(elecNum);error(' no jacksheet');end
elecList = unique(elecNam_stripped);


% read grid_structs_for_bipolar
[dim1,dim2,stripOrDepth,startVal] = readGridStructs(gridStructFileName,docDir);

%Write vox mother file
cd(talDir)
fid = fopen(voxMomFileName,'w');
for e = 1: size(voxElecs,1)
    ind = find(elecNum==voxElecs(e));
    groupInd = max(find(startVal<=elecNum(ind)));
    fprintf(fid,'%s\t%d %d %d %s\t %d %d',elecNam{ind},voxXYZ(e,1),voxXYZ(e,2),voxXYZ(e,3),stripOrDepth{groupInd},dim1(groupInd),dim2(groupInd))
    fprintf(fid,'\n')
end
fclose(fid)


%_________________Functions Below__________________
function [dim1,dim2,stripOrDepth,startVal] = readGridStructs(gridStructFileName,docDir);
cd(docDir);
grids_struct_for_bipolar

%get dimensions
dims = cat(1,grids_struct{:,1});
dim1 = dims(:,2);
dim2 = dims(:,3);
startVal = dims(:,1);

%get stripOrDepth
depthInd = strcmp({grids_struct{:,2}},'hipp'); %although, this will miss frontal deptht
stripOrDepth = cell(size(dims,1),1);
stripOrDepth(depthInd) = {'D'};
stripOrDepth(~depthInd) = {'S'};
gridInd = dim2>1;
stripOrDepth(gridInd) = {'G'};
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