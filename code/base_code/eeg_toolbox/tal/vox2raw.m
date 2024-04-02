function vox2raw(vox_coords,imageroot,ct_combined,ct2mr,mr_brain,mr2standard,reverseXcoord,fslPath);
%VOX2RAW - Convert VOX_coords.txt to RAW_coords.txt.mni
%
%  USAGE:
%   When both CT and MR images exist:
%    vox2raw(vox_coords,imageroot,ct_combined,ct2mr,mr_brain,mr2standard,reverseXcoord);
%
%   or
%
%   When only CT images exist:
%    vox2raw(vox_coords,imageroot,ct_combined,ct_registered,reverseXcoord);
%
%  VARIABLES:
%   vox_coords = path to voxel coordinates file (e.g., 'VOX_coords.txt')
%   imageroot = path to image files (e.g., 'images/combined')
%   ct_combined = CT_combined file (e.g., 'UP004_CT_combined')
%   ct_registered = CT_registered file (e.g., 'UP004_CT_registered')
%   ct2mr = CT2MR file (e.g., 'UP004_CT2MR')
%   mr_brain = MR_brain file (e.g., 'UP004_MR_brain')
%   mr2standard = MR2standard file (e.g., 'UP004_MR2standard')
%   reverseXcoord = Reverse the sign of the X-coordinate to flip
%                   the left and right sides of the coordinates
%                   (1 = yes, 0 = no).
%                   NB: Only reverses the sign in RAW_coords.txt.mni
%
%  EXAMPLE:
%   When both CT and MR images exist and the X-coordinate does not
%   need to be reversed:
%    vox2raw('VOX_coords.txt','images/combined','UP004_CT_combined','UP004_CT2MR','UP004_MR_brain','UP004_MR2standard',0);
%
%   or
%
%   When only CT images exist and the X-coordinate needs to be reversed:
%    vox2raw('VOX_coords.txt','images/combined','UP004_CT_combined','UP004_CT_registered',1);
%
%   NB: This function assumes that you are in the directory that
%       contains VOX_coords.txt, which should be SUBNUM/tal
%
%       Make sure FSL is installed correctly, as this script requres a
%       binary file in /Applications/FSL/fsl/bin
%

% Create VOX_coords_nochan.txt
VOX_coords = load(vox_coords);
[voxPathstr,voxName,voxExt] = fileparts(vox_coords);
% remove the channel column
VOX_coords_nochan = VOX_coords(:,2:4);
% save a VOX_coords_nochan.txt without the channel number column
fid = fopen([voxName,'_nochan.txt'],'w');
fprintf(fid,'%g %g %g\n',VOX_coords_nochan');
fclose(fid);

if ~exist('fslPath','var')||isempty(fslPath)
  fslPath = '/usr/global/fsl-2.0.10/fsl/bin/';
end

% Run img2stdcoord or img2talcoord to convert to MNI space
% If CT and MR images exist
if nargin >= 7
  
  % touch it to see if you have permissions
  if ~exist('MM_coords_MR.txt','file')
    touchedIt = system(['touch  MM_coords_MR.txt']);
    if touchedIt~=0
      fprintf('\n\ndo not have permisions to write at %s\n\n',pwd)
      return
    end
  end
    
  % convert voxel space to mm
  if exist(fullfile(fslPath,'img2stdcoord'),'file')
    unix(['img2stdcoord -img ',fullfile(imageroot,ct_combined),' -std ',fullfile(imageroot,ct2mr),' -xfm ',fullfile(imageroot,ct2mr),'.mat ',voxName,'_nochan.txt > MM_coords_MR.txt']);
  else
    fprintf('\n\n\n\n')
    error(sprintf('vox2raw:fileNotFound: Cannot find %s. Check FSL installation.',fullfile(fslPath,'img2stdcoord')));
  end
      
  % pause for a second because annoying.
  pause(3);
  
  % trim the last line because img2stdcoord/img2talcoord prints the last
  % line twice
  MM_coords_MR = load('MM_coords_MR.txt');
  MM_coords_MR = MM_coords_MR(1:end-1,:);
  % save
  fid = fopen('MM_coords_MR.txt','w');
  fprintf(fid,'%g %g %g\n',MM_coords_MR');
  fclose(fid);
  
  if isempty(MM_coords_MR)
     fprintf('ERROR With img2stdcoords')
     delete('MM_coords_MR.txt')
     return
  end
  
  % touch it to see if you have permissions
  if ~exist('RAW_coords.txt.mni','file')
    touchedIt = system(['touch  RAW_coords.txt.mni']);
    if touchedIt~=0
      fprintf('\n\ndo not have permisions to write RAW_coords.txt.mni at %s\n\n',pwd)
      return
    end
  end
  
  % convert mm to MNI space
  if exist(fullfile(fslPath,'img2stdcoord'),'file')
    unix(['img2stdcoord -mm -img ',fullfile(imageroot,mr_brain),' -std ',fullfile(imageroot,mr2standard),' -xfm ',fullfile(imageroot,mr2standard),'.mat MM_coords_MR.txt > RAW_coords.txt.mni']);
  else
    fprintf('\n\n\n\n')
    error(sprintf('vox2raw:fileNotFound: Cannot find %s. Check FSL installation.',fullfile(fslPath,'img2stdcoord')));
  end
    
  % trim the last line because img2stdcoord/img2talcoord prints the last
  % line twice
  RAW_coords = load('RAW_coords.txt.mni');
  RAW_coords = RAW_coords(1:end-1,:);
  if reverseXcoord == 1
    RAW_coords(:,1) = -1 * RAW_coords(:,1);
  end
  % save
  fid = fopen('RAW_coords.txt.mni','w');
  fprintf(fid,'%g %g %g\n',RAW_coords');
  fclose(fid);
  
elseif nargin == 5
  % If only CT images exist
  % Set variables correctly (inititally set incorrectly due to fewer varargin)
  ct_registered = ct2mr;
  reverseXcoord = mr_brain;
  % convert voxel space to MNI space
  if exist([fsl_dir,'/img2stdcoord'],'file')
    unix(['img2stdcoord -mm -img ',fullfile(imageroot,ct_combined),' -std ',fullfile(imageroot,ct_registered),' -xfm ',fullfile(imageroot,ct_registered),'.mat ',voxName,'_nochan.txt > RAW_coords.txt.mni']);
  elseif exist([fsl_dir,'/img2talcoord'],'file')
    unix(['img2talcoord -img ',fullfile(imageroot,ct_combined),' -tal ',fullfile(imageroot,ct_registered),' -xfm ',fullfile(imageroot,ct_registered),'.mat ',voxName,'_nochan.txt > RAW_coords.txt.mni']);
  else
    error('vox2raw:fileNotFound','Cannot find /Applications/FSL/fsl/bin/img2stdcoord or /Applications/FSL/fsl/bin/img2talcoord. Check FSL installation.');
  end
  % trim the last line because img2stdcoord/img2talcoord prints the last
  % line twice
  RAW_coords = load('RAW_coords.txt.mni');
  RAW_coords = RAW_coords(1:end-1,:);
  if reverseXcoord == 1
    RAW_coords(:,1) = -1 * RAW_coords(:,1);
  end
  fid = fopen('RAW_coords.txt.mni','w');
  fprintf(fid,'%g %g %g\n',RAW_coords');
  fclose(fid);
else
  error('vox2raw:wrongNumOfVars','Did not enter a supported number of variables.');
end

% Create RAW_coords.txt.mni with a channel number column
VOX_coords = load(vox_coords);
% grab the channel number column from VOX_coords
chans = VOX_coords(:,1);
RAW_coords = load('RAW_coords.txt.mni');
% add the channel number column to RAW_coords
RAW_coords = [chans,RAW_coords];
% write out a new RAW_coords.txt.mni file with channel numbers
fid = fopen('RAW_coords.txt.mni','w');
fprintf(fid,'%d %g %g %g\n',RAW_coords');
fclose(fid);
