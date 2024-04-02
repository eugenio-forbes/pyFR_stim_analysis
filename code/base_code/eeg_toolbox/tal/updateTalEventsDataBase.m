function updateTalEventsDataBase(subj, subjDir)
%
% FUNCTION:
%   
% DESCRIPTION:
%
% INPUT:
%
% OUTPUT:
%
% NOTES:
%   (1) written by jfburke [''] (john.fred.burke@gmail.com)
%
%

fprintf('\n\n')

% ste the file names
talLocsFileName = 'RAW_coords.txt';
mniLocsFileName = 'RAW_coords.txt.mni';
talEvdbFileName = sprintf('%s_talLocs_database.mat',subj);

% am I mounting RHINO?
% if ismac
%   MOUNT_DIR = '/Volumes/RHINO_root';
% else
%   MOUNT_DIR = '';
% end

% subjDir    = fullfile(MOUNT_DIR,'/data/eeg/',subj);
subjTalDir = fullfile(subjDir,'tal');
% subjDocDir = fullfile(subjDir,'docs');

% make the talairach coords if needed
talLocsFile = fullfile(subjTalDir,talLocsFileName);
% mniLocsFile = fullfile(subjTalDir,mniLocsFileName); 
talEvdbFile = fullfile(subjTalDir,talEvdbFileName); 
if exist(talEvdbFile,'file')
  fprintf('\tElectrodes Tal structure already exists.... EXITING\n\n')  
  ev = load(talEvdbFile);
  plotBipolElecMap_local(ev.events);
  return
end

% make the talaiarch coords file if it doesn't exist
if ~exist(talLocsFile,'file')
  fprintf('\nMaking talairach file...\n')
  origDir = pwd;
  cd(subjTalDir)
  convertFile_mni2tal(mniLocsFileName);
  cd(origDir);
  pause(.1)
  fprintf('done\n\n')
end

% read in the talaiarch coords
tal_contents = get_this_file_local(subjTalDir,talLocsFileName,'%d%f%f%f');
elec_inds    = tal_contents{1};
elec_xyz     = [tal_contents{2} tal_contents{3} tal_contents{4}];

% make a bunch of files if they are not made
fid_inf = open_this_file_local(subjTalDir,'inf.montage');
fid_rig = open_this_file_local(subjTalDir,'right.montage');
fid_lef = open_this_file_local(subjTalDir,'left.montage');
fid_bad = open_this_file_local(subjTalDir,'bad_leads.txt');

% make the events structure
events = struct('subject',[],'channel',[],'x',[],'y',[],'z',[],...
         'Loc1',[],'Loc2',[],'Loc3',[],'Loc4',[],'Loc5',[],...
         'Loc6',[],'isGood',[],'montage',[]);

fprintf('\n')
fprintf('Getting the Talaiarch information for each electrode: ')
ticker   = 0; 
tick_inc = 10;
for e = 1:size(elec_inds,1)
  if e./size(elec_inds,1)*100>=ticker
    fprintf(' %d%%',ticker)
    ticker = ticker + tick_inc;
  end
  events(e).subject = subj;
  events(e).channel = elec_inds(e);
  events(e).x       = elec_xyz(e,1);
  events(e).y       = elec_xyz(e,2);
  events(e).z       = elec_xyz(e,3);  
  events(e).isGood  = true;   
  events(e)=tal2Region(events(e),false);
  
  % take a guess at the montage
  switch upper(events(e).Loc2)
   case {'*','LIMBIC LOBE','SUB-LOBAR'}
    if ~isempty(fid_inf);fprintf(fid_inf,'%d\n',elec_inds(e));end
    events(e).montage = 'inf';
   otherwise
    if events(e).x<0
      if ~isempty(fid_inf);fprintf(fid_lef,'%d\n',elec_inds(e));end	
      events(e).montage = 'lsag';
    else
      if ~isempty(fid_rig);fprintf(fid_rig,'%d\n',elec_inds(e));end   
      events(e).montage = 'rsag';
    end    
  end
    
end
fprintf('\n\n')

% close a bunch of files (if they are open)
if ~isempty(fid_lef);fclose(fid_lef);end
if ~isempty(fid_rig);fclose(fid_rig);end
if ~isempty(fid_inf);fclose(fid_inf);end
if ~isempty(fid_bad);fclose(fid_bad);end
  
% save the file in the tal directory
save(talEvdbFile,'events');

% plot it 
plotBipolElecMap_local(events)

% -----------------------------------------------------------------
function fid = open_this_file_local(d,f)
  fil = fullfile(d,f);
  if exist(fil,'file')
    fprintf('%s already exists in %s\n',f,d)
    fid=[];
    return
  end
  
  fprintf('Making %s in %s\n',f,d)
  fid = fopen(fil,'w');  
  if fid==-1
   fprintf('\n\nERROR: cannot open %s in %s\n\n',f,d)
   error('done')
  end  

% -----------------------------------------------------------------
function out = get_this_file_local(d,f,thing)
  fil = fullfile(d,f);
  fid = fopen(fil);  
  if fid==-1
   fprintf('\n\nERROR: %s does not exist in %s\n\n',f,d)
   error('done')
  end  
  out = textscan(fid,thing);
  fclose(fid);

% -----------------------------------------------------------------  
function plotBipolElecMap_local(ev)
%
% plots elecs on a brain
%

global V F VS FS CS

subjCHANNELS  = [ev.channel];
x             = [ev.x]';
y             = [ev.y]';
z             = [ev.z]';
subjMont      = {ev.montage};

% load the brain pic
picpath = fileparts(which('tal3d'));
picfile1 = fullfile(picpath,'mni_cortical_surface.mat');
picfile2 = fullfile(picpath,'mni_depth_slice.mat');
load(picfile1);  
V = v;
F = f;
clear v f vdist vloc 
load(picfile2);  
VS = vs{6};
FS = fs{6};
CS = cs{6};
clear vs fs c

% set some variable
defCol = [.1 .1 .1];
highlightCol = [1 0 0];
outlineCol = [.01 .01 .01];
markerSize = 15; %set this to 150 and make figures []
bringToSurf = 1;
Marker = 'o';
lineWidth = 2;

[f1 hs(1)] = mkThisFig(1,[90 0]);hold on
[f2 hs(2)] = mkThisFig(2,'slice');hold on;xlim([-140 140])
[f3 hs(3)] = mkThisFig(3,[180 -90]);hold on;xlim([-140 140])
[f4 hs(4)] = mkThisFig(4,[-90 0]);hold on  

hasElecs = false(1,4);

% organize the elecs
for k=1:length(subjCHANNELS)  
  thisE = subjCHANNELS(k);
  thisP = [x(k) y(k) z(k)];  
  thisMont = subjMont{k};  
  
  switch upper(thisMont)
   case 'HIPP'
    thisP(3) = thisP(3) + 100;
    thisF = 2;
    hasElecs(2)=true;
   case 'INF'
    thisP(3) = thisP(3) - 100;
    thisF = 3;
    hasElecs(3)=true;
   case 'LSAG'
    thisP(1) = thisP(1) - 100;
    thisF = 4;
    hasElecs(4)=true;
   case 'RSAG' 
    thisP(1) = thisP(1) + 100;
    thisF = 1;
    hasElecs(1)=true;
   otherwise
    fprintf('bad mont 2\n\n')
    return
  end
  
  hanTemp = get(thisF,'Children');
  thisAx  = findobj(hanTemp,'type','axes');
  dots =plot3(thisP(1),thisP(2),thisP(3),'Parent',thisAx,...
	      'Marker',Marker,'visible','on',...
	      'MarkerSize',markerSize,'MarkerEdgeColor',outlineCol,...
	      'MarkerFaceColor',defCol,'LineWidth',lineWidth);
  
  T=text(thisP(1),thisP(2),thisP(3),num2str(thisE));
  set(T,'HorizontalAlignment','center','Parent',thisAx,...
	'Color',[.99 .99 .99]);
  
end

for k=1:4
  figure(k)
  if k~=2
    hLight = camlight('headlight');
    set(hLight,'Color',[1 1 1],'Style','infinite');
    lighting phong
    setBrainProps(hs(k),hasElecs(k));
  end
end

function [f hs] = mkThisFig(k,viewAZEL)
global V F VS FS CS
f = figure(k);
clf
if ~ischar(viewAZEL)
  view(viewAZEL)
  FVCD = repmat([.55 .55 .55],size(V,1),1);
  thisF = F;
  thisV = V;
else
  CS      = round(CS);
  BW_cMap = gray(125);
  FVCD = BW_cMap(CS,:);
  colormap(gray(128));
  view([0 90])
  thisF = FS;
  thisV = VS;
end

hs = patch('faces',thisF,'vertices',thisV,'edgecolor','none','FaceColor',...
	   'interp','FaceVertexCData',FVCD);
setBrainProps(hs);  