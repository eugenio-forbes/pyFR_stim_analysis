function subjbipolElecMap(subj,showNum,bpOnly)
%
% plots elecs on a brain
%

global V F VS FS CS

if ~exist('showNum','var')||isempty(showNum)
  showNum=false;
end
if ~exist('bpOnly','var')||isempty(bpOnly)
  bpOnly=false;
end

if ~ismac
  resDataRoot = '/data4/real_time/subjDataBase/';
  homeDir = '/home1/jfburke/matlab/burkeSVN/';
  talDataRoot = '/data/eeg/tal';
else
  resDataRoot = '/Volumes/subjDataBase/';
  homeDir = '/Volumes/RHINO/matlab/burkeSVN/';
  talDataRoot = '/Volumes/subjDataBase/';
end

% run the config and get the tal events
talDataBaseFile  = 'allTalLocs_GM.mat';
TALEVENTS     =load(fullfile(talDataRoot,talDataBaseFile));
allTalEvents  =TALEVENTS.events;
subjTalEvents = allTalEvents(strcmp({allTalEvents.subject},subj));
subjCHANNELS  = [subjTalEvents.channel];
x = [subjTalEvents.x]';
y = [subjTalEvents.y]';
z = [subjTalEvents.z]';
clear TALEVENTS

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
clear vs fs cs
[subjElec subjMont]=getManualBipolElecs(subj);

% set some variable
defCol = [.1 .1 .1];
highlightCol = [1 0 0];
outlineCol = [.01 .01 .01];
markerSize = 15; %set this to 150 and make figures []
bipolmarkerSize = 10;
bringToSurf = 1;
Marker = 'o';
lineWidth = 2;
bpCol     = 'r';

[f1 hs(1)] = mkThisFig(1,[90 0]);hold on
[f2 hs(2)] = mkThisFig(2,'slice');hold on;xlim([-140 140])
[f3 hs(3)] = mkThisFig(3,[180 -90]);hold on;xlim([-140 140])
[f4 hs(4)] = mkThisFig(4,[-90 0]);hold on  

hasElecs = false(1,4);

% organize the elecs
for k=1:size(subjElec,1)  
  e1 = subjElec(k,1);
  e2 = subjElec(k,2);
  e1Ind = find(subjCHANNELS==e1);
  e2Ind = find(subjCHANNELS==e2);        
  if isempty(e1Ind) || isempty(e2Ind)
    fprintf('!!!WARNING!!!.  \n\n')
    keyboard
    return
  end
  thisP1 = [x(e1Ind) y(e1Ind) z(e1Ind)];
  thisP2 = [x(e2Ind) y(e2Ind) z(e2Ind)];
  thisBP = bipolarHalfDistance(thisP1,thisP2);
  
  thisMont = subjMont{k};  
  switch upper(thisMont)
   case 'HIPP'
    thisP1(3) = thisP1(3) + 100;
    thisP2(3) = thisP2(3) + 100;
    thisBP(3) = thisBP(3) + 100;
    thisF = 2;
    hasElecs(2)=true;
   case 'INF'
    thisP1(3) = thisP1(3) - 100;
    thisP2(3) = thisP2(3) - 100;
    thisBP(3) = thisBP(3) - 100;
    thisF = 3;
    hasElecs(3)=true;
   case 'LEFT'
    thisP1(1) = thisP1(1) - 100;
    thisP2(1) = thisP2(1) - 100;
    thisBP(1) = thisBP(1) - 100;
    thisF = 4;
    hasElecs(4)=true;
   case 'RIGHT' 
    thisP1(1) = thisP1(1) + 100;
    thisP2(1) = thisP2(1) + 100;
    thisBP(1) = thisBP(1) + 100;
    thisF = 1;
    hasElecs(1)=true;
   otherwise
    fprintf('bad mont 2\n\n')
    keyboard
    return
  end
  
  hanTemp = get(thisF,'Children');
  thisAx  = findobj(hanTemp,'type','axes');
  dot=plot3(thisBP(1),thisBP(2),thisBP(3),'Parent',thisAx,'Marker',Marker,...
	    'MarkerSize',bipolmarkerSize,'MarkerEdgeColor',outlineCol,...
	    'MarkerFaceColor',bpCol,'LineWidth',lineWidth);

  if ~bpOnly
    dots =plot3(thisP1(1),thisP1(2),thisP1(3),'Parent',thisAx,...
		'Marker',Marker,'visible','on',...
		'MarkerSize',markerSize,'MarkerEdgeColor',outlineCol,...
		'MarkerFaceColor',defCol,'LineWidth',lineWidth);
    if showNum
      T=text(thisP1(1),thisP1(2),thisP1(3),num2str(e1));
      set(T,'HorizontalAlignment','center','Parent',thisAx,...
	    'Color',[.99 .99 .99]);
    end
    dots =plot3(thisP2(1),thisP2(2),thisP2(3),'Parent',thisAx,...
		 'Marker',Marker,'visible','on',...
		'MarkerSize',markerSize,'MarkerEdgeColor',outlineCol,...
		'MarkerFaceColor',defCol,'LineWidth',lineWidth);
    if showNum
      T=text(thisP2(1),thisP2(2),thisP2(3),num2str(e2));
      set(T,'HorizontalAlignment','center','Parent',thisAx,...
	    'Color',[.99 .99 .99]);
    end
  end
  
end

for k=1:4
  figure(k+33)
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