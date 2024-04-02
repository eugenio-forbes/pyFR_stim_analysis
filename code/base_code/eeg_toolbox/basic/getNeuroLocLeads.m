function [elsOUT,errFlag,jacOUT,disOUT] = getNeuroLocLeads(subj,loc1,loc2,distVec)
% This function interfaces with neurorad loc.txt
% Inputs:
% subj .... 'TJ039'
% loc1 .... mtl region of interest (see below)
% loc2 .... which 3 x 3 quadrant of the mtl (see an example sheet)
% distOUT .... distance from the hippocampal head 
%
%
% Notes:
% "the convention is right"
%
%

if ismac
  mounDir = '/Volumes/RHINO_root';
else
  mounDir = '';
end

% defaults
errFlag = 0;

if ~exist('loc1','var')||isempty(loc1)
  fprintf('\n\nplease enter a loc1 field\n\n')
  return
end
if ~exist('loc2','var')||isempty(loc2)
  skipCompass = true;
else
  skipCompass = false;
end
if ~exist('distVec','var')||isempty(distVec)
  distVec = [0 Inf];
end


% directories
dataDir = '/data/eeg';
dataDir = fullfile(mounDir,dataDir);
subjDir = fullfile(dataDir,subj);
docsDir = fullfile(subjDir,'docs');
talaDir = fullfile(subjDir,'tal');

% files
locName = 'neurorad_localization.txt';
locFile = fullfile(talaDir,locName);
if ~exist(locFile,'file')
  %fprintf('\n\n\t file does not exist\n\n')
  elsOUT=[]; jacOUT=[]; disOUT=[];
  errFlag = 1;
  return
end
jacName = 'jacksheet.txt';
jacFile = fullfile(docsDir,jacName);
if ~exist(jacFile,'file')
  elsOUT=[]; jacOUT=[]; disOUT=[];
  fprintf('\n\n\t jacksheet does not exist\n\n')
  errFlag = 1;
  return
end

% get the contents of the localization
fid = fopen(locFile,'r');
Loc_hem   = {};
Loc_dis   = []; 
Loc_nam   = {};
Loc_num   = []; 
Loc_roi   = {};
Loc_com   = []; 
Loc_nbs   = {};
nextLineNewElectrode = false;
gotTheNewElectrode   = false;

while true

  X = fgetl(fid);
  if ~ischar(X);break;end
  if ischar(X) && strcmp(X,'**')
    nextLineNewElectrode = true;
    continue
  end
  if isempty(X)
    gotTheNewElectrode = false;
    continue
  end
 
  
  % this parses the first line of the code
  if nextLineNewElectrode
    NEWELEC = textscan(X,'%s%s%s','delimiter','\t');
    nextLineNewElectrode = false;
    gotTheNewElectrode   = true;
    
    
    %special case when reach the sphenoid to hippocampus distance (mbm addition 1.27.14)
    if isempty(NEWELEC{3})
     
      while true
	doubleCheck=regexp(NEWELEC{1},'_','split');
	if false(doubleCheck{1}{2}=='sphen2hipp')
	  fprintf('\n\n\something wrong with sphen2hipp data uptake\n\n')
	  errFlag=1
	  return
	end
	distWithUnits = NEWELEC{2}{1};
	dist          = (regexp(distWithUnits,' ','split'));
	dist          = str2num(dist{1});%
	if doubleCheck{1}{1}=='R'
	  rSphenDist = dist;
	else
	  lSphenDist = dist;
	end

	%ugly code to go to next line and check if second (ie other side) measurement
	X = fgetl(fid);
	if ~ischar(X) || isempty(X);break;end %move on, equivalent to break above if there are no distance measurements at end of .txt file
	%grab next line if exists
	NEWELEC = textscan(X,'%s%s%s','delimiter','\t');
      end
   
    end
   
    this_hem = NEWELEC{1};
    this_dis = textscan(NEWELEC{2}{1},'%f%s');
    this_dis = this_dis{1};
    this_nam = NEWELEC{3};  
    continue
  end
  
  %this parses the meat of the localization
  if gotTheNewElectrode
    NEWLINE  = textscan(X,'%d%s%s%s','delimiter','\t');
    Loc_hem  = cat(1,Loc_hem,this_hem);
    Loc_dis  = cat(1,Loc_dis,this_dis);
    Loc_nam  = cat(1,Loc_nam,this_nam);
    Loc_num  = cat(1,Loc_num,NEWLINE{1});
    Loc_roi  = cat(1,Loc_roi,NEWLINE{2}{1});
    if strcmp(NEWLINE{3}{1},'--')
      com_tmp = [-999 -999];
    else
      this_com = textscan(NEWLINE{3}{1},'%d%d','delimiter',',');    
      com_tmp  = [this_com{1} this_com{2}];
    end    
    Loc_com  = cat(1,Loc_com,com_tmp);
    if ~isempty(NEWLINE{4})
      Loc_nbs = cat(1,Loc_nbs,NEWLINE{4}{1});
    else
      Loc_nbs = cat(1,Loc_nbs,{'---'});
    end    
  end  
end
fclose(fid);

%subtract sphenoid-hipp distance from sphen-elec distance to calculate distance of electrode from head of hippocampus; mbm addition 1.27.14

if exist('rSphenDist','var') || exist('lSphenDist','var')
    for d=1:length(Loc_dis)
    if Loc_hem{d}=='R'
      toSubtract=rSphenDist;
    else
      toSubtract=lSphenDist; 
    end
    Loc_dis(d)=Loc_dis(d)-toSubtract;
  end
else
  Loc_dis(:)=[];
end



%  [Loc_hem num2cell(Loc_dis) Loc_nam num2cell(Loc_num) Loc_roi ...
%   num2cell(Loc_com) Loc_nbs]
  
% get contents the jacksheet 
fid = fopen(jacFile,'r');
JAC=textscan(fid,'%d%s','delimiter','\t');
fclose(fid);
jac_nam = regexprep(JAC{2},'\d','');
jac_nam = regexprep(jac_nam,'-REF','');
jac_nam = regexprep(jac_nam,'\W','');
jac_num = str2double(regexprep(JAC{2},'\D',''));
els_num = JAC{1};


% get the region of interest indices

switch upper(loc1)
 case {'HIPP',...
       'H','HIP','HIPPO','HIPPOCAMPUS','GET_MONEY','GET MONEY'}
  locInd = strcmp(upper(Loc_roi),'HIPP');
 case {'HS',...
       'HIPP SULCUS', 'HIPP_SULCUS','HIPPSULCUS',...
       'HIPPOCAMPAL SULCUS','HIPPOCAMPAL_SULCUS','HIPPOCAMPALSULCUS',...
       'HIPP SULC','HIPP_SULC','HIPPSULC'}
  locInd = strcmp(upper(Loc_roi),'HS');
 case {'PHC',...
      'P','PARAHIPP CORTEX','PARAHIPP_CORTEX','PARAHIPPCORTEX',...
      'PARAHIPPOCAMPAL CORTEX','PARAHIPPOCAMPAL_CORTEX','PARAHIPPOCAMPALCORTEX'}
  locInd = strcmp(upper(Loc_roi),'PHC');
 case {'CS',...
       'COLL SULCUS', 'COLL_SULCUS','COLLSULCUS',...
       'COLLATERAL SULCUS','COLLATERAL_SULCUS','COLLATERALSULCUS',...
       'COLL SULC','COLL_SULC','COLLSULC'}  
  locInd = strcmp(upper(Loc_roi),'CS');
 case {'FG',...
      'F','FUSIFORM','FUSIFORM GYRUS','FUSIFORM_GYRUS','FUSIFORMGYRUS'}
  locInd = strcmp(upper(Loc_roi),'FG');
 case {'ITG',...
      'INF TEMP GYRUS','INF_TEMP_GYRUS','INFTEMPGYRUS',...
      'INF TEMPORAL GYRUS','INF_TEMPORAL_GYRUS','INFTEMPORALGYRUS',...
      'INFERIOR TEMPORAL GYRUS','INFERIOR_TEMPORAL_GYRUS','INFERIORTEMPORALGYRUS'}
  locInd = strcmp(upper(Loc_roi),'ITG');
 case {'MTG',...
      'MID TEMP GYRUS','MID_TEMP_GYRUS','MIDTEMPGYRUS',...
      'MID TEMPORAL GYRUS','MID_TEMPORAL_GYRUS','MIDTEMPORALGYRUS',...
      'MIDDLE TEMPORAL GYRUS','MIDDLE_TEMPORAL_GYRUS','MIDDLETEMPORALGYRUS'}
  locInd = strcmp(upper(Loc_roi),'MTG');
 case {'STG',...
      'SUP TEMP GYRUS','SUP_TEMP_GYRUS','SUPTEMPGYRUS',...
      'SUP TEMPORAL GYRUS','SUP_TEMPORAL_GYRUS','SUPTEMPORALGYRUS',...
      'SUPERIOR TEMPORAL GYRUS','SUPERIOR_TEMPORAL_GYRUS','SUPERIORTEMPORALGYRUS'}
  locInd = strcmp(upper(Loc_roi),'STG');
 case {'MTL'}
  locInd = strcmp(upper(Loc_roi),'MTL'); %unspecified MTL (see TJ070)
 case {'AMY'}
  locInd = strcmp(upper(Loc_roi),'AMY'); 
  Loc_dis = []; % this is to ensure that amyg elecs dont get thrown out because they are ant. to hipp head. agr 07-30-14
 case {'EC','ENTORHINALCORTEX', 'ENT CTX', 'ENTORHINAL CORTEX', 'ENTORHINAL_CTX', 'ENTORHINAL_CORTEX'}
  locInd = strcmp(upper(Loc_roi),'PHC'); 
  locInd2 = strcmp(upper(Loc_roi),'AMY');
  allLocs={locInd,locInd2};
  distVec = [-inf 1.5];
  %this distance vector and compass localizations override any inputed into function
  loc2_1    =   [1 2;1 3;2 2;2 3;3 2;3 3]; %phc
  loc2_2    =   [2 3;3 3;3 2];%amyg
  allCompassGroups={loc2_1,loc2_2};
  skipCompass = false;
 case {'PRC', 'PERIRHINAL CTX', 'PERIRHINAL CORTEX', 'PERIRHINALCORTEX', 'PERIRHINALCTX', 'PRCTX', 'PR_CTX'}
  locInd = strcmp(upper(Loc_roi),'PHC'); 
  locInd2 = strcmp(upper(Loc_roi),'FG');
  allLocs={locInd,locInd2};
  %this distance vector and compass localizations override any inputed into function
  distVec = [-inf 1.5];
  loc2_1    =   [1 1; 1 2; 1 3]; %phc
  loc2_2    =   [1 2; 1 3; 2 2; 2 3; 3 2; 3 3];%fusiform
  allCompassGroups={loc2_1,loc2_2};
  skipCompass = false;
 otherwise
  fprintf('\n\nbad loc1 field\n')
  return
end

% get the distance indices
if ~isempty(Loc_dis)
  disInd = Loc_dis >= distVec(1) & Loc_dis <= distVec(2);
else
  disInd = true(size(locInd));
end

% get the compass indices
% first correct for left/right, the convention is right
if ~skipCompass
  comInd = false(length(Loc_hem),1);
  allCompassInds=repmat({comInd},1,length(allLocs));
  for k=1:length(Loc_hem)
    thisComp = Loc_com(k,:);
    if strcmp(upper(Loc_hem(k)),'L')
      if thisComp(2)==1
	thisComp(2)=3;
      elseif thisComp(2)==3
	thisComp(2)=1;
      end
    end
    
    %if one compass index
    if size(loc2,1)==1
      comInd(k) = ismember(thisComp,loc2,'rows');   
    %if multiple compass indices and/or multipe locations
    %generates booleans for each location (e.g. PHC)  and assoc'd compass locations
    else
      
      for l=1:length(allLocs)
	for j=1:size(allCompassGroups{l})
	  hold = ismember(thisComp,allCompassGroups{l}(j,:),'rows');
	  if hold==1
	    allCompassInds{l}(k)=true; %triggers 'yes' in given compass group if one of compass locations matches
	    break
	  end
	end
      end
    end
  end
else
  comInd = true(length(Loc_hem),1);
end


% final indices
if ~exist('allLocs','var')
  finalInd = locInd & disInd & comInd;
else
  compileBooleans=false(length(Loc_hem),length(allLocs));
  for u=1:length(allLocs)
    compileBooleans(:,u)=allLocs{u}&disInd&allCompassInds{u};
  end
  for r=1:length(Loc_hem)
    if sum(compileBooleans(r,:))>=1
      finalInd(r)=true;
    else
      finalInd(r)=false;
    end
  end
end


if ~isempty(Loc_dis)
  disOUT   = Loc_dis(finalInd);
else
  disOUT   = [];
end

elsOUT=[];
jacOUT={};
for k=1:length(finalInd)
  if ~finalInd(k);continue;end
  thisJacNam = Loc_nam{k};
  thisJacNum = Loc_num(k);
    
  thisNamInd = strcmp(jac_nam,thisJacNam);
  thisNumInd = jac_num==thisJacNum;
  thisInd    = find(thisNamInd & thisNumInd);
  if length(thisInd)~=1
    disp((sprintf('%s%d not found',thisJacNam,thisJacNum)))
  end
  elsOUT = cat(1,elsOUT,els_num(thisInd));  
  jacOUT = cat(1,jacOUT,JAC{2}(thisInd));
end


  
function checkLOCS_local(LOC)
  unLOC = upper(unique(LOC));
  possPositions      = {'HIPP','HS','PHC','CS','FG','ITG','MTG','STG'};
  is_a_possible_name = ismember(unLOC,possPositions);
  if sum(~is_a_possible_name)>0
    error('there is an impossible name')
  end