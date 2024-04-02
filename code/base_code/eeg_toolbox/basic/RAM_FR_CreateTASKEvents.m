function events=RAM_FR_CreateTASKEvents(subject,expDir,session,forceSESSION)
%
% FUNCTION:
%   events=RAM_FR_CreateTASKEvents(subject,expDir,session,forceSESSION)
% 
% DECRIPTION:
%   extracts the events associated with pyFR.
%
% INPUTS:
%   subject......... 'FR200'
%   expDir.......... '/data/eeg/FR200/behavioral/pyFR/'
%   session......... 0 = searches 'session_0' in expDir
%   forceSESSION.... [optional] 1 = forces session to this number
%
% OUTPUTS:
%   events....... the events structure
%
% NOTES:
%   (1) written by jfburke on 10/2011 (john.fred.burke@gmail.com)
%   (2) intrusion code
%         vocalization = -999;
%         XLI          = -1;
%         PLI          = integer of number of lists back;
%         correct      = 0;
%

clear global
global SUBJECT SESSION LIST events stimParams currParamSet versionNum

SUBJECT = subject;
SESSION = session;
LIST    = -999;
currParamSet = 0;
versionNum = ''; 

thisSessDir = sprintf('session_%d',SESSION);
sessFile    = fullfile(expDir,thisSessDir,'session.log');

%^ use the latin1 coding scheme to handle freiburg patients, who
%have umlauts (I think that you need latin1 to take care of
%umlauts) john.fred.burke@gmail.com: June 19, 2013
fid = fopen(sessFile,'r','native','latin1');
if fid==-1
  fprintf('session %d..no session.log file found.\n',SESSION);   
  fprintf('EXITING\n\n');
  return
end

% you can change the session
if exist('forceSESSION','var') && ~isempty(forceSESSION)
  SESSION=forceSESSION;
end

% get experimental variables
LL        = 12;%getExpInfo_local(expDir,'listLen');
numTrials = 25;%getExpInfo_local(expDir,'numTrials');
numSess   = 18;%getExpInfo_local(expDir,'numSessions'); 
NOUNPOOL  = getNounPool_local(expDir);
evCounter = 0;
events    = [];

while true
  thisLine = fgetl(fid);
  if ~ischar(thisLine);return;end

  % get the third string before the underscore
  % Because I messed up session.log
  if thisLine(1)=='('
      xTOT = textscan(thisLine,'(%fL, %f)\t%*s\t%s','delimiter','\t');
  else
      xTOT=textscan(thisLine,'%f\t%f\t%s');
  end
  thisTYPE    = xTOT{3}{1};
  
  % based on the type write different fields for this event
  switch upper(thisTYPE)    
   %-------------------------------------------------------------------
   case {'B','E','SESS_END','ORIENT', 'STIM_ON', 'REC_END', 'DISTRACT_START', 'DISTRACT_END','SESSION_SKIPPED'}
    x=textscan(thisLine,'%f%f%s');
    evCounter = evCounter + 1;
    mkNewEvent_local(evCounter,x{1},x{2});  
    appendNewEvent_local(evCounter,'type',x{3}{1});      
   %-------------------------------------------------------------------
   case 'STIM_PARAMS'
    %keyboard
    x=textscan(thisLine,'%f%f%s%d%s%f%s%f%s%f');
    ParamSet = x{4};
    stimParams(ParamSet).Loc = [x{6};x{8}];
    %stimParams(ParamSet).Loc{1} = x{6}{1}; stimParams(ParamSet).Loc{2} = x{8}{1};
    stimParams(ParamSet).Amp = x{10};
   
   %-------------------------------------------------------------------
   case 'SESS_START'
    x=textscan(thisLine,'%f%f%s%d%s%s');
    thisSess  = x{4};
    if (thisSess-1)~=session
      error('sessions dont match');
    end
    if thisSess>numSess
      error('exceeds max session');
    end    
    versionNum = x{6}{1};
    evCounter = evCounter + 1;
    mkNewEvent_local(evCounter,x{1},x{2});  
    appendNewEvent_local(evCounter,'type',x{3}{1});       
   %-------------------------------------------------------------------
   case 'TRIAL'
       if thisLine(1)=='('
           x = textscan(thisLine,'(%fL, %f)\t%*s\t%s\t%d\t%s\t%s','delimiter','\t');
       else
           x=textscan(thisLine,'%f\t%f\t%s\t%d\t%s\t%s');
       end
    thisTRIAL = x{4};
    if thisTRIAL>numTrials
      error('exceeds max trials')
    end
    LIST = thisTRIAL;
    evCounter = evCounter + 1;
    if ~isempty(x{5})
        currParamSet = str2double(num2str(x{6}{1}));
    end
    mkNewEvent_local(evCounter,x{1},x{2});
    appendNewEvent_local(evCounter,'type',x{3}{1});

   %-------------------------------------------------------------------
   case 'WORD'
    x=textscan(thisLine,'%f%f%s%s%s%d%s');
    thisWORD = x{5}{1};
    thisSP   = x{6};
    if ~isempty(x{7})
        thisStim = x{7}{1};
    else
        thisStim = 'NOSTIM';
    end
    evCounter = evCounter + 1;
    itemNo   = find(strcmp(thisWORD,NOUNPOOL));
    if isempty(itemNo);
      error('word not found')
    end
    
    if strcmp(thisStim,'STIM')
        isStim = 1;
    else
        isStim = 0;
    end
    mkNewEvent_local(evCounter,x{1},x{2});  
    appendNewEvent_local(evCounter,'type',x{3}{1});
    appendNewEvent_local(evCounter,'serialpos',thisSP+1);
    appendNewEvent_local(evCounter,'item',thisWORD);
    appendNewEvent_local(evCounter,'itemno',itemNo); 
    appendNewEvent_local(evCounter,'recalled',0);
    
    appendNewEvent_local(evCounter,'isStim',isStim);

  %-------------------------------------------------------------------
   case 'REC_START'
    x=textscan(thisLine,'%f%f%s');
    evCounter = evCounter + 1;
    mkNewEvent_local(evCounter,x{1},-999);  
    appendNewEvent_local(evCounter,'type',x{3}{1});
       
    allWordsSoFarInd = strcmp({events.type},'WORD')';    
    allWordsSoFar    = {events(allWordsSoFarInd).item}';    
    allListsSoFar    = [events(allWordsSoFarInd).list]';
    allWordNumSoFar  = [events(allWordsSoFarInd).itemno]';
    wordsThisList    = allWordsSoFar(allListsSoFar==LIST);
    itemnoThisList   = allWordNumSoFar(allListsSoFar==LIST);
    
    % loop throught ht .ann file
    parseFileName = sprintf('%d',LIST-1);
    annFileName   = [parseFileName '.ann'];
    parFileName   = [parseFileName '.par']; 
    annFile       = fullfile(expDir,thisSessDir,annFileName);
    parFile       = fullfile(expDir,thisSessDir,parFileName);
    if ~exist(annFile)
      fid2=fopen(parFile,'r');
    else
      fid2=fopen(annFile,'r');
    end
    if fid2==-1;error('no .ann or .par files');end
    itemsIsaidThisList = [];
    while true
      tmpAnnLine=fgetl(fid2);
      if ~ischar(tmpAnnLine);break;end
      if numel(tmpAnnLine)==0;continue;end
      if strcmp(tmpAnnLine(1),'#');continue;end
      x2=textscan(tmpAnnLine,'%f%f%s');
      thisRT       = round(x2{1});
      thisWordNum  = x2{2};
      thisRecWord  = x2{3}{1};
      itemsIsaidThisList = cat(1,itemsIsaidThisList,[thisWordNum thisRT]);
      
      %---------------------------------------------------------------------
      % now see what kind of recall it is
      recalled = ismember(thisWordNum,itemnoThisList);
      thisRecallType = 'REC_WORD';      
      % is it correct?
      if recalled==1;
	intrusion=0;
	itemno=thisWordNum;
      else      
	% is it a vocalization?
	isVOC = strcmp(thisRecWord,'VV')||strcmp(thisRecWord,'<>')||strcmp(thisRecWord,'!');
	if isVOC
	  if thisWordNum~=-1;error('should never happen');end
	  intrusion = -999;
	  itemno=-999;
	  thisRecallType = 'REC_WORD_VV';
	else
	  % is it an XLI?
	  if thisWordNum==-1
	    intrusion = -1;
	    itemno=-1;
	  end	
	  % is it a PLI?
	  if thisWordNum>0
	    origInd  = find(strcmp(allWordsSoFar,thisRecWord));
	    if ~isempty(origInd)
	      origList = allListsSoFar(origInd);
	      intrusion = LIST - origList;
	      itemno    = thisWordNum;
	    else
	      %fprintf('LIST %2.0i: %-10.10s was not shown yet,  must be an XLI\n',LIST,thisRecWord)
	      intrusion = -1;
	      itemno=-1;
	    end
	  end % is a PLI 	  
	end % is a vocalization	
      end % is recalled    
            
      %---------------------------------------------------------------------
      % now add the event
      evCounter = evCounter + 1;
      mkNewEvent_local(evCounter,x{1}+thisRT,-999);  
      appendNewEvent_local(evCounter,'type',thisRecallType);
      appendNewEvent_local(evCounter,'item',thisRecWord);
      appendNewEvent_local(evCounter,'itemno',itemno);    
      appendNewEvent_local(evCounter,'intrusion',intrusion);
      appendNewEvent_local(evCounter,'rectime',thisRT);      
      %---------------------------------------------------------------------
      % now go back and add the rectime and recalled to the WORD
      % events (if correctly recalled)
      if recalled
	wordInd = [events.session]==SESSION & [events.list]==LIST ...
		  & strcmp({events.type},'WORD') & [events.itemno]==thisWordNum;
	
	% take first time they said it
	wordInd = find(wordInd);	
	if isempty(wordInd)
	  error('recalled word not on the list?')
	end
	% get attributes
	wordItemNo_tmp = events(wordInd).itemno;
	wordSess_tmp   = events(wordInd).session;
	wordList_tmp   = events(wordInd).list;	
	% check them
	if wordItemNo_tmp~=thisWordNum;error('wrong word');end
	if wordSess_tmp~=SESSION;error('wrong word');end
	if wordList_tmp~=LIST;error('wrong word');end
	% fill in the missing info	
	firstRTind = find(itemsIsaidThisList(:,1)==thisWordNum);
	if isempty(firstRTind);error('bad');end	  
	
	events(wordInd).recalled = true;
	events(wordInd).rectime = itemsIsaidThisList(firstRTind(1),2);

      end % was recalled      
    end % looping over recalls here   
   otherwise 
    error('type not recognized')
  end
     
end
keyboard
function mkNewEvent_local(evCounter,mstime,msoffset)
global SUBJECT SESSION LIST events stimParams currParamSet versionNum

events(evCounter).subject       = SUBJECT;
events(evCounter).session       = SESSION;
events(evCounter).list          = LIST;
events(evCounter).serialpos     = -999;
events(evCounter).type          = -999;
events(evCounter).item          = 'X';
events(evCounter).itemno        = -999;
events(evCounter).recalled      = -999;
events(evCounter).mstime        = mstime;
events(evCounter).msoffset      = msoffset;
events(evCounter).rectime       = -999;
events(evCounter).intrusion     = -999;
events(evCounter).isStim        = -999;
events(evCounter).expVersion    = versionNum;
events(evCounter).stimLoc       = 'X';
events(evCounter).stimAmp       = 'X';
if currParamSet==0
    events(evCounter).stimLoc=[nan; nan];
    events(evCounter).stimAmp=nan;
else
    events(evCounter).stimLoc = stimParams(currParamSet).Loc;
    events(evCounter).stimAmp = stimParams(currParamSet).Amp;
end

function appendNewEvent_local(evCounter,varargin)
global events
nVar = length(varargin)/2;
for v=1:nVar
  thisVarField = varargin{2*(v-1)+1};
  thisVarData  = varargin{2*(v-1)+2};
  events(evCounter)=setfield(events(evCounter),thisVarField,thisVarData);
end

function [out] = getExpInfo_local(expDir,str2get)
fid_foo1 = fopen(fullfile(expDir,'config.py'),'r');
while true
  thisLine = fgetl(fid_foo1);
  if ~ischar(thisLine);break;end
  if numel(thisLine)==0;continue;end
  if strcmp(thisLine(1),'#');continue;end
  possible_str=textscan(thisLine,'%s%f','Delimiter','=');
  X = regexprep(possible_str{1}{1},' ','');
  if strcmp(X,str2get)
    out=possible_str{2};
    break
  end
end
fclose (fid_foo1);

function [X] = getNounPool_local(expDir)
fid_foo1 = fopen(fullfile(expDir,'RAM_wordpool.txt'),'r','native','latin1');
if fid_foo1==-1;error('RAM_wordpool.txt not found');end
X = {};
while true
  tmp = fgetl(fid_foo1);
  if ~ischar(tmp);break;end
  X   = cat(1,X,{tmp});  
end
fclose (fid_foo1);
if ismember(upper('pass'),X)
  error('pass is in the nounpool!')
end