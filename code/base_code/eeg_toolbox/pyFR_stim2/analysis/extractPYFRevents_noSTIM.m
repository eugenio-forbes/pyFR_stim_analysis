function events=extractPYFRevents_noSTIM(subject,expDir,session,forceSESSION)
%
% FUNCTION:
%   events=extractPYFRevents(subject,expDir,session,forceSESSION)
% 
% DECRIPTION:
%   extracts the events associated with pyFR.
%
% INPUTS:
%   subject......... 'TJ055'
%   expDir.......... '/data/eeg/TJ055/behavioral/pyFR_stim2/TJ055_Maddox'
%   session......... 0 = searches 'session_0' in expDir
%   forceSESSION.... [optional] 1 = forces session to this number
%
% OUTPUTS:
%   events....... the events structure
%
% NOTES:
%   (1) written by jfburke on 02/2013 (john.fred.burke@gmail.com)
%   (2) intrusion code
%         vocalization = -999;
%         XLI          = -1;
%         PLI          = integer of number of lists back;
%         correct      = 0;
%

clear global
global SUBJECT SESSION LIST events IMAGERY

SUBJECT = subject;
SESSION = session;
LIST    = -999;
IMAGERY = 'X';

thisSessDir = sprintf('session_%d',SESSION);
sessFile    = fullfile(expDir,thisSessDir,'session.log');

fid = fopen(sessFile,'r');
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
numTrials = getExpInfo_local(expDir,'numTrials');
numSess   = getExpInfo_local(expDir,'numSessions'); 
NOUNPOOL  = getNounPool_local(expDir);

evCounter         = 0;
events            = [];
sessLineCounter   = 0; 

while true
  thisLine = fgetl(fid);
  if ~ischar(thisLine);return;end

  % get the third string before the underscore
  xTOT=textscan(thisLine,'%f%f%s');
  thisTYPE    = xTOT{3}{1};
  
  % based on the type write different fields for this event
  switch upper(thisTYPE)    
   %-------------------------------------------------------------------
   case {'B','E','SESS_END','ORIENT','REC_STUDY','REC_STOP',...
	 'ENODE_BLANK_START','ENODE_BLANK_END','START_MATH','END_MATH'}
    sessLineCounter = sessLineCounter + 1;
    x=textscan(thisLine,'%f%f%s');
    evCounter = evCounter + 1;
    mkNewEvent_local(evCounter,x{1},x{2});  
    appendNewEvent_local(evCounter,'type',x{3}{1});      
   %-------------------------------------------------------------------
   case 'SESS_START'
    sessLineCounter = sessLineCounter + 1;
    x=textscan(thisLine,'%f%f%s%d');
    thisSess  = x{4};
    if (thisSess-1)~=session
      error('sessions dont match');
    end
    if thisSess>numSess
      error('exceeds max session');
    end    
    evCounter = evCounter + 1;
    mkNewEvent_local(evCounter,x{1},x{2});  
    appendNewEvent_local(evCounter,'type',x{3}{1});       
   %-------------------------------------------------------------------
   case 'TRIAL'
    sessLineCounter = sessLineCounter + 1;
    x=textscan(thisLine,'%f%f%s%d%s%s');
    thisTRIAL = x{4};
    if thisTRIAL>numTrials
      error('exceeds max trials')
    end
    LIST = thisTRIAL;
    IMAGERY = x{5}{1};
    evCounter = evCounter + 1;
    mkNewEvent_local(evCounter,x{1},x{2});  
    appendNewEvent_local(evCounter,'type',x{3}{1});
   %-------------------------------------------------------------------
   case 'WORD'
    sessLineCounter = sessLineCounter + 1;
    x=textscan(thisLine,'%f%f%s%s%d');
    thisWORD = x{4}{1};
    thisSP   = x{5};
    evCounter = evCounter + 1;    
    itemNo   = find(strcmp(thisWORD,NOUNPOOL));
    mkNewEvent_local(evCounter,x{1},x{2});  
    appendNewEvent_local(evCounter,'type',x{3}{1});
    appendNewEvent_local(evCounter,'serialpos',thisSP+1);
    appendNewEvent_local(evCounter,'item',thisWORD);
    appendNewEvent_local(evCounter,'itemno',itemNo); 
    appendNewEvent_local(evCounter,'recalled',0); 
  %-------------------------------------------------------------------
   case 'REC_START'
    sessLineCounter = sessLineCounter + 1;
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

function mkNewEvent_local(evCounter,mstime,msoffset)
global SUBJECT SESSION LIST events IMAGERY

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
events(evCounter).imagery       = IMAGERY;

function appendNewEvent_local(evCounter,varargin)
global events
nVar = length(varargin)/2;
for v=1:nVar
  thisVarField = varargin{2*(v-1)+1};
  thisVarData  = varargin{2*(v-1)+2};
  events(evCounter)=setfield(events(evCounter),thisVarField,thisVarData);
end

function [out] = getExpInfo_local(expDir,str2get);
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

function [X] = getNounPool_local(expDir);
fid_foo1 = fopen(fullfile(expDir,'tor_words_only.txt'),'r');
if fid_foo1==-1;error('word-pool.txt not found');end
X = textscan(fid_foo1,'%s');
X = X{1};
%X=textread('word-pool.txt','%s');
%X=textread(fullfile(expDir,'word-pool.txt'),'%s');
fclose (fid_foo1);
