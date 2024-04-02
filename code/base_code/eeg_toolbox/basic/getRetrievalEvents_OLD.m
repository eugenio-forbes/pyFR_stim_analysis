function [recEv]= getRetrievalEvents_NEW(ev,isoTimeRecEv,recCntlIsoTime,rSeed,isScalp);
%
% this function puts some logic in to make sure that only recalled
% events are included if thye do not have anther recalled events
% 'minNumSec' to either side of them.
%
% It also gets a control condition for the recalled events.  The
% control condition is the same number of events per list except
% placed at random intervals throught the retrieval interval
%

if ~exist('isScalp','var')||isempty(isScalp)
  isScalp = false;
end

if ~exist('rSeed','var')||isempty(rSeed)
  thisIsTheSeed = sum(100*clock);
else
  thisIsTheSeed = rSeed;
end
s = RandStream.create('mt19937ar','seed',thisIsTheSeed);
RandStream.setGlobalStream(s);

% get the subject
thisSubj=unique({ev.subject});
useThisSubj=thisSubj{1};
if ~isScalp
  HospID = useThisSubj(1:2);
  ptNum  = str2double(useThisSubj(3:5));
else
  HospID = 'LTP';
  ptNum  = str2double(useThisSubj(4:5));
end

% These patients' rec_start field has -999 and not the list
foo1 = strcmp(useThisSubj(1:2),'TJ') & ptNum<=31;
foo2 = (strcmp(useThisSubj(1:2),'UP') & (ptNum>=11 & ptNum<=22));
has999InRecStartListField = foo1 | foo2;
sess   = [ev.session];
unSess = unique(sess);

recEv_out  = [];
recEv_cntl = [];

if ~isScalp
  rec_period_sec = 45; % length of retrieval period in seconds
else
  rec_period_sec = 75; % length of retrieval period in seconds
end

for s=1:length(unSess)
  thisSess       = unSess(s);
  evThisSess     = ev([ev.session]==thisSess);
  if ~isScalp
    listThisSess   = [evThisSess.list];
  else
    listThisSess   = [evThisSess.trial];
  end
  unListThisSess = unique(listThisSess);
  
  % a pace holder for the name of the eeg file for this session    
  thisEEGfile = [];    
      
  % in ideal world, we would just read the list from the rec_start
  % events, but many of the patients in the database have -999 in
  % the list field for the rec_start events.      
  if has999InRecStartListField
    listTmpCounter=0;
  end
  
  % loop thru all lists in this trial
  for l=1:length(unListThisSess)
    thisList          = unListThisSess(l);
    %fprintf(' %d %d\n',thisSess,thisList)
    
    % lists where they did not recall anything  
    % TJ003, sess 0, list 15
    
    % get all the events for this list.
    if ~isScalp
      evThisSessList    = evThisSess([evThisSess.list]==thisList);
    else
      evThisSessList    = evThisSess([evThisSess.trial]==thisList);
    end
    worEvThisSessListInd = strcmp({evThisSessList.type},'WORD');
    recEvThisSessListInd = strcmp({evThisSessList.type},'REC_WORD');
    vocEvThisSessListInd = strcmp({evThisSessList.type},'REC_WORD_VV');
    
    % get all the times that the patient opened his/her mouth and
    % something came out.
    worEvThisSessList = evThisSessList(worEvThisSessListInd);
    vocEvThisSessList = evThisSessList(vocEvThisSessListInd | recEvThisSessListInd);
  
    % see if there were any words for this list
    if isempty(worEvThisSessList)
      fprintf(' No words for session %d list %d\n',thisSess,thisList)
      continue
    end
    
    % get the start of the recall period.  If the REC_START field
    % is in there (because it has the correct list number in the
    % list field) then use that.  If this subjevty has -999 in the
    % list field for REC_START events, we have to find that differently.
    if ~has999InRecStartListField
      staEvThisSessListInd = strcmp({evThisSessList.type},'REC_START');
      staEvThisSessList    = evThisSessList(staEvThisSessListInd);
      if isempty(staEvThisSessList)
	fprintf('\n\n\tERROR!! I could not find start time\n\n')
	continue
	error('')
      end      
    else
      % there is no rec_start event in evThisSessList.  So you have
      % to find it.
      lastWordPresentedInd   = find([evThisSess.mstime]==worEvThisSessList(end).mstime);
      if length(lastWordPresentedInd)>1; error('should never happen');end
      potentialRecStartEvInd = lastWordPresentedInd  + 1;      
      if potentialRecStartEvInd>length(evThisSess)
	fprintf(' \n\n\tThere was no REC_START or REC events for sess %d list %d\n\n',...
		thisSess,thisList)
	continue
      end
      potentialRecStartEv    = evThisSess(potentialRecStartEvInd);
	
      % check that it is the start time foir this list
      if ~strcmp(potentialRecStartEv.type,'REC_START')
      	fprintf('\n\ntERROR!! I tried, but could not find start time\n\n')
	continue
	error(' ')
      end
      staEvThisSessList = potentialRecStartEv;       
    end
    staSamThisSessList = staEvThisSessList.eegoffset; 
    
        
  
    
    % get THE EEG FROM AN ARBITRARY ELECTRODE FOR THIS FILE
    if isempty(thisEEGfile) | ~strcmp(thisEEGfile,staEvThisSessList.eegfile);
      thisEEGfile = staEvThisSessList.eegfile;
      testElec    = 1;
      
      
      while true
	try
	  if testElec>=100
	    error(' this is jacked up')
	  end
	  EEG  = look(staEvThisSessList.eegfile,testElec,[],true); 
	  break
	catch
	  testElec = testElec + 1;
	  continue
	end
      end
    end
            
    % get the end point of this retrieval period and make sure that
    % it is within the EEG.
    thisRs             = GetRateAndFormat(thisEEGfile);
    rec_period_sam     = rec_period_sec*thisRs;  %(sec/recPeriod * sam/sec)  
    endSamThisSessList = staSamThisSessList + rec_period_sam;
    if endSamThisSessList>length(EEG)
      fprintf(' the eeg does not contain the end list %d, sess %d',thisList,thisSess)
      continue
    end        
        
    % now filter the events. get rid of all events that are not
    % properly isolated.
    endRecPeriodTime_ms = staEvThisSessList.mstime + rec_period_sec*1000;
    allTimes_ms         = [staEvThisSessList.mstime ... start of recall
		          [vocEvThisSessList.mstime] ... %rec/voc times
		          endRecPeriodTime_ms]; % end of recall		   
    timeBeforeEvent_ms  = allTimes_ms - [-Inf allTimes_ms(1:end-1)];
    timeAfterEvent_ms   = [allTimes_ms(2:end) Inf] - allTimes_ms;
    allVocIntervals_ms  = [timeBeforeEvent_ms(2:end-1)' timeAfterEvent_ms(2:end-1)'];
    wellIsolatedVocInd  = allVocIntervals_ms(:,1)>isoTimeRecEv(1) & ...
	                  allVocIntervals_ms(:,2)>isoTimeRecEv(2);
    wellIsolatedVocEv   = vocEvThisSessList(wellIsolatedVocInd);

    % now go back and clean the isolated events
    theseFinalRecEvents = wellIsolatedVocEv(strcmp({wellIsolatedVocEv.type},'REC_WORD'));            	    
    isVOCThisSessList   = strcmp({theseFinalRecEvents.item},'<>') | ...
                      	  strcmp({theseFinalRecEvents.item},'VV') | ...
	                  strcmp({theseFinalRecEvents.item},'!');    
    theseFinalRecEvents(isVOCThisSessList) = [];                        
    numFinalRecEvents   = length(theseFinalRecEvents);
    recEv_out           = cat(2,recEv_out,theseFinalRecEvents);   
        
    % now get the control events
    samNumbersInThisRetPeriod = staSamThisSessList:endSamThisSessList;
    msTimesInthisRetPeriod    = linspace(staEvThisSessList.mstime,...
					 endRecPeriodTime_ms,...
					 length(samNumbersInThisRetPeriod));
    dummyVariable             = zeros(1,length(samNumbersInThisRetPeriod));
    for foo=1:length(vocEvThisSessList)
      thisOff = vocEvThisSessList(foo).eegoffset;
      [~,ind]=min(abs(samNumbersInThisRetPeriod-thisOff));
      dummyVariable(ind)=1;
    end
    
    % now loop through and see if you can find some control events,  They should 
    % be isolated  by the isoTimes
    isoTimeSam_b    = floor(isoTimeRecEv(1)/1000*thisRs);
    isoTimeSam_f    = floor(isoTimeRecEv(2)/1000*thisRs);
    [~,randSamples] = sort(rand(1,length(dummyVariable)));

    % now get the events
    %clf;plot(dummyVariable);ax=gca;
    %dot = line('Parent',ax,'XData',0,'YData',0.5,'Marker','*','MarkerSize',10,'Color','r'); 
    %dummyVariable_orig = dummyVariable;hold on
    
    theseFinalRecControlEv = [];    
    numCntlEvThisList      = 0;
    for foo2=1:length(randSamples)
      thisSample = randSamples(foo2);
      if thisSample+isoTimeSam_f>length(dummyVariable);continue;end
      if thisSample-isoTimeSam_b<1;continue;end
	    
      indToTest    = thisSample-isoTimeSam_b:thisSample+isoTimeSam_f;         
      if sum(dummyVariable(indToTest))==0
	% marke that you got it	
	dummyVariable(indToTest)=1;	
	numCntlEvThisList = numCntlEvThisList + 1;
	
	%set(dot,'XData',thisSample);	
	%plot(dummyVariable,'Color','g','Marker','*')
	%plot(dummyVariable_orig,'Color','b','Marker','*')
	%fprintf('I got an event\n')	

	this_ms_time = floor(msTimesInthisRetPeriod(thisSample));
	
	theseFinalRecControlEv(end+1).subject=useThisSubj;
	theseFinalRecControlEv(end).session=thisSess;
	theseFinalRecControlEv(end).eegfile=thisEEGfile;
	theseFinalRecControlEv(end).eegoffset=samNumbersInThisRetPeriod(thisSample);
	theseFinalRecControlEv(end).list=thisList; 	
	theseFinalRecControlEv(end).mstime=this_ms_time; 	
	theseFinalRecControlEv(end).listStartTime=staEvThisSessList.mstime;
	theseFinalRecControlEv(end).rectime=this_ms_time-staEvThisSessList.mstime;
	%theseFinalRecControlEv(end).mstime=
      end  
      
    end
    %hold off;pause
    
    
    % re-oreder the rec_contl events
    if ~isempty(theseFinalRecControlEv)          
      times_temp_foo   = [theseFinalRecControlEv.eegoffset];
      [~,ind_temp_foo] = sort(times_temp_foo);
      theseFinalRecControlEv = theseFinalRecControlEv(ind_temp_foo);
    end
    
    % concatenate the cntl events
    recEv_cntl = cat(2,recEv_cntl,theseFinalRecControlEv);
    
    
  end % lists  
end % sessions

isZERO = [recEv_out.intrusion]==0;
isNEG1 = [recEv_out.intrusion]==-1;
isPLI  = [recEv_out.intrusion]>0;
isVOC  = strcmp({recEv_out.item},'<>') | ...
	 strcmp({recEv_out.item},'VV') | ...
	 strcmp({recEv_out.item},'!');
	 
recEv.rec     = recEv_out(isZERO);
recEv.pli     = recEv_out(isPLI);
recEv.xli     = recEv_out(isNEG1 & ~isVOC);
recEv.recCntl = recEv_cntl;

% filter by the closest recCntlEvent
recCtl_final = [];
neighborhood = 3000;
recCtl_copy  = recEv.recCntl;
maxnum       = 2;
for k=1:length(recEv.rec)
  timeDiff                 = abs([recCtl_copy.rectime]-recEv.rec(k).rectime);
  [sortedTimeDiff sortInd] = sort(timeDiff);
  indToKeep                = [];
  for k=1:min(maxnum,length(sortedTimeDiff))
    thisTimeDiff = sortedTimeDiff(k);
    thisOrigind  = sortInd(k);
    if thisTimeDiff<=neighborhood
      indToKeep = cat(2,indToKeep,thisOrigind);
    end    
  end              
  recCtl_final = cat(2,recCtl_final,recCtl_copy(indToKeep));
  recCtl_copy(indToKeep) = [];
  if isempty(recCtl_copy)
    break
  end
end

%indInNeighborhood_un = unique(indInNeighborhood_all);
%recCtl_final         = recCtl_copy(indInNeighborhood_un);

recEv.CntlOld = recEv.recCntl;
recEv.recCntl = recCtl_final;


%--------------------------------------------------------------------------------
function [sta_ms sta_off ] = getRetrievalStart_local(recStart,recEv);
  sta_ms_tmp  = recStart.mstime;
  sta_off_tmp = recStart.eegoffset;
  
  for k=1:length(recEv)
    checkThisTime = recEv(k).mstime-recEv(k).rectime;
    if ~isequal(sta_ms_tmp,checkThisTime)
      myErr('this is a double check that failed')
    end    
  end
  sta_ms=sta_ms_tmp;
  sta_off=sta_off_tmp;