function analyzePyfrStim2(subj,elStimNames,elStimSites)
%
%
%
%
%
%  elStimNames = {'Entorhinal Cortex','Hippocampus'};
%  elStimSites = {{'LSTA1,2 (EC)'},...
%	       {'pos=LAH2,LMH1 neg=LMH2,LPH1','pos=LPH1,2 neg=LMH1,2'}};
%
%

% default
if ~ismac
  rootDir='';
else
  rootDir='/Volumes/RHINO_root';
end

% subject
exp          = 'pyFR_stim2';
subjFileName =  sprintf('%s_events.mat',subj);
subjFile     = fullfile(rootDir,'/data/events',exp,subjFileName);
if ~exist(subjFile,'file')
  fprintf('\n\n\tEvents file is not made\n\n')
  return
end

% load events
ev = load(subjFile);
ev = ev.events;

% get unique conditions
evShamTrials    = ev([ev.stimcode]==5);
evBadTrials     = ev([ev.stimcode]==6);

% get the sham baseline
[frt_b numWords_b numCorrs_b] = getAllTimeToFirstRecall(evShamTrials); 

% loop thorugh by electrode
elecStimulated = unique({ev.stim_elec})';


for k=1:length(elStimNames)
  thisName   = elStimNames{k};
  thisRegion = elStimSites{k};
  if length(thisRegion)==1
    elEv = ev(strcmp({ev.stim_elec},thisRegion{1}));
  else
    elEv = ev(strcmp({ev.stim_elec},thisRegion{1})|...
	      strcmp({ev.stim_elec},thisRegion{2}));
  end

  evEncodeTrials  = elEv([elEv.stimcode]==1);
  evStorageTrials = elEv([elEv.stimcode]==2);
  evRetrievTrials = elEv([elEv.stimcode]==3);
  evMiddRetTrials = elEv([elEv.stimcode]==4);
  
  [frt_e numWords_e numCorrs_e] = getAllTimeToFirstRecall(evEncodeTrials); 
  [frt_s numWords_s numCorrs_s] = getAllTimeToFirstRecall(evStorageTrials); 
  [frt_r numWords_r numCorrs_r] = getAllTimeToFirstRecall(evRetrievTrials); 
  [frt_m numWords_m numCorrs_m] = getAllTimeToFirstRecall(evMiddRetTrials); 
   
  % plot the accuracies
  figure(1)
  barVect_1 = [numCorrs_b numCorrs_e numCorrs_s numCorrs_r];
  barVect_2 = [numWords_b numWords_e numWords_s numWords_r];
  barVect   = barVect_1./barVect_2;
  bar(1:length(barVect_1),barVect_1./barVect_2)
  for foo=1:length(barVect_1)
    text(foo-.3,barVect(foo)-.02*barVect(foo),...
	 sprintf('%d/%d',barVect_1(foo),barVect_2(foo)),...
	 'BackgroundColor',[.5 .5 .5],'FontWeight','Bold','FontSize',15,...
	 'FontName','courier','HorizontalAlignment','left',...
	 'VerticalAlignment','top')
  end
  set(gca,'XTickLabel',{'bas','enc','stor','ret'})
  ylim([.55 .8])
  xlim([.5 4.5])
  formatPic(gca)
  xlabel('stim condition')
  ylabel('Recall Accuracy')
  title(thisName)
  box off
  grid on
    
  % anova on the time  
  D = [frt_b; frt_e; frt_s; frt_r;];
  G = [repmat({'bas'},size(frt_b)); repmat({'enc'},size(frt_e));...
       repmat({'stor'},size(frt_s)); repmat({'ret'},size(frt_r))];
  p = anova1(D,G,'DISPLAYOPT','off');
  b = boxplot(D,G);
  formatPic(gca)
  set(gca,'fontSize',12)
  xlabel('stim condition')
  ylabel('Time to first recall (ms)')
  title(thisName)
  box off
  
  keyboard
  
end


function [frt numWords numCorrs] = getAllTimeToFirstRecall(E)
  unSess  = unique([E.session]);
  frt     = [];
  numTrial = 0;
  numWords = 0;
  numCorrs = 0;
  for s=1:length(unSess)
    E_sess = E([E.session]==unSess(s));
    
    unList = unique([E_sess.list]);
    for l=1:length(unList)
      numTrial      = numTrial + 1;
      E_sessList    = E_sess([E_sess.list]==unList(l));
      encEv         = E_sessList(strcmp({E_sessList.type},'WORD'));
      
      wordsThisList = {encEv.item};       
      corrThisList  = [encEv.recalled];
      numWords      = numWords + length(corrThisList);
      numCorrs      = numCorrs + sum(corrThisList);
      
      recStartEvInd = find(strcmp({E_sessList.type},'REC_START'));      
      recStartEv = E_sessList(recStartEvInd);
      
      for count=1:1000
	if count==1000; error('stuck in a loop');end
	nextWordEv = E_sessList(recStartEvInd+count);
	if strcmp(nextWordEv.type,'REC_STOP')
	  fprintf('no recalled words logged: sess %d, list %d\n',...
		  unSess(s),unList(l))
	  break
	end		 
	if ismember(nextWordEv.item,wordsThisList)
	  frt = cat(1,frt,nextWordEv.mstime-recStartEv.mstime);
	  break
	end	
      end % end for loop
      
    end        
  end