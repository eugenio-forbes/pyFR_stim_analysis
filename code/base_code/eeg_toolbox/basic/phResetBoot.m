 function phResetBoot(type,chan,dur,off,samprate,band,test,shuffles,compare)


    left_subj_array =  {'087','088','092','101','103','096'};
     
    %left_subj_array =  {'088'};
     
    right_subj_array = {'087','092','094','100','103','096'};
    
    %right_subj_array = {'094', '103'};


startStr = num2str(off);
endStr = num2str(off+dur);

for m = left_subj_array
        

    channel = chan;
    
    m = m{1};
    subj = sprintf('%s',m);
    side = 'L';

    tsubj = subj;
    
    lockDir = fullfile('~/DBS_events/kareem_analysis/sigtests',['lockDir-' side test band],compare);

    if ~exist(lockDir,'dir')
        mkdir(lockDir)
        
    else
        foo = 1;
    end
    
    
    
    lockFile = fullfile(lockDir,[tsubj '-LOCK']);
   if ~exist(lockFile,'file')
     system(['touch ' lockFile]);
   else
     fprintf('  Skipping %s: b/c it''s locked.\n', subj)
    continue;
   end
   
    subjdir = fullfile('~/DBS_data/kareem_analysis/sigtests/RedoSN/phaseReset',date,[side band test subj]);
    
    if ~exist(subjdir,'dir')
        mkdir(subjdir)
        
    else
        foo = 1;
    end
    
    saveFile = fullfile(subjdir,[startStr '-' endStr 'phReset31.mat']);
    
    beh_load_path = sprintf('/data/continuous/CABG/data/beh/CABG%sor%s%s.mat',subj,side,type);
    load(beh_load_path);
    
    beh_events = behData;
    
    events_load_path = sprintf('~/DBS_events/kareemCABG/CABG%sor%s%sevents.mat',subj,side,type);
    
    x = load(events_load_path);
    
    events = x.all_events;
    
     xx = events;
        
    Rs = GetRateAndFormat(xx(1));

    nyquist_freq = floor(Rs/2);

    freqs = eeganalparams('freqs');
    valid_freqs_ind = find(freqs<nyquist_freq);

    valid_freqs = freqs(valid_freqs_ind);
    
    lineNoise = [59 61];
    
      stim_events = filterStruct(events,'feedback==0');
    feed_events = filterStruct(events,'feedback==1');
 
    %only for the special phase reset area
    
   % feed_events = stim_events;
    
  corrVect = [feed_events.corr]==1;
  errVect = [feed_events.corr]==0;
  unWinVect = beh_events.unexpWin==1;
  unLossVect = beh_events.unexpLoss==1;
    expWinVect =  beh_events.expWin==1;
  expLossVect = beh_events.expLoss==1;
 
  unVect = unWinVect | unLossVect;
  expVect = expWinVect | expLossVect;
  
    alphaPhase = gethilbertphase(chan,feed_events,dur,off,1000,[8 14],60,samprate);
    thetaPhase = gethilbertphase(chan,feed_events,dur,off,1000,[4 8],60,samprate);
    betaPhase = gethilbertphase(chan,feed_events,dur,off,1000,[16 24],60,samprate);
    deltaPhase = gethilbertphase(chan,feed_events,dur,off,1000,[1.5 4],60,samprate);
    
    
       if strcmp(band,'alpha')
        testPhase = alphaPhase;
       elseif strcmp(band,'theta')
        testPhase = thetaPhase;
       elseif strcmp(band,'beta')
        testPhase = betaPhase;
       elseif strcmp(band,'delta')
         testPhase = deltaPhase;
        
       end
 
       if strcmp(test,'corr')
       
       testPhase = testPhase(corrVect,:);
       
       elseif strcmp(test,'err')
       
       testPhase = testPhase(errVect,:);
       
       elseif strcmp(test,'unCorr')
       
       testPhase = testPhase(unWinVect,:);
       
       elseif strcmp(test,'unErr')
       
       testPhase = testPhase(unLossVect,:);
       
       elseif strcmp(test,'expCorr')
       
       testPhase = testPhase(expWinVect,:);
       
       elseif strcmp(test,'expErr')
       
       testPhase = testPhase(expLossVect,:);
       
       elseif strcmp(test, 'allFB')
           
        foo = 1;
        
        elseif strcmp(test, 'allUn')
           
        testPhase = testPhase(unVect,:);
           
       elseif strcmp(test,'allExp')
           
         testPhase = testPhase(expVect,:);
        
        
       end
       
       
   psum = [];
   barSum = [];
   
   
   for n = 1:size(testPhase,2)
       
       [p1 rBar] = rayleigh(squeeze(testPhase(:,n)));
       
       psum(n) = p1;
       barSum(n) = rBar;
       
   end
       
   allpFoo = [];
   
   %barSumFoo = [];
   
   for r = 1:shuffles
       
       fooInd = [];
       
       fooInd = testPhase(randperm(size(testPhase,1)*size(testPhase,2)));
       
       
       fooInd = reshape(fooInd,size(testPhase,1),size(testPhase,2));
       
       psumFoo = [];
       
       for n = 1:size(fooInd,2)
       
       [pfoo, ~] = rayleigh(squeeze(fooInd(:,n)));
       
       psumFoo(n) = pfoo;
       
       end
       
       allpFoo(r,:) = psumFoo;
       
   end
   
   %save(saveFile,'allpFoo','psum');
   
  if exist(lockFile,'file')
     system(['rm ' lockFile]);
  end
 
  keyboard
  
  end


for mm = right_subj_array
           
        
    %subj_counter = subj_counter+1;
    
    side = 'R';
    
    channel = chan;
    
    mm = mm{1};
    subj = sprintf('%s',mm);

    tsubj = subj;
    
        lockDir = fullfile('~/DBS_events/kareem_analysis/sigtests',['lockDir-' side test band],compare);

    if ~exist(lockDir,'dir')
        mkdir(lockDir)
        
    else
        foo = 1;
    end
    
    
    
    lockFile = fullfile(lockDir,[tsubj '-LOCK']);
   if ~exist(lockFile,'file')
     system(['touch ' lockFile]);
   else
     fprintf('  Skipping %s: b/c it''s locked.\n', subj)
    continue;
   end
   
    subjdir = fullfile('~/DBS_data/kareem_analysis/sigtests/RedoSN/phaseReset',date,[side band test subj]);
    
    if ~exist(subjdir,'dir')
        mkdir(subjdir)
        
    else
        foo = 1;
    end
    
    saveFile = fullfile(subjdir,[startStr '-' endStr 'phReset31.mat']);
    
    beh_load_path = sprintf('/data/continuous/CABG/data/beh/CABG%sor%s%s.mat',subj,side,type);
    load(beh_load_path);
    
    beh_events = behData;
    
    events_load_path = sprintf('~/DBS_events/kareemCABG/CABG%sor%s%sevents.mat',subj,side,type);
    
    x = load(events_load_path);
    
    events = x.all_events;
    
     xx = events;
    
    
    Rs = GetRateAndFormat(xx(1));

    nyquist_freq = floor(Rs/2);

    freqs = eeganalparams('freqs');
    valid_freqs_ind = find(freqs<nyquist_freq);

    valid_freqs = freqs(valid_freqs_ind);
    
    lineNoise = [59 61];
    
      stim_events = filterStruct(events,'feedback==0');
    feed_events = filterStruct(events,'feedback==1');
    
    %only for the special phase reset area
    
   % feed_events = stim_events;
    
  corrVect = [feed_events.corr]==1;
  errVect = [feed_events.corr]==0;
  unWinVect = beh_events.unexpWin==1;
  unLossVect = beh_events.unexpLoss==1;
    expWinVect =  beh_events.expWin==1;
  expLossVect = beh_events.expLoss==1;
    
    unVect = unWinVect | unLossVect;
  expVect = expWinVect | expLossVect;
  
    alphaPhase = gethilbertphase(chan,feed_events,dur,off,1000,[8 12],60,samprate);
    thetaPhase = gethilbertphase(chan,feed_events,dur,off,1000,[4 8],60,samprate);
    betaPhase = gethilbertphase(chan,feed_events,dur,off,1000,[16 24],60,samprate);
    deltaPhase = gethilbertphase(chan,feed_events,dur,off,1000,[1.5 4],60,samprate);
    
       if strcmp(band,'alpha')
        testPhase = alphaPhase;
       elseif strcmp(band,'theta')
        testPhase = thetaPhase;
       elseif strcmp(band,'beta')
        testPhase = betaPhase;
       elseif strcmp(band,'delta')
         testPhase = deltaPhase;
        
       end
 
       if strcmp(test,'corr')
       
       testPhase = testPhase(corrVect,:);
       
       elseif strcmp(test,'err')
       
       testPhase = testPhase(errVect,:);
       
       elseif strcmp(test,'unCorr')
       
       testPhase = testPhase(unWinVect,:);
       
       elseif strcmp(test,'unErr')
       
       testPhase = testPhase(unLossVect,:);
       
       elseif strcmp(test,'expCorr')
       
       testPhase = testPhase(expWinVect,:);
       
       elseif strcmp(test,'expErr')
       
       testPhase = testPhase(expLossVect,:);
       
       elseif strcmp(test, 'allFB')
           
           foo = 1;
           
       elseif strcmp(test, 'allUn')
           
           testPhase = testPhase(unVect,:);
           
       elseif strcmp(test,'allExp')
           
           testPhase = testPhase(expVect,:);
           
        
        
       end
       
       
   psum = [];
   barSum = [];
   
   
   for n = 1:size(testPhase,2)
       
       [p1 rBar] = rayleigh(squeeze(testPhase(:,n)));
       
       psum(n) = p1;
       barSum(n) = rBar;
       
   end
       
   allpFoo = [];
   
   %barSumFoo = [];
   
   for r = 1:shuffles
       
       fooInd = [];
       
       fooInd = testPhase(randperm(size(testPhase,1)*size(testPhase,2)));
       
       
       fooInd = reshape(fooInd,size(testPhase,1),size(testPhase,2));
       
       psumFoo = [];
       
       for n = 1:size(fooInd,2)
       
       [pfoo, ~] = rayleigh(squeeze(fooInd(:,n)));
       
       psumFoo(n) = pfoo;
       
       end
       
       allpFoo(r,:) = psumFoo;
       
   end
   
   %save(saveFile,'allpFoo','psum');
   
  if exist(lockFile,'file')
     system(['rm ' lockFile]);
  end
 
  keyboard  
  
end


   
       
       
       
       
       
       
       