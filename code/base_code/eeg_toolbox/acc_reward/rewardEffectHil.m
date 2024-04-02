function rewardEffectHil(elecStruct,encodeDur,encodeOffset,binWidth,binNumber,kurtSkew,normVar,band)


%need to convert binWidth from msec into samples (assume rsamp of 250)
binWidth = floor(binWidth./4);

if strcmp(band,'low')
lowFreq=35; highFreq = 70;
elseif strcmp(band,'high')
lowFreq = 70; highFreq = 120;
elseif strcmp(band,'st')
    lowFreq = 2; highFreq = 4.5;
elseif strcmp(band,'nt')
    lowFreq = 4.5; highFreq = 9;
elseif strcmp(band,'alpha')
    lowFreq = 8;highFreq = 14;
elseif strcmp(band,'beta')
    lowFreq = 16; highFreq = 26;
end


subj_array = {elecStruct.subject};
subjCounter = 0;

for m = subj_array
    subjCounter = subjCounter+1;
    m = m{1};
    %m = m{1};
    %do this twice for odd reasons of the subject structure;
    subjname = sprintf('%s',m);
	
    
    subjname
    

    dataDir = fullfile('/Users/legab/iEEGDataBase/acc/Results/rewardHil',subjname);
 

		if ~exist(dataDir,'dir')
			mkdir(dataDir);
		else
			foo=1;
		end
	
    
    elecArray = elecStruct(subjCounter).electrodes;
    electvector = [];
    for n = 1:length(elecArray)
        electvector(n) = str2num(elecArray{n});
    end
    
      %break up the events according to type

    loadPath  = fullfile('/Volumes/SEAGATE/behData/acc_reward', [subjname '_events.mat']);
	load(loadPath,'events');
	
    %this gets rid of events for which there is no eegfile
    
    events = filterStruct(events,'~strcmp(eegfile,'''')');
    
    fbEvents = filterStruct(events,'strcmp(type,''FB'')');
    cueEvents = filterStruct(events,'strcmp(type,''CUE'')');
    
    
        %the usual steps for the sample rate
    
    lineNoise = getLineNoise(subjname);


    Rs = GetRateAndFormat(events(1));

    nyquist_freq = floor(Rs/2);

    freqs = eeganalparams('freqs');
    valid_freqs_ind = find(freqs<nyquist_freq);

    valid_freqs = freqs(valid_freqs_ind);
    freqs = valid_freqs;


	bufferMS = 500;
    
    elecCounter =  0;

    for l = electvector
	
    elecCounter = elecCounter +1;
    
    if kurtSkew == 0

			foo=1;


			elseif kurtSkew == 1
			

			d = gete_ms(l,fbEvents,encodeDur,encodeOffset,bufferMS);
			k = kurtosis(d,[],2); kThresh = 4;
			fbEvents = fbEvents(k<kThresh);
            
			d = gete_ms(l,cueEvents,encodeDur,encodeOffset,bufferMS);
			k = kurtosis(d,[],2); kThresh = 4;
			cueEvents = cueEvents(k<kThresh); 
            
    end
    
 
    posEvents = filterStruct(fbEvents,'feedback==1');
    negEvents = filterStruct(fbEvents,'feedback==0');
    neutEvents = filterStruct(fbEvents,'feedback==-1');
    
    posTopEvents = filterStruct(posEvents,'valence==2');
    posMidEvents = filterStruct(posEvents,'valence==1');
    posLowEvents = filterStruct(posEvents,'valence==0');
    
    cueTopEvents = filterStruct(cueEvents,'valence==2');
    cueMidEvents = filterStruct(cueEvents,'valence==1');
    cueLowEvents = filterStruct(cueEvents,'valence==0');
    cueNeutEvents = filterStruct(cueEvents,'valence==-1');
    
    %now get the power values for each event type
    
    
    [~,posPow] = gethilbertphase(l,posEvents,encodeDur,encodeOffset,bufferMS,[lowFreq highFreq],lineNoise,250);
    [~,negPow] = gethilbertphase(l,negEvents,encodeDur,encodeOffset,bufferMS,[lowFreq highFreq],lineNoise,250);
    [~,neutPow] = gethilbertphase(l,neutEvents,encodeDur,encodeOffset,bufferMS,[lowFreq highFreq],lineNoise,250);

    [~,posTop] = gethilbertphase(l,posTopEvents,encodeDur,encodeOffset,bufferMS,[lowFreq highFreq],lineNoise,250);
    [~,posMid] = gethilbertphase(l,posMidEvents,encodeDur,encodeOffset,bufferMS,[lowFreq highFreq],lineNoise,250);
    [~,posLow] = gethilbertphase(l,posLowEvents,encodeDur,encodeOffset,bufferMS,[lowFreq highFreq],lineNoise,250);

    [~,cueTop] = gethilbertphase(l,cueTopEvents,encodeDur,encodeOffset,bufferMS,[lowFreq highFreq],lineNoise,250);
    [~,cueMid] = gethilbertphase(l,cueMidEvents,encodeDur,encodeOffset,bufferMS,[lowFreq highFreq],lineNoise,250);
    [~,cueLow] = gethilbertphase(l,cueLowEvents,encodeDur,encodeOffset,bufferMS,[lowFreq highFreq],lineNoise,250);
    [~,cueNeut] = gethilbertphase(l,cueNeutEvents,encodeDur,encodeOffset,bufferMS,[lowFreq highFreq],lineNoise,250);

    %smooth all the rows
    
    posPow = smoothRows(posPow);
    negPow = smoothRows(negPow);
    neutPow = smoothRows(neutPow);
    
    posTop = smoothRows(posTop);
    posMid = smoothRows(posMid);
    posLow = smoothRows(posLow);
    
    cueTop = smoothRows(cueTop);
    cueMid = smoothRows(cueMid);
    cueLow = smoothRows(cueLow);
    cueNeut = smoothRows(cueNeut);
    
    
     %binWidth decides the size of the time windows to use for
             %the SME analysis binNumber is the number of step bins
             
             %convert binWidth to sample number based on resampled rate of
             %250
             binCenters = [];
          
             binCenters = floor(linspace(0,size(posPow,2),binNumber));
                          
             okBins= [];
             okBins = find(binCenters<(binCenters(end)-binWidth));
             binCenters = binCenters(1:okBins(end));
                 
             %convert binCenters so first element is 1
             binCenters(1) = 1;
             %convert binCenters to get rid of the last element
             binCenters = binCenters(1:end-1);
             
             
             posPowSmooth = [];
             negPowSmooth = [];
             neutPowSmooth = [];
             
             posTopSmooth = [];
             posMidSmooth = [];
             posLowSmooth = [];
             
             cueTopSmooth = [];
             cueMidSmooth = [];
             cueLowSmooth= [];
             cueNeutSmooth = [];
             
             stepCounter =  0;

             
             for q = binCenters
                 stepCounter = stepCounter+1;
                 
                 posPowSmooth(:,stepCounter) = squeeze(mean(posPow(:,q:q+binWidth),2));
                 negPowSmooth(:,stepCounter) = squeeze(mean(negPow(:,q:q+binWidth),2));
                 neutPowSmooth(:,stepCounter) = squeeze(mean(neutPow(:,q:q+binWidth),2));
                 
                 posTopSmooth(:,stepCounter) = squeeze(mean(posTop(:,q:q+binWidth),2));
                 posMidSmooth(:,stepCounter) = squeeze(mean(posMid(:,q:q+binWidth),2));
                 posLowSmooth(:,stepCounter) = squeeze(mean(posLow(:,q:q+binWidth),2));

                 cueTopSmooth(:,stepCounter) = squeeze(mean(cueTop(:,q:q+binWidth),2));
                 cueMidSmooth(:,stepCounter) = squeeze(mean(cueMid(:,q:q+binWidth),2));
                 cueLowSmooth(:,stepCounter) = squeeze(mean(cueLow(:,q:q+binWidth),2));
                 cueNeutSmooth(:,stepCounter) = squeeze(mean(cueNeut(:,q:q+binWidth),2));
                 
             end 
            
             
             
             subplot(2,1,1);
             hold off
             plot(mean(posTopSmooth,1));
             hold on
             plot(mean(posMidSmooth,1),'g');
             plot(mean(posLowSmooth,1),'y');
             plot(mean(negPowSmooth,1),'r');
             plot(mean(neutPowSmooth,1),'k');
             
             subplot(2,1,2);
             hold off
             plot(mean(cueTopSmooth,1));
             hold on
             plot(mean(cueMidSmooth,1),'g');
             plot(mean(cueLowSmooth,1),'y');
             plot(mean(cueNeutSmooth,1),'k');
             
             printPath = fullfile(dataDir,['rewardHil' band num2str(l) '.eps']);

             print(gcf,'-depsc2',printPath);

            
             
    end
end

    
    
