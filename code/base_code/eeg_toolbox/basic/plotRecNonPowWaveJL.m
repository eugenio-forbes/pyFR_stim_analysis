clc
clear

%Wrote plotRecNonowWave again so that it's no longer a function
encodeDur = 1800;
encodeOffset = 0;
kurtSkew = 1;
normVar = 1;

%loads single patient file
% load('UT002_events.mat')

loadPath  = fullfile('/Users/legab/ExperimentData/behData/pyFR_Stim/UT005a/session_0/eventsStim0.mat');
load(loadPath);

%this changes the directory of the events to for my mac, comment it out and
%it should work for you
%events = arrayfun(@(s) setfield(s,'eegfile','/Users/JL/Documents/MATLAB/eeg_data/UT002/eeg.reref/UT002_17Nov14_1808'),events);
% plotRecNonPowWave(2000,0,0,1);

recEvents = filterStruct(events,'recalled==1');
nonEvents = filterStruct(events,'recalled==0');
encodeEvents = [recEvents nonEvents];

lineNoise = [59 61]; %logspacing 2-128


    Rs = GetRateAndFormat(events(1));

    nyquist_freq = floor(Rs/2);

    freqs = (2.^((8:60)/8)); %new value from eeganalparams.txt
    valid_freqs_ind = find(freqs<nyquist_freq);

    valid_freqs = freqs(valid_freqs_ind);
    freqs = valid_freqs;


	bufferMS = 500;
    
%     elecCounter =  0;
%     allSig = [];allRec = [];allNon = [];

%     for l = electvector
	
%         elecStr = num2str(l);
%     elecCounter = elecCounter +1;

	%kurtSkew is a variable to decide whether or not to correct for kurtosis or skewed events

			if kurtSkew == 0

			foo=1;


			elseif kurtSkew == 1
			

			d = gete_ms(1,encodeEvents,encodeDur,encodeOffset,bufferMS);
			k = kurtosis(d,[],2); kThresh = 4;
			encodeEventsTemp = encodeEvents(k<kThresh);


			drec = gete_ms(1,recEvents,encodeDur,encodeOffset,bufferMS);
			krec = kurtosis(drec,[],2);
			recEventsTemp = recEvents(krec<kThresh);
			

			dnon = gete_ms(1,nonEvents,encodeDur,encodeOffset,bufferMS);
			knon = kurtosis(dnon,[],2);
			nonEventsTemp = nonEvents(knon<kThresh);
			
			elseif kurtSkew == 2
			
			foo = 1;
			
			%need to add the skew stuff here

				
            end
            
            [~,recPow] = getphasepow(21,recEventsTemp,encodeDur,encodeOffset,bufferMS,'filtfreq',[59 61],'resampledrate',250,'freqs',freqs); %the is the part with the problem
            [~,nonPow] = getphasepow(21,nonEventsTemp,encodeDur,encodeOffset,bufferMS,'filtfreq',lineNoise,'resampledrate',250,'freqs',freqs);
            [~,encodePow] = getphasepow(1,encodeEventsTemp,encodeDur,encodeOffset,bufferMS,'filtfreq',lineNoise,'resampledrate',250,'freqs',freqs);
             recPow = smoothRows(recPow);
             nonPow = smoothRows(nonPow);
            
             if normVar==1
            
            normPow = squeeze(mean(encodePow,3));
            normPowSTD = std(normPow,0,1);
            normPowMean = squeeze(mean(normPow,1));
            encodePow = (encodePow - repmat(normPowMean,[size(encodePow,1) 1 size(encodePow,3)]))./repmat(normPowSTD,[size(encodePow,1) 1 size(encodePow,3)]);

            recPow = (recPow - repmat(normPowMean,[size(recPow,1) 1 size(recPow,3)]))./repmat(normPowSTD,[size(recPow,1) 1 size(recPow,3)]);
            nonPow = (nonPow - repmat(normPowMean,[size(nonPow,1) 1 size(nonPow,3)]))./repmat(normPowSTD,[size(nonPow,1) 1 size(nonPow,3)]);
            
        else
            recPow = log10(recPow);
			nonPow = log10(nonPow);
			encodePow = log10(encodePow);
             end
   
             recPow= recPow(:,:,15:(end-15));
             nonPow = nonPow(:,:,15:(end-15));
             
             set(gcf,'name','Left Precuneus (Lead 21)','numbertitle','off')
             
             
             subplot(2,2,1)
             contourf(1:size(recPow,3),freqs,squeeze(mean(recPow,1)),40,'edgecolor','none');
             set(gca,'yscale','log','ytick',[2 4 8 16 32 64 128]);
             xxx =caxis;  
             title('Recalled');
             xlabel('time(ms)');
             ylabel('frequency(Hz)');
             colorbar;
             
             subplot(2,2,2)
             contourf(1:size(nonPow,3),freqs,squeeze(mean(nonPow,1)),40,'edgecolor','none');
             set(gca,'yscale','log','ytick',[2 4 8 16 32 64 128]);
             caxis([xxx]);  
             title('Nonrecalled');
             xlabel('time(ms)');
             ylabel('frequency(Hz)');
             colorbar;
             
%             printPath = fullfile(printDir,['recPowWave' subjname '_' elecStr '.eps']);
%             print(gcf,'-depsc2',printPath);
             
            
         
   