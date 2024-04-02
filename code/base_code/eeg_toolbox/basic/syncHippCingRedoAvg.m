

function syncAPHippM(events,encodeDur,encodeOff,kurtSkew)

  tic  
%     subjCounter = 0;
%     p_matrixRec = [];
%     p_matrixNon = [];

    
%     for nn = 1:length(events)
%debugging

        
        
            
    
        
%         subjCounter = subjCounter+1;
        
        
        
          %events = hippEvents(nn);
      
  subjname = events(1).subject;
  
%   elecStructIdx = find(strcmp(events(1).subject,{hippStruct.subject}));
%   
%   elecStructIdx = elecStructIdx(1);
%   
%   elecsHipp = hippStruct(elecStructIdx).electrodes;
        
       
%   subjNumer = strcmp(subjname,'CC002');
   
%   if subjNumer ~=1
 %      continue
  % else
 %      foo=1;
 %  end
        
        printPath = sprintf('/Users/brad/Documents/HippSync/%s/psync250',subjname);
        if ~exist(printPath);
            mkdir(printPath);
        else
            foo=1;
        end
        
        
        
        
        
%         tempStruct = cingStruct(strcmp({cingStruct.subject},subjname));
%         
%         elecsCing = [tempStruct.electrodes];
%    
        
	events = filterStruct(events,'~strcmp(eegfile,'''')');

	baseEvents = filterStruct(events,'strcmp(type,''ORIENT'')');

	recEvents = filterStruct(events,'recalled==1');
	nonEvents = filterStruct(events,'recalled==0');
    encodeEvents = [recEvents nonEvents];
    
	recallstr = sprintf('strcmp(type,''REC_WORD'')');
	x_voc  = filterStruct(events,recallstr);
	%this gives you all of the vocalizations but now you have to pull out only the correct 		%vocalizations, no intrusions, etc.

%     retrEvents =getRetrievalEvents_NEW(events,[1500 500],[1500 500]);
%     
%     retrEvents = retrEvents.rec;

lineNoise = [59 61];


Rs = GetRateAndFormat(events(1));

nyquist_freq = floor(Rs/2);

freqs = eeganalparams('freqs');

valid_freqs_ind = find(freqs<nyquist_freq);

valid_freqs = freqs(valid_freqs_ind);
freqs = valid_freqs;
freqs = freqs(1:45);

	bufferMS = 500;
	encodeOffset = encodeOff;
        
    antCounter = 0;
    postCounter = 0;
    pairCounter = 0;
    
    zMatrixRec = [];
    zMatrixNon =[];
    rBarMatrixRec = [];
    rBarMatrixNon =[];
    p_matrixRec = [];
    p_matrixNon = [];
        
        for l = [2]%elecsHipp
            
            

            antCounter = antCounter+1;

% tempSide = hippStruct(elecStructIdx).side;
% sideValue = tempSide(hippCounter);
% 
% elecsCingTemp = elecsCing([tempStruct.side]==sideValue);


            bufferMS = 500;
            encodeOffset = encodeOff;
	%time window for encoding is set as an input value

	%loop through all the relevant electrodes for this subject

	
	%kurtSkew is a variable to decide whether or not to correct for kurtosis or skewed events

    
   
    
			if kurtSkew == 0

			foo=1;


			elseif kurtSkew == 1
			

			d = gete_ms(l,encodeEvents,encodeDur,encodeOffset,bufferMS);
			k = kurtosis(d,[],2); kThresh = 4;
			encodeEventsTemp = encodeEvents(k<kThresh);


			drec = gete_ms(l,recEvents,encodeDur,encodeOffset,bufferMS);
			krec = kurtosis(drec,[],2);
			recEventsTemp = recEvents(krec<kThresh);
			

			dnon = gete_ms(l,nonEvents,encodeDur,encodeOffset,bufferMS);
			knon = kurtosis(dnon,[],2);
			nonEventsTemp = nonEvents(knon<kThresh);
            
%             numRec = length(recEventsTemp);
%             randInd = randperm(length(nonEventsTemp));
%             randInd = randInd(1:numRec);
%             nonEventsTemp = nonEventsTemp(randInd);
            
%             dretr = gete_ms(l,retrEvents,1400,-1400,bufferMS);
%             kretr = kurtosis(dretr,[],2);
%             retrEventsTemp = retrEvents(kretr<kThresh);
			
			
			elseif kurtSkew == 2
			
			foo = 1;
			
			%need to add the skew stuff here

				
            end
			
            %[encodePhHipp,~] = getphasepow(l,encodeEventsTemp,encodeDur,encodeOffset,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);
            [recAHipp,~] = getphasepow(l,recEventsTemp,encodeDur,encodeOffset,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);
            [nonAHipp,~] = getphasepow(l,nonEventsTemp,encodeDur,encodeOffset,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);
            %[retrPhHipp,~] = getphasepow(l,retrEventsTemp,1400,-1400,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);

            postCounter = 0;

                    for c = [5]%elecsCingTemp
    
                    postCounter = postCounter+1;
                    pairCounter = pairCounter+1;
    
                    %[encodePhCing,~] = getphasepow(c,encodeEventsTemp,encodeDur,encodeOffset,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);
                    [recPHipp,~] = getphasepow(c,recEventsTemp,encodeDur,encodeOffset,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);
                    [nonPHipp,~] = getphasepow(c,nonEventsTemp,encodeDur,encodeOffset,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);
                    %[retrPhCing,~] = getphasepow(c,retrEventsTemp,1400,-1400,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);


                        for ff = 1:size(recAHipp,2)
    
                                for tt = 1:size(recAHipp,3)
        
                                        ph1 =squeeze(recAHipp(:,ff,tt));
                                        ph2 = squeeze(recPHipp(:,ff,tt));

                                        nShuffles=250;

                                        [p_1,rBar]=rayleigh(ph1-ph2);



%to do the bootstrap from Lachaux et al. 1999
                                        s=size(ph1); 
                                        fakeP=repmat(nan,[nShuffles s(2)]);
                                                for n=1:nShuffles
    
                                                inds=randperm(s(1));
                                                d=ph1-ph2(inds,:);
                                                fakeP(n,:)=rayleigh(d);
                                                [foo1 tempR]=rayleigh(d);
                                                fakeP(n,:) = tempR;
                                                end
  
                                        p=mean(repmat(rBar,[nShuffles 1])<fakeP,1);
                
                                        p_matrixRec(antCounter,postCounter,ff,tt) = p;
                                        rBarMatrixRec(antCounter,postCounter,ff,tt) = rBar;

                                end
                        end
        



                        for ff = 1:size(nonAHipp,2)
    
                                for tt = 1:size(nonAHipp,3)
        
                                        ph1 =squeeze(nonAHipp(:,ff,tt));
                                        ph2 = squeeze(nonPHipp(:,ff,tt));

                                        nShuffles=250;

                                        [p_1,rBar]=rayleigh(ph1-ph2);



%to do the bootstrap from Lachaux et al. 1999
                                        s=size(ph1); 
                                        fakeP=repmat(nan,[nShuffles s(2)]);
                                                for n=1:nShuffles
    
                                                    inds=randperm(s(1));
                                                    d=ph1-ph2(inds,:);
                                                    fakeP(n,:)=rayleigh(d);
                                                    [foo1 tempR]=rayleigh(d);
                                                    fakeP(n,:) = tempR;
                                                end
  
                                        p=mean(repmat(rBar,[nShuffles 1])<fakeP,1);
                
                                        p_matrixNon(antCounter,postCounter,ff,tt) = p;
                                        rBarMatrixNon(antCounter,postCounter,ff,tt) = rBar;


                                end
                        end  
            rBarMatrixNon = log10(rBarMatrixNon);
            %   rBarMatrixRetr = squeeze(mean(rBarMatrixRetr,4));
            
            zMatrixNon = p_matrixNon;
            zMatrixNon(zMatrixNon<.0001) = .0001;
            zMatrixNon(zMatrixNon>.9999) = .9999;
            zMatrixNon = -norminv(zMatrixNon);
                
%           zMatrixRetr = squeeze(mean(zMatrixRetr,4));
            
            rBarMatrixRec = log10(rBarMatrixRec);
            % rBarMatrix = squeeze(mean(rBarMatrix,4));
            
            zMatrixRec = p_matrixRec;
            zMatrixRec(zMatrixRec<.0001) = .0001;
            zMatrixRec(zMatrixRec>.9999) = .9999;
            zMatrixRec = -norminv(zMatrixRec);
            
            timeVect = (encodeOffset+1):2:(encodeOffset+encodeDur);
        
           subplot(2,2,1)
           xxx = [];
           xxx = squeeze(zMatrixNon(antCounter,postCounter,:,:));
           xxx = xxx(1:45,:);
           contourf(timeVect,freqs,xxx,40,'edgecolor','none')
           set(gca,'yscale','log','ytick',[2.^(1:8)]);
           title_string = sprintf('AHip ch%d PHip ch%d NonRec', l,c);
           title(title_string);
           caxis([-2 4])
           colorbar
            
            
           subplot(2,2,2)
           yyy = [];
           yyy = squeeze(zMatrixRec(antCounter,postCounter,:,:));
           yyy = yyy(1:45,:);
           contourf(timeVect,freqs,yyy,40,'edgecolor','none')
           set(gca,'yscale','log','ytick',[2.^(1:8)]);
           title_string = sprintf('AHip ch%d PHip ch%d Rec', l,c);
           title(title_string);
           caxis([-2 4])
           colorbar
           
           strps = sprintf('APHippSynch Ch%d Ch%d',l,c);
           printFile = fullfile(printPath, strps);
           print(gcf,'-dpng','-append',printFile);
        
            toc 
            beep
                end
        end
   
            
            
            
            
%             zMatrix = squeeze(mean(zMatrix,4));      
        toc
        beep
    
        strmat = sprintf('APHippSynch.mat');
        saveFile = fullfile(printPath, strmat);
        save(saveFile,'zMatrixRec','zMatrixNon','rBarMatrixRec','rBarMatrixNon','pairCounter');

