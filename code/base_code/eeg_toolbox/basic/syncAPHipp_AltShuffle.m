

% function syncAPHippRedoAvg(events,encodeDur,encodeOff,kurtSkew)

encodeDur = 1800;
encodeOff = 0;
kurtSkew = 1;
subjnum = {'CC002' 'CC007' 'CC013' 'CC016' 'CC020' 'CC024' 'UT001' 'UT008' 'UT010' 'UT012'};
Aelec = {[73:75] [33:35] [61:64] [157:161] [97:100] [31:33] [11,23:24] [1:3] [1:3] [1:3]};
Pelec = {[53:56] [13:15] [28:30] [121:123] [83:84] [21:23] [19,21,22,25] [23:24] [21:24] [21:24]};
subj.subjnum = subjnum;
subj.Aelec = Aelec;
subj.Pelec = Pelec;



for ind = 1:length(subj.subjnum)
    
    tic
    
    %Make directory for saving pics and zMat
    subjname = subj.subjnum{ind};
    printPath = sprintf('/Users/brad/Documents/JL/HippSync/%s',subjname);
    if ~exist(printPath);
        mkdir(printPath);
    else
        foo=1;
    end
    
    
    
    %Make the lock file directory
    lockfileDir = sprintf('/Users/brad/Documents/JL/HippSync/%s/rankLockFiles/',subjname);
    if ~exist(lockfileDir,'dir')
        mkdir(lockfileDir)
    end
    
    this_subj_lock = strcat(subjname,'_LOCK');
    thisSubjLockPath = strcat(lockfileDir,this_subj_lock);
    this_subj_done = strcat(subjname,'_DONE');
    ThisSubjDonePath =  strcat(lockfileDir,this_subj_done);
    this_subj_error = strcat(subjname,'_ERROR');
    ThisSubjErrorPath =  strcat(lockfileDir,this_subj_error);
    
    
    
    
    
    
    if ~exist(thisSubjLockPath,'file')
        %Save the lockfile
        fid = fopen(thisSubjLockPath,'w');fclose('all');
        try
            
            Eventpath = sprintf('/Users/brad/ExperimentData/subjFiles/%s/behavioral',subjname);
            behContent = dir(Eventpath);
            filenames = {behContent.name};
            
            
            findEvent = strfind(filenames,'FR1');
            notMatch = cellfun(@isempty, findEvent);
            Matchname = filenames(~notMatch);
            if isempty(Matchname) == 0
                Eventpath = fullfile(Eventpath, 'FR1/events.mat');
            else
                findEvent = strfind(filenames,'pyFR');
                notMatch = cellfun(@isempty, findEvent);
                Matchname = filenames(~notMatch);
                if isempty(Matchname) == 0
                    Eventpath = fullfile(Eventpath, 'pyFR/events.mat');
                else
                    findEvent = strfind(filenames,'CatFR1');
                    notMatch = cellfun(@isempty, findEvent);
                    Matchname = filenames(~notMatch);
                    if isempty(Matchname) == 0
                        Eventpath = fullfile(Eventpath, 'CatFR1/events.mat');
                    else
                        error('No non-stim FR task exist for subject')
                    end
                end
            end
            
            %     for task = ['catFR' 'pyFR' 'FR1']
            %         findEvent = strfind(filenames,task);
            %         notMatch = cellfun(@isempty, findEvent);
            %         Matchname = filenames(~notMatch);
            %         if isempty(Matchname) == 0
            %             Eventpath = fullfile(Eventpath,task);
            %         else end
            %     end
            
            
            
            
            Eventpath = char(Eventpath);
            load(Eventpath);
            
            
            
            
            
            
            
            
            
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
            
            for A = subj.Aelec{ind}%elecsHipp
                
                
                
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
                    
                    
                    d = gete_ms(A,encodeEvents,encodeDur,encodeOffset,bufferMS);
                    k = kurtosis(d,[],2); kThresh = 4;
                    encodeEventsTemp = encodeEvents(k<kThresh);
                    
                    
                    drec = gete_ms(A,recEvents,encodeDur,encodeOffset,bufferMS);
                    krec = kurtosis(drec,[],2);
                    recEventsTemp = recEvents(krec<kThresh);
                    
                    
                    dnon = gete_ms(A,nonEvents,encodeDur,encodeOffset,bufferMS);
                    knon = kurtosis(dnon,[],2);
                    nonEventsTemp = nonEvents(knon<kThresh);
                    
                    numRec = length(recEventsTemp);
                    randInd = randperm(length(nonEventsTemp));
                    randInd = randInd(1:numRec);
                    nonEventsTemp = nonEventsTemp(randInd);
                    
                    %             dretr = gete_ms(A,retrEvents,1400,-1400,bufferMS);
                    %             kretr = kurtosis(dretr,[],2);
                    %             retrEventsTemp = retrEvents(kretr<kThresh);
                    
                    
                elseif kurtSkew == 2
                    
                    foo = 1;
                    
                    %need to add the skew stuff here
                    
                    
                end
                
                %[encodePhHipp,~] = getphasepow(A,encodeEventsTemp,encodeDur,encodeOffset,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);
                [recAHipp,~] = getphasepow(A,recEventsTemp,encodeDur,encodeOffset,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);
                [nonAHipp,~] = getphasepow(A,nonEventsTemp,encodeDur,encodeOffset,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);
                %[retrPhHipp,~] = getphasepow(A,retrEventsTemp,1400,-1400,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);
                
                postCounter = 0;
                
                for P = subj.Pelec{ind}%elecsCingTemp
                    
                    postCounter = postCounter+1;
                    pairCounter = pairCounter+1;
                    
                    %[encodePhCing,~] = getphasepow(P,encodeEventsTemp,encodeDur,encodeOffset,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);
                    [recPHipp,~] = getphasepow(P,recEventsTemp,encodeDur,encodeOffset,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);
                    [nonPHipp,~] = getphasepow(P,nonEventsTemp,encodeDur,encodeOffset,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);
                    %[retrPhCing,~] = getphasepow(P,retrEventsTemp,1400,-1400,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);
                    
                    
                    for ff = 1:size(recAHipp,2)
                        
                        for tt = 1:size(recAHipp,3)
                            
                              [reczVal recpVal]=ranksum_shuff(recAHipp(:,ff,tt),recPHipp(:,ff,tt),250);
%                             ph1 =squeeze(recAHipp(:,ff,tt));
%                             ph2 = squeeze(recPHipp(:,ff,tt));
%                             
%                             nShuffles=250;
%                             
%                             [p_1,rBar]=rayleigh(ph1-ph2);
%                             
%                             
%                             
%                             %to do the bootstrap from Lachaux et al. 1999
%                             s=size(ph1);
%                             fakeP=repmat(nan,[nShuffles s(2)]);
%                             for n=1:nShuffles
%                                 
%                                 inds=randperm(s(1));
%                                 d=ph1-ph2(inds,:);
%                                 fakeP(n,:)=rayleigh(d);
                                
%                             end
                              zMatrixRec(antCounter,postCounter,ff,tt) = reczVal; 
%                             p=mean(repmat(p_1,[nShuffles 1])>fakeP,1);
                            
%                             p_matrixRec(antCounter,postCounter,ff,tt) = p;
%                             rBarMatrixRec(antCounter,postCounter,ff,tt) = rBar;
                            
                        end
                    end
                    
                    
                    
                    
                    for ff = 1:size(nonAHipp,2)
                        
                        for tt = 1:size(nonAHipp,3)
                            
                              [nonzVal nonpVal]=ranksum_shuff(nonAHipp(:,ff,tt),nonPHipp(:,ff,tt),250);
%                             ph1 =squeeze(nonAHipp(:,ff,tt));
%                             ph2 = squeeze(nonPHipp(:,ff,tt));
%                             
%                             nShuffles=250;
%                             
%                             [p_1,rBar]=rayleigh(ph1-ph2);
%                             
%                             
%                             
%                             %to do the bootstrap from Lachaux et al. 1999
%                             s=size(ph1);
%                             fakeP=repmat(nan,[nShuffles s(2)]);
%                             for n=1:nShuffles
%                                 
%                                 inds=randperm(s(1));
%                                 d=ph1-ph2(inds,:);
%                                 fakeP(n,:)=rayleigh(d);
                                
%                         end
                             zMatrixNon(antCounter,postCounter,ff,tt) = nonzVal;
%                             p=mean(repmat(p_1,[nShuffles 1])>fakeP,1);
%                             
%                             p_matrixNon(antCounter,postCounter,ff,tt) = p;
%                             rBarMatrixNon(antCounter,postCounter,ff,tt) = rBar;
                            
                            
                        end
                    end
                    
                    
%                     rBarMatrixNon = log10(rBarMatrixNon);
%                     %   rBarMatrixRetr = squeeze(mean(rBarMatrixRetr,4));
%                     
%                     zMatrixNon = p_matrixNon;
%                     zMatrixNon(zMatrixNon<.0001) = .0001;
%                     zMatrixNon(zMatrixNon>.9999) = .9999;
%                     zMatrixNon = -norminv(zMatrixNon);
%                     
%                     %           zMatrixRetr = squeeze(mean(zMatrixRetr,4));
%                     
%                     rBarMatrixRec = log10(rBarMatrixRec);
%                     % rBarMatrix = squeeze(mean(rBarMatrix,4));
%                     
%                     zMatrixRec = p_matrixRec;
%                     zMatrixRec(zMatrixRec<.0001) = .0001;
%                     zMatrixRec(zMatrixRec>.9999) = .9999;
%                     zMatrixRec = -norminv(zMatrixRec);
                    
                    timeVect = (encodeOffset+1):2:(encodeOffset+encodeDur);
                    
                    subplot(1,2,1)
                    xxx = [];
                    xxx = squeeze(zMatrixNon(antCounter,postCounter,:,:));
                    xxx = xxx(1:45,:);
                    contourf(timeVect,freqs,xxx,40,'edgecolor','none')
                    set(gca,'yscale','log','ytick',[2.^(1:8)]);
                    title_string = sprintf('AHip ch%d PHip ch%d NonRec', A,P);
                    title(title_string);
                    caxis([-2 4])
                    colorbar
                    
                    
                    subplot(1,2,2)
                    yyy = [];
                    yyy = squeeze(zMatrixRec(antCounter,postCounter,:,:));
                    yyy = yyy(1:45,:);
                    contourf(timeVect,freqs,yyy,40,'edgecolor','none')
                    set(gca,'yscale','log','ytick',[2.^(1:8)]);
                    title_string = sprintf('AHip ch%d PHip ch%d Rec', A,P);
                    title(title_string);
                    caxis([-2 4])
                    colorbar
                    
                    strps = sprintf('APHippSynch Ch%d Ch%d',A,P);
                    printFile = fullfile(printPath, strps);
                    print(gcf,'-dpng','-append',printFile);
                    
                    %            strps = sprintf('APHippSynch ElectrodeAvg');
                    %            printFile = fullfile(printPath, strps);
                    %            print(gcf,'-dpng','-append',printFile);
                    
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
            
            errorFlag = 0;
            
        catch e
            errorFlag = 1;
            fid = fopen(ThisSubjErrorPath,'w');fclose('all');
            save(ThisSubjErrorPath,'e');
        end
        
        if errorFlag == 0
            fid = fopen(ThisSubjDonePath,'w');fclose('all');
        end
    else continue
        
    end
end

