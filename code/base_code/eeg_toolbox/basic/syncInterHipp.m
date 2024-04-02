

% function syncAPHippRedoAvg(events,encodeDur,encodeOff,kurtSkew)

encodeDur = 1800;
encodeOff = 0;
kurtSkew = 1;
subjnum = {'CC002' 'CC007' 'CC009' 'CC013' 'CC013' 'CC014' 'CC015' 'CC016' 'CC016' 'CC020' 'CC024' 'UT001' 'UT001' 'UT005' 'UT005' 'UT008' 'UT008' 'UT010' 'UT010' 'UT011' 'UT011' 'UT012' 'UT012' 'UT013' 'UT013'};
% Relec = {[61:64] [178:183] [1:2] [1:6] [101:104] [1:3] [71:74] [1:3] [121:122]};
% Lelec = {[153:154] [157:161] [11,23:24] [103:108] [1:3] [101:104] [1:4] [111:113] [1:2]};
Antelec = {[73:75] [33:35] [93:96] [61:64] [153:154] [79:82] [95:97] [157:161] [178:183] [97:100] [31:33] [1:2] [11,23:24] [1:6] [103:108] [101:104] [1:3] [1:3] [101:104] [71:74] [1:4] [1:3] [111:113] [121:122] [1:2]};
subj.subjnum = subjnum;
% subj.Relec = Relec;
% subj.Lelec = Lelec;
subj.Antelec = Antelec;



for ind = 1:length(subj.subjnum)
    
    tic
    
    %Make directory for saving pics and zMat
    subjname = subj.subjnum{ind};
    printPath = sprintf('/Users/brad/Documents/JL/InterHippSyncHil/%s',subjname);
    if ~exist(printPath);
        mkdir(printPath);
    else
        foo=1;
    end
    
    
    
    %Make the lock file directory
    lockfileDir = sprintf('/Users/brad/Documents/JL/InterHippSyncHil/%s/LockFiles/',subjname);
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
            bandPass = [2 4;4 9;9 14;14 30;30 70];
            
            antCounter = 0;
            postCounter = 0;
            pairCounter = 0;
            
            zMatrixRec = [];
            zMatrixNon =[];
            rBarMatrixRec = [];
            rBarMatrixNon =[];
            p_matrixRec = [];
            p_matrixNon = [];
          
           
            subjnck = nchoosek(subj.Antelec{ind},2);
            for p = 1:length(subjnck)
                
            for one = subjnck(p,1)%elecsHipp
                
                
                
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
                    
                    
                    d = gete_ms(one,encodeEvents,encodeDur,encodeOffset,bufferMS);
                    k = kurtosis(d,[],2); kThresh = 4;
                    encodeEventsTemp = encodeEvents(k<kThresh);
                    
                    
                    drec = gete_ms(one,recEvents,encodeDur,encodeOffset,bufferMS);
                    krec = kurtosis(drec,[],2);
                    recEventsTemp = recEvents(krec<kThresh);
                    
                    
                    dnon = gete_ms(one,nonEvents,encodeDur,encodeOffset,bufferMS);
                    knon = kurtosis(dnon,[],2);
                    nonEventsTemp = nonEvents(knon<kThresh);
                    
                    numRec = length(recEventsTemp);
                    randInd = randperm(length(nonEventsTemp));
                    randInd = randInd(1:numRec);
                    nonEventsTemp = nonEventsTemp(randInd);
                    
                    %             dretr = gete_ms(one,retrEvents,1400,-1400,bufferMS);
                    %             kretr = kurtosis(dretr,[],2);
                    %             retrEventsTemp = retrEvents(kretr<kThresh);
                    
                    
                elseif kurtSkew == 2
                    
                    foo = 1;
                    
                    %need to add the skew stuff here
                    
                    
                end
                
                %[encodePhHipp,~] = getphasepow(one,encodeEventsTemp,encodeDur,encodeOffset,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);
                [recAHipp,~] = getphasepow(one,recEventsTemp,encodeDur,encodeOffset,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);
                [nonAHipp,~] = getphasepow(one,nonEventsTemp,encodeDur,encodeOffset,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);
                %[retrPhHipp,~] = getphasepow(one,retrEventsTemp,1400,-1400,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);
                 %[recAHipp,~] = gethilbertphase2(one,recEventsTemp,encodeDur,encodeOffset,bufferMS,bandPass,60,500);
                 %[nonAHipp,~] = gethilbertphase2(one,recEventsTemp,encodeDur,encodeOffset,bufferMS,bandPass,60,500);

                
                postCounter = 0;
                
                for two = subjnck(p,2)%elecsCingTemp
                    
                    postCounter = postCounter+1;
                    pairCounter = pairCounter+1;
                    
                    %[encodePhCing,~] = getphasepow(two,encodeEventsTemp,encodeDur,encodeOffset,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);
                    [recPHipp,~] = getphasepow(two,recEventsTemp,encodeDur,encodeOffset,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);
                    [nonPHipp,~] = getphasepow(two,nonEventsTemp,encodeDur,encodeOffset,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);
                    %[retrPhCing,~] = getphasepow(two,retrEventsTemp,1400,-1400,bufferMS,'filtfreq',lineNoise,'resampledrate',500,'freqs',freqs);
                     %[recPHipp,~] = gethilbertphase2(two,recEventsTemp,encodeDur,encodeOffset,bufferMS,bandPass,60,500);
                     %[nonPHipp,~] = gethilbertphase2(two,recEventsTemp,encodeDur,encodeOffset,bufferMS,bandPass,60,500);

                    
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
                                
                            end
                            
                            p=mean(repmat(p_1,[nShuffles 1])>fakeP,1);
                            
                            p_matrixRec(pairCounter,ff,tt) = p;
                            rBarMatrixRec(pairCounter,ff,tt) = rBar;
                            
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
                                
                            end
                            
                            p=mean(repmat(p_1,[nShuffles 1])>fakeP,1);
                            
                            p_matrixNon(pairCounter,ff,tt) = p;
                            rBarMatrixNon(pairCounter,ff,tt) = rBar;
                            
                            
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
                    xxx = squeeze(zMatrixNon(pairCounter,:,:));
                    xxx = xxx(1:45,:);
                    contourf(timeVect,freqs,xxx,40,'edgecolor','none')
                    set(gca,'yscale','log','ytick',[2.^(1:8)]);
                    title_string = sprintf('1AHip ch%d 2AHip ch%d Non', one,two);
                    title(title_string);
                    colorbar
                    
                    
                    subplot(2,2,2)
                    yyy = [];
                    yyy = squeeze(zMatrixRec(pairCounter,:,:));
                    yyy = yyy(1:45,:);
                    contourf(timeVect,freqs,yyy,40,'edgecolor','none')
                    set(gca,'yscale','log','ytick',[2.^(1:8)]);
                    title_string = sprintf('1AHip ch%d 2AHip ch%d Rec', one,two);
                    title(title_string);
                    colorbar
                    
                    strps = sprintf('InterHippSynch Ch%d Ch%d',one,two);
                    printFile = fullfile(printPath, strps);
                    print(gcf,'-dpng','-append',printFile);
                    
                    %            strps = sprintf('APHippSynch ElectrodeAvg');
                    %            printFile = fullfile(printPath, strps);
                    %            print(gcf,'-dpng','-append',printFile);
                    
                    toc
                    beep
                end
            end
            
            end
            
            
            
            
            %             zMatrix = squeeze(mean(zMatrix,4));
            toc
            beep
            
            strmat = sprintf('InterHippSynch.mat');
            saveFile = fullfile(printPath, strmat);
            save(saveFile,'zMatrixRec','zMatrixNon','pairCounter');
            
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

