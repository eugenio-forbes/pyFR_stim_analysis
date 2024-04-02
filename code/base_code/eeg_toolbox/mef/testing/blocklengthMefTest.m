
clear;

origFilepath = '/data/eeg/TJ008/eeg.reref/TJ008_28Dec09_1053.080';
% origFilepath = '/data/eeg/TJ023/eeg.reref/TJ023_23Nov10_1447.113';
% origFilepath = '/data/eeg/TJ027/eeg.reref/TJ027_24Mar11_1031.093';
origParamDir = fileparts(origFilepath);

mefFilepath = 'completeTest.mef';

subjectEncodePass = 'a';
sessionEncodePass = 'b';
subjectDecodePass = 'a';
sessionDecodePass = 'b';

%% process original
[samplerate,nBytes,dataformat,gain] = GetRateAndFormat(origParamDir);

eegfile = fopen(origFilepath,'r','l'); % NOTE: the 'l' means that it came from a PC!

if eegfile==-1
    error('ERROR: EEG File did not open\n');
end

% read the eeg data
fseek(eegfile,0,-1);
origEEG2 = fread(eegfile,inf,'int16')';

% close the file
fclose(eegfile);

origEEG4 = int32(origEEG2);

%% iteration test
blockLengths = unique(round(logspace(log10(1),log10(20),10)));
nB = length(blockLengths);

tRats = zeros(nB,1);
sRats = zeros(nB,1);

Niter = 10;
N = 300;

for iB = 1:nB
    tic 
    
    %%% mef
    % save as mef
    raw2mef(origEEG4,samplerate,mefFilepath,subjectEncodePass,sessionEncodePass,blockLengths(iB));

    % reload mef
    mefHeader = read_mef_header(mefFilepath,sessionDecodePass);
    mefEEG = decomp_mef(mefFilepath,0,mefHeader.number_of_samples,sessionDecodePass)';

    % compare
    if( sum(mefEEG ~= origEEG2) )
        error('Decoded mef not equal to original');
    end

    %%% timing

    nSamples = length(origEEG4);

    tOrig = zeros(Niter,1);
    tMef = zeros(Niter,1);

    iDuration = 100;
    iStart = floor((nSamples - iDuration)*rand(Niter,N));

    for i = 1:Niter
        tic
        for j = 1:N
            eegfile = fopen(origFilepath,'r','l'); % NOTE: the 'l' means that it came from a PC!

            % read the eeg data
            fseek(eegfile,iStart(i,j),-1);
            origEEGpart = fread(eegfile,iDuration,'int16')';

            % close the file
            fclose(eegfile);
        end
        tOrig(i) = toc;

        tic
        for j = 1:N
            mefEEGpart = decomp_mef(mefFilepath,iStart(i,j),iStart(i,j)+iDuration,sessionDecodePass)';
        end
        tMef(i) = toc;
    end

    tOrigAvg = mean(tOrig);
    tMefAvg = mean(tMef);

    tRats(iB) = tOrigAvg / tMefAvg;

    sOrig = dir(origFilepath);
    sMef = dir(mefFilepath);
    sRats(iB) = sOrig.bytes / sMef.bytes;
    
    fprintf('Done block length %d s\n', blockLengths(iB));
end

figure(1); clf;
subplot(1,2,1);
plot(blockLengths,tRats);
xlabel('block lengths [s]'); ylabel('time ratio (orig/mef)');

subplot(1,2,2);
plot(blockLengths,sRats);
xlabel('block lengths [s]'); ylabel('size ratio (orig/mef)');
