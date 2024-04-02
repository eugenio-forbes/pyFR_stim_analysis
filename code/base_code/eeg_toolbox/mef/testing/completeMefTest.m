
clear;

% origFilepath = '/data/eeg/TJ003/eeg.reref/TJ003_30Mar09_1121.044';
% origFilepath = '/data/eeg/TJ005_1/eeg.noreref/TJ005_1_18Aug09_1008.091';
origFilepath = '/data/eeg/TJ008/eeg.reref/TJ008_28Dec09_1053.080';
% origFilepath = '/data/eeg/TJ023/eeg.reref/TJ023_23Nov10_1447.113'; % this is strange (no high frequency)
% origFilepath = '/data/eeg/TJ027/eeg.reref/TJ027_24Mar11_1031.093'; % this is strange (no high frequency)
origParamDir = fileparts(origFilepath);

mefFilepath = 'completeTest.mef';

subjectEncodePassword = 'a';
sessionEncodePassword = 'b';
subjectDecodePassword = 'a';
sessionDecodePassword = 'b';

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

%% filtering [optional]
FILTER = false;
if( ~FILTER )
    % no filtering
    origEEG4filt = origEEG4;
else
    nyq = samplerate/2;
    
    % plot original spectrum
    eegFFT = fft(double(origEEG4));
    freqs = linspace(0,samplerate,length(eegFFT));
    
    figure(1); clf;
    plot( freqs, abs(eegFFT) );
    
    origEEG4filt = double(origEEG4);
    
    % low-pass filter to remove high frequency noise
    origEEG4filt = buttfilt( origEEG4filt, min(350,nyq-50), samplerate, 'low' );

    % notch filter at multiples of 60 Hz
    notchFreqs = bsxfun(@plus, 60 * (1:floor(nyq/60))', 3*[-1 1]);
    origEEG4filt = buttfilt( origEEG4filt, notchFreqs, samplerate, 'stop', 1 );
    
    origEEG4filt = int32( origEEG4filt );
    
    % plot spectrum
    eegFFT = fft(double(origEEG4filt));
    freqs = linspace(0,samplerate,length(eegFFT));
    
    figure(2); clf;
    plot( freqs, abs(eegFFT) );
    
    % compare filtered with unfiltered over window
    wind = 5000:10000;
    T = wind/samplerate;
    
    figure(3); clf;
    plot( T, origEEG4(wind), 'b', T, origEEG4filt(wind), 'r' );
end

%% mef
% save as mef
blockLength = 5;
raw2mef(origEEG4filt,samplerate,mefFilepath,subjectEncodePassword,sessionEncodePassword,blockLength);

% reload mef
mefHeader = read_mef_header(mefFilepath,sessionDecodePassword);
mefEEG = decomp_mef(mefFilepath,0,mefHeader.number_of_samples,sessionDecodePassword)';

% compare
if( sum(mefEEG ~= origEEG4filt) )
    error('Decoded mef not equal to original');
end

%% timing
nSamples = length(origEEG4);

Niter = 20;
N = 1000;
% Niter = 1;
% N = 20;

tOrig = zeros(Niter,1);
tMef = zeros(Niter,1);

iDuration = 1000;
% iStart = floor((nSamples - iDuration)*rand(Niter,N));
iStart = floor((nSamples)*rand(Niter,N));
% iStart = 3*ones(Niter,N);

KEEP_OPEN_ORIG = true;
KEEP_OPEN_MEF = true;

for i = 1:Niter
    origEEGparts = zeros(N,iDuration,'int16');
    mefEEGparts = zeros(N,iDuration,'int16');

    if( ~KEEP_OPEN_ORIG )
        tic
        for j = 1:N
            eegfile = fopen(origFilepath,'r','l'); % NOTE: the 'l' means that it came from a PC!

            % read the eeg data
            fseek(eegfile,2*iStart(i,j),-1);
            readEEG = fread(eegfile,iDuration,'int16')';
            origEEGparts(j,1:length(readEEG)) = readEEG;

            % close the file
            fclose(eegfile);
        end
        tOrig(i) = toc;
    else
        tic
        eegfile = fopen(origFilepath,'r','l');
        for j = 1:N
            fseek(eegfile,2*iStart(i,j),-1);
            readEEG = fread(eegfile,iDuration,'int16')';
            origEEGparts(j,1:length(readEEG)) = readEEG;
        end
        fclose(eegfile);
        tOrig(i) = toc;
    end
        
    if( ~KEEP_OPEN_MEF )    
        tic
        for j = 1:N
            mefEEGparts(j,:) = decomp_mef(mefFilepath,iStart(i,j)+1,iStart(i,j)+iDuration,sessionDecodePassword)';
        end
        tMef(i) = toc;
    else
        tic
        [mefEEGparts,mefEEGpartLengths] = decomp_mef_events(mefFilepath,iStart(i,:),iDuration,sessionDecodePassword);
        mefEEGparts = int16(mefEEGparts');
        tMef(i) = toc;
    end
    
    if( sum(sum(origEEGparts ~= mefEEGparts )) )
        error('Decoded MEF not equal to original')
    end
    if( sum( mefEEGpartLengths ~= iDuration ) )
        warning('Some parts not full duration')
%         keyboard;
    end
end

% figure(4); clf;
% j = 1;
% plot(1:iDuration,origEEGparts(j,:),1:iDuration,mefEEGparts(j,:));

figure(5); clf; hold on;
bar( [mean(tOrig) mean(tMef)] );
errorbar( [mean(tOrig) mean(tMef)], 1.96/sqrt(Niter)*[std(tOrig) std(tMef)], 'r.' );

tOrigAvg = mean(tOrig);
tMefAvg = mean(tMef);

tRat = tOrigAvg / tMefAvg

sOrig = dir(origFilepath);
sMef = dir(mefFilepath);
sRat = sOrig.bytes / sMef.bytes

