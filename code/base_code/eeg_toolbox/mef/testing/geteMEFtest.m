
% before running this script, change the mefFileroot to the directory where
% you want the MEF files to be placed. Ensure this directory contains an
% 'eeg.reref' folder (for the MEF files) and an 'eeg.noreref' folder (for
% the params.txt files)

%%% MEF file parameters
mefFileroot = '/home1/erichuns/mef/';
subjectEncodePassword = 'XJnfklsuy332';
sessionEncodePassword = '';
blockLength = 5;

%%% subject parameters and events
sub='TJ005';
lead=32;
duration=1000;offset=-200;buffer=200; resamp=1000;
times=offset:(1000/resamp):offset+duration-1;

[events,eventDetails]=get_sub_events('pymms',sub,0);

% targs = filterStruct(events,'strcmp(type,''CUE'')&istarget==1');
% lures=filterStruct(events,'strcmp(type,''CUE'')&istarget==0');
targs = events;
targsMEF = targs;

%%% create MEF files
channelSuffix = sprintf('.%03i',lead);

targFileList = {targs.eegfile};
targFileSet = unique(targFileList);
nF = length(targFileSet);

sizeRats = zeros(nF,1);

for i = 1:nF
    eegfname = [targFileSet{i}, channelSuffix];
    fMask = strcmp(targFileSet{i}, targFileList);
    fInd = find(fMask);
    fEvents = targs(fMask);
    
    [samplerate,nBytes,dataformat,gain] = GetRateAndFormat(fEvents(1));
    
    eegfile = fopen(eegfname,'r','l');

    if eegfile==-1
        error('ERROR: EEG File did not open\n');
    end

    % read the eeg data
    fseek(eegfile,0,-1);
    origEEG4 = int32( fread(eegfile,inf,'int16')' );

    % close the file
    fclose(eegfile);
    
    [~,mefFilename] = fileparts(targFileSet{i});
    mefFilepath = [mefFileroot 'eeg.reref/' mefFilename];
    mefFilepathExt = [mefFilepath channelSuffix];
    
    raw2mef(origEEG4,samplerate,mefFilepathExt,subjectEncodePassword,sessionEncodePassword,blockLength);
    
    mefFilepathParams = [mefFileroot 'eeg.noreref/' mefFilename '.params.txt'];
    paramsFile = fopen(mefFilepathParams,'w+');
    if paramsFile == -1
        error('ERROR opening params.txt file');
    end
    fwrite(paramsFile,sprintf('samplerate %f\n', samplerate));
    fwrite(paramsFile,sprintf('dataformat ''%s''\n', dataformat));
    fwrite(paramsFile,sprintf('gain %f\n', gain));
    fclose(paramsFile);
    
    for j = fInd
        targsMEF(j).eegfile = mefFilepath;
    end
    
    sOrig = dir(eegfname);
    sMef = dir(mefFilepathExt);
    sizeRats(i) = sOrig.bytes / sMef.bytes;
end
    
tic
targ_d = gete_ms(lead,targs,duration,offset,buffer,[59 ],'low',1,resamp,[offset 0]);
toc

global IS_MEF
IS_MEF = false;

tic
targ_d3 = gete_ms2(lead,targs,duration,offset,buffer,[59 ],'low',1,resamp,[offset 0]);
toc

IS_MEF = true;
tic
targ_dM = gete_ms2(lead,targsMEF,duration,offset,buffer,[59 ],'low',1,resamp,[offset 0]);
toc

% fprintf('gete_ms2 errors: %d\n', sum(targ_d3(:) ~= targ_d(:)) );
fprintf('gete_msMEF errors: %d\n', sum(targ_dM(:) ~= targ_d(:)) );
fprintf('average compression ratio: %0.2f\n', mean(sizeRats));
