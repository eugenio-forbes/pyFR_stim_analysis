function [ev_out]=get_baseline_random(ev_in,secBetweenEv,secRandJitter,rSeed,recType)
%
% FUNCTION:
%   [ev_out]=getECoG_baseline_randTime(ev_in,secBetweenEv,secRandJitter,rSeed)
%
% DESCRIPTION:
%   gets randomly placed events in the eegfiles
%
% INPUTS:
%   ev_in............ events
%   secBetweenEv..... seconds between events (60)  
%   secRandJitter.... jitter around events times (10)
%   randSeed......... ranodm seed [optional]
%   recType..........'eeg' or other ('lfp' etc.) (optional; def = 'eeg')
%
% OUTPUTS:
%   ev_out........... baseline events
%
% NOTES:
%   (1) written by jfburke 10/2011 (john.fred.burke@gmail.com)
%   (2) AGR (ashwinramayya@gmail.com 04/2013).
%           -now works when session labels have more than one character
%           -additional input added: recording field (default 'eegfile'),
%           but can enter 'lfpfile' etc. if needed.

basEv = [];
counter=0;
%set the recType
if ~exist('recType','var')||isempty(recType)
    recType = 'eeg';
end

% set the random seed
if ~exist('rSeed','var')||isempty(rSeed)
  s = RandStream.create('mt19937ar','seed',sum(100*clock));
else
  s = RandStream.create('mt19937ar','seed',rSeed);
end
RandStream.setDefaultStream(s);

if isa(ev_in(1).session,'numeric')
    str2numFlag = 1;
    for i = 1:length(ev_in)
      ev_in(i).session = num2str(ev_in(i).session);
    end
else
    str2numFlag = 0;
end

uniSess    = unique({ev_in.session});
numSess    = length(uniSess);
thisSubj=unique({ev_in.subject});
useThisSubj=thisSubj{1};

for s=1:numSess  
  thisSess=uniSess(s);  
  if length({ev_in.session})~=length(ev_in);
    error('very bad');
  end  
  
  evThisSess = ev_in(strcmp({ev_in.session},thisSess));        
  eegFilesThisSess=unique({evThisSess.([recType 'file'])})';
  nFiles = length(eegFilesThisSess);
  
  % loop through events
  for f=1:nFiles
    thisEEGfile = eegFilesThisSess{f};
    if isempty(thisEEGfile);continue;end
  
    % get ino for this session
    [sampleFreq,nBytes,dataformat,gain] = GetRateAndFormat(thisEEGfile);
    sampleFreq = round(sampleFreq);
    if isempty(sampleFreq);error('NO PARAMS FILE');end
    
    % get the data type
    for k=1:1000
      filestem=sprintf('%s.%.3i',thisEEGfile, k);
      if ~exist(filestem,'file'); %try a different format
          filestem = sprintf('%s.%.i',thisEEGfile, k);
      end
      if ~exist(filestem,'file'); continue; end
      fileFeatures=dir(filestem);
      switch dataformat
       case {'int16','short'}
	numSamples=floor(fileFeatures.bytes/2);
       otherwise
	error('Are all dataformats int16 or short?')
      end
      break
    end
    
    intervalSam = floor(secBetweenEv*sampleFreq);    
    randomSam   = floor(secRandJitter*sampleFreq);    
    startSam    = intervalSam*2;
    endSam      = numSamples-intervalSam*2;  
    
    basSamplesThisFile_raw      = startSam:intervalSam:endSam;    
    nBasSamplesThisFile         = length(basSamplesThisFile_raw);  
    basSamplesThisFile_jittered = basSamplesThisFile_raw + ...	
	round((rand(1,nBasSamplesThisFile)-.5)*randomSam);    
    
    % see if you had any overflow    
    basSamplesThisFile_jittered(basSamplesThisFile_jittered<startSam)=startSam;    
    basSamplesThisFile_jittered(basSamplesThisFile_jittered>endSam)=endSam;    
    
    for k=1:nBasSamplesThisFile
      counter=counter+1;
      basEv(counter).subject=useThisSubj;
      basEv(counter).session=thisSess;
      basEv(counter).([recType 'file'])=thisEEGfile;
      basEv(counter).([recType 'offset'])=basSamplesThisFile_jittered(k);
    end    
    
  end
end
ev_out = basEv;

if str2numFlag == 1
   for i = 1:length(ev_out)
      ev_out(i).session = str2double(ev_out(i).session); 
   end
end
