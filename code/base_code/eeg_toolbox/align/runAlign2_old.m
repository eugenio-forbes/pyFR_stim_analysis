function runAlign2(beh_file,eeg_file,chan_files,log_files,ms_field,isfrei,isegi)
%RUNALIGN2 - Wrapper for the pulsealign and logalign functions.
%
% This is JJ's modified version of runAlign.  I changed it to call a
%different set of functions that I've begun working on to do a
%better job of aligning.  Critically, these functions will use
%*all* of the pulses that are manually selected by the alignTool,
%rather than only a subset, which is what the current tools do (12/8/2009).
%
% This is a relatively specialized wrapper for the EEG toolbox
% alignment functions.  You provide a list of behavioral sync pulse
% files, a list of eeg sync pulse files, a list of matching eeg
% channel files, and a list of log files or files containg events
% structures to process.
%
% This function will modify the files passed into the log_files
% parameters, saving original copies.
%
% FUNCTION:
%   runAlign(beh_file,eeg_file,chan_files,log_files,isfrei)
%
% INPUT ARGS:
%   beh_file = {'behfile1','behfile2'};
%   eeg_file = {'eegfile1','eegfile2'};
%   chan_files = {'/data/eeg/file1.001','/data/eeg/file2.001'};
%   log_files = {'events.mat'};
%   ms_field = 'mstime';
%   isfrei = 1;       % Read in freiburg format
%
%
%
%

if length(eeg_file)>1,
  error('Josh has only designed this code to work with a single EEG file; if you have multiple files, then you should use the standard runAlign.m function. sorry...');
end
         
if ~exist('ms_field','var')
  ms_field = 'mstime';
end

if ~exist('isfrei','var')
  isfrei = 1;
end

if ~exist('isegi','var')
  isegi = 0;
end

threshMS = 10;
window = 100;

% load in the beh_ms and the pulses
beh_ms = [];
for f = 1:length(beh_file) 
  % read in free recall data
  beh_ms = [beh_ms ; textread(beh_file{f},'%n%*[^\n]','delimiter','\t')];
end

% sort it
beh_ms = sort(beh_ms);

% loop over pulses and run pulsealign to get sets of matched inputs
beh_in = cell(1,length(eeg_file));
eeg_in = cell(1,length(eeg_file));
for f = 1:length(eeg_file)
    
  % load eeg pulses
  if isfrei
    [s,pulse_str] = system(['grep -i SYNC1 ' eeg_file{f} ' | cut -d, -f 1']);
    pulses = strread(pulse_str,'%u');
  elseif isegi
    % open DIN file 
    eegsyncID = fopen(eeg_file{f});
    eegsync = fread(eegsyncID, inf, 'int8');
    fclose(eegsyncID);
    pulses = find(eegsync>0);
  else
    [s,pulse_str] = system(['cut -f 1 ' eeg_file{f}]);
    pulses = strread(pulse_str,'%u');
  end

  %JJ: I removed this code below, because I cannot think of a reason that it
  %needs to be here.  this also allowed me to remove the samplerate
  %parameter, which I think is delightful
  
% $$$   pulse_ms = pulses*1000/samplerate;
% $$$ 
% $$$   % remove all pulses under 100ms (Part of Start and End pulses)
% $$$   dp = diff(pulse_ms);
% $$$   yp = find(dp < 100);
% $$$   pulse_ms(yp+1) = [];
% $$$   pulses(yp+1) = [];

% $$$   % run pulsealign
% $$$   mywin = min(round(length(pulses)/2),window);
% $$$   if mywin < 5
% $$$     mywin = 5;
% $$$   end
  
%  [beh_in{f},eeg_in{f}] = pulsealign(beh_ms,pulses,samplerate,threshMS,mywin);

beh_in{f}=beh_ms;
eeg_in{f}=pulses;

end

% run logalign with the beh_ms and eeg_offset info
superaligner(beh_in{1},eeg_in{1},chan_files{1},log_files{1},ms_field);


function superaligner(beh_ms,eeg_offset,eeg_files,log_files,ms_field)
%this is a modified version of logalign.  the big difference is that it does *not* require that beh_ms and eeg_offset
%are matched.  the main requirement is that eeg_offset is a subset of beh_ms.

%%%%%%%%%%%%%%%%%
%first lets make a *rough* estimate of the sampling rate
%this will be the starting point for the optimization, bitches.
experimentDurationSec=range(beh_ms)/1000; %in seconds
[path,efile,ext] = fileparts(eeg_files);   %get length of recording
chan = str2num(ext(2:end));
event = struct('eegfile',  fullfile(path,efile));
eeg = gete(chan,event,0);
durationSamples = length(eeg{1});
approxSampRate=durationSamples/experimentDurationSec;
%%%%%%%%%%%%%%%%%

eeg_offset_s=eeg_offset./approxSampRate;
beh_s=beh_ms./1000;

windSize=30;

eegBlockInds=1:windSize:length(eeg_offset_s);

%go in blocks of windsize, looking for significant matches in the inter-pulse intervals
beh_d_s=diff(beh_s);
for b=1:length(eegBlockInds)-1
  eeg_d_s=diff(eeg_offset_s(eegBlockInds(b):eegBlockInds(b+1)));
  clear r;
  for i=1:(length(beh_d_s)-length(eeg_d_s))  
    r(i)=corr(eeg_d_s,beh_d_s(i+[0:windSize-1]));
  end
  
  [blockR(b),beginWindowOffset(b)]=max(r);

end


keyboard



recordingStart_s=beh_s(1);  %starting offset for regression
beh_s=beh_s-recordingStart_s;

%srs=approxSampRate+[-10:.05:1];
srs=512+[-1:.0005:1];


%cost=@(x) (findMaxDev(beh_s,eeg_offset./x(1)+x(2)));
%cost=@(x) (a=eeg_offset./x(1);findMaxDev(beh_s,a-a(1)+beh_s(x(2))));

%keyboard

%eeg_offset=eeg_offset(10:end);

maxDev=inf+zeros(length(srs),length(beh_s));
for s=1:length(srs)
  curSR=srs(s);
  disp([curSR min(maxDev(:))]);
  eeg_offset_s=eeg_offset./curSR;
  eeg_offset_s2=eeg_offset_s-eeg_offset_s(1);  %now these start at 0
%    for t=77:83
for t=1:length(beh_s)      
      eeg_offset_s3=eeg_offset_s2+beh_s(t);
      [maxDev(s,t),error]=findMaxDev(beh_s,eeg_offset_s3);
    end
end

%cost=@(x) (findMaxDev(beh_s,eeg_offset./x(1)+x(2)));

keyboard


function [maxDev_s,error]=findMaxDev(actualPulse_s,predPulse_s)
error=zeros(size(predPulse_s));
%yes, I know this could be vectorized. but i wanted to save memory in case of lots of pulses
for i=1:length(predPulse_s) 
%  [error(i),ind(i)]=min(abs(predPulse_s(i)-actualPulse_s));
  [error(i)]=min(abs(predPulse_s(i)-actualPulse_s));
end;

%closest=actualPulse_s(ind);
%maxDev_s=max(error);
maxDev_s=median(error);
%maxDev_s=sum(error);



