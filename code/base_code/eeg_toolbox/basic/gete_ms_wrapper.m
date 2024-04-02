function[eeg] = gete_ms_wrapper(elecNum,events,durationMS,offsetMS,bufferMS,filtfreq,filttype,filtorder,resampleFreq)
%This wrapper returns raw eeg traces for a set of events. 
%Inputs
%elecNum        %electrode number to load. Can be a single channel (eg,[1]) or a set
                %of channels (for bipolar, eg,[1 2]). see getBipolarSubjElecs in
                %eeg_toolbox
%events         events structure (with 'eegfile', and 'eegoffset' fields).
                %see ps2_makeEvents and alignTool
%optional:
%durationMS     %how long to clip
%offsetMS       %where to start relative to mstime
%bufferMS       %buffer to add to beginning and end of clips
%filtfreq       %filters out particular frequencies (e.g., [59.5 60.5]);
                %see buttfilts
%filtorder      %e.g, 4 (see buttfult)

%Note:
%resampleFreq (last arg of gete) is set to 1000 so that gete returns eeg in ms.
                
%Written by AGR 10-31-12

%Parse Inputs
%
if ~exist('durationMS','var') || isempty(durationMS)
    durationMS = 1500;
end
if ~exist('offsetMS','var') || isempty(offsetMS)
    offsetMS = -500;
end    
if ~exist('bufferMS','var') || isempty(bufferMS)
    bufferMS = 0;
end
if ~exist('filtfreq','var') || isempty(filtfreq)
    filtfreq = [];
end 
if ~exist('filttype','var') || isempty(filttype)
    filttype = [];
end
if ~exist('filtorder','var') || isempty(filtorder)
    filtorder = [];
end
if ~exist('resampleFreq','var') || isempty(resampleFreq)
    resampleFreq = [];
end

% Update durationMS and offsetMS with buffer
%(Check this: buffer input to gete_ms doesnt seem to work when not
%using a filter)

durationMS = durationMS + 2*bufferMS; % 01-09-18 
offsetMS = offsetMS-bufferMS;

%check whether it is a bipolar montage or not
if size(elecNum,2) == 1 % not bipolar, proceed normally
   [eeg] = gete_ms(elecNum, events,durationMS,offsetMS,bufferMS,...
            filtfreq,filttype,filtorder,resampleFreq);
    
elseif size(elecNum,2) == 2 % then it is bipolar
    %First, rename events.eegfile to point to ,noreref
    for k=1:length(events)
         events(k).eegfile=regexprep(events(k).eegfile,'eeg.reref','eeg.noreref');
    end
    
    %Second, gete from both channels and subtract them from each other 
    [eeg1] = gete_ms(elecNum(1), events,durationMS,offsetMS,0,...
        filtfreq,filttype,filtorder,resampleFreq);
    [eeg2] = gete_ms(elecNum(2), events,durationMS,offsetMS,0,...
        filtfreq,filttype,filtorder,resampleFreq);
    
    eeg = eeg1 - eeg2;
end