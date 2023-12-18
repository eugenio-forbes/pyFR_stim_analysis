function [samplerate,nBytes,dataformat,gain] = GetRateAndFormat(event)
%GETRATEANDFORMAT - Get the samplerate, gain, and format of eeg data.
%
% function [samplerate,nBytes,dataformat,gain] = GetRateAndFormat(event)
%

if ischar(event)
    % event is actually a path to an EEG file
    path = event;
else
    path = event.eegfile;
end

% Look in "noreref" directory for session-specific parameters file
fs = filesep;
sepInds = strfind(path,fs);

if isempty(sepInds) % relative path (current folder)
    dir = '';
    parentDir = ['..' fs]; % parent dir is up one level
    name = path;
elseif length(sepInds) == 1 % relative path (child folder)
    dir = path(1:sepInds(end));
    parentDir = ''; % parent dir is current dir
    name = path(sepInds(end)+1:end);
else
    dir = path(1:sepInds(end));
    parentDir = path(1:sepInds(end-1));
    name = path(sepInds(end)+1:end);
end

nameDotInds = strfind(name,'.');
if ~isempty(nameDotInds) && nameDotInds(1) ~= 1
    baseName = name(1:nameDotInds(1)-1); % strip all extensions
else
    baseName = name;
end

% Modified 09/10/18 Sarah Seger
% Added line to account for bipolar eegs --> will check eeg.bipolar path 
if any(strfind(dir, 'eeg.bipolar'))
    sessionParamsPath = [parentDir 'eeg.bipolar' fs baseName '.params.txt'];
else
    sessionParamsPath = [parentDir 'eeg.noreref' fs baseName '.params.txt'];
end

%sessionParamsPath = '/Volumes/project/TIBIR/Lega_lab/shared/lega_ansir/subjFiles/UT068/eeg.noreref/UT068_03Oct17_1100.params.txt';
% EH - trying to open the file is a quicker way to check for its existence
% than using exist(sessionParamsPath,'file')
file = fopen(sessionParamsPath,'rt'); 
if( file == -1 )
    sessionParamsPath = [dir 'params.txt'];
    file = fopen(sessionParamsPath,'rt'); 
end

params = eegparams({'samplerate','gain','dataformat'},file);

samplerate = params{1};
if( isempty(samplerate) )
    % EH - can't do anything if the sample rate isn't present (no default)
    gain = [];
    dataformat = '';
    nBytes = [];
    return;
end

if( ~isempty(params{2}) )
    gain = params{2};
else
    gain = 1.0;
end

if( ~isempty(params{3}) )
    dataformat = params{3};
else
    dataformat = 'short';
end

switch dataformat
    case {'short','int16'}
        nBytes = 2;
    case {'single','int32'}
        nBytes = 4;
    case {'double','int64'}
        nBytes = 8;
    otherwise
        error('BAD DATA FORMAT!');
end
