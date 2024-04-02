function [eog,perif] = getcapinfo(captype)
%GETCAPINFO   Get electrode information for a type of scalp EEG cap.
%   [EOG,PERIF] = GETCAPINFO(CAPTYPE) gets information about the electrodes
%   of CAPTYPE caps. EOG is a cell array, where EOG{1} gives the
%   right electrooculogram channels, and EOG{2} gives the left. These
%   channels are used for blink detection.
%
%   PERIF gives a list of peripheral electrodes that should be excluded
%   when rereferencing EEG data, as they tend to have poor or unusual signal.

if ~exist('captype','var')
  % for backwards compatibility, default is older cap
	captype = 'GSN200';
end

switch captype
	case 'GSN200'
	% EGI Geodesic Sensor Net, assumed to be 129 channels
	eog = {[26 127], [8 126]};
	% peripheral electrodes tend to flip up during the session, and get poor
	% signal; always exclude them from rereferencing
	perif = [127 126 17 128 125 120 44 49 56 63 69 74 82 89 95 100 108 114];
	
	case 'HCGSN'
	% EGI Hydrocel Geodesic Sensor Net, assumed to be 129 channels
	eog = {[25 127], [8 126]};
	% peripheral electrodes get better signal compared to GSN200, so don't
	% exclude them; just exclude the lower EOG channels since they have 
	% opposite polarity during eye-movement artifacts
	perif = [126 127]; 
end
