function ev_OUT = convertEventsToEegNoreref(ev,subj)
%
% FUNCTION:
%   convertEventsToEegNoreref
%
% DESCRIPTION:
%   This kjsdhfjkhsgd
%
% INPUT:
%   l;ksjdflksjdf
%
% OUTPUT:
%   lksdjflkjshd
%
% NOTES:
%   (1) written by jfburke [''] (john.fred.burke@gmail.com)
%
%

thisSubjHospID = subj(1:2);

switch thisSubjHospID
 case {'TJ','UP'}
  thisEEGnorerefCorrection = 'eeg.noreref';
 case {'CH'}
  error('RAFI please fill this in')
 case {'BW'}
  error('RAFI please fill this in')
 otherwise
  error('did not recognize subject hospital ID')
end

ev_OUT = ev;
for k=1:length(ev)  
  thisEEGfile_orig      = ev(k).eegfile;  
  thisEEGfile_corrected = regexprep(thisEEGfile_orig,'eeg.reref',thisEEGnorerefCorrection);  
  ev_OUT(k).eegfile     = thisEEGfile_corrected;
end