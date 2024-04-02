function extractPYFR_STIM_wrapper(subject,expDir,session,forceSession,isStim)
%
% FUNCTION:
%   extractPYFR_STIMevents(SUBJECT,EXPDIR,SESSION,FORCESESSION)
% 
% DESCRIPTION:
%   wrapper functions that makes and saves free recall stim and
%   math events
%
% INPUTS:
%   SUBJECT.........'TJ055'
%   EXPDIR..........path to 'session_['SESSION']' directory.  Examples:
%                   '/data/eeg/TJ055/behavioral/pyFR_stim2/TJ055_Maddox'
%   SESSION.........0 = looks in 'session_0' in EXPDIR
%   FORCESESSION....1 = [optional] sets session to 1 (despite the
%                       fact that behavioral data are in session_0)
%                       Leave blank or empty if session number is
%                       same as SESSION
%
% OUTPUTS:
%   Makes and save two events:
%     (1) events.mat......... contains 'events'
%     (2) MATH_events.mat.... contains 'events' and 'MATHcfg
%
% NOTES:
%   (1) written by jfburke on 2013-02-12 (john.fred.burke@gmail.com)
%



% check tp see if subject names match
fprintf('\n')
if isempty(regexp(expDir,subject))
  fprintf('  WARNNG: %s not found in %s\n',upper(subject),upper(expDir))
  fprintf('          you might be making an error.\n')
  fprintf('          please check this before making events. EXITING\n\n')
  fprintf('               !!! NO EVENTS SAVED !!! \n\n')
  return
end

% set defaults
if ~exist('forceSession','var')
  forceSession = [];
end

% get the directories
thisSessDirNAME = sprintf('session_%d',session);
thisSessDir     = fullfile(expDir,thisSessDirNAME);
 evFile         = fullfile(thisSessDir,'events.mat');
mevFile         = fullfile(thisSessDir,'MATH_events.mat');

% print opening line
fprintf('  Making FREE RECALL STIM and MATH EVENTS for ')
fprintf('%s, session %d: \n',subject,session)

%--------------------------------------------
fprintf('    %-15.15s','FREE RECALL: ')
if ~exist(evFile,'file')
  if ~isStim
    events=extractPYFRevents_noSTIM(subject,expDir,session,forceSession);
  else
    events=extractPYFRevents_STIM(subject,expDir,session,forceSession);
  end
  if isempty(events)
    return
  end
  save(evFile,'events');
  clear events
  fprintf('DONE.\n')
else
  fprintf('SKIPPING (events exist).\n')
end

%--------------------------------------------
fprintf('    %-15.15s','MATH:')
if ~exist(mevFile,'file')
  [events MATHcfg]=extractPYFR_MATHevents_STIM(subject,expDir,session,forceSession);
  save(mevFile,'events','MATHcfg');
  clear events
  fprintf('DONE.\n')  
else
  fprintf('SKIPPING (events exist).\n')  
end

fprintf('\n\n')