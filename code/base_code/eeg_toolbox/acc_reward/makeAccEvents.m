function events = makeAccEvents(behPath,subject,session)


%This file will create an events structure for the reward task.  The
%events = makeAccEvents(behPath,catPath,subject, session)
%
% behPath is the path to the actual log file 
% catVect is an array of the category image preferences in order from most
% desirable to least desirable, as exported by the session_0 logFile  


behlogFile = fullfile(behPath,'session.log');

%read in the session.log file into a series of arrays, one for each column
[c1,c2,c3,c4,c5,c6,c7,c8]=textread(behlogFile,'%n%n%s%s%s%s%s%s','delimiter','\t');

%get the number of log file entries
logFileLength = size(c1,1);

%make a vector that is the length of all entries in the logFile
idx_trial_cue = 1:logFileLength;

%this will create a vector that tells you each instance of a cue in the log
%file indexed by logFile entry number
idx_trial_cue = idx_trial_cue(strcmp(c3,'CUE'));
idx_trial_cue = idx_trial_cue';

%now create a vector that assigns a category valence to each trial, i.e.
%what level of reward was on offer for each trial  valence of -1 means that
%it was a neutral feedback trial; 2 = most, 1 = mid, 0 = least preferred

trial_idx_valence = c7(idx_trial_cue);

%convert these to doubles

for t = 1: length(trial_idx_valence)
    if strcmp(trial_idx_valence(t),'0')
        temp(t) = 0;
    else
        if strcmp(trial_idx_valence(t),'-1')
        temp(t) = -1;
        else
            if strcmp(trial_idx_valence(t),'1')
                temp(t) = 1;
            else
                if strcmp(trial_idx_valence(t),'2')
                    temp(t) = 2;
                end
            end
        end
    end
end 
        
trial_idx_valence = temp;       

%now create a vector the same length as the logFIle that assigns a trial
%number to each logFile entry

numCueEvents = size(idx_trial_cue,1);
listTrialNums = 1:numCueEvents;
trialNum = zeros(logFileLength,1);
valenceNum = zeros(logFileLength,1);

trialNum(trialNum==0) = 999;
valenceNum(valenceNum==0) = 999;

for n = 1:(numCueEvents-1)

    trialNum(idx_trial_cue(n):(idx_trial_cue(n+1)-1)) = listTrialNums(n);
    valenceNum(idx_trial_cue(n):(idx_trial_cue(n+1)-1)) = trial_idx_valence(n);
end

trialNum(idx_trial_cue(end):(end-2)) = listTrialNums(end);
valenceNum(idx_trial_cue(end):(end-2)) = trial_idx_valence(end);

%trialNum(end-1:end) = 999;
%valenceNum(end-1:end) = 999;


%change the block_cue events to 999 for trial and valence

trialNum(strcmp(c3,'BLOCK_CUE'))= 999;
valenceNum(strcmp(c3,'BLOCK_CUE')) = 999;

%create a vector the length of the logFile that has the feedback type for
%each feedback event.  The fb events themselves are already coded by the
%feedback type but the correct/incorrect/neutral needs to be coded as 1 for
%positive feedback, 0 for negative feedback, and -1 for neutral feedback
%all non-feedback events are coded as 999

fbVect = zeros(logFileLength,1);
fbVect(fbVect==0) = 999;
fbVect(strcmp(c4,'You got the shot!'))=1;
fbVect(strcmp(c4,'Too fast.'))=0;
fbVect(strcmp(c4,'Too slow.'))=0;
fbVect(strcmp(c4,'Camera reloaded.'))=-1;

%now create a vector that tells you which side the stimulus was located on
%during the game play

sideVect = cell(logFileLength,1);
[sideVect{strcmp(c4,'right')}] = deal('right');
[sideVect{strcmp(c4,'left')}] = deal('left');
bothIdx = strcmp(c4,'right') | strcmp(c4,'left');
[sideVect{~bothIdx}] = deal('999');


event = struct('subject',[], 'session', [], 'type',{},...
    'trial', [], 'valence', [], 'feedback',[],...
    'side',{});

events = repmat(event,1,logFileLength);


for e=1:logFileLength
    
    events(e).subject = subject;
    events(e).session = session;
    events(e).type = c3{e};
    events(e).trial = trialNum(e);
    events(e).valence = valenceNum(e);
    events(e).feedback = fbVect(e);
    events(e).side = sideVect{e};
    events(e).mstime = c1(e);
    events(e).msoffset = 0;
end


%pare down the events, such that we get rid of all the move events except
%for the one in the middle

events = events((strcmp(c3,'MOVE') & strcmp(c4,'2')) | ~strcmp(c3,'MOVE'));

eventRoot = fullfile(behPath,'/events.mat');
saveEvents(events,eventRoot);


