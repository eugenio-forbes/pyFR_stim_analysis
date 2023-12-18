function n09_check_events(varargin)
if isempty(varargin)
    %%% Directory information
    root_directory = '/directory/with/analysis_folder';
    analysis_folder_name = 'pyFR_stim_analysis';
else
    root_directory = varargin{1};
    analysis_folder_name = varargin{2};
end

%%% List directories
subjects_directory = fullfile(root_directory,'subject_files');
analysis_directory = fullfile(root_directory,analysis_folder_name);
list_directory = fullfile(analysis_directory,'lists');
data_directory = fullfile(analysis_directory,'data');
statistics_directory = fullfile(analysis_directory,'statistics');
exclusion_directory = fullfile(analysis_directory,'exclusion_lists');

%%% Load subject list and excluded subjects so far
load(fullfile(list_directory,'subject_list.mat'),'subject_list');
load(fullfile(exclusion_directory,'excluded_subjects.mat'),'excluded_subjects');

%%% Load session list and excluded sessions so far
load(fullfile(list_directory,'session_list.mat'),'session_list');
load(fullfile(exclusion_directory,'excluded_sessions.mat'),'excluded_sessions');

%%% Load electrode list and excluded electrodes so far
load(fullfile(list_directory,'electrode_list.mat'),'electrode_list');
load(fullfile(exclusion_directory,'excluded_electrodes.mat'),'excluded_electrodes');

%%% Loop through sessions to load events files, change EEG file to current
%%% noreref folder path, and correct any errors in the files. Initialize
%%% table columns for events file information to be able to identify events
%%% files that might need to be realigned, and to exclude subjects with
%%% incomplete sessions.

n_sessions = height(session_list);
bad_session = false(n_sessions,1);
word_duration = NaN(n_sessions,1);
retrieval_duration = NaN(n_sessions,1);
encoding_duration = NaN(n_sessions,1);
math_duration = NaN(n_sessions,1);
lists = NaN(n_sessions,1);
subject = cell(n_sessions,1);
task = cell(n_sessions,1);
session = cell(n_sessions,1);
subject_ID = NaN(n_sessions,1);
session_ID = NaN(n_sessions,1);
n_conditions = NaN(n_sessions,1);

for idx = 1:n_sessions
    %%% Get session information
    this_subject = session_list.subject{idx};
    subject{idx} = this_subject;
    this_task = session_list.task{idx};
    task{idx} = this_task;
    this_session = session_list.session{idx};
    session{idx} = this_session;
    subject_ID(idx) = session_list.subject_ID(idx);
    session_ID(idx) = session_list.session_ID(idx);
    n_conditions(idx) = session_list.n_conditions(idx);
    
    %%% Path to original events file in subject's directory
    events_file_path = fullfile(subjects_directory,this_subject,'behavioral',this_task,this_session,'events.mat');

    %%% Load events and convert to table
    load(events_file_path,'events')
    events = struct2table(events);
    
    %%% 'B' type event marks the beginning of the session, however some 
    %%% sessions could have been restarted a few times and there could be
    %%% more than one row with 'B' and no events in between. So deleting
    %%% rows up to the last 'B' that is list number 1 so that no actual
    %%% events are deleted.
    if sum(strcmp(events.type,'B'))>1
        last_B = find(strcmp(events.type,'B') & events.list <= 1,1,'last');
        events(1:last_B-1,:) = [];
    end
    
    %%% Serial position should be only one value and output as a numerical
    %%% array from table, but a few subjects had events with errors that
    %%% led to there being more than 1 number in serial position and as a
    %%% result the column is cell. Correcting with true serial position

    if iscell(events.serialpos) && any(cellfun(@length,events.serialpos)>1)
        bad_serialpos = find(cellfun(@length,events.serialpos)>1);
        for jdx = 1:length(bad_serialpos)
            this_bad_serialpos_idx = bad_serialpos(jdx);
            bad_serialpos_list = events.list(this_bad_serialpos_idx);
            bad_serialpos_type = events.type(this_bad_serialpos_idx);
            if strcmp(bad_serialpos_type,'WORD')
                event_indices = find(events.list == bad_serialpos_list & strcmp(events.type,bad_serialpos_type));
                true_serialpos = find(event_indices == this_bad_serialpos_idx);
                events.serialpos(this_bad_serialpos_idx) = {int32(true_serialpos)};
                clear true_serialpos event_indices
            end
            clear this_bad_serialpos_idx bad_serialpos_list bad_serialpos_type
        end
        events.serialpos = cellfun(@(x) int32(x),events.serialpos);
    end
    
    %%% For some events, a few words appeared more than once throughout
    %%% encoding lists. Extra-list intrusions are marked with -1 in
    %%% intrusion column. Prior list intrusions have a value corresponding
    %%% to the difference of the list number in which the intrusion
    %%% happened and the list in which the word was encoded. Since some
    %%% words could have appeared more than once in encoding, some prior
    %%% list intrusion could have more than one value, thus resulting in
    %%% column of type cell. Some values could represent number of lists in
    %%% which the word reappeared after the intrusion was made. Correcting
    %%% with list number difference of closest list prior to when intrusion
    %%% was made.
    
    if iscell(events.intrusion) && any(cellfun(@length,events.intrusion)>1)
        bad_intrusions = find(cellfun(@length,events.intrusion)>1);
        for jdx = 1:length(bad_intrusions)
            this_bad_intrusion_idx = bad_intrusions(jdx);
            bad_intrusion_item = events.item(this_bad_intrusion_idx);
            bad_intrusion_list = events.list(this_bad_intrusion_idx);
            if any(strcmp(events.item,bad_intrusion_item) & strcmp(events.type,'WORD'))
                word_idx = strcmp(events.item,bad_intrusion_item) & strcmp(events.type,'WORD');
                word_list = events.list(word_idx);
                word_list = word_list(1);
                if word_list < bad_intrusion_list
                    events.intrusion{this_bad_intrusion_idx} = bad_intrusion_list - word_list;
                else
                    events.intrusion{this_bad_intrusion_idx} = -1;
                end
                clear word_idx word_list
            else
                events.intrusion{this_bad_intrusion_idx} = -1;
            end
            clear this_bad_intrusion_idx bad_intrusion_item bad_intrusion_list
        end
        events.intrusion = cellfun(@(x) int32(x),events.intrusion);
    end
    
    %%%Looping through lists to make sure that words were correctly marked
    %%%as recalled and assigning recall time to those that were incorrectly
    %%%marked as recall. While looping through lists gathering information
    %%%about event durations.
    
    unique_lists = unique(events.list);
    unique_lists = unique_lists(unique_lists >= 0);
    n_lists = length(unique_lists);
    lists(idx) = n_lists;
    temp_word_duration = NaN(1,n_lists);
    temp_retrieval_duration = NaN(1,n_lists);
    temp_encoding_duration = NaN(1,n_lists);
    temp_math_duration = NaN(1,n_lists);
    for jdx = 1:length(unique_lists)
        this_list = unique_lists(jdx);
        retrieved_indices = strcmp(events.type,'REC_WORD') & events.list == this_list;
        retrieved_words = events.item(retrieved_indices);
        if any(ismember(events.item,retrieved_words) & strcmp(events.type,'WORD') & events.list == this_list & events.recalled == 0)
            bad_word_idx = find(ismember(events.item,retrieved_words) & strcmp(events.type,'WORD') & events.list == this_list & events.recalled == 0);
            for kdx = 1:length(bad_word_idx)
                this_bad_word_idx = bad_word_idx(kdx);
                bad_word = events.item{this_bad_word_idx};
                retrieved_word_idx = strcmpi(events.item,bad_word) & strcmp(events.type,'REC_WORD') & events.list == this_list;
                recall_time = events.rectime(retrieved_word_idx);
                events.recall(this_bad_word_idx) = 1;
                events.rectime(this_bad_word_idx) = recall_time;
                clear this_bad_word_idx bad_word retrieved_word_idx recall_time
            end
        end
        orient_time = events.eegoffset(strcmp(events.type,'ORIENT') & events.list == this_list);
        word_times = events.eegoffset(strcmp(events.type,'WORD')  & events.list == this_list);
        retrieval_start = events.eegoffset(strcmp(events.type,'REC_START') & events.list == this_list);
        retrieval_stop = events.eegoffset(ismember(events.type,{'REC_STOP','E'}) & events.list == this_list);
        
        if length(word_times) < 10 || isempty(retrieval_start)
            lists(idx) = lists(idx)-1;
            events(events.list == this_list,:) = [];
            continue
        else
            retrieval_stop = retrieval_stop(1);
            temp_word_duration(jdx) = mean(diff(word_times),'omitnan');
            temp_retrieval_duration(jdx) = retrieval_stop-retrieval_start;
            temp_encoding_duration(jdx) = word_times(end) - orient_time;
            temp_math_duration(jdx) = retrieval_start - word_times(end);
        end
    end
    word_duration(idx) = mean(temp_word_duration,'omitnan')/1000;
    retrieval_duration(idx) = mean(temp_retrieval_duration,'omitnan')/1000;
    encoding_duration(idx) = mean(temp_encoding_duration,'omitnan')/1000;
    math_duration(idx) = mean(temp_math_duration,'omitnan')/1000;
    
    %%% Spanish session item numbers correspond to the items in English.
    %%% However there are could have been a few more words in Spanish that
    %%% would have a number greater than the total number of words in the
    %%% English word pool. So that latent semantic analysis can be
    %%% performed for all subject, replacing these item numbers with a
    %%% random number between 1 and 308.
    
    if any(events.itemno > 308)
        bad_itemno_idx = find(events.itemno > 308);
        taken_care_of = false(length(bad_itemno_idx),1);
        for jdx = 1:length(bad_itemno_idx)
            this_bad_itemno_idx = bad_itemno_idx(jdx);
            bad_itemno_type = events.type{this_bad_itemno_idx};
            bad_itemno_word = events.item{this_bad_itemno_idx};
            bad_itemno_list = events.list(this_bad_itemno_idx);
            new_itemno = randi([1,308],1,1);
            if strcmp(bad_itemno_type,'WORD')
                events.itemno(this_bad_itemno_idx) = new_itemno;
                was_recalled = find(strcmpi(events.item,bad_itemno_word) & strcmp(events.type,'REC_WORD')  & events.list == bad_itemno_list);
                if ~isempty(was_recalled)
                    events.itemno(was_recalled) = new_itemno;
                    taken_care_of(bad_itemno_idx == was_recalled) = true;
                end
            else
                if ~taken_care_of(bad_itemno_idx == this_bad_itemno_idx)
                    events.itemno(this_bad_itemno_idx) = new_itemno;
                end
            end
        end
    end
    
    %%% All possible changes needed to be made to eegfile paths listed
    %%% below.
    
    if any(cellfun(@isempty,events.eegfile))
        empty_cells = cellfun(@isempty,events.eegfile);
        events.eegfile(empty_cells) = repelem({''},sum(empty_cells),1);
    end
    if any(contains(events.eegfile,'/Volumes'))
        events.eegfile = strrep(events.eegfile,'/Volumes','');
    end
    if any(contains(events.eegfile,'/some/bad/pattern'))
        events.eegfile = strrep(events.eegfile,'/some/bad/pattern','/fixed/pattern');
    end
    
    %%% Saving a copy of events to analysis directory with the EEG file
    %%% path switched from reref (common average reference) to noreref so
    %%% that own referencing method can be used.
    
    if any(contains(events.eegfile,{'eeg.reref','eeg.bipolar'}))
        events.eegfile = strrep(events.eegfile,'eeg.bipolar','eeg.noreref');
        events.eegfile = strrep(events.eegfile,'eeg.reref','eeg.noreref');
    end
    
    %%% For excluding incomplete sessions
    
    if lists(idx) >= n_conditions(idx)*10
        save_path = fullfile(data_directory,this_subject,this_task,this_session);
        if ~isfolder(save_path)
            mkdir(save_path);
        end
        save(fullfile(save_path,'events.mat'),'events');
    else
        bad_session(idx) = true;
    end
end

%%% Add excluded sessions to exclusion list. Identify possibly subjects and
%%% electrodes that could have been excluded as a result of this exclusion.
%%% Reduce subject, session, and electrode list. 
new_excluded_sessions = session_list(bad_session,{'subject','task','session','subject_ID','session_ID'});
reason_for_exclusion = repelem({'session incomplete'},height(new_excluded_sessions),1);
new_excluded_sessions.reason_for_exclusion = reason_for_exclusion;
excluded_sessions = [excluded_sessions;new_excluded_sessions];
session_list = session_list(~ismember(session_list.session_ID,excluded_sessions.session_ID),:);

new_excluded_subjects = subject_list(~ismember(subject_list.subject_ID,session_list.subject_ID),{'subject','subject_ID'});
reason_for_exclusion = repelem({'session incomplete'},height(new_excluded_subjects),1);
new_excluded_subjects.reason_for_exclusion = reason_for_exclusion;
excluded_subjects = [excluded_subjects;new_excluded_subjects];
subject_list = subject_list(~ismember(subject_list.subject_ID,excluded_subjects.subject_ID),:);

new_excluded_electrodes = electrode_list(ismember(electrode_list.session_ID,new_excluded_sessions.session_ID),{'subject','task','session','channel_number','subject_ID','session_ID','electrode_ID'});
reason_for_exclusion = repelem({'session incomplete'},height(new_excluded_electrodes),1);
new_excluded_electrodes.reason_for_exclusion = reason_for_exclusion;
electrode_list = electrode_list(~ismember(electrode_list.session_ID,excluded_sessions.session_ID),:);
excluded_electrodes = [excluded_electrodes;new_excluded_electrodes]; 

%%% Save events info
events_info = table(subject,task,session,subject_ID,session_ID,word_duration,encoding_duration,math_duration,retrieval_duration,lists,n_conditions);
events_info(bad_session,:) = [];
save(fullfile(statistics_directory,'events_info.mat'),'events_info');

%%%Save all lists
save(fullfile(list_directory,'subject_list.mat'),'subject_list');
save(fullfile(list_directory,'session_list.mat'),'session_list');
save(fullfile(list_directory,'electrode_list.mat'),'electrode_list');
save(fullfile(exclusion_directory,'excluded_subjects.mat'),'excluded_subjects');
save(fullfile(exclusion_directory,'excluded_sessions.mat'),'excluded_sessions');
save(fullfile(exclusion_directory,'excluded_electrodes.mat'),'excluded_electrodes');

%%% Loop through subjects that had stimulation free session to gather
%%% copies of events for coherence and stimulation site analyses. Gather
%%% information of event durations to identify misaligned events files.
%%% Being more liberal about these events because only measuring overall coherence
%%% of hippocampal contacts and stimulation side during encoding and retrieval events
%%% Change path of EEG file from reref to noreref to use own referencing method
%%% and exclude events that were not within recordings.

subject_list = subject_list(subject_list.has_stim_free_session,:);

n_subjects = height(subject_list);
word_duration = NaN(n_subjects,1);
retrieval_duration = NaN(n_subjects,1);
encoding_duration = NaN(n_subjects,1);
math_duration = NaN(n_subjects,1);
lists = NaN(n_subjects,1);
subject = cell(n_subjects,1);
subject_ID = NaN(n_subjects,1);

for idx = 1:n_subjects
    %%% Get session information
    this_subject = subject_list.subject{idx};
    subject{idx} = this_subject;
    subject_ID(idx) = subject_list.subject_ID(idx);
    
    %%% Path to original events file in subject's directory
    events_file_path = fullfile(subjects_directory,this_subject,'behavioral/FR1/session_0/events.mat');
    
    %%% Getting single most complete session. Some subject had incomplete sessions.
    switch this_subject
        case {'subject_with_better_session_1'}
            events_file_path = strrep(events_file_path,'session_0','session_1');
        case {'subject_with_better_session_2'}
            events_file_path = strrep(events_file_path,'session_0','session_2');
        case {'subject_with_different_folder'}
            events_file_path = strrep(events_file_path,'bad_folder','good_folder');
    end
    
    %%% Load events and convert to table
    load(events_file_path,'events')
    events = struct2table(events);
    
        %%% All possible changes needed to be made to eegfile paths listed
    %%% below.    
    if any(cellfun(@isempty,events.eegfile))
        empty_cells = cellfun(@isempty,events.eegfile);
        events.eegfile(empty_cells) = repelem({''},sum(empty_cells),1);
    end
    if any(contains(events.eegfile,'/Volumes'))
        events.eegfile = strrep(events.eegfile,'/Volumes','');
    end
    if any(contains(events.eegfile,'/some/bad/pattern'))
        events.eegfile = strrep(events.eegfile,'/some/bad/pattern','/fixed/pattern');
    end
    if any(contains(events.eegfile,'.h5'))
        events.eegfile = strrep(events.eegfile,'.h5','');
    end
    
    no_offset = strcmp(events.eegfile,'') | events.eegoffset == 0;
    events(no_offset,:) = [];
    
    %%%Looping through lists to gather information about event durations.
    unique_lists = unique(events.list);
    unique_lists = unique_lists(unique_lists >= 0);
    n_lists = length(unique_lists);
    lists(idx) = n_lists;
    temp_word_duration = NaN(1,n_lists);
    temp_retrieval_duration = NaN(1,n_lists);
    temp_encoding_duration = NaN(1,n_lists);
    temp_math_duration = NaN(1,n_lists);
    
    for jdx = 1:n_lists
        this_list = unique_lists(jdx);
        orient_time = events.eegoffset(ismember(events.type,{'ORIENT','ORIENT_START'}) & events.list == this_list);
        word_times = events.eegoffset(strcmp(events.type,'WORD')  & events.list == this_list);
        retrieval_start = events.eegoffset(strcmp(events.type,'REC_START') & events.list == this_list);
        retrieval_stop = events.eegoffset(ismember(events.type,{'REC_END','E'}) & events.list == this_list);
        
        if length(word_times) < 10 || isempty(retrieval_start)
            lists(idx) = lists(idx)-1;
            events(events.list == this_list,:) = [];
            continue
        else
            retrieval_stop = retrieval_stop(1);
            temp_word_duration(jdx) = mean(diff(word_times),'omitnan');
            temp_retrieval_duration(jdx) = retrieval_stop-retrieval_start;
            if ~isempty(orient_time)
                temp_encoding_duration(jdx) = word_times(end) - orient_time;
            end
            temp_math_duration(jdx) = retrieval_start - word_times(end);
        end
    end
    word_duration(idx) = mean(temp_word_duration,'omitnan')/1000;
    retrieval_duration(idx) = mean(temp_retrieval_duration,'omitnan')/1000;
    encoding_duration(idx) = mean(temp_encoding_duration,'omitnan')/1000;
    math_duration(idx) = mean(temp_math_duration,'omitnan')/1000;
    
    %%% Saving a copy of events to analysis directory with the EEG file
    %%% path switched from reref (common average reference) to noreref so
    %%% that own referencing method can be used.
    if any(contains(events.eegfile,{'eeg.reref','eeg.bipolar'}))
        events.eegfile = strrep(events.eegfile,'eeg.bipolar','eeg.noreref');
        events.eegfile = strrep(events.eegfile,'eeg.reref','eeg.noreref');
    end
    
    %%% Save events copy
    save_path = fullfile(data_directory,this_subject,'FR1');
    if ~isfolder(save_path)
        mkdir(save_path);
    end
    save(fullfile(save_path,'events.mat'),'events');
    
end

%%% Save events info
stim_free_events_info = table(subject,subject_ID,word_duration,encoding_duration,math_duration,retrieval_duration,lists);
save(fullfile(statistics_directory,'stim_free_events_info.mat'),'stim_free_events_info');
end