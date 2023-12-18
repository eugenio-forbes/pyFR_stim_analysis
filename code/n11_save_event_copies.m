function n11_save_event_copies(varargin)
if isempty(varargin)
    %%% Directory information
    root_directory = '/directory/with/analysis_folder';
    analysis_folder_name = 'pyFR_stim_analysis';
else
    root_directory = varargin{1};
    analysis_folder_name = varargin{2};
end

%%% Declare directories
analysis_directory = fullfile(root_directory,analysis_folder_name);
list_directory = fullfile(analysis_directory,'lists');
data_directory = fullfile(analysis_directory,'data');

%%% Load subject list
load(fullfile(list_directory,'subject_list.mat'),'subject_list');

%%% Load session list
load(fullfile(list_directory,'session_list.mat'),'session_list');
n_sessions = height(session_list);

%%% Load electrode list
load(fullfile(list_directory,'electrode_list.mat'),'electrode_list');

%%% Loop through sessions to save pyFR_stim events copies for each
%%% electrode to be able to process data with parfor without risking
%%% conflicts from access to files

for idx = 1:n_sessions
    %%% Get information to load events fo
    subject = session_list.subject{idx};
    task = session_list.task{idx};
    session = session_list.session{idx};
    session_ID = session_list.session_ID(idx);
    
    %%% Get hipppocampal channel numbers for this session
    session_electrodes = electrode_list.session_ID == session_ID;
    channel_numbers = electrode_list.channel_number(session_electrodes);

    %%% Load events. These save events copies have already been converted to table
    session_directory = fullfile(data_directory,subject,task,session);
    events_file = fullfile(session_directory,'events.mat');
    load(events_file,'events')
    
    %%% Remove events with empty eegoffset or offset less than 0. Not done
    %%% previously to conserve all behavioral data
    no_file_or_offset = strcmp(events.eegfile,'') | isempty(events.eegoffset) | events.eegoffset <0;
    events(no_file_or_offset,:) = [];
    
    %%% Loop through channels to save events
    par_save(session_directory,channel_numbers,events);
end

%%% Loop through subject list of subjects with stim free session to save a
%%% copy of events for hippocampal and stimulation electrode channel to be
%%% able to process data with parfor loop
subject_list = subject_list(subject_list.has_stim_free_session,:);
n_subjects = height(subject_list);

for idx = 1:n_subjects
    %%% Get information to load stim free events
    subject = subject_list.subject{idx};
    subject_ID = subject_list.subject_ID(idx);
    stimulation_electrode_channel = subject_list.stimulation_electrode_channel(idx);
    
    %%% Get hipppocampal channel numbers for this subject
    subject_electrodes = electrode_list.subject_ID == subject_ID;
    channel_numbers = electrode_list.channel_number(subject_electrodes);
    channel_numbers = unique(channel_numbers); %%% Since some subjects have more than one pyFR_stim session
    all_channels = [channel_numbers;stimulation_electrode_channel];
    
    %%% Load events. These save events copies have already been converted to table
    session_directory = fullfile(data_directory,subject,'FR1');
    events_file = fullfile(session_directory,'events.mat');
    load(events_file,'events')
    
    %%% Remove events with empty eegoffset or offset less than 0.
    no_file_or_offset = strcmp(events.eegfile,'') | isempty(events.eegoffset) | events.eegoffset <0;
    events(no_file_or_offset,:) = [];
    
    %%% Loop through channels to save events
    par_save(session_directory,all_channels,events);
end

end

function par_save(session_directory,channel_numbers,events)
n_channels = length(channel_numbers);
for idx = 1:n_channels
    this_channel = channel_numbers(idx);
    save_path = fullfile(session_directory,num2str(this_channel));
    if ~isfolder(save_path)
        mkdir(save_path);
    end
    file_name = fullfile(save_path,'events.mat');
    save(file_name,'events')

end
end