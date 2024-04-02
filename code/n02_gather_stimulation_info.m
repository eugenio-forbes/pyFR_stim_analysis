function n02_gather_stimulation_info(varargin)
if isempty(varargin)
    %%% Directory information
    root_directory = '/directory/to/pyFR_stim_analysis';
else
    root_directory = varargin{1};
end

%%% List directories
subjects_directory = fullfile(root_directory,'subject_files');
list_directory = fullfile(root_directory,'lists');

%%% Load session list
load(fullfile(list_directory,'session_list.mat'),'session_list');

%%% Loop through sessions and load events to retrieve:
%%%  -label of stimulated electrode (field: stim_elec)
%%%  -amplitude of stimulation pulses in mA (field: stim_mamp)
%%%  -frequency of stimulation pulses in Hz (field: stim_hz)
%%%  -total number of lists completed in the session
%%%  -number of unique conditions in session (field: stimcode)
%%%  -whether each of these experimental conditions were present in the events:
%%%     -stimulation in encoding (ES): value of 1
%%%     -stimulation in retrieval (RS): value of 3
%%%     -no stimulation (NS): value of 5

n_sessions = height(session_list);

stimulation_electrode_label = cell(n_sessions,1);
stimulation_pulse_amplitude = cell(n_sessions,1);
stimulation_pulse_frequency = cell(n_sessions,1);
n_lists_completed = NaN(n_sessions,1);
n_conditions = NaN(n_sessions,1);
has_encoding_stimulation = false(n_sessions,1);
has_retrieval_stimulation = false(n_sessions,1);
has_no_stimulation = false(n_sessions,1);

for idx = 1:n_sessions
    %%% Get subject, task, and session names for each session
    this_subject = session_list.subject{idx};
    this_task = session_list.task{idx};
    this_session = session_list.session{idx};
    
    %%% To get events file path
    events_file = fullfile(subjects_directory,this_subject,'behavioral',this_task,this_session,'events.mat');
    
    %%% Load events (Matlab structure)
    load(events_file,'events');
    
    %%% Get stimulation electrode label
    stim_elec = {events.stim_elec};
    stim_elec = stim_elec(~cellfun(@isempty,stim_elec) & ~strcmp(stim_elec,'X')); %Remove empty rows or rows containing only 'X'.
    stim_elec = cellfun(@char,stim_elec,'UniformOutput',false); %A few were saved as strings instead of character arrays.
    label = upper(unique(stim_elec)); %A few were in lower case.
    if ~isempty(label)
        label = label{:};
        if contains(label,'-') %A few had bipolar notation. Stim electrode is first one
            hyphen_location = strfind(label,'-');
            label = label(1:hyphen_location-1);
        end
        stimulation_electrode_label{idx} = label;
    end
    
    %%% Get pulse amplitude
    stim_mamp = [events.stim_mamp];
    stim_mamp = stim_mamp(stim_mamp > 0);
    stimulation_pulse_amplitude{idx} = unique(stim_mamp);
    
    %%% Get pulse frequency
    stim_hz = [events.stim_hz];
    stim_hz = stim_hz(stim_hz > 0);
    stimulation_pulse_frequency{idx} = unique(stim_hz);    
    
    %%% Get number of lists completed
    encoding_events = events(strcmp({events.type},'WORD'));
    recall_stop_events = events(strcmp({events.type},'REC_STOP'));
    first_viewed_list = encoding_events(1).list;
    last_completed_list = recall_stop_events(end).list;
    n_lists_completed(idx) = last_completed_list - first_viewed_list + 1;
    
    %%% Get unique experimental conditions in session
    stim_codes = [events.stimcode];
    stim_codes = stim_codes(stim_codes>0);
    stim_codes = unique(stim_codes);
    
    n_conditions(idx) = length(stim_codes);
    has_encoding_stimulation(idx) = ismember(1,stim_codes);
    has_retrieval_stimulation(idx) = ismember(3,stim_codes);
    has_no_stimulation(idx) = ismember(5,stim_codes);    
end

%%% Add information to session list
session_list.stimulation_electrode_label = stimulation_electrode_label;
session_list.stimulation_pulse_amplitude = stimulation_pulse_amplitude;
session_list.stimulation_pulse_frequency = stimulation_pulse_frequency;
session_list.n_lists_completed = n_lists_completed;
session_list.n_conditions = n_conditions;
session_list.has_encoding_stimulation = has_encoding_stimulation;
session_list.has_retrieval_stimulation = has_retrieval_stimulation;
session_list.has_no_stimulation = has_no_stimulation;

%%% Some empty fields need to be manually added by looking at configuration
%%% files for session.

%%% Save updated session list
save(fullfile(list_directory,'session_list.mat'),'session_list');
end
