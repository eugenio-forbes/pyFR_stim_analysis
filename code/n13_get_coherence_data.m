function n13_get_coherence_data(varargin)
if isempty(varargin)
    %%% Directory information
    root_directory = '/directory/to/pyFR_stim_analysis';
    parpool_n = 16;
else
    root_directory = varargin{1};
    parpool_n = varargin{2};
end

%%% List directories
list_directory = fullfile(root_directory,'lists');
data_directory = fullfile(root_directory,'data');
table_directory = fullfile(root_directory,'tables');
error_directory = fullfile(root_directory,'error_logs');
if ~isfolder(table_directory)
    mkdir(table_directory);
end
if ~isfolder(error_directory)
    mkdir(error_directory);
end

%%% Define event types and stimcodes
encoding_type = 'WORD';
retrieval_type = 'REC_WORD';

%%% Load subject list
load(fullfile(list_directory,'subject_list.mat'),'subject_list');

%%% Load electrode list
load(fullfile(list_directory,'electrode_list.mat'),'electrode_list');

%%% Reduce subject and electrode list to those that had stimulation free
%%% session
has_stim_free_session = subject_list.has_stim_free_session;
subject_list = subject_list(has_stim_free_session,:);
has_stim_free_session = ismember(electrode_list.subject_ID,subject_list.subject_ID);
electrode_list = electrode_list(has_stim_free_session,:);

%%% Get counts for reduced subject and electrode list
n_subjects = height(subject_list);
n_electrodes = height(electrode_list);

%%% Initialize cell arrays to hold table data that will ultimately be
%%% merged into a single table for encoding, retrieval, encoding vs retrieval
%%% and anatomical variation LME analyses and correlation analyses
encoding_phase_clustering_data = cell(n_electrodes,1);
retrieval_phase_clustering_data = cell(n_electrodes,1);

%%%Transform electrode list columns into logicals in order to make tables
%%%less heavy
electrode_list.stimulation_lateral = strcmp(electrode_list.stimulation_group,'IPL'); %%% Positive t-statistics for stimulation group terms would indicate a greater effect for IPL
electrode_list.stimulation_left = strcmp(electrode_list.stimulation_hemisphere,'left'); %%% Used for random effect of hemisphere nested in stimulation_group
electrode_list.anterior = strcmp(electrode_list.longitudinal_location,'anterior'); %%%Used for separating tables and also positive t-statistics comparing hippocampal region effects would indicate a greater effect for anterior hippocampus.
electrode_list.left = strcmp(electrode_list.hemisphere,'left'); %%% Used for separating table and also positive t-statistics comparing hippocampal region effects would indicate a greater effect for left hippocampus.

%%%Transform subject list columns as well
subject_list.stimulation_lateral = strcmp(subject_list.stimulation_group,'IPL');
subject_list.stimulation_left = strcmp(subject_list.stimulation_hemisphere,'left');

%%%Initialize parpool if it hasn't been initialized
pool_object = gcp('nocreate');
if isempty(pool_object)
    parpool(parpool_n);
end

encoding_phase_clustering_exclusions = false(n_electrodes,1);
retrieval_phase_clustering_exclusions = false(n_electrodes,1);

parfor idx = 1:n_electrodes
    
    %%% Get information about electrode used for accessing data and processing
    subject = electrode_list.subject{idx};
    channel_number = electrode_list.channel_number(idx);
    reference_channels = electrode_list.reference_channels{idx};
    stimulation_electrode_channel = electrode_list.stimulation_electrode_channel(idx);
    stimulation_references = electrode_list.stimulation_references{idx};
    
    %%% Get information about electrode
    electrode_info = electrode_list(idx,{'subject_ID','session_ID','electrode_ID','anterior','left','stimulation_lateral','stimulation_left'});
    
    %%% Each electrode has a copy of events file to be able to load through
    %%% parfor with no risk of conflict between other electrodes of the
    %%% same sesion
    electrode_directory = fullfile(data_directory,subject,'FR1',num2str(channel_number));
    events_file = fullfile(electrode_directory,'events.mat');
    
    %%%Load events (saved as table)
    events = load_events_copy(events_file);
    
    %%% Get encoding phase clustering
    is_encoding = strcmp(events.type,encoding_type);
    encoding_events = events(is_encoding,:);
    [encoding_phase_clustering_table,error_flag,error_message] = get_phase_clustering_tables('encoding',encoding_events,channel_number,reference_channels,stimulation_electrode_channel,stimulation_references,electrode_info);
    if ~error_flag
        encoding_phase_clustering_data{idx} = encoding_phase_clustering_table;
    else
        encoding_phase_clustering_exclusions(idx) = true;
        error_file = fullfile(error_directory,sprintf('%s_%d_encoding_phase_error.txt',subject,channel_number));
        fid = fopen(error_file,'w');
        fprintf(fid,error_message); fclose(fid);
    end
    
    %%% Get retrieval phase clustering
    is_retrieval = strcmp(events.type,retrieval_type) & ~strcmp(events.item,'<>'); %%%excluding vocalizations
    retrieval_events = events(is_retrieval,:);
    [retrieval_phase_clustering_table,error_flag,error_message] = get_phase_clustering_tables('retrieval',retrieval_events,channel_number,reference_channels,stimulation_electrode_channel,stimulation_references,electrode_info);
    if ~error_flag
        retrieval_phase_clustering_data{idx} = retrieval_phase_clustering_table;
    else
        retrieval_phase_clustering_exclusions(idx) = true;
        error_file = fullfile(error_directory,sprintf('%s_%d_retrieval_phase_error.txt',subject,channel_number));
        fid = fopen(error_file,'w');
        fprintf(fid,error_message); fclose(fid);
    end
    
end

%%% Initialize cell arrays to identify electrodes of bad subjects and to
%%% hold data of good subjects
bad_electrode_cells = false(n_electrodes,1);
bad_subject_cells = false(n_subjects,1);
averaged_encoding_phase_data = cell(n_subjects,1);
averaged_retrieval_phase_data = cell(n_subjects,1);

for idx = 1:n_subjects
    %%% Get information about subject for accessing electrode data
    subject_ID = subject_list.subject_ID(idx);

    subject_info = subject_list(idx,{'subject_ID','stimulation_lateral','stimulation_left'});
    
    %%% Get subject electrode indices to retrieve data
    subject_electrodes = electrode_list.subject_ID == subject_ID;    
    n_subject_electrodes = sum(subject_electrodes);
    
    subject_encoding_data = encoding_phase_clustering_data(subject_electrodes);
    empty_encoding_cells = cellfun(@isempty,subject_encoding_data);
    
    %%% Exclude subjects and respective electrodes if more than half of the electrodes were excluded
    if sum(empty_encoding_cells) > floor(n_subject_electrodes/2)
        bad_electrode_cells(subject_electrodes) = true;
        bad_subject_cells(idx) = true;
    else
        subject_encoding_data(empty_encoding_cells) = [];
        averaged_encoding_phase_data{idx} = get_subject_average_data(subject_encoding_data,subject_info);
    end
    
    %%% Get retrieval data of subject's respective electrodes and cound
    %%% empty cells (excluded electrodes because of error or more than
    %%% half of the events were thrown out).
    subject_retrieval_data = retrieval_phase_clustering_data(subject_electrodes);
    empty_retrieval_cells = cellfun(@isempty,subject_retrieval_data);
     
    %%% Exclude subjects and respective electrodes if more than half of the electrodes were excluded
    if sum(empty_retrieval_cells) > floor(n_subject_electrodes/2)
        bad_electrode_cells(subject_electrodes) = true;
        bad_subject_cells(idx) = true;
    else
        subject_retrieval_data(empty_retrieval_cells) = [];
        averaged_retrieval_phase_data{idx} = get_subject_average_data(subject_retrieval_data,subject_info);
    end      
        
end

%%% Delete bad electrode cells and bad subject cells from data
encoding_phase_clustering_data(bad_electrode_cells) = [];
retrieval_phase_clustering_data(bad_electrode_cells) = [];
averaged_encoding_phase_data(bad_subject_cells) = [];
averaged_retrieval_phase_data(bad_subject_cells) = [];

%%% Concatenate all data to make tables ready for analysis
encoding_band_phase_clustering_table = vertcat(encoding_phase_clustering_data{:}); clear encoding_band_data
retrieval_band_phase_clustering_table = vertcat(retrieval_phase_clustering_data{:}); clear retrieval_band_data

encoding_subject_phase_clustering_table = vertcat(averaged_encoding_phase_data{:});
retrieval_subject_phase_clustering_table = vertcat(averaged_retrieval_phase_data{:});

%%% Save all tables
save(fullfile(table_directory,'encoding_band_phase_clustering_table.mat'),'encoding_band_phase_clustering_table');
save(fullfile(table_directory,'retrieval_band_phase_clustering_table.mat'),'retrieval_band_phase_clustering_table');
save(fullfile(table_directory,'encoding_subject_phase_clustering_table.mat'),'encoding_subject_phase_clustering_table');
save(fullfile(table_directory,'retrieval_subject_phase_clustering_table.mat'),'retrieval_subject_phase_clustering_table');

%%% Save exclusions to list_directory
save(fullfile(list_directory,'phase_exclusions.mat'),'bad_electrode_cells','bad_subject_cells',...
    'encoding_phase_clustering_exclusions','retrieval_phase_clustering_exclusions');
end

function events = load_events_copy(events_file)
load(events_file,'events');
end

function [phase_clustering_table,error_flag,error_message] = get_phase_clustering_tables(event_type,events,channel_number,reference_channels,stim_site_channel,stim_site_references,electrode_info)
%%%Initialize variables to return;
error_flag = false;
error_message = [];
phase_clustering_table=[];

try
    %%% Get power normalized to the average baseline of no stimulation condition
    [hippocampus_phase,frequencies,hippocampal_good_indices] = phase_WM_reference(event_type,events,channel_number,reference_channels);
    [stimulation_site_phase,~,stim_site_good_indices] = phase_WM_reference(event_type,events,stim_site_channel,stim_site_references);

    %%% Define conditions and frequency bands
    slow_theta_frequencies = frequencies >= 2 & frequencies <= 4;
    fast_theta_frequencies = frequencies >  4 & frequencies <= 8;
    
    %%% Get original number of events
    n_events = height(events);
    n_hippocampus = length(hippocampal_good_indices);
    n_stim_site = length(stim_site_good_indices);
    
    if ~isempty(hippocampus_phase) && ~isempty(stimulation_site_phase)
        
        %%% Get shared event indices between hippocampus and stimulation
        %%% site to match both phase matrices
        shared_indices_hippocampus = ismember(hippocampal_good_indices,stim_site_good_indices);
        shared_indices_stim_site = ismember(stim_site_good_indices,hippocampal_good_indices);
        
        n_shared = sum(shared_indices_hippocampus);
        
        if n_shared > floor(n_events/2)
            %%% If the amount of shared events is more than half the number
            %%% of events calculate intertrial phase clustering
            hippocampus_pase = hippocampus_phase(shared_indices_hippocampus,:,:);
            stimulation_site_phase = stimulation_site_phase(shared_indices_stim_site,:,:);
            
            phase_differences = hippocampus_pase - stimulation_site_phase;
            intertrial_phase_clustering = squeeze(abs(mean(exp(1i*phase_differences),1)));
            
            %%% Get mean intertrial phase clustering within band and across
            %%% event time
            slow_theta_phase_clustering = squeeze(mean(intertrial_phase_clustering(slow_theta_frequencies,:),'all'));
            fast_theta_phase_clustering = squeeze(mean(intertrial_phase_clustering(fast_theta_frequencies,:),'all'));
            
            %%% Add information about band to analysis band table
            phase_clustering_table = electrode_info;
            phase_clustering_table.slow_theta_phase_clustering = slow_theta_phase_clustering;
            phase_clustering_table.fast_theta_phase_clustering = fast_theta_phase_clustering;
            
        else
            percent_shared = n_shared/n_events;
            error('Less than half (%.2f) of good indices between hippocampus and stimulation site shared.',percent_shared);
        end
    else
        %%% Error related to events being thrown out. Use data to describe.
        if isempty(hippocampal_good_indices) && isempty(stim_site_good_indices)
            error('All events were thrown out for both sites.');
        else
            if isempty(hippocampal_good_indices)
                percent_hippocampus = 1;
            else
                percent_hippocampus = n_hippocampus/n_events;
            end
            
            if isempty(hippocampal_good_indices)
                percent_stim_site = 1;
            else
                percent_stim_site = n_stim_site/n_events;
            end
            
            error('%.2f percent events hippocampus remaining. %.2f percent stimulation site remaining.',percent_hippocampus,percent_stim_site);
        end
    end
catch this_error %%% Error could also be due to power extraction error
    error_flag = true;
    error_message = getReport(this_error, 'extended', 'hyperlinks', 'off');
end

end

function averaged_band_table = get_subject_average_data(subject_data,subject_info)
subject_ID = subject_info.subject_ID;
stimulation_lateral = subject_info.stimulation_lateral;
stimulation_left = subject_info.stimulation_left;

electrode_tables = vertcat(subject_data{:});

n_anatomical_locations = 4;
is_anterior_left = electrode_tables.anterior & electrode_tables.left;
is_anterior_right = electrode_tables.anterior & ~electrode_tables.left;
is_posterior_left = ~electrode_tables.anterior & electrode_tables.left;
is_posterior_right = ~electrode_tables.anterior & ~electrode_tables.left;

anatomical_indices = [{is_anterior_left};{is_anterior_right};{is_posterior_left};{is_posterior_right}];

averaged_data = cell(n_anatomical_locations,1);

for idx = 1:n_anatomical_locations
    these_indices = anatomical_indices{idx};
    if sum(these_indices) > 0
        these_electrodes = electrode_tables(these_indices,:);
        anterior = these_electrodes.anterior(1);
        left = these_electrodes.left(1);
        slow_theta_phase_clustering = mean(these_electrodes.slow_theta_phase_clustering);
        fast_theta_phase_clustering = mean(these_electrodes.fast_theta_phase_clustering);
        averaged_data{idx} = table(subject_ID,stimulation_lateral,stimulation_left,anterior,left,...
            slow_theta_phase_clustering,fast_theta_phase_clustering);
    end        
end

averaged_band_table = vertcat(averaged_data{:});

end

function [phase,frequencies,good_indices] = phase_WM_reference(event_type,events,channel_number,reference_channels)
%%% Signal parameters
%%% signal acquired at 1000 Hz, not resampling
sampling_rate = 1000;
nyquist = sampling_rate/2;
buffer = 500; %%% in ms, flanking periods of segment

switch event_type
    case 'encoding'
        duration = 1600; %%% in ms, the time of the baseline period (500ms) + 1500ms of word presentation + 100ms shortly after disappearance
        offset = 0; %%% in ms, prior to word display to acquire singal
    case 'retrieval'
        duration = 1000; %%% time for baseline period (500ms) + 1000 ms prior to onset of verbalization of item in recall
        offset = -1000; %%% time prior to onset of verbalization 
end

%%% Parameters for exclusion of events
kurtosis_threshold = 4;

%%% Filtering parameters to, highpass (0.5Hz) bandstop line noise (60Hz) 
%%% and harmonics.
filter_order = 3;
highpass_frequency = 0.5;
notch_frequencies = [(55:60:475)',(65:60:485)'];

%%% Parameters for wavelet convolution
frequencies = 2.^((8:64)/8); %%% Logarithmically spaced powers of 2 between 2 and 128

width = 6; % width of the wavelet for wavelet convolution

%%% Convert events to struct if needed
if istable(events)
    events = table2struct(events);
end

n_events = length(events);

%%%% Make array of indices to keep track of which events remained after
%%%% processing
good_indices = 1:n_events;
bad_eeg_offset = [events.eegoffset] + offset - buffer < 0;
good_indices(bad_eeg_offset) = [];
events(bad_eeg_offset) = [];

%%% Remove any offset and drift from signals with highpass over 0.5Hz prior to referencing
[highpass_b,highpass_a] = butter(filter_order,highpass_frequency/nyquist,'high');

%%% Get stim site eeg
stim_site_eeg = gete_ms(channel_number,events,duration+(buffer*2),offset-buffer,0);
stim_site_eeg = filtfilt(highpass_b,highpass_a,stim_site_eeg')';

%%% Gather eegs for reference channels
reference_eegs = NaN(length(events),duration+(buffer*2),length(reference_channels));
for idx = 1:length(reference_channels)
    temp_reference_eeg = gete_ms(reference_channels(idx),events,duration+buffer*2,offset,0);
    temp_reference_eeg = filtfilt(highpass_b,highpass_a,temp_reference_eeg')';
    reference_eegs(:,:,idx) = temp_reference_eeg; clear temp_reference_eeg;
end

%%% Subtract mean of reference channels to get referenced site eeg
stim_site_eeg = stim_site_eeg - squeeze(mean(reference_eegs,3));

%%% Notch filter site EEG signal
for idx = 1:size(notch_frequencies,1)
    filter_frequencies = notch_frequencies(idx,:);
    [notch_b,notch_a] = butter(filter_order,filter_frequencies/nyquist,'stop');
    stim_site_eeg = filtfilt(notch_b,notch_a,stim_site_eeg')';
end

%%% Identify events that exceed kurtosis threshold
kurtosis_values = kurtosis(stim_site_eeg(:,buffer+1:end-buffer)');
nan_values = isnan(kurtosis_values);
if all(nan_values)
    phase = [];
    good_indices = [];
    return; %%% This could mean the signal was flat across events. Stop processing.
end

bad_kurtosis_indices = isnan(kurtosis_values) | kurtosis_values > kurtosis_threshold;
good_indices(bad_kurtosis_indices) = [];
stim_site_eeg(bad_kurtosis_indices,:) = [];

if length(good_indices) < ceil(n_events/2)
    phase = [];
    return; %%% No need to continue since more than half of the events have been thrown out;
end

root_median_square = sqrt(sum((stim_site_eeg - median(stim_site_eeg,2)).^2,2));
artifactual = root_median_square >= prctile(root_median_square,95);

good_indices(artifactual) = [];
stim_site_eeg(artifactual,:) = [];

if length(good_indices) < ceil(n_events/2)
    phase = [];
    return; %%% No need to continue since more than half of the event were thrown out
end

%%% Get the raw power
[phase,~] = multiphasevec3(frequencies,stim_site_eeg,sampling_rate,width);

%%% Remove the buffer
phase = phase(:,:,buffer+1:end-buffer);
end
