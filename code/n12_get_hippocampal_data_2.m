function n12_get_hippocampal_data_2(varargin)
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
table_directory = fullfile(root_directory,'tables_time');
error_directory = fullfile(root_directory,'error_logs');
if ~isfolder(table_directory)
    mkdir(table_directory);
end
if ~isfolder(error_directory)
    mkdir(error_directory);
end

%%% Define event types and stimcodes
encoding_type = 'WORD';
retrieval_type = 'REC_START';
encoding_stimulation_code = 1;
retrieval_stimulation_code = 3;

%%% Load subject list
load(fullfile(list_directory,'subject_list.mat'),'subject_list');
n_subjects = height(subject_list);

%%% Load electrode list
load(fullfile(list_directory,'electrode_list.mat'),'electrode_list');
n_electrodes = height(electrode_list);

%%% Initialize cell arrays to hold table data that will ultimately be
%%% merged into a single table for encoding, retrieval, encoding vs retrieval
%%% and anatomical variation LME analyses and correlation analyses
encoding_data = cell(n_electrodes,1);
retrieval_data = cell(n_electrodes,1);

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

encoding_analysis_exclusions = false(n_electrodes,1);
retrieval_analysis_exclusions = false(n_electrodes,1);

for idx = 1:n_electrodes
    
    %%% Get information about electrode used for accessing data and processing
    subject = electrode_list.subject{idx};
    task = electrode_list.task{idx};
    session = electrode_list.session{idx};
    channel_number = electrode_list.channel_number(idx);
    reference_channels = electrode_list.reference_channels{idx};
    has_encoding_stimulation = electrode_list.has_encoding_stimulation(idx);
    has_retrieval_stimulation = electrode_list.has_retrieval_stimulation(idx);
    has_no_stimulation = electrode_list.has_no_stimulation(idx);
    
    %%% Get information about electrode that will go into LME analysis tables
    electrode_info = electrode_list(idx,{'subject_ID','session_ID','electrode_ID','anterior','left','stimulation_lateral','stimulation_left'});
    
    %%% Each electrode has a copy of events file to be able to load through
    %%% parfor with no risk of conflict between other electrodes of the
    %%% same sesion
    electrode_directory = fullfile(data_directory,subject,task,session,num2str(channel_number));
    events_file = fullfile(electrode_directory,'events.mat');
    
    %%%Load events (saved as table)
    events = load_events_copy(events_file);
    rec_start = events.eegoffset(strcmp(events.type,'REC_START'));
    rec_stop = events.eegoffset(strcmp(events.type,'REC_STOP'));
    retrieval_duration = mean(rec_stop-rec_start)/1000;
    if retrieval_duration < 17500
        retrieval_duration = 14500;
    elseif retrieval_duration < 25000
        retrieval_duration = 19500;
    else
        retrieval_duration = 29500;
    end
    electrode_info.retrieval_duration = retrieval_duration;
    
    %%% Only sessions that had both encoding stimulation and control (no stimualtion) included in encoding analysis
    if has_encoding_stimulation && has_no_stimulation
        is_encoding = strcmp(events.type,encoding_type) & events.serialpos == 1;
        encoding_events = events(is_encoding,:);
        encoding_conditions = encoding_events.stimcode == encoding_stimulation_code; %%% Positive t-statistics would indicate greater power for stimulation condition
        [encoding_power_table,error_flag,error_message] = get_power_analysis_tables('encoding',encoding_events,encoding_conditions,channel_number,reference_channels,electrode_info);
        if ~error_flag
            encoding_data{idx} = encoding_power_table;
        else
            encoding_analysis_exclusions(idx) = true;
            error_file = fullfile(error_directory,sprintf('%s_%s_%s_%d_encoding_power_error.txt',subject,task,session,channel_number));
            fid = fopen(error_file,'w');
            fprintf(fid,error_message); fclose(fid);
        end 
    end
    
    %%% Only sessions that had both retrieval stimulation and control included in retrieval analysis
    if has_retrieval_stimulation && has_no_stimulation
        is_retrieval = strcmp(events.type,retrieval_type);
        retrieval_events = events(is_retrieval,:);
        retrieval_conditions = retrieval_events.stimcode == retrieval_stimulation_code; %%% Positive t-statistics would indicate greater power for stimulation condition
        [retrieval_power_table,error_flag,error_message] = get_power_analysis_tables('retrieval',retrieval_events,retrieval_conditions,channel_number,reference_channels,electrode_info);
        if ~error_flag
            retrieval_data{idx} = retrieval_power_table;
        else
            retrieval_analysis_exclusions(idx) = true;
            error_file = fullfile(error_directory,sprintf('%s_%s_%s_%d_retrieval_power_error.txt',subject,task,session,channel_number));
            fid = fopen(error_file,'w');
            fprintf(fid,error_message); fclose(fid);
        end 
    end
   
end

%%% Initialize cell arrays to identify electrodes of bad subjects and to
%%% hold data of good subjects
bad_electrode_cells = false(n_electrodes,1);

for idx = 1:n_subjects
    %%% Get information about subject for accessing electrode data
    subject_ID = subject_list.subject_ID(idx);

    subject_info = subject_list(idx,{'subject_ID','stimulation_lateral','stimulation_left'});
    
    %%% Get subject electrode indices to retrieve data
    subject_electrodes = electrode_list.subject_ID == subject_ID;    
    n_subject_electrodes = sum(subject_electrodes);
    
    sample_electrode = find(subject_electrodes,1,'first');
    
    has_encoding_stimulation = electrode_list.has_encoding_stimulation(sample_electrode);
    has_retrieval_stimulation = electrode_list.has_retrieval_stimulation(sample_electrode);
    has_no_stimulation = electrode_list.has_no_stimulation(sample_electrode);
    
    if has_encoding_stimulation && has_no_stimulation
        %%% Get retrieval data of subject's respective electrodes and cound
        %%% empty cells (excluded electrodes because of error or more than
        %%% half of the events were thrown out).
        subject_encoding_data = encoding_data(subject_electrodes);
        empty_encoding_cells = cellfun(@isempty,subject_encoding_data);
        
        %%% Exclude subjects and respective electrodes if more than half of the electrodes were excluded
        if sum(empty_encoding_cells) > floor(n_subject_electrodes/2)
            bad_electrode_cells(subject_electrodes) = true;
        end        
    end
    
    if has_retrieval_stimulation && has_no_stimulation
        %%% Get retrieval data of subject's respective electrodes and cound
        %%% empty cells (excluded electrodes because of error or more than
        %%% half of the events were thrown out).
        subject_retrieval_data = retrieval_data(subject_electrodes); 
        empty_retrieval_cells = cellfun(@isempty,subject_retrieval_data);
        
        %%% Exclude subjects and respective electrodes if more than half of the electrodes were excluded
        if sum(empty_retrieval_cells) > floor(n_subject_electrodes/2)
            bad_electrode_cells(subject_electrodes) = true;
        end      
    end
    
end

%%% Delete bad electrode cells and bad subject cells from data
encoding_data(bad_electrode_cells) = [];
retrieval_data(bad_electrode_cells) = [];

%%% Concatenate all data to make tables ready for analysis
encoding_time_table = vertcat(encoding_data{:}); clear encoding_data
retrieval_time_table = vertcat(retrieval_data{:}); clear retrieval_data

%%% Save all tables
save(fullfile(table_directory,'encoding_time_table.mat'),'encoding_time_table','-v7.3');
save(fullfile(table_directory,'retrieval_time_table.mat'),'retrieval_time_table','-v7.3');
end

function events = load_events_copy(events_file)
load(events_file,'events');
end

function [analysis_table,error_flag,error_message] = get_power_analysis_tables(event_type,events,active_stimulation,channel_number,reference_channels,electrode_info)
%%%Initialize variables to return;
error_flag = false;
error_message = [];
analysis_table = [];
retrieval_duration = electrode_info.retrieval_duration;

try
    %%% Get power normalized to the average baseline of no stimulation condition
    [power,frequencies,good_indices] = power_WM_reference(event_type,events,active_stimulation,channel_number,reference_channels,retrieval_duration);
    
    if ~isempty(power)
        n_frequencies = length(frequencies);
        
        total_time = size(power,2);
        first_segment_indices = 1:ceil(total_time/5);
        last_segment_indices = ceil(total_time-total_time/5):total_time;
        
        first_segments = squeeze(mean(power(:,:,first_segment_indices),3));
        last_segments = squeeze(mean(power(:,:,last_segment_indices),3));
        
        %%% Define conditions and frequency bands
        power = [first_segments(:);last_segments(:)];
        frequency = repmat(1:n_frequencies,length(good_indices),1);
        frequency = frequency(:);
        frequency = repmat(frequency,2,1);
        condition = repmat(active_stimulation(good_indices),n_frequencies*2,1);
        trial = repmat((1:length(good_indices))',n_frequencies*2,1);
        time = [zeros(length(good_indices)*n_frequencies,1);ones(length(good_indices)*n_frequencies,1)];
       
        %%% Multiply by 1000 to conserve one decimal of percent change (or 3 decimals of tstats)
        %%% and convert to int16 (so possible range of -3276.7% to 3276.7% change in power
        %%% and .1% least significant bit of change in power) (for tstats possible range of
        %%% -32.767 to 32.767 and least significant bit of .001)
        
        %%% Repeat electrode info n_frequencies*n_samples number of rows to make power table.
        analysis_table = repmat(electrode_info,n_frequencies*length(good_indices)*2,1);
        analysis_table.power = int16(power*1000);
        analysis_table.condition = condition;
        analysis_table.frequency = frequency;
        analysis_table.time = time;
        analysis_table.trial = trial;
        
    else
        %%% Error related to events being thrown out. Use data to describe.
        if isempty(good_indices)
            error('All events were thrown out');
        else
            remaining_conditions = active_stimulation(good_indices);
            percent_events = length(remaining_conditions)/length(active_stimulation);
            percent_stimulation = sum(remaining_conditions)/sum(active_stimulation);
            percent_control = sum(~remaining_conditions)/sum(~active_stimulation);
            error('%.2f of events remaining. %.2f stimulation remaining, %.2f control remaining.\n',percent_events,percent_stimulation,percent_control)
        end
    end
catch this_error %%% Error could also be due to power extraction error
    error_flag = true;
    error_message = getReport(this_error, 'extended', 'hyperlinks', 'off');
end

end

function [power,frequencies,good_indices] = power_WM_reference(event_type,events,active_stimulation,channel_number,reference_channels,retrieval_duration)
%%% Signal parameters
%%% signal acquired at 1000 Hz, not resampling
sampling_rate = 1000;
resampled_rate = 500;
nyquist = resampled_rate/2;
buffer = 500; %%% in ms, flanking periods of segment
baseline_duration = 300;

switch event_type
    case 'encoding'
        event_duration = 24000; %%% in ms, the time of the baseline period (500ms) + 1500ms of word presentation + 100ms shortly after disappearance
        event_offset = -500; %%% in ms, prior to word display to acquire singal
        baseline_offset = -1300; %%% samples to be used for calculating baseline power
    case 'retrieval'
        event_duration = retrieval_duration; %%% time for baseline period (500ms) + 1000 ms prior to onset of verbalization of item in recall
        event_offset = 500; %%% time prior to onset of verbalization 
        baseline_offset = 0;
end

%%% Parameters for exclusion of events
kurtosis_threshold = 4;

%%% Filtering parameters to, highpass (0.5Hz) bandstop line noise (60Hz) 
%%% and stimulation noise (100Hz) and harmonics.
filter_order = 1;
notch_frequencies = [58 62;98 102;118 122];

%%% Parameters for wavelet convolution
frequencies = 2.^((8:56)/8); %%% Logarithmically spaced powers of 2 between 2 and 128

width = 6; % width of the wavelet for wavelet convolution

%%% Convert events to struct if needed
if istable(events)
    events = table2struct(events);
end

n_events = length(events);
n_stimulation = sum(active_stimulation);
n_control = sum(~active_stimulation);

%%%% Make array of indices to keep track of which events remained after
%%%% processing
good_indices = 1:n_events;
bad_eeg_offset = [events.eegoffset] + event_offset - buffer < 0;
good_indices(bad_eeg_offset) = [];
events(bad_eeg_offset) = [];
active_stimulation(bad_eeg_offset) = [];    

%%% Get hippocampal eeg
hippocampal_eeg = gete_ms(channel_number,events,event_duration+(buffer*2),event_offset-buffer,0);
hippocampal_eeg = filtfilt(highpass_b,highpass_a,hippocampal_eeg')';

%%% Gather eegs for reference channels
reference_eegs = NaN(length(events),event_duration+(buffer*2),length(reference_channels));
for idx = 1:length(reference_channels)
    temp_reference_eeg = gete_ms(reference_channels(idx),events,event_duration+buffer*2,event_offset,0);
    temp_reference_eeg = filtfilt(highpass_b,highpass_a,temp_reference_eeg')';
    reference_eegs(:,:,idx) = temp_reference_eeg; clear temp_reference_eeg;
end

%%% Subtract mean of reference channels to get reference hippocampal eeg
hippocampal_eeg = hippocampal_eeg - squeeze(mean(reference_eegs,3));

%%% Notch filter hippocampal EEG signal
for idx = 1:size(notch_frequencies,1)
    filter_frequencies = notch_frequencies(idx,:);
    [notch_b,notch_a] = butter(filter_order,filter_frequencies/nyquist,'stop');
    hippocampal_eeg = filtfilt(notch_b,notch_a,hippocampal_eeg')';
end

%%% Identify events that exceed kurtosis threshold
kurtosis_values = kurtosis(hippocampal_eeg(:,buffer+1:end-buffer)');
nan_values = isnan(kurtosis_values);
if all(nan_values)
    power = [];
    good_indices = [];
    return; %%% This could mean the signal was flat across events. Stop processing.
end

bad_kurtosis_indices = isnan(kurtosis_values) | kurtosis_values > kurtosis_threshold;
good_indices(bad_kurtosis_indices) = [];
hippocampal_eeg(bad_kurtosis_indices,:) = [];
active_stimulation(bad_kurtosis_indices,:) = [];

if isempty(good_indices)
    power = [];
    return; %%% No need to continue since more than half of the events have been thrown out;
end

root_median_square = sqrt(sum((hippocampal_eeg - median(hippocampal_eeg,2)).^2,2));
artifactual = root_median_square >= prctile(root_median_square,95);

good_indices(artifactual) = [];
hippocampal_eeg(artifactual,:) = [];
active_stimulation(artifactual,:) = [];

if isempty(good_indices)
    power = [];
    return; %%% No need to continue since more than half of the event were thrown out
end

%%% Get the raw power
[~,raw_power] = multiphasevec3(frequencies,hippocampal_eeg,sampling_rate,width);

%%% Remove the buffer
raw_power = raw_power(:,:,buffer+1:end-buffer);

%%% Normalize the power
baseline_power = raw_power(:,:,baseline_indices);
baseline_power(active_stimulation,:,:) = [];
baseline_power = mean(baseline_power,[1 3]);
power = raw_power(:,:,event_indices); clear raw_power
power = (power - baseline_power)./baseline_power;
end
