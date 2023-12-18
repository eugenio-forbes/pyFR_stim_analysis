function n14_get_stimulation_site_data(varargin)
if isempty(varargin)
    %%% Directory information
    root_directory = '/directory/with/analysis_folder';
    analysis_folder_name = 'pyFR_stim_analysis';
else
    root_directory = varargin{1};
    analysis_folder_name = varargin{2};
end

%%% List directories
analysis_directory = fullfile(root_directory,analysis_folder_name);
list_directory = fullfile(analysis_directory,'lists');
data_directory = fullfile(analysis_directory,'data');
table_directory = fullfile(analysis_directory,'tables');
error_directory = fullfile(analysis_directory,'error_logs');
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

%%% Reduce subject and electrode list to those that had stimulation free
%%% session
has_stim_free_session = subject_list.has_stim_free_session;
subject_list = subject_list(has_stim_free_session,:);
n_subjects = height(subject_list);

%%% Initialize cell arrays to hold table data that will ultimately be
%%% merged into a single table for encoding and retrieval plots of
%%% stimulation site in stimulation free session of free recall
encoding_data = cell(n_subjects,1);
retrieval_data = cell(n_subjects,1);

%%% Initialize cell arrays to hold directories of electrodes included in
%%% fooof analysis
% encoding_fooof_directories_stim_site = cell(n_subjects,1);
% retrieval_fooof_directories_stim_site = cell(n_subjects,1);

%%%Transform subject list columns as well
subject_list.stimulation_lateral = strcmp(subject_list.stimulation_group,'IPL');
subject_list.stimulation_left = strcmp(subject_list.stimulation_hemisphere,'left');

%%%Initialize parpool if it hasn't been initialized
pool_object = gcp('nocreate');
if isempty(pool_object)
    parpool(24);
end

encoding_stim_site_exclusions = false(n_subjects,1);
retrieval_stim_site_exclusions = false(n_subjects,1);
% fooof_bosc_stim_site_exclusions = false(n_subjects,1);

parfor idx = 1:n_subjects
    
    %%% Get information about subject and stimulation site used for accessing data and processing
    subject = subject_list.subject{idx};
    stimulation_electrode_channel = subject_list.stimulation_electrode_channel(idx);
    stimulation_references = subject_list.stimulation_references{idx};
    
    %%% Get information about electrode that will go into LME analysis tables
    subject_info = subject_list(idx,{'subject_ID','stimulation_lateral','stimulation_left','stimulation_coordinates'});

    %%% For stimulation site get events from subject´s directory and save
    %%% data to stimulation electrode directory
    subject_directory = fullfile(data_directory,subject,'FR1')
    electrode_directory = fullfile(subject_directory,num2str(stimulation_electrode_channel));
    events_file = fullfile(subject_directory,'events.mat');
    
    %%%Load events (saved as table)
    events = load_events_copy(events_file);
    
    %%% Get encoding data
    is_encoding = strcmp(events.type,encoding_type);
    encoding_events = events(is_encoding,:);
    [encoding_power_table,error_flag,error_message] = get_stim_site_tables('encoding',encoding_events,stimulation_electrode_channel,stimulation_references,subject_info);
    if ~error_flag
        encoding_data{idx} = encoding_power_table;
    else
        encoding_stim_site_exclusions(idx) = true;
        error_file = fullfile(error_directory,sprintf('%s_%d_encoding_stim_site_error.txt',subject,stimulation_electrode_channel));
        fid = fopen(error_file,'w');
        fprintf(fid,error_message); fclose(fid);
    end
    
    %%% Get retrieval power data
    is_retrieval = strcmp(events.type,retrieval_type) & ~strcmp(events.item,'<>'); %%%excluding vocalizations
    retrieval_events = events(is_retrieval,:);
    [retrieval_power_table,error_flag,error_message] = get_stim_site_tables('retrieval',retrieval_events,stimulation_electrode_channel,stimulation_references,subject_info);
    if ~error_flag
        retrieval_data{idx} = retrieval_power_table;
    else
        retrieval_stim_site_exclusions(idx) = true;
        error_file = fullfile(error_directory,sprintf('%s_%d_retrieval_stim_site_error.txt',subject,stimulation_electrode_channel));
        fid = fopen(error_file,'w');
        fprintf(fid,error_message); fclose(fid);
    end
       
    [error_flag,error_message] = get_fooof_and_bosc_data_stim_site('encoding',encoding_events,electrode_directory,stimulation_electrode_channel,stimulation_references);
    if ~error_flag
        encoding_fooof_directories_stim_site{idx} = electrode_directory;
    else
        error_file = fullfile(error_directory,sprintf('%s_%d_stim_site_encoding_pre-fBOSC_error.txt',subject,stimulation_electrode_channel));
        fid = fopen(error_file,'w');
        fprintf(fid,error_message); fclose(fid);
    end
    
    [error_flag,error_message] = get_fooof_and_bosc_data_stim_site('retrieval',retrieval_events,electrode_directory,stimulation_electrode_channel,stimulation_references);
    if ~error_flag
        retrieval_fooof_directories_stim_site{idx} = electrode_directory;
    else
        fooof_bosc_stim_site_exclusions(idx) = true;
        error_file = fullfile(error_directory,sprintf('%s_%d_stim_site_retrieval_pre-fBOSC_error.txt',subject,stimulation_electrode_channel));
        fid = fopen(error_file,'w');
        fprintf(fid,error_message); fclose(fid);
    end
end

%%% Concatenate all data to make tables ready for analysis
encoding_stim_site_table = vertcat(encoding_data{:}); clear encoding_data
retrieval_stim_site_table = vertcat(retrieval_data{:}); clear retrieval_data

%%% Save all tables
save(fullfile(table_directory,'encoding_stim_site_table.mat'),'encoding_stim_site_table','-v7.3');
save(fullfile(table_directory,'retrieval_stim_site_table.mat'),'retrieval_stim_site_table','-v7.3');

%%% Print non-empty fooof data directories to text file in list directory
empty_cells = cellfun(@isempty,encoding_fooof_directories_stim_site);
encoding_fooof_directories_stim_site(empty_cells) = [];
writecell(encoding_fooof_directories_stim_site,fullfile(list_directory,'encoding_fooof_directories_stim_site.txt'))

empty_cells = cellfun(@isempty,retrieval_fooof_directories_stim_site);
retrieval_fooof_directories_stim_site(empty_cells) = [];
writecell(retrieval_fooof_directories_stim_site,fullfile(list_directory,'retrieval_fooof_directories_stim_site.txt'))

%%% Save exclusions to list_directory
save(fullfile(list_directory,'stim_site_exclusions.mat'),'encoding_stim_site_exclusions','retrieval_stim_site_exclusions','fooof_bosc_stim_site_exclusions');
end

function events = load_events_copy(events_file)
load(events_file,'events');
end

function [stim_site_table,error_flag,error_message] = get_stim_site_tables(event_type,events,stimulation_electrode_channel,stimulation_references,subject_info)
%%%Initialize variables to return;
error_flag = false;
error_message = [];
stim_site_table=[];

try
    %%% Get power that is zscore normalized to the baseline period
    [power,frequencies,good_indices] = tstat_power_WM_reference(event_type,events,stimulation_electrode_channel,stimulation_references);
    
    if ~isempty(power)
        
        %%% Define frequency bands
        slow_theta_frequencies = frequencies >= 2 & frequencies <= 4;
        fast_theta_frequencies = frequencies >  4 & frequencies <= 8;
        
        %%% Get median across events of mean within band/across time for slow and fast theta power
        slow_theta_power = mean(power(slow_theta_frequencies));
        fast_theta_power = mean(power(fast_theta_frequencies));
        
        %%% Add information about band to analysis band table
        stim_site_table = subject_info;
        stim_site_table.slow_theta_power = int16(slow_theta_power*1000);
        stim_site_table.fast_theta_power = int16(fast_theta_power*1000);
        
    else
        %%% Error related to events being thrown out. Use data to describe.
        if isempty(good_indices)
            error('All events were thrown out');
        else
            percent_events = length(good_indices)/height(events);
            error('%.2f of events remaining.\n',percent_events)
        end
    end
catch this_error %%% Error could also be due to power extraction error
    error_flag = true;
    error_message = getReport(this_error, 'extended', 'hyperlinks', 'off');
end

end

function [power,frequencies,good_indices] = tstat_power_WM_reference(event_type,events,stimulation_electrode_channel,stimulation_references)
%%% Signal parameters
%%% signal acquired at 1000 Hz, not resampling
sampling_rate = 1000;
nyquist = sampling_rate/2;
buffer = 500; %%% in ms, flanking periods of segment

switch event_type
    case 'encoding'
        duration = 2100; %%% in ms, the time of the baseline period (500ms) + 1500ms of word presentation + 100ms shortly after disappearance
        offset = -500; %%% in ms, prior to word display to acquire singal
        baseline_indices = 1:300; %%% samples to be used for calculating baseline power
        event_indices = 501:2100; %%% sample to be used for analysis of ecoding;
    case 'retrieval'
        duration = 1500; %%% time for baseline period (500ms) + 1000 ms prior to onset of verbalization of item in recall
        offset = -1500; %%% time prior to onset of verbalization 
        baseline_indices = 1:300;
        event_indices = 501:1500; 
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

%%% Get hippocampal eeg
stim_site_eeg = gete_ms(stimulation_electrode_channel,events,duration+(buffer*2),offset-buffer,0);
stim_site_eeg = filtfilt(highpass_b,highpass_a,stim_site_eeg')';

%%% Gather eegs for reference channels
reference_eegs = NaN(length(events),duration+(buffer*2),length(stimulation_references));
for idx = 1:length(stimulation_references)
    temp_reference_eeg = gete_ms(stimulation_references(idx),events,duration+buffer*2,offset,0);
    temp_reference_eeg = filtfilt(highpass_b,highpass_a,temp_reference_eeg')';
    reference_eegs(:,:,idx) = temp_reference_eeg; clear temp_reference_eeg;
end

%%% Subtract mean of reference channels to get reference hippocampal eeg
stim_site_eeg = stim_site_eeg - squeeze(mean(reference_eegs,3));

%%% Notch filter hippocampal EEG signal
for idx = 1:size(notch_frequencies,1)
    filter_frequencies = notch_frequencies(idx,:);
    [notch_b,notch_a] = butter(filter_order,filter_frequencies/nyquist,'stop');
    stim_site_eeg = filtfilt(notch_b,notch_a,stim_site_eeg')';
end

%%% Identify events that exceed kurtosis threshold
kurtosis_values = kurtosis(stim_site_eeg(:,buffer+1:end-buffer)');
nan_values = isnan(kurtosis_values);
if all(nan_values)
    power = [];
    good_indices = [];
    return; %%% This could mean the signal was flat across events. Stop processing.
end

bad_kurtosis_indices = isnan(kurtosis_values) | kurtosis_values > kurtosis_threshold;
good_indices(bad_kurtosis_indices) = [];
stim_site_eeg(bad_kurtosis_indices,:) = [];

if length(good_indices) < ceil(n_events/2)
    power = [];
    return; %%% No need to continue since more than half of the events have been thrown out;
end

root_median_square = sqrt(sum((stim_site_eeg - median(stim_site_eeg,2)).^2,2));
artifactual = root_median_square >= prctile(root_median_square,95);

good_indices(artifactual) = [];
stim_site_eeg(artifactual,:) = [];

if length(good_indices) < ceil(n_events/2)
    power = [];
    return; %%% No need to continue since more than half of the event were thrown out
end

%%% Get the raw power
[~,raw_power] = multiphasevec3(frequencies,stim_site_eeg,sampling_rate,width);

%%% Remove the buffer
raw_power = raw_power(:,:,buffer+1:end-buffer);

%%% Normalize the power
baseline_power = raw_power(:,:,baseline_indices);
event_power = raw_power(:,:,event_indices); clear raw_power
mean_baseline_power = squeeze(mean(baseline_power,3));
mean_event_power = squeeze(mean(event_power,3));
[~,~,~,stats] = ttest(mean_event_power,mean_baseline_power,'Dim',1,'tail','both');
power = stats.tstat;
end

function [error_flag,error_message] = get_fooof_and_bosc_data_stim_site(event_type,events,saving_directory,stimulation_electrode_channel,stimulation_references)
%%% Initialize output variables
error_flag = false;
error_message = [];

try
    %%% Signal parameters
    %%% signal acquired at 1000 Hz, not resampling
    sampling_rate = 1000;
    nyquist = sampling_rate/2;
    padding = 2000; %%% in ms, flanking periods of segment (BOSC recommends 2 sets of 1000 pads for their algorithm)
    
    switch event_type
        case 'encoding'
            duration = 1600; %%% in ms, the time of the baseline period (500ms) + 1500ms of word presentation + 100ms shortly after disappearance
            offset = 0; %%% in ms, prior to word display to acquire singal
            event_indices = [1+padding,duration+padding];
        case 'retrieval'
            duration = 1000; %%% time for baseline period (500ms) + 1000 ms prior to onset of verbalization of item in recall
            offset = -1000; %%% time prior to onset of verbalization
            event_indices = [1+padding,duration+padding];
    end
    
    %%% Parameters for exclusion of events
    kurtosis_threshold = 4;
    
    %%% Filtering parameters to highpass (0.5Hz)
    %%% Only high pass filtering to perform rereferncing. Avoiding other
    %%% filters to minimize distortions of the signal and the spectrum.
    %%% Aiming to reduce noise with rereferencing. Only fitting frequencies
    %%% from 2 to
    filter_order = 3;
    highpass_frequency = 0.5;

    %%% Parameters for wavelet convolution
    %%%Since power was calculated for 5000ms or more, conserving a
    %%%frequency resolution of 0.25 for generating power spectra. Fitting
    %%% 1 to 40 to minimize influence of line noise power in fitting of
    %%% aperiodic component. fooof requires linearly increasing frequencies
    frequencies = 1:0.25:45;
    width = 6; % width of the wavelet for wavelet convolution
    
    %%% Convert events to struct if needed
    if istable(events)
        events = table2struct(events);
    end
    
    n_events = length(events);
    
    %%%% Make array of indices to keep track of which events remained after
    %%%% processing
    good_indices = 1:n_events;
    bad_eeg_offset = [events.eegoffset] + offset - padding < 0;
    good_indices(bad_eeg_offset) = [];
    events(bad_eeg_offset) = [];
    
    %%% Remove any offset and drift from signals with highpass over 0.5Hz prior to referencing
    [highpass_b,highpass_a] = butter(filter_order,highpass_frequency/nyquist,'high');
    
    %%% Get hippocampal eeg
    stim_site_eeg = gete_ms(stimulation_electrode_channel,events,duration+(padding*2),offset-padding,0);
    stim_site_eeg = filtfilt(highpass_b,highpass_a,stim_site_eeg')';
    
    %%% Gather eegs for reference channels
    reference_eegs = NaN(length(events),duration+(padding*2),length(stimulation_references));
    for idx = 1:length(stimulation_references)
        temp_reference_eeg = gete_ms(stimulation_references(idx),events,duration+(padding*2),offset-padding,0);
        temp_reference_eeg = filtfilt(highpass_b,highpass_a,temp_reference_eeg')';
        reference_eegs(:,:,idx) = temp_reference_eeg; clear temp_reference_eeg;
    end
    
    %%% Subtract mean of reference channels to get referenced stim_site eeg
    stim_site_eeg = stim_site_eeg - squeeze(mean(reference_eegs,3));
    
    %%% Identify events that exceed kurtosis threshold
    kurtosis_values = kurtosis(stim_site_eeg(:,padding+1:end-padding)');
    nan_values = isnan(kurtosis_values);
    if all(nan_values)
        error('All kurtosis values NaN.');
    end
    
    bad_kurtosis_indices = isnan(kurtosis_values) | kurtosis_values > kurtosis_threshold;
    good_indices(bad_kurtosis_indices) = [];
    stim_site_eeg(bad_kurtosis_indices,:) = [];
    
    if length(good_indices) < ceil(n_events/2)
        error('More than half of the events were thrown out.');
    end
    
    root_median_square = sqrt(sum((stim_site_eeg - median(stim_site_eeg,2)).^2,2));
    artifactual = root_median_square >= prctile(root_median_square,95);
    
    good_indices(artifactual) = [];
    stim_site_eeg(artifactual,:) = [];
    
    if length(good_indices) < ceil(n_events/2)
        error('More than half of the events were thrown out.');
    end
    
    %%% Get the raw power
    [~,raw_power] = multiphasevec3(frequencies,stim_site_eeg,sampling_rate,width);
    
    %%% Get arithmetic mean across time to get power spectra for fooof
    power_spectra = squeeze(mean(raw_power(:,:,event_indices(1):event_indices(end)),3));
    
    %%% Get geometric mean for fitting aperiodic component to calculate
    %%% power threshold for fBOSC
    mean_power_spectrum = 10.^(squeeze(mean(log10(raw_power(:,:,event_indices(1):event_indices(end))),[1 3])));
    
    %%% Save data for fooof
    fooof_directory = fullfile(saving_directory,sprintf('%s_fooof',event_type));
    if ~isfolder(fooof_directory)
        mkdir(fooof_directory);
    end
    save(fullfile(fooof_directory,'data.mat'),'power_spectra','frequencies','-v7.3');

    %%% Save data for fBOSC
    bosc_directory = fullfile(saving_directory,sprintf('%s_bosc',event_type));
    if ~isfolder(bosc_directory)
        mkdir(bosc_directory);
    end
    save(fullfile(bosc_directory,'data.mat'),'mean_power_spectrum','frequencies','good_indices','padding','event_indices','-v7.3');
    
catch this_error %%% Error could also result from power extraction
    error_flag = true;
    error_message = getReport(this_error, 'extended', 'hyperlinks', 'off');
end
end