function n12_get_hippocampal_data(varargin)
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
encoding_vs_retrieval_data = cell(n_electrodes,1);
encoding_band_data = cell(n_electrodes,1);
retrieval_band_data = cell(n_electrodes,1);

%%% Initialize cell arrays to hold directories of elctrodes included in
%%% fooof analysis
encoding_fooof_directories_hippocampus = cell(n_electrodes,1);
retrieval_fooof_directories_hippocampus = cell(n_electrodes,1);

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
encoding_vs_retrieval_analysis_exclusions = false(n_electrodes,1);
fooof_bosc_analysis_exclusions = false(n_electrodes,1);

parfor idx = 1:n_electrodes
    
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
    
    %%% Only sessions that had both encoding stimulation and control (no stimualtion) included in encoding analysis
    if has_encoding_stimulation && has_no_stimulation
        is_encoding = strcmp(events.type,encoding_type);
        encoding_events = events(is_encoding,:);
        encoding_conditions = encoding_events.stimcode == encoding_stimulation_code; %%% Positive t-statistics would indicate greater power for stimulation condition
        [encoding_power_table,encoding_band_table,error_flag,error_message] = get_power_analysis_tables('encoding',encoding_events,encoding_conditions,channel_number,reference_channels,electrode_info);
        if ~error_flag
            encoding_data{idx} = encoding_power_table;
            encoding_band_data{idx} = encoding_band_table;
        else
            encoding_analysis_exclusions(idx) = true;
            error_file = fullfile(error_directory,sprintf('%s_%s_%s_%d_encoding_power_error.txt',subject,task,session,channel_number));
            fid = fopen(error_file,'w');
            fprintf(fid,error_message); fclose(fid);
        end 
    end
    
    %%% Only sessions that had both retrieval stimulation and control included in retrieval analysis
    if has_retrieval_stimulation && has_no_stimulation
        is_retrieval = strcmp(events.type,retrieval_type) & ~strcmp(events.item,'<>'); %%%excluding vocalizations
        retrieval_events = events(is_retrieval,:);
        retrieval_conditions = retrieval_events.stimcode == retrieval_stimulation_code; %%% Positive t-statistics would indicate greater power for stimulation condition
        [retrieval_power_table,retrieval_band_table,error_flag,error_message] = get_power_analysis_tables('retrieval',retrieval_events,retrieval_conditions,channel_number,reference_channels,electrode_info);
        if ~error_flag
            retrieval_data{idx} = retrieval_power_table;
            retrieval_band_data{idx} = retrieval_band_table;
        else
            retrieval_analysis_exclusions(idx) = true;
            error_file = fullfile(error_directory,sprintf('%s_%s_%s_%d_retrieval_power_error.txt',subject,task,session,channel_number));
            fid = fopen(error_file,'w');
            fprintf(fid,error_message); fclose(fid);
        end 
    end
    
    %%% All sessions with encoding and retrieval stimulation included in
    %%% encoding vs retrieval analysis
    if has_encoding_stimulation && has_retrieval_stimulation
        is_encoding_or_retrieval = ismember(events.type,{encoding_type,retrieval_type}) & ~strcmp(events.item,'<>');
        all_events = events(is_encoding_or_retrieval,:);
        is_retrieval = strcmp(all_events.type,retrieval_type);
        all_events_conditions = (~is_retrieval & all_events.stimcode == encoding_stimulation_code) | (is_retrieval & all_events.stimcode == retrieval_stimulation_code); %%% Positive t-statistics would indicate greater power for stimulation condition
        all_events_types = is_retrieval; %%% Positive t-statistic would indicate a greater effect for retrieval type
        [encoding_vs_retrieval_table,error_flag,error_message] = get_contrast_table(all_events,all_events_conditions,all_events_types,channel_number,reference_channels,electrode_info);
        if ~error_flag
            encoding_vs_retrieval_data{idx} = encoding_vs_retrieval_table;
        else
            encoding_vs_retrieval_analysis_exclusions(idx) = true;
            error_file = fullfile(error_directory,sprintf('%s_%s_%s_%d_encoding_vs_retrieval_power_error.txt',subject,task,session,channel_number));
            fid = fopen(error_file,'w');
            fprintf(fid,error_message); fclose(fid);
        end 
    end
    
    %%% For all sessions get power spectra and time frequency matrices for
    %%% fooof and BOSC analysis
    is_encoding = strcmp(events.type,encoding_type);
    is_retrieval = strcmp(events.type,retrieval_type) & ~strcmp(events.item,'<>'); %%%excluding vocalizations

    encoding_events = events(is_encoding,:);
    encoding_conditions = encoding_events.stimcode == encoding_stimulation_code;
    [error_flag,error_message] = get_fooof_and_bosc_data('encoding',encoding_events,encoding_conditions,electrode_directory,channel_number,reference_channels);
    if ~error_flag
        encoding_fooof_directories_hippocampus{idx} = electrode_directory;
    else
        error_file = fullfile(error_directory,sprintf('%s_%s_%s_%d_encoding_pre-fBOSC_error.txt',subject,task,session,channel_number));
        fid = fopen(error_file,'w');
        fprintf(fid,error_message); fclose(fid);
    end
    
    retrieval_events = events(is_retrieval,:);
    retrieval_conditions = retrieval_events.stimcode == retrieval_stimulation_code;
    [error_flag,error_message] = get_fooof_and_bosc_data('retrieval',retrieval_events,retrieval_conditions,electrode_directory,channel_number,reference_channels);
    if ~error_flag
        retrieval_fooof_directories_hippocampus{idx} = electrode_directory;
    else
        fooof_bosc_analysis_exclusions(idx) = true;
        error_file = fullfile(error_directory,sprintf('%s_%s_%s_%d_retrieval_pre-fBOSC_error.txt',subject,task,session,channel_number));
        fid = fopen(error_file,'w');
        fprintf(fid,error_message); fclose(fid);
    end
end

%%% Initialize cell arrays to identify electrodes of bad subjects and to
%%% hold data of good subjects
bad_electrode_cells = false(n_electrodes,1);
bad_subject_cells = false(n_subjects,1);
averaged_encoding_data = cell(n_subjects,1);
averaged_retrieval_data = cell(n_subjects,1);

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
        subject_encoding_data = encoding_band_data(subject_electrodes);
        empty_encoding_cells = cellfun(@isempty,subject_encoding_data);
        
        %%% Exclude subjects and respective electrodes if more than half of the electrodes were excluded
        if sum(empty_encoding_cells) > floor(n_subject_electrodes/2)
            bad_electrode_cells(subject_electrodes) = true;
            bad_subject_cells(idx) = true;
        else
            subject_encoding_data(empty_encoding_cells) = [];
            averaged_encoding_data{idx} = get_subject_average_data(subject_encoding_data,subject_info);
        end        
    end
    
    if has_retrieval_stimulation && has_no_stimulation
        %%% Get retrieval data of subject's respective electrodes and cound
        %%% empty cells (excluded electrodes because of error or more than
        %%% half of the events were thrown out).
        subject_retrieval_data = retrieval_band_data(subject_electrodes);
        empty_retrieval_cells = cellfun(@isempty,subject_retrieval_data);
        
        %%% Exclude subjects and respective electrodes if more than half of the electrodes were excluded
        if sum(empty_retrieval_cells) > floor(n_subject_electrodes/2)
            bad_electrode_cells(subject_electrodes) = true;
            bad_subject_cells(idx) = true;
        else
            subject_retrieval_data(empty_retrieval_cells) = [];
            averaged_retrieval_data{idx} = get_subject_average_data(subject_retrieval_data,subject_info);
        end      
    end
    
end

%%% Delete bad electrode cells and bad subject cells from data
encoding_data(bad_electrode_cells) = [];
retrieval_data(bad_electrode_cells) = [];
encoding_band_data(bad_electrode_cells) = [];
retrieval_band_data(bad_electrode_cells) = [];
encoding_vs_retrieval_data(bad_electrode_cells) = [];
encoding_fooof_directories_hippocampus(bad_electrode_cells) = [];
retrieval_fooof_directories_hippocampus(bad_electrode_cells) = [];
averaged_encoding_data(bad_subject_cells) = [];
averaged_retrieval_data(bad_subject_cells) = [];

%%% Concatenate all data to make tables ready for analysis
encoding_analysis_table = vertcat(encoding_data{:}); clear encoding_data
encoding_band_table = vertcat(encoding_band_data{:}); clear encoding_band_data
retrieval_analysis_table = vertcat(retrieval_data{:}); clear retrieval_data
retrieval_band_table = vertcat(retrieval_band_data{:}); clear retrieval_band_data
encoding_vs_retrieval_analysis_table = vertcat(encoding_vs_retrieval_data{:}); clear encoding_vs_retrieval_data

encoding_subject_average_table = vertcat(averaged_encoding_data{:});
retrieval_subject_average_table = vertcat(averaged_retrieval_data{:});

%%% Save all tables
save(fullfile(table_directory,'encoding_analysis_table.mat'),'encoding_analysis_table','-v7.3');
save(fullfile(table_directory,'encoding_band_table.mat'),'encoding_band_table','-v7.3');
save(fullfile(table_directory,'retrieval_analysis_table.mat'),'retrieval_analysis_table','-v7.3');
save(fullfile(table_directory,'retrieval_band_table.mat'),'retrieval_band_table','-v7.3');
save(fullfile(table_directory,'encoding_vs_retrieval_analysis_table.mat'),'encoding_vs_retrieval_analysis_table','-v7.3');
save(fullfile(table_directory,'encoding_subject_average_table.mat'),'encoding_subject_average_table','-v7.3');
save(fullfile(table_directory,'retrieval_subject_average_table.mat'),'retrieval_subject_average_table','-v7.3');

%%% Print non-empty fooof data directories to text file in list directory
empty_cells = cellfun(@isempty,encoding_fooof_directories_hippocampus);
encoding_fooof_directories_hippocampus(empty_cells) = [];
writecell(encoding_fooof_directories_hippocampus,fullfile(list_directory,'encoding_fooof_directories_hippocampus.txt'))

empty_cells = cellfun(@isempty,retrieval_fooof_directories_hippocampus);
retrieval_fooof_directories_hippocampus(empty_cells) = [];
writecell(retrieval_fooof_directories_hippocampus,fullfile(list_directory,'retrieval_fooof_directories_hippocampus.txt'))

%%% Save exclusions to list_directory
save(fullfile(list_directory,'power_exclusions.mat'),'bad_electrode_cells','bad_subject_cells',...
    'encoding_analysis_exclusions','retrieval_analysis_exclusions',...
    'fooof_bosc_analysis_exclusions','encoding_vs_retrieval_analysis_exclusions');
end

function events = load_events_copy(events_file)
load(events_file,'events');
end

function [analysis_table,analysis_band_table,error_flag,error_message] = get_power_analysis_tables(event_type,events,active_stimulation,channel_number,reference_channels,electrode_info)
%%%Initialize variables to return;
error_flag = false;
error_message = [];
analysis_table = [];
analysis_band_table=[];

try
    %%% Get power normalized to the average baseline of no stimulation condition
    [power,frequencies,good_indices] = power_WM_reference(event_type,events,active_stimulation,channel_number,reference_channels);
    
    if ~isempty(power)
        
        %%% Define conditions and frequency bands
        active_stimulation = active_stimulation(good_indices);
        slow_theta_frequencies = frequencies >= 2 & frequencies <= 4;
        fast_theta_frequencies = frequencies >  4 & frequencies <= 8;
        
        %%% Get two-tailed two-sample t-test between stim and no stim
        stimulation_power = power(active_stimulation,:,:);
        control_power = power(~active_stimulation,:,:);
        clear power
        [~,~,~,statistics] = ttest2(stimulation_power,control_power,'dim',1,'tail','both','vartype','unequal');
        tstats = squeeze(statistics.tstat);
        
        %%% Get mean t-statistic within band for slow and fast theta
        slow_theta_tstats = mean(tstats(slow_theta_frequencies,:),1);
        fast_theta_tstats = mean(tstats(fast_theta_frequencies,:),1);
        
        %%% Get median across events of mean within band for slow and fast theta power
        slow_theta_stimulation_power = squeeze(median(mean(stimulation_power(:,slow_theta_frequencies,:),2),1))';
        slow_theta_control_power = squeeze(median(mean(control_power(:,slow_theta_frequencies,:),2),1))';
        fast_theta_stimulation_power = squeeze(median(mean(stimulation_power(:,fast_theta_frequencies,:),2),1))';
        fast_theta_control_power = squeeze(median(mean(control_power(:,fast_theta_frequencies,:),2),1))';
        
        %%% Get median stimulation power across events for all frequencies and
        %%% times and turn into column vector to make table
        stimulation_power = squeeze(median(stimulation_power,1));
        control_power = squeeze(median(control_power,1));
        [n_frequencies,n_samples] = size(control_power);
        power = [stimulation_power(:);control_power(:)]; clear stimulation_power control_power

        %%% Convert conditions, frequencies and time to logical, int8, and
        %%% int16 columns; respectively.
        
        condition_stimulation = true(n_frequencies*n_samples,1);
        condition_control = false(n_frequencies*n_samples,1);
        condition = [condition_stimulation;condition_control];
        
        frequencies = int8(1:n_frequencies)';
        frequency = repmat(frequencies,1,n_samples);
        frequency = repmat(frequency(:),2,1);
  
        times = int16(1:n_samples);
        time = repmat(times,n_frequencies,1);
        time = repmat(time(:),2,1);
        
        
        %%% Multiply by 1000 to conserve one decimal of percent change (or 3 decimals of tstats)
        %%% and convert to int16 (so possible range of -3276.7% to 3276.7% change in power
        %%% and .1% least significant bit of change in power) (for tstats possible range of
        %%% -32.767 to 32.767 and least significant bit of .001)
        
        %%% Repeat electrode info n_frequencies*n_samples number of rows to make power table.
        analysis_table = repmat(electrode_info,n_frequencies*n_samples*2,1);
        analysis_table.power = int16(power*1000);
        analysis_table.condition = condition;
        analysis_table.frequency = frequency;
        analysis_table.time = time;
        
        %%% Add information about band to analysis band table
        analysis_band_table = electrode_info;
        analysis_band_table.slow_theta_stimulation_power = {int16(slow_theta_stimulation_power*1000)};
        analysis_band_table.slow_theta_control_power = {int16(slow_theta_control_power*1000)};
        analysis_band_table.fast_theta_stimulation_power = {int16(fast_theta_stimulation_power*1000)};
        analysis_band_table.fast_theta_control_power = {int16(fast_theta_control_power*1000)};
        analysis_band_table.slow_theta_tstats = {int16(slow_theta_tstats*1000)};
        analysis_band_table.fast_theta_tstats = {int16(fast_theta_tstats*1000)};
        
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

function [contrast_table,error_flag,error_message] = get_contrast_table(events,conditions,is_retrieval,channel_number,reference_channels,electrode_info)
%%%Initialize variables to return;
error_flag = false;
error_message = [];
contrast_table=[];

try
    %%% Get power normalized to the average baseline of no stimulation condition
    encoding_events = events(~is_retrieval,:);
    encoding_conditions = conditions(~is_retrieval);
    [encoding_power,~,encoding_good_indices] = power_WM_reference('encoding',encoding_events,encoding_conditions,channel_number,reference_channels);
    
    retrieval_events = events(is_retrieval,:);
    retrieval_conditions = conditions(is_retrieval);
    [retrieval_power,~,retrieval_good_indices] = power_WM_reference('retrieval',retrieval_events,retrieval_conditions,channel_number,reference_channels);

    if ~isempty(encoding_power) && ~isempty(retrieval_power)
        %%% Filter conditions for good indices
        encoding_conditions = encoding_conditions(encoding_good_indices);
        retrieval_conditions = retrieval_conditions(retrieval_good_indices);
        
        %%% Get mean power across time for all events and frequencies and
        %%% turn into column vector to make table
        encoding_power = squeeze(mean(encoding_power,3));
        retrieval_power = squeeze(mean(retrieval_power,3));
        
        [n_encoding_events,n_frequencies] = size(encoding_power);
        [n_retrieval_events,~] = size(retrieval_power);
        n_events = n_encoding_events + n_retrieval_events;
       
        power = int16([encoding_power(:);retrieval_power(:)]*1000); clear encoding_power retrieval_power

        %%% Convert conditions, frequencies and event_type int8, and
        %%% int16 columns; respectively.
        encoding_conditions = repmat(encoding_conditions,1,n_frequencies);
        retrieval_conditions = repmat(retrieval_conditions,1,n_frequencies);
        condition = [encoding_conditions(:);retrieval_conditions(:)];
        
        frequencies = int8(1:n_frequencies)';
        encoding_frequencies = repmat(frequencies,n_encoding_events,1);
        retrieval_frequencies = repmat(frequencies,n_retrieval_events,1);
        frequency = [encoding_frequencies;retrieval_frequencies];
        
        encoding_type = false(n_encoding_events*n_frequencies,1);
        retrieval_type = true(n_retrieval_events*n_frequencies,1);
        event_type = [encoding_type;retrieval_type];
        
        %%% Multiply by 1000 to conserve one decimal of percent change (or 3 decimals of tstats)
        %%% and convert to int16 (so possible range of -3276.7% to 3276.7% change in power
        %%% and .1% least significant bit of change in power)
        
        %%% Repeat electrode info n_frequencies*n_samples number of rows to make power table.
        contrast_table = repmat(electrode_info,n_events*n_frequencies,1);
        contrast_table.power = power;
        contrast_table.condition = condition;
        contrast_table.frequency = frequency;
        contrast_table.event_type = event_type;        
    else
        %%% Error related to events being thrown out. Use data to describe.
        if isempty(encoding_good_indices) && isempty(retrieval_good_indices)
            error('All events were thrown out');
        else
            if isempty(encoding_good_indices)
                percent_encoding_events = 1;
                percent_encoding_stimulation = 0;
                percent_encoding_control = 0;
            else
                remaining_encoding_conditions = encoding_conditions(encoding_good_indices);
                percent_encoding_events = length(remaining_encoding_conditions)/length(encoding_conditions);
                percent_encoding_stimulation = sum(remaining_encoding_conditions)/sum(encoding_conditions);
                percent_encoding_control = sum(~remaining_encoding_conditions)/sum(~encoding_conditions);
            end
            
            if isempty(retrieval_good_indices)
                percent_retrieval_events = 1;
                percent_retrieval_stimulation = 0;
                percent_retrieval_control = 0;
            else
                remaining_retrieval_conditions = retrieval_conditions(retrieval_good_indices);
                percent_retrieval_events = length(remaining_retrieval_conditions)/length(retrieval_conditions);
                percent_retrieval_stimulation = sum(remaining_retrieval_conditions)/sum(retrieval_conditions);
                percent_retrieval_control = sum(~remaining_retrieval_conditions)/sum(~retrieval_conditions);
            end
            error(['%.2f encoding remaining. %.2f encoding stimulation remaining, %.2f encoding control remaining.\n',...
                '%.2f retrieval remaining. %.2f retrieval stimulation remaining, %.2f retrieval control remaining.\n '],...
                percent_encoding_events,percent_encoding_stimulation,percent_encoding_control,...
                percent_retrieval_events,percent_retrieval_stimulation,percent_retrieval_control);
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
        slow_theta_stimulation_power = {mean(vertcat(these_electrodes.slow_theta_stimulation_power{:}),1)};
        slow_theta_control_power = {mean(vertcat(these_electrodes.slow_theta_control_power{:}),1)};
        slow_theta_tstats = {mean(vertcat(these_electrodes.slow_theta_tstats{:}),1)};
        fast_theta_stimulation_power = {mean(vertcat(these_electrodes.fast_theta_stimulation_power{:}),1)};
        fast_theta_control_power = {mean(vertcat(these_electrodes.fast_theta_control_power{:}),1)};
        fast_theta_tstats = {mean(vertcat(these_electrodes.fast_theta_tstats{:}),1)};
        averaged_data{idx} = table(subject_ID,stimulation_lateral,stimulation_left,anterior,left,...
            slow_theta_stimulation_power,slow_theta_control_power,slow_theta_tstats,...
            fast_theta_stimulation_power,fast_theta_control_power,fast_theta_tstats);
    end        
end

averaged_band_table = vertcat(averaged_data{:});

end

function [power,frequencies,good_indices] = power_WM_reference(event_type,events,active_stimulation,channel_number,reference_channels)
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
%%% and stimulation noise (100Hz) and harmonics.
filter_order = 3;
highpass_frequency = 0.5;
notch_frequencies = [[(55:60:475)',(65:60:485)'];[(95:100:395)',(105:100:405)']];
notch_frequencies(11,:) = []; %%% To avoid repeat notch filtering

%%% Parameters for wavelet convolution
frequencies = 2.^((8:64)/8); %%% Logarithmically spaced powers of 2 between 2 and 128

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
bad_eeg_offset = [events.eegoffset] + offset - buffer < 0;
good_indices(bad_eeg_offset) = [];
events(bad_eeg_offset) = [];
active_stimulation(bad_eeg_offset) = [];    

%%% Remove any offset and drift from signals with highpass over 0.5Hz prior to referencing
[highpass_b,highpass_a] = butter(filter_order,highpass_frequency/nyquist,'high');

%%% Get hippocampal eeg
hippocampal_eeg = gete_ms(channel_number,events,duration+(buffer*2),offset-buffer,0);
hippocampal_eeg = filtfilt(highpass_b,highpass_a,hippocampal_eeg')';

%%% Gather eegs for reference channels
reference_eegs = NaN(length(events),duration+(buffer*2),length(reference_channels));
for idx = 1:length(reference_channels)
    temp_reference_eeg = gete_ms(reference_channels(idx),events,duration+buffer*2,offset,0);
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

if length(good_indices) < ceil(n_events/2) || sum(active_stimulation) < ceil(n_stimulation/2) || sum(~active_stimulation) < ceil(n_control/2)
    power = [];
    return; %%% No need to continue since more than half of the events have been thrown out;
end

root_median_square = sqrt(sum((hippocampal_eeg - median(hippocampal_eeg,2)).^2,2));
artifactual = root_median_square >= prctile(root_median_square,95);

good_indices(artifactual) = [];
hippocampal_eeg(artifactual,:) = [];
active_stimulation(artifactual,:) = [];

if length(good_indices) < ceil(n_events/2) || sum(active_stimulation) < ceil(n_stimulation/2) || sum(~active_stimulation) < ceil(n_control/2)
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

function [error_flag,error_message] = get_fooof_and_bosc_data(event_type,events,conditions,saving_directory,channel_number,reference_channels)
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
            event_indices = [1+padding,1600+padding];
        case 'retrieval'
            duration = 1000; %%% time for baseline period (500ms) + 1000 ms prior to onset of verbalization of item in recall
            offset = -1000; %%% time prior to onset of verbalization
            event_indices = [1+padding,1000+padding];
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
    n_stimulation = sum(conditions);
    n_control = sum(~conditions);
    
    %%%% Make array of indices to keep track of which events remained after
    %%%% processing
    good_indices = 1:n_events;
    bad_eeg_offset = [events.eegoffset] + offset - padding < 0;
    good_indices(bad_eeg_offset) = [];
    events(bad_eeg_offset) = [];
    conditions(bad_eeg_offset) = [];    
    
    %%% Remove any offset and drift from signals with highpass over 0.5Hz prior to referencing
    [highpass_b,highpass_a] = butter(filter_order,highpass_frequency/nyquist,'high');
    
    %%% Get hippocampal eeg
    hippocampal_eeg = gete_ms(channel_number,events,duration+(padding*2),offset-padding,0);
    hippocampal_eeg = filtfilt(highpass_b,highpass_a,hippocampal_eeg')';
    
    %%% Gather eegs for reference channels
    reference_eegs = NaN(length(events),duration+(padding*2),length(reference_channels));
    for idx = 1:length(reference_channels)
        temp_reference_eeg = gete_ms(reference_channels(idx),events,duration+(padding*2),offset-padding,0);
        temp_reference_eeg = filtfilt(highpass_b,highpass_a,temp_reference_eeg')';
        reference_eegs(:,:,idx) = temp_reference_eeg; clear temp_reference_eeg;
    end
    
    %%% Subtract mean of reference channels to get reference hippocampal eeg
    hippocampal_eeg = hippocampal_eeg - squeeze(mean(reference_eegs,3));
    
    %%% Identify events that exceed kurtosis threshold
    
    kurtosis_values = kurtosis(hippocampal_eeg(:,padding+1:end-padding)');
    nan_values = isnan(kurtosis_values);
    if all(nan_values)
        error('All kurtosis values NaN.');
    end
    
    bad_kurtosis_indices = isnan(kurtosis_values) | kurtosis_values > kurtosis_threshold;
    good_indices(bad_kurtosis_indices) = [];
    hippocampal_eeg(bad_kurtosis_indices,:) = [];
    conditions(bad_kurtosis_indices,:) = [];
    
    if length(good_indices) < ceil(n_events/2) || sum(conditions) < ceil(n_stimulation/2) || sum(~conditions) < ceil(n_control/2)
        error('More than half of the events were thrown out.');
    end
    
    root_median_square = sqrt(sum((hippocampal_eeg - median(hippocampal_eeg,2)).^2,2));
    artifactual = root_median_square >= prctile(root_median_square,95);
    
    good_indices(artifactual) = [];
    hippocampal_eeg(artifactual,:) = [];
    conditions(artifactual,:) = [];
    
    if length(good_indices) < ceil(n_events/2) || sum(conditions) < ceil(n_stimulation/2) || sum(~conditions) < ceil(n_control/2)
        error('More than half of the events were thrown out.');
    end
    
    %%% Get the raw power
    [~,raw_power] = multiphasevec3(frequencies,hippocampal_eeg,sampling_rate,width);
    
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
    save(fullfile(fooof_directory,'data.mat'),'power_spectra','frequencies','conditions','-v7.3');

    %%% Save data for fBOSC
    bosc_directory = fullfile(saving_directory,sprintf('%s_bosc',event_type));
    if ~isfolder(bosc_directory)
        mkdir(bosc_directory);
    end
    save(fullfile(bosc_directory,'data.mat'),'mean_power_spectrum','frequencies','conditions','good_indices','padding','event_indices','-v7.3');
    
catch this_error %%% Error could also result from power extraction
    error_flag = true;
    error_message = getReport(this_error, 'extended', 'hyperlinks', 'off');
end
end
