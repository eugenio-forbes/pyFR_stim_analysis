function n17_fooof_bosc_stim_site(varargin)
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
table_directory = fullfile(root_directory,'tables');
plots_directory = fullfile(root_directory,'plots');

%%% Load subject list to associate fooof and bosc data with stimulation site information
load(fullfile(list_directory,'subject_list.mat'),'subject_list');
has_stim_free_session = subject_list.has_stim_free_session;
subject_list = subject_list(has_stim_free_session,:);

%%% Declare regions of interest and periods of interest
periods = {'encoding';'retrieval'};

%%%Initialize parpool
pool_object = gcp('nocreate');
if isempty(pool_object)
    parpool(parpool_n)
end

for idx = 1:length(periods)
    period = periods{idx};
    
    electrode_directories_file = fullfile(list_directory,sprintf('%s_fooof_directories_stim_site.txt',period));
    electrode_directories = readcell(electrode_directories_file,'Delimiter',','); %%%Need to use a delimeter not present in directories so that they are not split into columns
    
    directory_info = readcell(electrode_directories_file,'Delimiter','/'); %%%This will split directories into columns to match info.
    electrode_info = match_directories_stim_site(directory_info,subject_list);
    
    [fooof_table,frequencies] = make_fooof_table(electrode_info,electrode_directories,period);
    save(fullfile(table_directory,sprintf('%s_fooof_table_stim_site',period)),'fooof_table');
    
    bosc_table = make_bosc_table(electrode_info,electrode_directories,period);
    save(fullfile(table_directory,sprintf('%s_bosc_table_stim_site',period)),'bosc_table');

    make_fooof_histogram_plots(plots_directory,fooof_table,period,frequencies);
    make_bosc_p_episode_plots(plots_directory,bosc_table,period,frequencies);
    
    clean_up(electrode_directories,period);
end

end
function electrode_info = match_directories_stim_site(directory_info,electrode_list)
%%% Initialize variable for extracting electrode info that matches
%%% order in listed directories
variables_to_keep = {'subject_ID','stimulation_group','stimulation_hemisphere',...
    'stimulation_electrode_channel','stimulation_references'};

%%% Get electrode info
subjects = electrode_list.subject;
indices = cellfun(@(x) find(strcmp(subjects,x)),directory_info(:,end-2));
electrode_info = electrode_list(indices,variables_to_keep);

end
function [fooof_table,frequencies] = make_fooof_table(electrode_info,electrode_directories,period)
n_electrodes = length(electrode_directories);

%%% Initialize cell array to hold histograms for every electrode
electrode_histograms = cell(n_electrodes,1);

%%% Load first directory data to retrieve frequencies that will be used to
%%% make histograms
folder = sprintf('%s_fooof',period);
sample_data_file = fullfile(electrode_directories{1},folder,'data.mat');
load(sample_data_file,'frequencies')

%%% Bin edges of hiscounts will be frequencies
edges = frequencies;

parfor idx = 1:n_electrodes
    this_directory = fullfile(electrode_directories{idx},folder);
    fooof_results = fooof_par_load(this_directory);
    trial_histcounts = histcounts(vertcat(fooof_results.center_frequencies{:}),'BinEdges',edges);
    electrode_histograms{idx} = trial_histcounts/height(fooof_results);
end

fooof_table = electrode_info;
fooof_table.histograms = electrode_histograms;
end
function fooof_results = fooof_par_load(data_directory)

fooof_results_file = fullfile(data_directory,'fooof_results.json');
undecoded_data = importdata(fooof_results_file);

n_trials = length(undecoded_data);
r_squared = NaN(n_trials,1);
center_frequencies = cell(n_trials,1);
power = cell(n_trials,1);

for idx = 1:n_trials
    temp_results = jsondecode(undecoded_data{idx});
    peak_params = temp_results.peak_params_;
    r_squared(idx) = temp_results.r_squared_;
    if ~isempty(peak_params)
        center_frequencies{idx} = peak_params(:,1);
        power{idx} = peak_params(:,2);
    end
end

fooof_results = table;
fooof_results.r_squared = r_squared;
fooof_results.center_frequencies = center_frequencies;
fooof_results.power = power;

end

function bosc_table = make_bosc_table(electrode_info,electrode_directories,period)
n_electrodes = length(electrode_directories);

%%% Initialize cell array to hold p-episode values for every frequency
%%% for every electrode
p_episodes = cell(n_electrodes,1);

parfor idx = 1:n_electrodes
    this_directory = electrode_directories{idx};
    this_info = electrode_info(idx,:);

    temp_p_episodes = bosc_par_load(this_directory,this_info,period);
    p_episodes{idx} = mean(vertcat(temp_p_episodes{:}),1);
end

bosc_table = electrode_info;
bosc_table.p_episodes = p_episodes;

end

function p_episodes = bosc_par_load(electrode_directory,electrode_info,period)
sampling_rate = 1000;
channel_number = electrode_info.stimulation_electrode_channel;
references = electrode_info.stimulation_references{:};

data_directory = fullfile(electrode_directory,sprintf('%s_bosc',period));

data_file = fullfile(data_directory,'data.mat');
load(data_file,'frequencies')
load(data_file,'event_indices')
load(data_file,'good_indices')
load(data_file,'padding')

wavelet_width = 6;
threshold_percentile = 0.95;
threshold_n_cycles = 3;
threshold_n_samples = (threshold_n_cycles./frequencies)*sampling_rate;

events_file = fullfile(electrode_directory,'events.mat');
load(events_file,'events');

switch period
    case 'encoding'
        encoding_type = strcmp(events.type,'WORD');
        events = events(encoding_type,:);
    case 'retrieval'
        retrieval_type = strcmp(events.type,'REC_WORD') & ~strcmp(events.item,'<>');
        events = events(retrieval_type,:);
end

events = events(good_indices,:);

fooof_results_file = fullfile(data_directory,'fooof_results.json');
undecoded_data = importdata(fooof_results_file);
fooof_results = jsondecode(undecoded_data{:});
aperiodic_params = fooof_results.aperiodic_params_;
offset = aperiodic_params(1);
knee = aperiodic_params(2);
exponent = aperiodic_params(3); clear aperiodic_params

log_power = offset - log10(abs(knee) + frequencies.^(exponent));
aperiodic_fit_power = 10.^(log_power);
threshold_power = chi2inv(threshold_percentile,2)*aperiodic_fit_power/2;

pad1 = ceil(padding/2);
pad2 = floor(padding/2);
raw_power = get_bosc_data(period,events,event_indices,padding,channel_number,references,frequencies,wavelet_width);
raw_power = raw_power(:,:,pad1+1:end-pad1);
[n_events,n_frequencies,n_samples] = size(raw_power);
p_episodes = cell(n_events,1);

for idx = 1:n_events
    event_power = squeeze(raw_power(idx,:,:));
    oscillations_detected = zeros(n_frequencies,n_samples);
    for fdx = 1:length(frequencies)
        oscillations_detected(fdx,:) = BOSC_detect(event_power(fdx,:),threshold_power(fdx),threshold_n_samples(fdx),sampling_rate);
    end
    oscillations_detected = oscillations_detected(:,pad2+1:end-pad2);
    p_episodes{idx} = mean(oscillations_detected,1);
end

end

function raw_power = get_bosc_data(event_type,events,event_indices,padding,channel_number,references,frequencies,wavelet_width)

%%% Signal parameters
%%% signal acquired at 1000 Hz, not resampling
sampling_rate = 1000;
nyquist = sampling_rate/2;

switch event_type
    case 'encoding'
        offset = 0; %%% in ms, prior to word display to acquire singal
    case 'retrieval'
        offset = -1000; %%% time prior to onset of verbalization
end

duration = event_indices(end);

%%% Filtering parameters to highpass (0.5Hz)
%%% Only high pass filtering to perform rereferncing. Avoiding other
%%% filters to minimize distortions of the signal and the spectrum.
%%% Aiming to reduce noise with rereferencing. Only fitting frequencies
%%% from 2 to
filter_order = 3;
highpass_frequency = 0.5;

%%% Convert events to struct if needed
if istable(events)
    events = table2struct(events);
end

%%% Remove any offset and drift from signals with highpass over 0.5Hz prior to referencing
[highpass_b,highpass_a] = butter(filter_order,highpass_frequency/nyquist,'high');

%%% Get hippocampal eeg
stim_site_eeg = gete_ms(channel_number,events,duration+(padding*2),offset-padding,0);
stim_site_eeg = filtfilt(highpass_b,highpass_a,stim_site_eeg')';

%%% Gather eegs for reference channels
reference_eegs = NaN(length(events),duration+(padding*2),length(references));
for idx = 1:length(references)
    temp_reference_eeg = gete_ms(references(idx),events,duration+(padding*2),offset-padding,0);
    temp_reference_eeg = filtfilt(highpass_b,highpass_a,temp_reference_eeg')';
    reference_eegs(:,:,idx) = temp_reference_eeg; clear temp_reference_eeg;
end

%%% Subtract mean of reference channels to get reference hippocampal eeg
stim_site_eeg = stim_site_eeg - squeeze(mean(reference_eegs,3));

%%% Get the raw power
[~,raw_power] = multiphasevec3(frequencies,stim_site_eeg,sampling_rate,wavelet_width);

end

function make_fooof_histogram_plots(plots_directory,fooof_table,period,frequencies) 

%%% Plot parameteres
figure_width = 250;
figure_height = 150;
bar_width = min(diff(frequencies));
bar_frequencies = frequencies(1:end-1) + bar_width/2;
n_bars = length(bar_frequencies);
only_color = [205 234 192]/255;
x_limits = [frequencies(1),frequencies(end)];
tick_frequencies = mod(frequencies,4) == 0;
x_ticks = unique([frequencies(1),frequencies(tick_frequencies)]);
x_tick_labels = arrayfun(@(x) num2str(x),x_ticks,'UniformOutput',false);

%%% Make a unique plot directory for the histogram plots of this period
this_plot_directory = fullfile(plots_directory,sprintf('%s_fooof_histograms_stim_site',period));
if ~isfolder(this_plot_directory)
    mkdir(this_plot_directory);
end

%%% Make combinations of anatomical regions and stimulation groups to
%%% acquire histograms plots for every combination
anatomical_regions = {'left','right','both'};
stimulation_groups = {'IPL','retrosplenial'};
[Ax,Bx] = ndgrid(1:numel(anatomical_regions),1:numel(stimulation_groups));
n_combos = length(Ax(:));
anatomical_regions = anatomical_regions(Ax(:));
stimulation_groups = stimulation_groups(Bx(:));

%%% Get information about anatomic location to filter table results for
%%% every combination
is_left = strcmp(fooof_table.stimulation_hemisphere,'left');

for idx = 1:n_combos
    anatomical_region = anatomical_regions{idx};
    stimulation_group = stimulation_groups{idx};
    this_plot_file_name = fullfile(this_plot_directory,sprintf('%s_%s_stim_site',stimulation_group,anatomical_region));
    
    switch anatomical_region
        case 'left'
            temp_table = fooof_table(is_left,:);
        case 'right'
            temp_table = fooof_table(~is_left,:);
        case 'both'
            temp_table = fooof_table;
    end
    
    %%% Get electrode results specific to stimulation group
    same_group = strcmp(temp_table.stimulation_group,stimulation_group);
    temp_table = temp_table(same_group,:);
    
    if ~isempty(temp_table)
        
        %%% Get averages of electrode results for stimulation and control weighted by the number of
        %%% electrodes each subject has
        electrode_weights = arrayfun(@(x) 1/sum(temp_table.subject_ID == x),temp_table.subject_ID);
        n_subjects = length(unique(temp_table.subject_ID));
        histograms = vertcat(temp_table.histograms{:});
        plot_histogram = sum(histograms.*electrode_weights,1)/n_subjects;
        
        %%% Define y limits based on maxima of results
        y_limits = [0, max(plot_histogram)*1.1];
        y_ticks = y_limits(1):diff(y_limits)/5:y_limits(end);
        y_tick_labels = repelem({''},1,length(y_ticks));
        y_tick_labels{y_ticks == y_limits(end)} = sprintf('%.3f',y_limits(end));
        y_tick_labels{y_ticks == 0} = '0';
        
        %%% Initialize array to hold handles for bar objects
        bar_handles = cell(n_bars,1);
        
        %%% Create a new figure
        figure('Units','pixels','Visible','off','Position',[0 0 figure_width figure_height]);
        
        hold on

        %%% Bars of the histograms of the condition overlap up to lower
        %%% value, which is plotted in mixed color, one of the two
        %%% conditions has a greater value than the other (value
        %%% difference) which is stacked to overlapping bar with condition
        %%% specific color (residual_color).
        
        for fdx = 1:n_bars
            bar_handles{fdx} = bar(bar_frequencies(fdx),plot_histogram(fdx),'FaceColor','flat','EdgeColor','none','BarWidth',bar_width);
            these_bars = bar_handles{fdx};
            these_bars(1).CData = only_color;
        end
        
        ylim(y_limits);yticks(y_ticks); yticklabels(y_tick_labels);
        xlim(x_limits);xticks(x_ticks);xticklabels(x_tick_labels);
        
        hold off
        pause(1)
        
        print(this_plot_file_name,'-dsvg');
        close all
    end
    
end
end

function make_bosc_p_episode_plots(plots_directory,bosc_table,period,frequencies)

%%% Plot parameteres
figure_width = 250;
figure_height = 150;
bar_width = min(diff(frequencies));
n_bars = length(frequencies);
only_color = [205 234 192]/255;
x_limits = [frequencies(1),frequencies(end)];
tick_frequencies = mod(frequencies,4) == 0;
x_ticks = unique([frequencies(1),frequencies(tick_frequencies)]);
x_tick_labels = arrayfun(@(x) num2str(x),x_ticks,'UniformOutput',false);

%%% Make a unique plot directory for the histogram plots of this period
this_plot_directory = fullfile(plots_directory,sprintf('%s_bosc_pepisodes_stim_site',period));
if ~isfolder(this_plot_directory)
    mkdir(this_plot_directory);
end

%%% Make combinations of anatomical regions and stimulation groups to
%%% acquire histograms plots for every combination
anatomical_regions = {'left','right','both'};
stimulation_groups = {'IPL','retrosplenial'};
[Ax,Bx] = ndgrid(1:numel(anatomical_regions),1:numel(stimulation_groups));
n_combos = length(Ax(:));
anatomical_regions = anatomical_regions(Ax(:));
stimulation_groups = stimulation_groups(Bx(:));

%%% Get information about anatomic location to filter table results for
%%% every combination
is_left = strcmp(bosc_table.stimulation_hemisphere,'left');

for idx = 1:n_combos
    anatomical_region = anatomical_regions{idx};
    stimulation_group = stimulation_groups{idx};
    this_plot_file_name = fullfile(this_plot_directory,sprintf('%s_%s_stim_site',stimulation_group,anatomical_region));
    
    switch anatomical_region
        case 'left'
            temp_table = bosc_table(is_left,:);
        case 'right'
            temp_table = bosc_table(~is_left,:);
        case 'both'
            temp_table = bosc_table;
    end
    
    %%% Get electrode results specific to stimulation group
    same_group = strcmp(temp_table.stimulation_group,stimulation_group);
    temp_table = temp_table(same_group,:);
    empty_cells = cellfun(@isempty,temp_table.p_episodes);
    temp_table(empty_cells,:) = [];
    
    if ~isempty(temp_table)
        
        %%% Get averages of electrode results for stimulation and control weighted by the number of
        %%% electrodes each subject has
        electrode_weights = arrayfun(@(x) 1/sum(temp_table.subject_ID == x),temp_table.subject_ID);
        n_subjects = length(unique(temp_table.subject_ID));
        p_episodes = vertcat(temp_table.p_episodes{:});
        plot_p_episode = sum(p_episodes.*electrode_weights,1)/n_subjects;
        
        %%% Define y limits based on maxima of results
        y_limits = [0, max(plot_p_episode)*1.1];
        y_ticks = y_limits(1):diff(y_limits)/5:y_limits(end);
        y_tick_labels = repelem({''},1,length(y_ticks));
        y_tick_labels{y_ticks == y_limits(end)} = sprintf('%.3f',y_limits(end));
        y_tick_labels{y_ticks == 0} = '0';
        
        %%% Initialize array to hold handles for bar objects
        bar_handles = cell(n_bars,1);
        
        %%% Create a new figure
        figure('Units','pixels','Visible','off','Position',[0 0 figure_width figure_height]);
        
        hold on

        %%% Bars of the p episodes of the condition overlap up to lower
        %%% value, which is plotted in mixed color, one of the two
        %%% conditions has a greater value than the other (value
        %%% difference) which is stacked to overlapping bar with condition
        %%% specific color (residual_color).
        
        for fdx = 1:n_bars
            bar_handles{fdx} = bar(frequencies(fdx),p_episodes(fdx),'FaceColor','flat','EdgeColor','none','BarWidth',bar_width);
            these_bars = bar_handles{fdx};
            these_bars(1).CData = only_color;
        end
        
        ylim(y_limits);yticks(y_ticks); yticklabels(y_tick_labels);
        xlim(x_limits);xticks(x_ticks);xticklabels(x_tick_labels);
        
        hold off
        pause(1)
        
        print(this_plot_file_name,'-dsvg');
        close all
    end
    
end
end

function clean_up(electrode_directories,period)
n_directories = length(electrode_directories);
parfor idx = 1:n_directories
    electrode_directory = electrode_directories{idx};
    fooof_directory = fullfile(electrode_directory,sprintf('%s_fooof',period));
    bosc_directory = fullfile(electrode_directory,sprintf('%s_bosc',period));
    if isfolder(fooof_directory)
        rmdir(fooof_directory,'s');
    end
    if isfolder(bosc_directory)
        rmdir(bosc_directory,'s');
    end
end
end
