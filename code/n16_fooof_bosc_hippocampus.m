function n16_fooof_bosc_hippocampus(varargin)
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

%%% Load electrode list to associate fooof and bosc data with electrode information
load(fullfile(list_directory,'electrode_list.mat'),'electrode_list');

%%% Declare regions of interest and periods of interest
periods = {'encoding';'retrieval'};

%%%Initialize parpool
pool_object = gcp('nocreate');
if isempty(pool_object)
    parpool(parpool_n)
end

for idx = 1:length(periods)
    period = periods{idx};
    
    electrode_directories_file = fullfile(list_directory,sprintf('%s_fooof_directories_hippocampus.txt',period));
    electrode_directories = readcell(electrode_directories_file,'Delimiter',','); %%%Need to use a delimeter not present in directories so that they are not split into columns
    
    directory_info = readcell(electrode_directories_file,'Delimiter','/'); %%%This will split directories into columns to match info.
    electrode_info = match_directories_and_list(directory_info,electrode_list);
    
    [fooof_table,aperiodic_table,frequencies] = make_fooof_table(electrode_info,electrode_directories,period);
    save(fullfile(table_directory,sprintf('%s_fooof_table_hippocampus',period)),'fooof_table');
    save(fullfile(table_directory,sprintf('%s_aperiodic_table_hippocampus',period)),'aperiodic_table');
    
    bosc_table = make_bosc_table(electrode_info,electrode_directories,period);
    save(fullfile(table_directory,sprintf('%s_bosc_table_hippocampus',period)),'bosc_table');

    make_fooof_histogram_plots(plots_directory,fooof_table,period,frequencies);
    make_aperiodic_plot(plots_directory,aperiodic_table,period);
    make_bosc_p_episode_plots(plots_directory,bosc_table,period,frequencies);
    
    clean_up(electrode_directories,period);
end

end
function electrode_info = match_directories_and_list(directory_info,electrode_list)
%%% Initialize variable for extracting electrode info that matches
%%% order in listed directories
variables_to_keep = {'subject_ID','electrode_ID','stimulation_group',...
    'stimulation_hemisphere','longitudinal_location','hemisphere',...
    'channel_number','reference_channels'};
n_directories = size(directory_info,1);
indices = NaN(n_directories,1);

%%% Get electrode info
subjects = electrode_list.subject;
tasks = electrode_list.task;
sessions = electrode_list.session;
channel_numbers = electrode_list.channel_number;

for idx = 1:n_directories
    %%% Get directory info
    subject = directory_info(idx,end-3);
    task = directory_info(idx,end-2);
    session = directory_info(idx,end-1);
    channel_number = directory_info{idx,end};
    
    %%% Compare electrode info to directory info
    same_subject = strcmp(subjects,subject);
    same_task = strcmp(tasks,task);
    same_session = strcmp(sessions,session);
    same_channel = channel_numbers == channel_number;
   
    %Find matching index for this electrode
    indices(idx) = find(same_subject & same_task & same_session & same_channel);
end

electrode_info = electrode_list(indices,variables_to_keep);

end
function [fooof_table,aperiodic_table,frequencies] = make_fooof_table(electrode_info,electrode_directories,period)
n_electrodes = length(electrode_directories);

%%% Initialize cell array to hold histograms for every electrode
stimulation_histograms = cell(n_electrodes,1);
control_histograms = cell(n_electrodes,1);
aperiodic_data = cell(n_electrodes,1);

%%% Load first directory data to retrieve frequencies that will be used to
%%% make histograms
folder = sprintf('%s_fooof',period);
sample_data_file = fullfile(electrode_directories{1},folder,'data.mat');
load(sample_data_file,'frequencies')

%%% Bin edges of hiscounts will be frequencies
edges = frequencies;

parfor idx = 1:n_electrodes
    this_directory = fullfile(electrode_directories{idx},folder);
    this_electrode_info = electrode_info(idx,:);
    [active_stimulation,fooof_results] = fooof_par_load(this_directory);
    
    %%% Can initially reject trials based on r-squared and then reject
    %%% peaks based on power
    temp_aperiodic_table = repmat(this_electrode_info,length(active_stimulation),1);
    temp_aperiodic_table.condition = active_stimulation;
    temp_aperiodic_table.exponent = fooof_results.exponent;
    temp_aperiodic_table.r_squared = fooof_results.r_squared;
    aperiodic_data{idx} = temp_aperiodic_table;
    
    stimulation_trials = fooof_results(active_stimulation,:);
    control_trials = fooof_results(~active_stimulation,:);
    stimulation_histcounts = histcounts(vertcat(stimulation_trials.center_frequencies{:}),'BinEdges',edges);
    stimulation_histograms{idx} = stimulation_histcounts/sum(active_stimulation);
    control_histcounts = histcounts(vertcat(control_trials.center_frequencies{:}),'BinEdges',edges);
    control_histograms{idx} = control_histcounts/sum(~active_stimulation);
    
end

fooof_table = electrode_info;
fooof_table.stimulation_histograms = stimulation_histograms;
fooof_table.control_histograms = control_histograms;
aperiodic_table = vertcat(aperiodic_data{:});

end
function [conditions,fooof_results] = fooof_par_load(data_directory)

data_file = fullfile(data_directory,'data.mat');
load(data_file,'conditions')

fooof_results_file = fullfile(data_directory,'fooof_results.json');
undecoded_data = importdata(fooof_results_file);

n_trials = length(undecoded_data);
exponent = NaN(n_trials,1);
r_squared = NaN(n_trials,1);
center_frequencies = cell(n_trials,1);
power = cell(n_trials,1);

for idx = 1:n_trials
    temp_results = jsondecode(undecoded_data{idx});
    aperiodic_params = temp_results.aperiodic_params_;
    peak_params = temp_results.peak_params_;
    r_squared(idx) = temp_results.r_squared_;
    exponent(idx) = aperiodic_params(3);
    if ~isempty(peak_params)
        center_frequencies{idx} = peak_params(:,1);
        power{idx} = peak_params(:,2);
    end
end

fooof_results = table;
fooof_results.exponent = exponent;
fooof_results.r_squared = r_squared;
fooof_results.center_frequencies = center_frequencies;
fooof_results.power = power;

end

function bosc_table = make_bosc_table(electrode_info,electrode_directories,period)
n_electrodes = length(electrode_directories);

%%% Initialize cell array to hold p-episode values for every frequency
%%% for every electrode
stimulation_p_episodes = cell(n_electrodes,1);
control_p_episodes = cell(n_electrodes,1);

parfor idx = 1:n_electrodes
    this_directory = electrode_directories{idx};
    this_info = electrode_info(idx,:);

    [active_stimulation,p_episodes] = bosc_par_load(this_directory,this_info,period);
    
    stimulation_trials = p_episodes(active_stimulation,:);
    control_trials = p_episodes(~active_stimulation,:);
    stimulation_p_episodes{idx} = mean(vertcat(stimulation_trials{:}),1);
    control_p_episodes{idx} = mean(vertcat(control_trials{:}),1);
end

bosc_table = electrode_info;
bosc_table.stimulation_p_episodes = stimulation_p_episodes;
bosc_table.control_p_episodes = control_p_episodes;

end

function [conditions,p_episodes] = bosc_par_load(electrode_directory,electrode_info,period)
sampling_rate = 1000;
channel_number = electrode_info.channel_number;
references = electrode_info.reference_channels{:};

data_directory = fullfile(electrode_directory,sprintf('%s_bosc',period));

data_file = fullfile(data_directory,'data.mat');
load(data_file,'frequencies')
load(data_file,'event_indices')
load(data_file,'good_indices')
load(data_file,'conditions')
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
hippocampal_eeg = gete_ms(channel_number,events,duration+(padding*2),offset-padding,0);
hippocampal_eeg = filtfilt(highpass_b,highpass_a,hippocampal_eeg')';

%%% Gather eegs for reference channels
reference_eegs = NaN(length(events),duration+(padding*2),length(references));
for idx = 1:length(references)
    temp_reference_eeg = gete_ms(references(idx),events,duration+(padding*2),offset-padding,0);
    temp_reference_eeg = filtfilt(highpass_b,highpass_a,temp_reference_eeg')';
    reference_eegs(:,:,idx) = temp_reference_eeg; clear temp_reference_eeg;
end

%%% Subtract mean of reference channels to get reference hippocampal eeg
hippocampal_eeg = hippocampal_eeg - squeeze(mean(reference_eegs,3));

%%% Get the raw power
[~,raw_power] = multiphasevec3(frequencies,hippocampal_eeg,sampling_rate,wavelet_width);

end

function make_fooof_histogram_plots(plots_directory,fooof_table,period,frequencies) 

%%% Plot parameteres
figure_width = 250;
figure_height = 150;
bar_width = min(diff(frequencies));
bar_frequencies = frequencies(1:end-1) + bar_width/2;
n_bars = length(bar_frequencies);
mixed_color = [205 234 192]/255;
color_stimulation = [255 179 15]/255;
color_control = [39 154 241]/255;
x_limits = [frequencies(1),frequencies(end)];
tick_frequencies = mod(frequencies,4) == 0;
x_ticks = unique([frequencies(1),frequencies(tick_frequencies)]);
x_tick_labels = arrayfun(@(x) num2str(x),x_ticks,'UniformOutput',false);

%%% Make a unique plot directory for the histogram plots of this period
this_plot_directory = fullfile(plots_directory,sprintf('%s_fooof_histograms_hippocampus',period));
if ~isfolder(this_plot_directory)
    mkdir(this_plot_directory);
end

%%% Make combinations of anatomical regions and stimulation groups to
%%% acquire histograms plots for every combination
anatomical_regions = {'anterior','posterior','anterior-left','anterior-right','posterior-left','posterior-right'};
stimulation_groups = {'IPL','retrosplenial'};
[Ax,Bx] = ndgrid(1:numel(anatomical_regions),1:numel(stimulation_groups));
n_combos = length(Ax(:));
anatomical_regions = anatomical_regions(Ax(:));
stimulation_groups = stimulation_groups(Bx(:));

%%% Get information about anatomic location to filter table results for
%%% every combination
is_anterior = strcmp(fooof_table.longitudinal_location,'anterior');
is_left = strcmp(fooof_table.hemisphere,'left');

for idx = 1:n_combos
    anatomical_region = anatomical_regions{idx};
    stimulation_group = stimulation_groups{idx};
    this_plot_file_name = fullfile(this_plot_directory,sprintf('%s_%s_hippocampus',stimulation_group,anatomical_region));
    
    switch anatomical_region
        case 'anterior'
            temp_table = fooof_table(is_anterior,:);
        case 'posterior'
            temp_table = fooof_table(~is_anterior,:);
        case 'anterior-left'
            temp_table = fooof_table(is_anterior & is_left,:);
        case 'anterior-right'
            temp_table = fooof_table(is_anterior & ~is_left,:);
        case 'posterior-left'
            temp_table = fooof_table(~is_anterior & is_left,:);
        case 'posterior-right'
            temp_table = fooof_table(~is_anterior & ~is_left,:);
    end
    
    %%% Get electrode results specific to stimulation group
    same_group = strcmp(temp_table.stimulation_group,stimulation_group);
    temp_table = temp_table(same_group,:);
    
    if ~isempty(temp_table)
        
        %%% Get averages of electrode results for stimulation and control weighted by the number of
        %%% electrodes each subject has
        electrode_weights = arrayfun(@(x) 1/sum(temp_table.subject_ID == x),temp_table.subject_ID);
        n_subjects = length(unique(temp_table.subject_ID));
        stimulation_histograms = vertcat(temp_table.stimulation_histograms{:});
        control_histograms = vertcat(temp_table.control_histograms{:});
        stimulation_histogram = sum(stimulation_histograms.*electrode_weights,1)/n_subjects;
        control_histogram = sum(control_histograms.*electrode_weights,1)/n_subjects;
        
        %%% Define y limits based on maxima of results
        y_limits = [0, max([stimulation_histogram,control_histogram])*1.1];
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
            value_stimulation = stimulation_histogram(fdx);
            value_control = control_histogram(fdx);
            lower_value = min([value_stimulation;value_control]);
            value_difference = abs(value_stimulation-value_control);
            if value_stimulation > value_control
                residual_color = color_stimulation;
            else
                residual_color = color_control;
            end
            bar_handles{fdx} = bar(bar_frequencies(fdx),[lower_value value_difference],'stacked','FaceColor','flat','EdgeColor','none','BarWidth',bar_width);
            these_bars = bar_handles{fdx};
            these_bars(1).CData = mixed_color;
            these_bars(2).CData = residual_color;
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

function make_aperiodic_plot(plots_directory,aperiodic_table,period)
%%% LME formulas for stimulation effect and interaction with stimulation group
random_effect = '1 + condition';
stimulation_effect_formula = sprintf('exponent ~ condition + (%s|subject_ID) + (%s|subject_ID:electrode_ID) + (%s|subject_ID:stimulation_hemisphere)',...
    random_effect,random_effect,random_effect);
stimulation_effect_index = 2;
interaction_formula = sprintf('exponent ~ condition*stimulation_group + (%s|subject_ID) + (%s|subject_ID:electrode_ID) + (%s|subject_ID:stimulation_hemisphere) + (%s|stimulation_group:stimulation_hemisphere)',...
    random_effect,random_effect,random_effect,random_effect);
interaction_index = 4;


%%% Plot parameteres
figure_width = 250;
figure_height = 150;
color_stimulation = [255 179 15]/255;
color_control = [39 154 241]/255;
colors_boxes = [[0 0 0];[0 0 0]];
x_limits = [0.5,2.5];
x_ticks = [1,2];

this_plot_directory = fullfile(plots_directory,sprintf('%s_aperiodic_plots_hippocampus',period));
if ~isfolder(this_plot_directory)
    mkdir(this_plot_directory);
end

anatomical_regions = {'anterior','posterior','anterior-left','anterior-right','posterior-left','posterior-right'};
stimulation_groups = {'IPL','retrosplenial'};
[Ax,Bx] = ndgrid(1:numel(anatomical_regions),1:numel(stimulation_groups));
n_combos = length(Ax(:));
anatomical_region = anatomical_regions(Ax(:));
stimulation_group = stimulation_groups(Bx(:));

t_statistic = NaN(n_combos,1);
p_value = NaN(n_combos,1);
df = NaN(n_combos,1);
n_subjects = NaN(n_combos,1);
n_subjects_IPL = NaN(n_combos,1);
n_subjects_RS = NaN(n_combos,1);
interaction_t = NaN(n_combos,1);
interaction_p = NaN(n_combos,1);
interaction_df = NaN(n_combos,1);

is_anterior = strcmp(aperiodic_table.longitudinal_location,'anterior');
is_left = strcmp(aperiodic_table.hemisphere,'left');

for idx = 1:n_combos
    this_anatomical_region = anatomical_region{idx};
    this_stimulation_group = stimulation_group{idx};
    this_plot_file_name = fullfile(this_plot_directory,sprintf('%s_%s_hippocampus',this_stimulation_group,this_anatomical_region));
    
    switch this_anatomical_region
        case 'anterior'
            temp_table = aperiodic_table(is_anterior,:);
        case 'posterior'
            temp_table = aperiodic_table(~is_anterior,:);
        case 'anterior-left'
            temp_table = aperiodic_table(is_anterior & is_left,:);
        case 'anterior-right'
            temp_table = aperiodic_table(is_anterior & ~is_left,:);
        case 'posterior-left'
            temp_table = aperiodic_table(~is_anterior & is_left,:);
        case 'posterior-right'
            temp_table = aperiodic_table(~is_anterior & ~is_left,:);
    end
    
    if ~isempty(temp_table)
        
        same_group = strcmp(temp_table.stimulation_group,this_stimulation_group);
        is_IPL = strcmp(temp_table.stimulation_group,'IPL');
        n_subjects_IPL(idx) = length(unique(temp_table.subject_ID(is_IPL)));
        n_subjects_RS(idx) = length(unique(temp_table.subject_ID(~is_IPL)));
        
        temp_table.stimulation_group = is_IPL;
        temp_table.condition = logical(temp_table.condition);
        temp_table.subject_ID = categorical(temp_table.subject_ID);
        temp_table.electrode_ID = categorical(temp_table.electrode_ID);
        temp_table.stimulation_hemisphere = categorical(temp_table.stimulation_hemisphere);
        
        interaction_lme = fitlme(temp_table,interaction_formula);
        interaction_t(idx) = interaction_lme.Coefficients.tStat(interaction_index);
        interaction_p(idx) = interaction_lme.Coefficients.pValue(interaction_index);
        interaction_df(idx) = interaction_lme.Coefficients.DF(interaction_index);
        
        temp_table = temp_table(same_group,:);
        
        if ~isempty(temp_table)
            temp_table.subject_ID = double(temp_table.subject_ID);
            unique_subjects = unique(temp_table.subject_ID);
            n_subjects(idx) = length(unique_subjects);
            
            exponent_boxes = NaN(n_subjects(idx),2);
            
            for sdx = 1:n_subjects(idx)
                subject_indices = temp_table.subject_ID == unique_subjects(sdx);
                subject_electrodes = temp_table(subject_indices,:);
                active_stimulation = subject_electrodes.condition;
                stimulation_exponents = subject_electrodes.exponent(active_stimulation);
                control_exponents = subject_electrodes.exponent(~active_stimulation);
                exponent_boxes(sdx,1) = mean(stimulation_exponents);
                exponent_boxes(sdx,2) = mean(control_exponents);
            end
            
            temp_table.subject_ID = categorical(temp_table.subject_ID);
            
            stimulation_effect_lme = fitlme(temp_table,stimulation_effect_formula);
            t_statistic(idx) = stimulation_effect_lme.Coefficients.tStat(stimulation_effect_index);
            p_value(idx) = stimulation_effect_lme.Coefficients.pValue(stimulation_effect_index);
            df(idx) = stimulation_effect_lme.Coefficients.DF(stimulation_effect_index);
        
            min_exponent = min(temp_table.exponent);
%             mean_exponent = mean(temp_table.exponent);
%             std_exponent = std(temp_table.exponent);
%             max_exponent = mean_exponent + (3*std_exponent);
            max_exponent = max(temp_table.exponent);
            diff_exponent = diff([min_exponent,max_exponent]);
            y_limits = [min_exponent-(0.1*diff_exponent),max_exponent+(0.1*diff_exponent)];
            y_ticks = min_exponent:(diff_exponent/6):max_exponent;
            y_tick_labels = repelem({''},1,length(y_ticks));
            y_tick_labels{y_ticks == min_exponent} = sprintf('%.2f',min_exponent);
            y_tick_labels{y_ticks == max_exponent} = sprintf('%.2f',max_exponent);
            
            
            figure('Units','pixels','Visible','off','Position',[0 0 figure_width figure_height]);
            hold on
            scatter(ones(n_subjects(idx),1),exponent_boxes(:,1),[],color_stimulation);
            scatter(repmat(2,n_subjects(idx),1),exponent_boxes(:,2),[],color_control);
            boxplot_handle = boxplot(exponent_boxes,[1,2],'Colors',colors_boxes,'Symbol','');
            set(boxplot_handle,{'linew'},{2});
            ylim(y_limits);yticks(y_ticks);yticklabels(y_tick_labels);
            xlim(x_limits);xticks(x_ticks);xticklabels([]);
            hold off
            
            print(this_plot_file_name,'-dsvg');
            close all
        end
    end
    
end

end

function make_bosc_p_episode_plots(plots_directory,bosc_table,period,frequencies)

%%% Plot parameteres
figure_width = 250;
figure_height = 150;
bar_width = min(diff(frequencies));
n_bars = length(frequencies);
mixed_color = [205 234 192]/255;
color_stimulation = [255 179 15]/255;
color_control = [39 154 241]/255;
x_limits = [frequencies(1),frequencies(end)];
tick_frequencies = mod(frequencies,4) == 0;
x_ticks = unique([frequencies(1),frequencies(tick_frequencies)]);
x_tick_labels = arrayfun(@(x) num2str(x),x_ticks,'UniformOutput',false);

%%% Make a unique plot directory for the histogram plots of this period
this_plot_directory = fullfile(plots_directory,sprintf('%s_bosc_pepisodes_hippocampus',period));
if ~isfolder(this_plot_directory)
    mkdir(this_plot_directory);
end

%%% Make combinations of anatomical regions and stimulation groups to
%%% acquire histograms plots for every combination
anatomical_regions = {'anterior','posterior','anterior-left','anterior-right','posterior-left','posterior-right'};
stimulation_groups = {'IPL','retrosplenial'};
[Ax,Bx] = ndgrid(1:numel(anatomical_regions),1:numel(stimulation_groups));
n_combos = length(Ax(:));
anatomical_regions = anatomical_regions(Ax(:));
stimulation_groups = stimulation_groups(Bx(:));

%%% Get information about anatomic location to filter table results for
%%% every combination
is_anterior = strcmp(bosc_table.longitudinal_location,'anterior');
is_left = strcmp(bosc_table.hemisphere,'left');

for idx = 1:n_combos
    anatomical_region = anatomical_regions{idx};
    stimulation_group = stimulation_groups{idx};
    this_plot_file_name = fullfile(this_plot_directory,sprintf('%s_%s_hippocampus',stimulation_group,anatomical_region));
    
    switch anatomical_region
        case 'anterior'
            temp_table = bosc_table(is_anterior,:);
        case 'posterior'
            temp_table = bosc_table(~is_anterior,:);
        case 'anterior-left'
            temp_table = bosc_table(is_anterior & is_left,:);
        case 'anterior-right'
            temp_table = bosc_table(is_anterior & ~is_left,:);
        case 'posterior-left'
            temp_table = bosc_table(~is_anterior & is_left,:);
        case 'posterior-right'
            temp_table = bosc_table(~is_anterior & ~is_left,:);
    end
    
    %%% Get electrode results specific to stimulation group
    same_group = strcmp(temp_table.stimulation_group,stimulation_group);
    temp_table = temp_table(same_group,:);
    empty_cells = cellfun(@isempty,temp_table.stimulation_p_episodes) | cellfun(@isempty,temp_table.control_p_episodes);
    temp_table(empty_cells,:) = [];
    
    if ~isempty(temp_table)
        
        %%% Get averages of electrode results for stimulation and control weighted by the number of
        %%% electrodes each subject has
        electrode_weights = arrayfun(@(x) 1/sum(temp_table.subject_ID == x),temp_table.subject_ID);
        n_subjects = length(unique(temp_table.subject_ID));
        stimulation_p_episodes = vertcat(temp_table.stimulation_p_episodes{:});
        control_p_episodes = vertcat(temp_table.control_p_episodes{:});
        stimulation_p_episode = sum(stimulation_p_episodes.*electrode_weights,1)/n_subjects;
        control_p_episode = sum(control_p_episodes.*electrode_weights,1)/n_subjects;
        
        %%% Define y limits based on maxima of results
        y_limits = [0, max([stimulation_p_episode,control_p_episode])*1.1];
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
            value_stimulation = stimulation_p_episode(fdx);
            value_control = control_p_episode(fdx);
            lower_value = min([value_stimulation;value_control]);
            value_difference = abs(value_stimulation-value_control);
            if value_stimulation > value_control
                residual_color = color_stimulation;
            else
                residual_color = color_control;
            end
            bar_handles{fdx} = bar(frequencies(fdx),[lower_value value_difference],'stacked','FaceColor','flat','EdgeColor','none','BarWidth',bar_width);
            these_bars = bar_handles{fdx};
            these_bars(1).CData = mixed_color;
            these_bars(2).CData = residual_color;
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
