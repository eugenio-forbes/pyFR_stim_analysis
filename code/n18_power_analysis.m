function n18_power_analysis(varargin)
if isempty(varargin)
    %%% Directory information
    root_directory = '/directory/to/pyFR_stim_analysis';
    parpool_n = 16;
else
    root_directory = varargin{1};
    parpool_n = varargin{2};
end

%%% List directories
table_directory = fullfile(root_directory,'tables');
plots_directory = fullfile(root_directory,'plots_revision4');
lock_directory = fullfile(root_directory,'locks');
lme_directory = fullfile(root_directory,'lme_results_revision4');

%%% Declare regions of interest and periods of interest
anatomical_regions = {'anterior'}; %{'anterior','posterior'};
hemispheres = {'ipsilateral','contralateral'}; %{'left','right','ipsilateral','contralateral'};
stimulation_groups = {'IPL','retrosplenial','both'};
periods = {'retrieval';'encoding'};
[Ax,Bx,Cx,Dx] = ndgrid(1:numel(periods),1:numel(stimulation_groups),1:numel(anatomical_regions),1:numel(hemispheres));
n_combos = length(Ax(:));
hemispheres = hemispheres(Dx(:));
anatomical_regions = anatomical_regions(Cx(:));
stimulation_groups = stimulation_groups(Bx(:));
periods = periods(Ax(:));
frequencies = 2.^((8:64)/8);
slow_theta_frequencies = frequencies >= 2 & frequencies <=4;

pause(rand*20)

%%%Initialize parpool
pool_object = gcp('nocreate');
if isempty(pool_object)
    parpool(parpool_n)
end

for idx = 1:n_combos
    period = periods{idx};
    anatomical_region = anatomical_regions{idx};
    hemisphere = hemispheres{idx};
    stimulation_group = stimulation_groups{idx};
    
    folder_name = sprintf('%s_power_analysis',period);
    this_lme_directory = fullfile(lme_directory,folder_name);
    this_lock_directory = fullfile(lock_directory,folder_name);
    this_plot_directory = fullfile(plots_directory,folder_name);
    if ~isfolder(this_lme_directory)
        mkdir(this_lme_directory);
    end
    if ~isfolder(this_lock_directory)
        mkdir(this_lock_directory);
    end
    if ~isfolder(this_plot_directory)
        mkdir(this_plot_directory);
    end
    
    file_name = sprintf('%s_%s_%s',stimulation_group,anatomical_region,hemisphere);
    lock_file = fullfile(this_lock_directory,[file_name '_lock.txt']);
    done_file = fullfile(this_lock_directory,[file_name '_done.txt']);
    error_file = fullfile(this_lock_directory,[file_name '_error.txt']);
    mat_file = fullfile(this_lme_directory,[file_name '.mat']);
    plot_file = fullfile(this_plot_directory,file_name);
    
    pause(rand*5)
    
    if ~isfile(lock_file)
        fid = fopen(lock_file,'w'); fclose(fid);
        try
            analysis_table_file = sprintf('%s_analysis_table',period);
            analysis_table = load(fullfile(table_directory,[analysis_table_file '.mat']));
            analysis_table = analysis_table.(analysis_table_file);
            band_table_file = sprintf('%s_band_table',period);
            band_table = load(fullfile(table_directory,[band_table_file '.mat']));
            band_table = band_table.(band_table_file);
            subject_average_table_file = sprintf('%s_subject_average_table',period);
            subject_average_table = load(fullfile(table_directory,[subject_average_table_file '.mat']));
            subject_average_table = subject_average_table.(subject_average_table_file);
           
            [analysis_table,band_table,subject_average_table] = perform_exclusions(root_directory,analysis_table,band_table,subject_average_table);
            
            analysis_table.anterior = logical(analysis_table.anterior);
            analysis_table.left = logical(analysis_table.left);
            analysis_table.stimulation_left = logical(analysis_table.stimulation_left);

            band_table.anterior = logical(band_table.anterior);
            band_table.left = logical(band_table.left);
            band_table.stimulation_left = logical(band_table.stimulation_left);
            
            subject_average_table.anterior = logical(subject_average_table.anterior);
            subject_average_table.left = logical(subject_average_table.left);
            subject_average_table.stimulation_left = logical(subject_average_table.stimulation_left);
            
            switch anatomical_region
                case 'anterior'
                    analysis_table = analysis_table(analysis_table.anterior,:);
                    band_table = band_table(band_table.anterior,:);
                    subject_average_table = subject_average_table(subject_average_table.anterior,:);
                case 'posterior_right'
                    analysis_table = analysis_table(~analysis_table.anterior,:);
                    band_table = band_table(~band_table.anterior,:);
                    subject_average_table = subject_average_table(~subject_average_table.anterior,:);
            end
            
            switch hemisphere
                case 'left'
                    analysis_table = analysis_table(analysis_table.left,:);
                    band_table = band_table(band_table.left,:);
                    subject_average_table = subject_average_table(subject_average_table.left,:);
                case 'right'
                    analysis_table = analysis_table(~analysis_table.left,:);
                    band_table = band_table(~band_table.left,:);
                    subject_average_table = subject_average_table(~subject_average_table.left,:);
                case 'ipsilateral'
                    is_left = analysis_table.left;
                    stimulation_left = analysis_table.stimulation_left;
                    ipsilateral = (is_left & stimulation_left) | (~is_left & ~stimulation_left);
                    analysis_table = analysis_table(ipsilateral,:);

                    is_left = band_table.left;
                    stimulation_left = band_table.stimulation_left;
                    ipsilateral = (is_left & stimulation_left) | (~is_left & ~stimulation_left);
                    band_table = band_table(ipsilateral,:);

                    is_left = subject_average_table.left;
                    stimulation_left = subject_average_table.stimulation_left;
                    ipsilateral = (is_left & stimulation_left) | (~is_left & ~stimulation_left);
                    subject_average_table = subject_average_table(ipsilateral,:);
                case 'contralateral'
                    is_left = analysis_table.left;
                    stimulation_left = analysis_table.stimulation_left;
                    ipsilateral = (is_left & stimulation_left) | (~is_left & ~stimulation_left);
                    analysis_table = analysis_table(~ipsilateral,:);

                    is_left = band_table.left;
                    stimulation_left = band_table.stimulation_left;
                    ipsilateral = (is_left & stimulation_left) | (~is_left & ~stimulation_left);
                    band_table = band_table(~ipsilateral,:);

                    is_left = subject_average_table.left;
                    stimulation_left = subject_average_table.stimulation_left;
                    ipsilateral = (is_left & stimulation_left) | (~is_left & ~stimulation_left);
                    subject_average_table = subject_average_table(~ipsilateral,:);

            end
            
            switch stimulation_group
                case {'IPL','retrosplenial'}
                    if strcmp(stimulation_group,'IPL')
                        same_group = logical(analysis_table.stimulation_lateral);
                    else
                        same_group = ~logical(analysis_table.stimulation_lateral);
                    end
                    analysis_table = analysis_table(logical(same_group),:);
                    is_IPL = logical(analysis_table.stimulation_lateral);
                    analysis_table.stimulation_group = is_IPL;
                    n_subjects = length(unique(analysis_table.subject_ID));
                    n_electrodes = length(unique(analysis_table.electrode_ID));
                    value_bundle = get_effect('stimulation',analysis_table);
                    t_statistics = value_bundle.t_statistics;
                    p_values = value_bundle.p_values;
                    [mask,perimeter] = get_benjamini_hochberg_correction(p_values,frequencies);
                    save(mat_file,'mask','perimeter','n_subjects','n_electrodes','value_bundle');
%                     load(mat_file,'t_statistics','mask','perimeter');
                    plot_effect(plot_file,'stimulation',t_statistics);
                    plot_mask(plot_file,mask,perimeter);
                case {'both'}
                    is_IPL = logical(analysis_table.stimulation_lateral);
                    n_subjects_IPL = length(unique(analysis_table.subject_ID(is_IPL)));
                    n_subjects_RS = length(unique(analysis_table.subject_ID(~is_IPL)));
                    n_electrodes_IPL = length(unique(analysis_table.electrode_ID(is_IPL)));
                    n_electrodes_RS = length(unique(analysis_table.electrode_ID(~is_IPL)));
                    analysis_table.stimulation_group = is_IPL;
                    value_bundle = get_effect('interaction',analysis_table);
                    t_statistics = value_bundle.t_statistics;
                    p_values = value_bundle.p_values;
                    [mask,perimeter] = get_benjamini_hochberg_correction(p_values,frequencies);
                    save(mat_file,'mask','perimeter','n_subjects_IPL','n_subjects_RS',...
                        'n_electrodes_IPL','n_electrodes_RS','value_bundle');
%                     load(mat_file,'t_statistics','mask','perimeter');
                    plot_effect(plot_file,'interaction',t_statistics);
                    plot_mask(plot_file,mask,perimeter);
                    plot_band_average(plot_file,'slow_theta','electrodes',band_table,mask,t_statistics(slow_theta_frequencies,:));
                    plot_band_average(plot_file,'slow_theta','subjects',subject_average_table,mask,t_statistics(slow_theta_frequencies,:));
            end
            fid = fopen(done_file,'w'); fclose(fid);
        catch this_error
            error_message = getReport(this_error, 'extended', 'hyperlinks', 'off');
            fid = fopen(error_file,'w');
            fprintf(fid,error_message);
            fclose(fid);
        end
    end
end

end

function [analysis_table,band_table,subject_average_table] = perform_exclusions(root_directory,analysis_table,band_table,subject_average_table)
bad_sesscodes = [7,8,40,33,35];
bad_electrodes = [322];
take_out = ismember(analysis_table.electrode_ID,bad_electrodes) | ismember(analysis_table.session_ID,bad_sesscodes);
analysis_table(take_out,:) = [];
take_out = ismember(band_table.electrode_ID,bad_electrodes) | ismember(band_table.session_ID,bad_sesscodes);
band_table(take_out,:) = [];

% %%% List directories
% list_directory = fullfile(root_directory,'lists');
% exclusion_directory = fullfile(root_directory,'exclusion_lists');
% 
% %%% Load subject, session, electrode lists and exclusion lists
% load(fullfile(list_directory,'electrode_list.mat'),'electrode_list')
% load(fullfile(list_directory,'session_list.mat'),'session_list')
% load(fullfile(list_directory,'subject_list.mat'),'subject_list')
% load(fullfile(exclusion_directory,'excluded_electrodes.mat'),'excluded_electrodes')
% load(fullfile(exclusion_directory,'excluded_sessions.mat'),'excluded_sessions')
% load(fullfile(exclusion_directory,'excluded_subjects.mat'),'excluded_subjects')
% 
% %%% Exclusions based on outlying power values
% power = double(analysis_table.power);
% 
% outliers = ((power-mean(power))./std(power))>4;
% outlier_electrodes = unique(analysis_table.electrode_ID(outliers));
% outlier_sessions = unique(analysis_table.session_ID(outliers));
% % outlier_subjects = unique(analysis_table.subject_ID(outliers));
% 
% bad_rows = ismember(analysis_table.electrode_ID,outlier_electrodes);
% analysis_table(bad_rows,:) = [];
% bad_rows = ismember(band_table.electrode_ID,outlier_electrodes);
% band_table(bad_rows,:) = [];
% % excluded = ismember(electrode_list.electrode_ID,outlier_electrodes);
% % try
% %     new_excluded_electrodes = electrode_list(excluded,{'subject','task','session','channel_number','subject_ID','session_ID','electrode_ID'});
% %     new_excluded_electrodes.reason_for_exclusion = repelem({'bad_power'},sum(excluded),1);
% %     excluded_electrodes = [excluded_electrodes;new_excluded_electrodes];
% % catch
% %         new_excluded_electrodes = electrode_list(excluded,{'subject','task','session','subject_ID','session_ID','electrode_ID'});
% %         new_excluded_electrodes.reason_for_exclusion = repelem({'bad_power'},sum(excluded),1);
% %         excluded_electrodes = [excluded_electrodes;new_excluded_electrodes];
% % end
% % 
% % n_available = arrayfun(@(x) sum(electrode_list.session_ID == x),outlier_sessions);
% % n_excluded = arrayfun(@(x) sum(new_excluded_electrodes.session_ID == x),outlier_sessions);
% % too_many_exclusions = n_excluded./n_available >= 0.5;
% % bad_sessions = outlier_sessions(too_many_exclusions);
% % bad_rows = ismember(analysis_table.session_ID,bad_sessions);
% % bad_session_electrodes = unique(analysis_table.electrode_ID(bad_rows));
% % analysis_table(bad_rows,:) = [];
% % bad_rows = ismember(band_table.session_ID,bad_sessions);
% % band_table(bad_rows,:) = [];
% % excluded = ismember(session_list.session_ID,bad_sessions);
% % new_excluded_sessions = session_list(excluded,{'subject','task','session','subject_ID','session_ID'});
% % new_excluded_sessions.reason_for_exclusion = repelem({'bad_power'},sum(excluded),1);
% % excluded_sessions = [excluded_sessions;new_excluded_sessions];
% % 
% % n_available = arrayfun(@(x) sum(electrode_list.subject_ID == x),outlier_subjects);
% % n_excluded = arrayfun(@(x) sum(new_excluded_electrodes.subject_ID == x),outlier_subjects);
% % too_many_exclusions = n_excluded./n_available >= 0.5;
% % bad_subjects = outlier_subjects(too_many_exclusions);
% % bad_rows = ismember(subject_average_table.subject_ID,bad_subjects);
% % subject_average_table(bad_rows,:) = [];
% % bad_rows = ismember(analysis_table.subject_ID,bad_subjects);
% % bad_subject_electrodes = unique(analysis_table.electrode_ID(bad_rows));
% % analysis_table(bad_rows,:) = [];
% % bad_rows = ismember(band_table.subject_ID,bad_subjects);
% % band_table(bad_rows,:) = [];
% % excluded = ismember(subject_list.subject_ID,bad_subjects);
% % new_excluded_subjects = subject_list(excluded,{'subject','subject_ID'});
% % new_excluded_subjects.reason_for_exclusion = repelem({'bad_power'},sum(excluded),1);
% % excluded_subjects = [excluded_subjects;new_excluded_subjects];
% % 
% % excluded = ismember(electrode_list.electrode_ID,[bad_session_electrodes;bad_subject_electrodes]);
% % try
% %     new_excluded_electrodes = electrode_list(excluded,{'subject','task','session','channel_number','subject_ID','session_ID','electrode_ID'});
% %     new_excluded_electrodes.reason_for_exclusion = repelem({'bad_power_subject_session'},sum(excluded),1);
% %     excluded_electrodes = [excluded_electrodes;new_excluded_electrodes];
% % catch
% %     new_excluded_electrodes = electrode_list(excluded,{'subject','task','session','subject_ID','session_ID','electrode_ID'});
% %     new_excluded_electrodes.reason_for_exclusion = repelem({'bad_power_subject_session'},sum(excluded),1);
% %     excluded_electrodes = [excluded_electrodes;new_excluded_electrodes];
% % end
% % 
% % save(fullfile(exclusion_directory,'excluded_subjects.mat'),'excluded_subjects');
% % save(fullfile(exclusion_directory,'excluded_sessions.mat'),'excluded_sessions');
% % save(fullfile(exclusion_directory,'excluded_electrodes.mat'),'excluded_electrodes');

end

function value_bundle = get_effect(effect_type,analysis_table)
value_bundle = struct;
analysis_table.subject_ID = categorical(analysis_table.subject_ID);
analysis_table.electrode_ID = categorical(analysis_table.electrode_ID);
analysis_table.stimulation_left = logical(analysis_table.stimulation_left);
analysis_table.power = double(analysis_table.power);

random_effect = '1 + condition';
switch effect_type
    case 'stimulation'
        formula = sprintf('power ~ condition + (%s|subject_ID) + (%s|subject_ID:electrode_ID) + (%s|subject_ID:stimulation_left) + (%s|stimulation_group:stimulation_left)',...
            random_effect,random_effect,random_effect,random_effect);
        formula = sprintf('power ~ condition + (%s|subject_ID) + (%s|subject_ID:electrode_ID)',...
            random_effect,random_effect);
        statistics_index = 2;
    case 'interaction'
        formula = sprintf('power ~ condition*stimulation_group + (%s|subject_ID) + (%s|subject_ID:electrode_ID) + (%s|subject_ID:stimulation_left) + (%s|stimulation_group:stimulation_left)',...
            random_effect,random_effect,random_effect,random_effect);
        statistics_index = 4;
end

n_frequencies = length(unique(analysis_table.frequency));
n_samples = length(unique(analysis_table.time));

t_statistics = NaN(n_frequencies,n_samples);
p_values = NaN(n_frequencies,n_samples);
degrees_of_freedom = NaN(n_frequencies,n_samples);
F_stat = NaN(n_frequencies,n_samples);
p_F_stat = NaN(n_frequencies,n_samples);
DF_numerator = NaN(n_frequencies,n_samples);
DF_denominator = NaN(n_frequencies,n_samples);


parfor tdx = 1:n_samples
    this_time = analysis_table.time == tdx;
    temp_table = analysis_table(this_time,:);
    frequency_tables = arrayfun(@(x) temp_table(temp_table.frequency==x,:),1:n_frequencies,'UniformOutput',false);
    for fdx = 1:n_frequencies
        frequency_table = frequency_tables{fdx};
        model = fitlme(frequency_table,formula);
        [betas,beta_names,stats] = fixedEffects(model,'DFMethod','satterthwaite','alpha',0.05);
        H = ones(1,height(beta_names));
        H(1) = 0;
        C = H*betas;
        t_statistics(fdx,tdx) = stats.tStat(statistics_index);
        p_values(fdx,tdx) = stats.pValue(statistics_index);
        degrees_of_freedom(fdx,tdx) = stats.DF(statistics_index);
        [p_F_stat(fdx,tdx),F_stat(fdx,tdx),DF_numerator(fdx,tdx),DF_denominator(fdx,tdx)] = coefTest(model,H,C,'DFMethod','satterthwaite');
    end
end
value_bundle.t_statistics = t_statistics;
value_bundle.p_values = p_values;
value_bundle.degrees_of_freedom = degrees_of_freedom;
value_bundle.p_F_stat = p_F_stat;
value_bundle.F_stat = F_stat;
value_bundle.DF_numerator = DF_numerator;
value_bundle.DF_denominator = DF_denominator;

end

function [mask,perimeter] = get_benjamini_hochberg_correction(p_values,frequencies)
alpha = 0.05;
method = 'pdep';
verbose = 'no';

[n_frequencies,n_samples] = size(p_values);
theta_frequencies = frequencies >= 2 & frequencies <= 8;

adjusted_p_values = NaN(n_frequencies,n_samples);

for fdx = 1:n_frequencies
[~,~,~,adjusted_p_values(fdx,:)] = fdr_bh(p_values(fdx,:),alpha,method,verbose);
end

mask = adjusted_p_values < alpha;
mask(~theta_frequencies,:) = false;
perimeter = bwboundaries(mask);
end

function plot_effect(plot_file,plot_type,t_statistics)
%%% Plot parameters
switch plot_type
    case 'stimulation'
        mapx = makecolormap_EF('sigmoid3');
    case 'interaction'
        mapx = makecolormap_EF('uniform3');
end
[n_frequencies,n_samples] = size(t_statistics);
y_ticks = 1:8:n_frequencies;
y_limits = [0.5 n_frequencies+0.5];
x_ticks = 0:200:n_samples;
x_limits = [0.5 n_samples+0.5];
figure_width = 80*(n_samples/1000);
figure_height = 80;

%%% Plotting
fig = figure('Units','pixels','Visible','off','Position', [0 0 figure_width figure_height]);
axes('Parent',fig,'Units', 'pixels','Position', [0,0,figure_width,figure_height]);
imagesc(t_statistics);
set(gca,'YDir','normal');
xlim(x_limits);xticks(x_ticks);xticklabels([]);
ylim(y_limits);yticks(y_ticks);yticklabels([]);
colormap(mapx);
caxis([-3 3]);

%%% Saving
print(plot_file,'-dsvg')
close all
end

function plot_mask(plot_file,mask,perimeter)
plot_file = [plot_file '_mask'];

%%% Plot parameters
[n_frequencies,n_samples] = size(mask);
y_ticks = 1:8:n_frequencies;
y_limits = [0.5 n_frequencies+0.5];
x_ticks = 0:200:n_samples;
x_limits = [0.5 n_samples+0.5];
figure_width = 80*(n_samples/1000);
figure_height = 80;

%%% Plotting
fig = figure('Units','pixels','Visible','off','Position', [0 0 figure_width figure_height]);
axes('Parent',fig,'Units', 'pixels','Position', [0,0,figure_width,figure_height]);
hold on
for bdx = 1:length(perimeter)
    boundary = perimeter{bdx};
    plot(boundary(:,2),boundary(:,1),'k','LineWidth',1.5)
end
xlim(x_limits);xticks(x_ticks);xticklabels([]);
ylim(y_limits);yticks(y_ticks);yticklabels([]);
hold off

%%% Saving
print(plot_file,'-dsvg')
close all
end

function plot_band_average(plot_file,band,table_type,band_table,mask,t_statistics)
plot_file = [plot_file sprintf('_%s_%s_band',table_type,band)];

%%% Plot parameters
[~,n_samples] = size(mask);
y_ticks = -2:0.5:2;
y_limits = [-2.25 2.25];
x_ticks = 0:200:n_samples;
x_limits = [0.5 n_samples+0.5];
figure_width = 80*(n_samples/1000);
figure_height = 70;
mapx = makecolormap_EF('uniform3');
color_IPL = mapx(end,:);
color_RS = mapx(1,:);

mask = sum(mask,1);
mask = mask>0;
t_statistics = mean(t_statistics,1);

%%% Data
switch band
    case 'slow_theta'
        band_data = double(vertcat(band_table.slow_theta_tstats{:}))/1000;
    case 'fast_theta'
        band_data = double(vertcat(band_table.fast_theta_tstats{:}))/1000;
end

is_IPL = logical(band_table.stimulation_lateral);
% n_IPL = sum(is_IPL);
% n_RS = sum(~is_IPL);
IPL_band_data = band_data(is_IPL,:);
RS_band_data = band_data(~is_IPL,:);

% IPL_mean = mean(IPL_band_data,1);
% IPL_sem = std(IPL_band_data,[],1)/sqrt(n_IPL);
% RS_mean = mean(RS_band_data,1);
% RS_sem = std(RS_band_data,[],1)/sqrt(n_RS);

%%% Plotting
fig = figure('Units','pixels','Visible','off','Position', [0 0 figure_width figure_height]);
axes('Parent',fig,'Units', 'pixels','Position', [0,0,figure_width,figure_height]);
hold on

for mdx = 1:n_samples
    t_statistic = t_statistics(mdx);
    mask_value = mask(mdx);
    if mask_value
        if t_statistic > 0
            patch_color = color_IPL;
        else
            patch_color = color_RS;
        end
        patch([mdx-1 mdx-1 mdx mdx mdx-1],[-2.5 2.5 2.5 -2.5 -2.5],patch_color,'FaceAlpha',0.2,'LineStyle','none')
    end
end
stdshade(IPL_band_data,0.4,color_IPL);
stdshade(RS_band_data,0.4,color_RS);
plot([0 n_samples],[0 0],'-k','LineWidth',1);

xlim(x_limits);xticks(x_ticks);xticklabels([]);
ylim(y_limits);yticks(y_ticks);yticklabels([]);

hold off

%%% Saving
print(plot_file,'-dsvg')
close all
end
