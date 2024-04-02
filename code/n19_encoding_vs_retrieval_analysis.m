function n19_encoding_vs_retrieval_analysis(varargin)
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
plots_directory = fullfile(root_directory,'plots');
lock_directory = fullfile(root_directory,'locks');
lme_directory = fullfile(root_directory,'lme_results');

%%% Declare regions of interest and periods of interest
anatomical_regions = {'anterior','posterior'};
hemispheres = {'ipsilateral','contralateral','left','right'};
stimulation_groups = {'IPL','retrosplenial'};
[Ax,Bx,Cx] = ndgrid(1:numel(stimulation_groups),1:numel(anatomical_regions),1:numel(hemispheres));
n_combos = length(Ax(:));
hemispheres = hemispheres(Cx(:));
anatomical_regions = anatomical_regions(Bx(:));
stimulation_groups = stimulation_groups(Ax(:));

analysis_table_file = 'encoding_vs_retrieval_analysis_table';
analysis_table = load(fullfile(table_directory,[analysis_table_file '.mat']));
analysis_table = analysis_table.(analysis_table_file);
analysis_table = check_previous_exclusions(root_directory,analysis_table);
is_left = logical(analysis_table.left);
is_anterior = logical(analysis_table.anterior);
stimulation_left = logical(analysis_table.stimulation_left);
ipsilateral = (stimulation_left & is_left) | (~stimulation_left & ~is_left);

%%%Initialize parpool
pool_object = gcp('nocreate');
if isempty(pool_object)
    parpool(parpool_n)
end

for idx = 1:n_combos
    anatomical_region = anatomical_regions{idx};
    stimulation_group = stimulation_groups{idx};
    hemisphere = hemispheres{idx};
    
    folder_name = 'encoding_vs_retrieval_power_analysis';
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
    
    if isfile(lock_file)
        fid = fopen(lock_file,'w'); fclose(fid);
        try
            
            switch anatomical_region
                case 'anterior'
                    this_analysis_table = analysis_table(is_anterior,:);
                case 'posterior'
                    this_analysis_table = analysis_table(~is_anterior,:);
            end
            
            switch hemisphere
                case 'left'
                    this_analysis_table = analysis_table(is_left,:);
                case 'right'
                    this_analysis_table = analysis_table(~is_left,:);
                case 'ipsilateral'
                    this_analysis_table = analysis_table(ipsilateral,:);
                case 'contralateral'
                    this_analysis_table = analysis_table(~ipsilateral,:);
            end
            
            is_IPL = logical(this_analysis_table.stimulation_lateral);
            this_analysis_table.stimulation_group = is_IPL;
            if strcmp(stimulation_group,'IPL')
                same_group = is_IPL;
            else
                same_group = ~is_IPL;
            end
            this_analysis_table = this_analysis_table(logical(same_group),:);
            n_subjects = length(unique(this_analysis_table.subject_ID));
            n_electrodes = length(unique(this_analysis_table.electrode_ID));
            [t_statistics,p_values,degrees_of_freedom] = encoding_vs_retrieval_lme(this_analysis_table);
            save(mat_file,'t_statistics','p_values','n_subjects','n_electrodes','degrees_of_freedom');
            plot_encoding_vs_retrieval(plot_file,t_statistics,p_values);

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

function analysis_table = check_previous_exclusions(root_directory,analysis_table)
exclusion_directory = fullfile(root_directory,'exclusion_lists');

%%% Load subject, session, electrode lists and exclusion lists
load(fullfile(exclusion_directory,'excluded_electrodes.mat'),'excluded_electrodes')

bad_electrodes = excluded_electrodes.electrode_ID;
exclude = ismember(analysis_table.electrode_ID,bad_electrodes);
analysis_table(exclude,:) = [];

end

function [t_statistics,p_values,degrees_of_freedom] = encoding_vs_retrieval_lme(analysis_table)
analysis_table.subject_ID = categorical(analysis_table.subject_ID);
analysis_table.electrode_ID = categorical(analysis_table.electrode_ID);
analysis_table.stimulation_left = logical(analysis_table.stimulation_left);
analysis_table.power = double(analysis_table.power);

random_effect = '1 + condition + event_type + condition*event_type';
formula = sprintf('power ~ condition*event_type + (%s|subject_ID) + (%s|subject_ID:electrode_ID) + (%s|subject_ID:stimulation_left) + (%s|stimulation_group:stimulation_left)',...
    random_effect,random_effect,random_effect,random_effect);
statistics_index = 4;

n_frequencies = length(unique(analysis_table.frequency));

t_statistics = NaN(n_frequencies,1);
p_values = NaN(n_frequencies,1);
degrees_of_freedom = NaN(n_frequencies,1);

frequency_tables = arrayfun(@(x) analysis_table(analysis_table.frequency==x,:),1:n_frequencies,'UniformOutput',false);
for fdx = 1:n_frequencies
    frequency_table = frequency_tables{fdx};
    model = fitlme(frequency_table,formula);
    t_statistics(fdx) = model.Coefficients.tStat(statistics_index);
    p_values(fdx) = model.Coefficients.pValue(statistics_index);
    degrees_of_freedom(fdx) = model.Coefficients.DF(statistics_index);
end

degrees_of_freedom = degrees_of_freedom(1);

end

function plot_encoding_vs_retrieval(plot_file,t_statistics,p_values)
%%% Plot parameters
mapx = makecolormap_EF('uniform2');
n_frequencies = length(t_statistics);
y_ticks = 1:8:n_frequencies;
y_limits = [0.5 n_frequencies+0.5];
x_ticks = -4:1:4;
x_limits = [-4 4];
figure_width = 60;
figure_height = 80;

color_encoding = mapx(1,:);
color_retrieval = mapx(end,:);

%%% Plotting
fig = figure('Units','pixels','Visible','off','Position', [0 0 figure_width figure_height]);
axes('Parent',fig,'Units', 'pixels','Position', [2,2,figure_width-4,figure_height-4]);
hold on
for fdx = 1:n_frequencies
    
    t_statistic = t_statistics(fdx);
    p_value = p_values(fdx);
    
    bar_y_position = [fdx-0.5,fdx+0.5,fdx+0.5,fdx-0.5,fdx-0.5];
    
    color_index = round((t_statistic * 1000) + 3000);
    if color_index < 1
        color_index = 1;
    elseif color_index > 6000
        color_index = 6000;
    end
    this_bar_color = mapx(color_index,:);
    
    if p_value < 0.05
        if t_statistic < 0
            this_significant_color = color_encoding;
        else
            this_significant_color = color_retrieval;
        end
        patch([-1,-1,1,1,-1]*x_limits(2),bar_y_position,this_significant_color,'FaceAlpha',0.6,'LineStyle','none');
    end
    
    patch([0,0,1,1,0]*t_statistic,bar_y_position,this_bar_color);
    
end
plot([0 0],y_limits,'-k');
xlim(x_limits);xticks(x_ticks);xticklabels([]);
ylim(y_limits);yticks(y_ticks);yticklabels([]);
hold off

%%% Saving
print(plot_file,'-dsvg')
close all
end
