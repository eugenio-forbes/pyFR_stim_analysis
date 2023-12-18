function n22_correlation_analyses(varargin)
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
table_directory = fullfile(analysis_directory,'tables');
plots_directory = fullfile(analysis_directory,'plots');
lock_directory = fullfile(analysis_directory,'locks');
statistics_directory = fullfile(analysis_directory,'statistics');

%%% Load subject behavioral, power, phase clustering, tables
load(fullfile(statistics_directory,'subject_behavioral_table.mat'),'subject_behavioral_table');
subject_behavioral_table = check_bad_subjects(analysis_directory,subject_behavioral_table);
behavioral_subjects = subject_behavioral_table.subject_ID;

%%% Declare regions of interest and periods of interest
correlation_metrics = {'overall_percent_recall';'tstat_recalls';'zvalue_recalls';...
    'total_PLIs';'total_XLIs';'total_intrusions';'TCF';'SCF'};
frequency_bands = {'slow_theta';'fast_theta'};
anatomical_regions = {'anterior';'posterior'};
hemispheres = {'left';'right';'ipsilateral';'contralateral'};
stimulation_hemispheres = {'left';'right';'ipsilateral';'contralateral'};
stimulation_groups = {'IPL';'retrosplenial'};
periods = {'encoding';'retrieval'};

%%% Make combinations of parameters to loop through
[Ax,Bx,Cx,Dx,Ex,Fx] = ndgrid(1:numel(periods),1:numel(stimulation_groups),1:numel(anatomical_regions),...
    1:numel(hemispheres),1:numel(frequency_bands),1:numel(correlation_metrics));
n_combos = length(Ax(:));

correlation_metric = correlation_metrics(Fx(:));
frequency_band = frequency_bands(Ex(:));
hemisphere = hemispheres(Dx(:));
anatomical_region = anatomical_regions(Cx(:));
stimulation_group = stimulation_groups(Bx(:));
period = periods(Ax(:));

%%% Intialize arrarys to hold correlation results and data for each
%%% combination of parameters
spearman_rho = NaN(n_combos,1);
r_squared = NaN(n_combos,1);
p_value = NaN(n_combos,1);
n_subjects = NaN(n_combos,1);

%%% Declare conditions for each period to extract information from table
encoding_conditions = [1,5];
retrieval_conditions = [3,5];
control_condition = 5;

%%%Initialize parpool
pool_object = gcp('nocreate');
if isempty(pool_object)
    parpool(24)
end

%%% Loop through combinations in parallel
parfor idx = 1:n_combos
    %%% Combination info
    this_period = period{idx};
    this_anatomical_region = anatomical_region{idx};
    this_stimulation_group = stimulation_group{idx};
    this_hemisphere = hemisphere{idx};
    this_frequency_band = frequency_band{idx};
    this_correlation_metric = correlation_metric{idx};
    
    %%% Use values to get file names for results, lock, done, and error
    %%% files
    folder_name = sprintf('%s_correlation_analyses/%s/%s',this_period,this_frequency_band,this_correlation_metric);
    this_lock_directory = fullfile(lock_directory,folder_name);
    this_plot_directory = fullfile(plots_directory,folder_name);
    if ~isfolder(this_lock_directory)
        mkdir(this_lock_directory);
    end
    if ~isfolder(this_plot_directory)
        mkdir(this_plot_directory);
    end
    
    file_name = sprintf('%s_%s_%s',this_stimulation_group,this_anatomical_region,this_hemisphere);
    lock_file = fullfile(this_lock_directory,[file_name '_lock.txt']);
    done_file = fullfile(this_lock_directory,[file_name '_done.txt']);
    error_file = fullfile(this_lock_directory,[file_name '_error.txt']);
    plot_file = fullfile(this_plot_directory,file_name);
    
    %%% Random length pause to separate nodes executing code simultaneously
    %%% Before they check for the existance of lock file
    pause(rand)
    if ~isfile(lock_file)
        fid = fopen(lock_file,'w'); fclose(fid);
        
        try 
            %%% Load table with power averaged by subject
            power_file = sprintf('%s_subject_average_table',this_period);
            power_table = load(fullfile(table_directory,[power_file '.mat']));
            power_table = power_table.(power_file);
                        
            %%% Flter table based on anatomical region
            is_anterior = power_table.anterior;
            switch this_anatomical_region
                case 'anterior'
                    power_table = power_table(is_anterior,:);
                case 'posterior'
                    power_table = power_table(~is_anterior,:);
            end
            
            %%% Filter table based on hemisphere
            is_left = power_table.left;
            stimulation_left = power_table.stimulation_left;
            is_ipsilateral = (is_left & stimulation_left) | (~is_left & ~stimulation_left);
            switch this_hemisphere
                case 'left'
                    power_table = power_table(is_left,:);
                case 'right'
                    power_table = power_table(~is_left,:);
                case 'ipsilateral'
                    power_table = power_table(is_ipsilateral,:);
                case 'contralateral'
                    power_table = power_table(~is_ipsilateral,:);
            end
        
            %%% Filter table based on group
            if strcmp(this_stimulation_group,'IPL')
                same_group = power_table.stimulation_lateral;
            else
                same_group = ~power_table.stimulation_lateral;
            end
            power_table = power_table(same_group,:);
            subject_IDs = power_table.subject_ID;
            
            %%% Identify behavioral table rows that have the same subjects
            %%% as filtered power table
            has_subjects = ismember(behavioral_subjects,subject_IDs);
            behavioral_table = subject_behavioral_table(has_subjects,:);
            
            %%% Filter behavioral table based on period relevant coditions
            switch this_period
                case 'encoding'
                    has_conditions = ismember(behavioral_table.condition,encoding_conditions);
                case 'retrieval'
                    has_conditions = ismember(behavioral_table.condition,retrieval_conditions);
            end
            behavioral_table = behavioral_table(has_conditions,:);
            
            is_control = behavioral_table.condition == control_condition;
            is_stimulation = ~is_control;
            
            %%% Obtain correlation metric from behavioral table based on
            %%% metric for metrics other than t-stat and z-value the
            %%% difference between experimental condition and control needs
            %%% to be calculated
            switch this_correlation_metric
                case {'tstat_recalls','zvalue_recalls'}
                    correlation_table = behavioral_table(is_stimulation,{'subject_ID',this_correlation_metric});
                    correlation_table.metric = correlation_table.(this_correlation_metric);
                otherwise
                    stimulation_values = behavioral_table(is_stimulation,{'subject_ID',this_correlation_metric});
                    control_values = behavioral_table(is_control,{'subject_ID',this_correlation_metric});
                    stimulation_IDs = stimulation_values.subject_ID;
                    control_IDs = control_values.subject_ID;
                    match_control = ismember(stimulation_IDs,control_IDs);
                    match_stimulation = ismember(control_IDs,stimulation_IDs);
                    stimulation_values = stimulation_values(match_control,:);
                    control_values = control_values(match_stimulation,:);
                    stimulation_IDs = stimulation_IDs(match_control);
                    control_IDs = control_IDs(match_stimulation);
                    new_order = arrayfun(@(x) find(control_IDs == x),stimulation_IDs);
                    control_values = control_values(new_order,:);
                    correlation_table = stimulation_values(:,{'subject_ID'});
                    correlation_table.metric = stimulation_values.(this_correlation_metric) - control_values.(this_correlation_metric);
            end
            
            %%% Exclude bad observations
            bad_values = isnan(correlation_table.metric);
            correlation_table(bad_values,:) = [];
            
            %%% Count subjects
            n_subjects(idx) = height(correlation_table);
            subject_IDs = correlation_table.subject_ID;
            
            %%% Get subject mean power values of relevant frequency band
            subject_indices = arrayfun(@(x) find(power_table.subject_ID == x),subject_IDs);
            power_table = power_table(subject_indices,:);
            switch this_frequency_band
                case 'slow_theta'
                    correlation_table.power = mean(double(vertcat(power_table.slow_theta_tstats{:}))/1000,2);
                case 'fast_theta'
                    correlation_table.power = mean(double(vertcat(power_table.fast_theta_tstats{:}))/1000,2);
            end
            
            %%% Count number of rows before performing Spearman Rank
            %%% Correlation
            n_rows = height(correlation_table)
            if n_rows > 3
                [spearman_rho(idx),r_squared(idx),p_value(idx)] = do_correlation('behavioral',plot_file,correlation_table);
            end
            
            %%% Succesful execution of code. Open and close done file
            fid = fopen(done_file,'w'); fclose(fid);
        catch this_error
            %%% Unsuccesful execution of code. Print error message in error
            %%% file
            error_message = getReport(this_error, 'extended', 'hyperlinks', 'off');
            fid = fopen(error_file,'w');
            fprintf(fid,error_message);
            fclose(fid);
        end
    end
end

%%% Create and save results table with pre-initialized and filled out
%%% columns
behavioral_correlation_table = table(correlation_metric,frequency_band,...
    hemisphere,anatomical_region,stimulation_group,period,...
    spearman_rho,r_squared,p_value,n_subjects);
save(fullfile(statistics_directory,'behavioral_correlation_table.mat'),'behavioral_correlation_table');

%%% Make cominations of periods, groups, anatomical regions, stimulation
%%% hemispheres, and frequency bands to make correlations of phase
%%% clustering and power changes and basiline stim site power
[Ax,Bx,Cx,Dx,Ex] = ndgrid(1:numel(periods),1:numel(stimulation_groups),1:numel(anatomical_regions),...
    1:numel(hemispheres),1:numel(frequency_bands));
n_combos = length(Ax(:));

correlation_metric = repelem({'phase_clustering'},n_combos,1);
frequency_band = frequency_bands(Ex(:));
hemisphere = stimulation_hemispheres(Dx(:));
anatomical_region = anatomical_regions(Cx(:));
stimulation_group = stimulation_groups(Bx(:));
period = periods(Ax(:));

%%% Initialize columns to hold correlation results for all combinations
spearman_rho = repelem({NaN(1,2)},n_combos,1);
r_squared = repelem({NaN(1,2)},n_combos,1);
p_value = repelem({NaN(1,2)},n_combos,1);
n_subjects = NaN(n_combos,1);

%%% Loop through combinations
for idx = 1:n_combos
    %%% Combo info
    this_period = period{idx};
    this_anatomical_region = anatomical_region{idx};
    this_stimulation_group = stimulation_group{idx};
    this_hemisphere = hemisphere{idx};
    this_frequency_band = frequency_band{idx};
    
    %%% Use combo info to get file names
    folder_name = sprintf('%s_correlation_analyses/%s/phase_clustering',this_period,this_frequency_band);
    this_lock_directory = fullfile(lock_directory,folder_name);
    this_plot_directory = fullfile(plots_directory,folder_name);
    if ~isfolder(this_lock_directory)
        mkdir(this_lock_directory);
    end
    if ~isfolder(this_plot_directory)
        mkdir(this_plot_directory);
    end
    
    file_name = sprintf('%s_%s_%s',this_stimulation_group,this_anatomical_region,this_hemisphere);
    lock_file = fullfile(this_lock_directory,[file_name '_lock.txt']);
    done_file = fullfile(this_lock_directory,[file_name '_done.txt']);
    error_file = fullfile(this_lock_directory,[file_name '_error.txt']);
    plot_file = fullfile(this_plot_directory,file_name);
    
    %%% Random pause to separate nodes with simultaneous execution
    pause(rand*5)
    if ~isfile(lock_file)
        fid = fopen(lock_file,'w'); fclose(fid);
        try
            %%% Get subject average power change, phase clustering, and
            %%% stim site base line power
            stim_site_file = sprintf('%s_stim_site_table',this_period);
            stim_site_table = load(fullfile(table_directory,[stim_site_file '.mat']));
            stim_site_table = stim_site_table.(stim_site_file);
            stim_site_table = check_bad_subjects(analysis_directory,stim_site_table);
            
            power_file = sprintf('%s_subject_average_table',this_period);
            power_table = load(fullfile(table_directory,[power_file '.mat']));
            power_table = power_table.(power_file);
            power_table = check_bad_subjects(analysis_directory,power_table);
            
            phase_file = sprintf('%s_subject_phase_clustering_table',this_period);
            phase_table = load(fullfile(table_directory,[phase_file '.mat']));
            phase_table = phase_table.(phase_file);
            phase_table = check_bad_subjects(analysis_directory,phase_table);
            
            %%% Filter tables based on combination
            phase_is_anterior = logical(phase_table.anterior);
            power_is_anterior = logical(power_table.anterior);
            
            switch this_anatomical_region
                case 'anterior'
                    phase_table = phase_table(phase_is_anterior,:);
                    power_table = power_table(power_is_anterior,:);
                case 'posterior'
                    phase_table = phase_table(~phase_is_anterior,:);
                    power_table = power_table(~power_is_anterior,:);
            end
            
            phase_is_left = logical(phase_table.left);
            power_is_left = logical(power_table.left);
            phase_stimulation_left = logical(phase_table.stimulation_left);
            power_stimulation_left = logical(power_table.stimulation_left);
            phase_is_ipsilateral = (phase_is_left & phase_stimulation_left) | (~phase_is_left & ~phase_stimulation_left);
            power_is_ipsilateral = (power_is_left & power_stimulation_left) | (~power_is_left & ~power_stimulation_left);
            
            switch this_hemisphere
                case 'left'
                    phase_table = phase_table(phase_is_left,:);
                    power_table = power_table(power_is_left,:);
                case 'right'
                    phase_table = phase_table(~phase_is_left,:);
                    power_table = power_table(~power_is_left,:);
                case 'ipsilateral'
                    phase_table = phase_table(phase_is_ipsilateral,:);
                    power_table = power_table(power_is_ipsilateral,:);
                case 'contralateral'
                    phase_table = phase_table(~phase_is_ipsilateral,:);
                    power_table = power_table(~power_is_ipsilateral,:);
            end
        
            if strcmp(this_stimulation_group,'IPL')
                same_group = logical(phase_table.stimulation_lateral);
            else
                same_group = ~logical(phase_table.stimulation_lateral);
            end
            phase_table = phase_table(same_group,:);
            
            has_subjects = ismember(power_table.subject_ID,phase_table.subject_ID);
            power_table = power_table(has_subjects,:);
            has_subjects = ismember(phase_table.subject_ID,power_table.subject_ID);
            phase_table = phase_table(has_subjects,:);
            has_subjects = ismember(stim_site_table.subject_ID,phase_table.subject_ID);
            stim_site_table = stim_site_table(has_subjects,:);
            
            subject_IDs = phase_table.subject_ID;
            
            subject_indices = arrayfun(@(x) find(power_table.subject_ID == x),subject_IDs);
            power_table = power_table(subject_indices,:);
            subject_indices = arrayfun(@(x) find(stim_site_table.subject_ID == x),subject_IDs);
            stim_site_table = stim_site_table(subject_indices,:);
            
            correlation_table = power_table(:,{'subject_ID'});
            
            switch this_frequency_band
                case 'slow_theta'
                    correlation_table.power = mean(double(vertcat(power_table.slow_theta_tstats{:}))/1000,2);
                    correlation_table.phase = phase_table.slow_theta_phase_clustering;
                    correlation_table.stim_site_power = double(stim_site_table.slow_theta_power)/1000;
                case 'fast_theta'
                    correlation_table.power = mean(double(vertcat(power_table.fast_theta_tstats{:}))/1000,2);
                    correlation_table.phase = phase_table.fast_theta_phase_clustering;
                    correlation_table.stim_site_power = double(stim_site_table.fast_theta_power)/1000;
            end
            
            
            %%% Check number of rows before performing Spearman Correlation
            n_rows = height(correlation_table);
            
            if n_rows > 3
                [spearman_rho{idx}(1),r_squared{idx}(1),p_value{idx}(1)] = do_correlation('phase',plot_file,correlation_table);
                [spearman_rho{idx}(2),r_squared{idx}(2),p_value{idx}(2)] = do_correlation('stim_site',plot_file,correlation_table);
            end
            
            %%% Open donce file
            fid = fopen(done_file,'w'); fclose(fid);
        catch this_error
            %%% Print error message
            error_message = getReport(this_error, 'extended', 'hyperlinks', 'off');
            fid = fopen(error_file,'w');
            fprintf(fid,error_message);
            fclose(fid);
        end
    end
end

phase_clustering_correlation_table = table(correlation_metric,frequency_band,...
    hemisphere,anatomical_region,stimulation_group,period,...
    spearman_rho,r_squared,p_value,n_subjects);

save(fullfile(statistics_directory,'phase_clustering_correlation_table.mat'),'phase_clustering_correlation_table');

end

function analysis_table = check_bad_subjects(analysis_directory,analysis_table)
exclusion_directory = fullfile(analysis_directory,'exclusion_lists');

%%% Load subject, session, electrode lists and exclusion lists
load(fullfile(exclusion_directory,'excluded_subjects.mat'),'excluded_subjects')

bad_subjects = excluded_subjects.subject_ID;
exclude = ismember(analysis_table.subject_ID,bad_subjects);
analysis_table(exclude,:) = [];

end

function [spearman_rho,r_squared,p_value] = do_correlation(effect_type,plot_file,correlation_table)

n_subjects = height(correlation_table);
colors = zeros(n_subjects,3);

switch effect_type
    case 'behavioral'
        dependent_variable = correlation_table.metric;
        independent_variable = correlation_table.power;
        enhanced_behavior = dependent_variable > 0;
        increased_power = independent_variable > 0;
        scenario1 = increased_power & enhanced_behavior;
        scenario2 = increased_power & ~enhanced_behavior;
        scenario3 = ~increased_power & enhanced_behavior;
        scenario4 = ~increased_power & ~enhanced_behavior;
        colors(scenario1,:) = repmat([230 30 100]/255,sum(scenario1),1);
        colors(scenario2,:) = repmat([30 100 230]/255,sum(scenario2),1);
        colors(scenario3,:) = repmat([190 30 190]/255,sum(scenario3),1);
        colors(scenario4,:) = repmat([190 190 30]/255,sum(scenario4),1);
        plot_y_axis = true;
        plot_x_axis = true;
        x_limits = [-3,3];
        x_ticks = -3:1:3;
        max_absolute = max(abs(dependent_variable));
        y_limits = [-1,1]*max_absolute*1.1;
        y_ticks = (-1:0.25:1)*max_absolute*1.1;
    case 'phase'
        plot_file = [plot_file '_phase'];
        dependent_variable = correlation_table.power;
        independent_variable = correlation_table.phase;
        increased_power = dependent_variable > 0;
        scenario1 = increased_power;
        scenario2 = ~increased_power;
        colors(scenario1,:) = repmat([30 230 90]/255,sum(scenario1),1);
        colors(scenario2,:) = repmat([230 130 30]/255,sum(scenario2),1);
        plot_y_axis = false;
        plot_x_axis = true;
        x_limits = [0,0.5];
        x_ticks = 0:0.1:0.5;
        y_limits = [-3,3];
        y_ticks = -3:1:3;
    case 'stim_site'
        plot_file = [plot_file '_stim_site'];
        dependent_variable = correlation_table.phase;
        independent_variable = correlation_table.stim_site_power;
        mapx = makecolormap_EF('uniform1');
        color_indices = round((independent_variable*2000)+3000);
        under1 = color_indices < 1;
        over6000 = color_indices > 6000;
        color_indices(under1) = ones(sum(under1),1);
        color_indices(over6000) = ones(sum(over6000),1)*6000;
        colors = cell2mat(arrayfun(@(x) mapx(x,:),color_indices,'UniformOutput',false));
        plot_y_axis = true;
        plot_x_axis = false;
        x_limits = [-3,3];
        x_ticks = -3:1:3;
        max_y = max(dependent_variable);
        y_limits = [0,1]*max_y*1.1;
        y_ticks = (0:0.25:1)*max_y*1.1;   
end

%%% Plot parameters
figure_width = 137;
figure_height = 85;
text_position = [min(x_limits)+(diff(x_limits)*0.1),min(y_limits)+(diff(y_limits)*0.05)];

%%% Plotting
fig = figure('Units','pixels','Visible','off','Position', [0 0 figure_width figure_height]);
axes('Parent',fig,'Units', 'pixels','Position', [0,0,figure_width,figure_height]);

hold on

scatter(independent_variable,dependent_variable,[],colors,'filled');

if plot_y_axis
    plot([0 0],y_limits,'-k')
end
if plot_x_axis
    plot(x_limits,[0 0],'-k')
end

xlim(x_limits);xticks(x_ticks);xticklabels([]);
ylim(y_limits);yticks(y_ticks);yticklabels([]);

try
    [spearman_rho,p_value] = corr(independent_variable,dependent_variable,'Type','Spearman');
    r_squared = spearman_rho.^2;
    text(text_position(1),text_position(2),sprintf('rho:%.2f r2:%.2f p:%.2f',spearman_rho,r_squared,p_value));
catch
    spearman_rho = NaN;
    r_squared = NaN;
    p_value = NaN;
end

hold off

%%% Saving
print(plot_file,'-dsvg')
close all
end