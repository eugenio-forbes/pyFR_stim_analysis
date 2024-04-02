function n20_anatomical_variation_analysis(varargin)
if isempty(varargin)
    %%% Directory information
    root_directory = '/directory/to/pyFR_stim_analysis';
    parpool_n = 16;
else
    root_directory = varargin{1};
    parpool_n = varargin{2};
end

%%% List directories
analysis_directory = fullfile(root_directory,username,analysis_folder_name);
table_directory = fullfile(analysis_directory,'tables');
list_directory = fullfile(analysis_directory,'lists');
plots_directory = fullfile(analysis_directory,'plots');
lock_directory = fullfile(analysis_directory,'locks');
lme_directory = fullfile(analysis_directory,'lme_results');
resources_directory = fullfile(analysis_directory,'resources');

%%% Load electrode list to get coordinates
load(fullfile(list_directory,'electrode_list.mat'),'electrode_list');

%%% Declare regions of interest and periods of interest
stimulation_groups = {'IPL','retrosplenial'};
periods = {'encoding';'retrieval'};
[Ax,Bx] = ndgrid(1:numel(periods),1:numel(stimulation_groups));
n_combos = length(Ax(:));
stimulation_groups = stimulation_groups(Bx(:));
periods = periods(Ax(:));

%%%Initialize parpool
pool_object = gcp('nocreate');
if isempty(pool_object)
    parpool(parpool_n)
end

for idx = 1:n_combos
    period = periods{idx};
    stimulation_group = stimulation_groups{idx};
    
    folder_name = 'anatomical_variation_analysis_hippocampus3';
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
    
    file_name = sprintf('%s_%s',period,stimulation_group);
    lock_file = fullfile(this_lock_directory,[file_name '_lock.txt']);
    done_file = fullfile(this_lock_directory,[file_name '_done.txt']);
    error_file = fullfile(this_lock_directory,[file_name '_error.txt']);
    mat_file = fullfile(this_lme_directory,[file_name '.mat']);
    plot_file = fullfile(this_plot_directory,file_name);
    
    pause(rand*5)
    
    if ~isfile(lock_file)
        fid = fopen(lock_file,'w'); fclose(fid);
        try
            electrode_table_file = sprintf('%s_band_table',period);
            electrode_table = load(fullfile(table_directory,[electrode_table_file '.mat']));
            electrode_table = electrode_table.(electrode_table_file);
            
%             electrode_table = check_previous_exclusions(analysis_directory,electrode_table);
            is_IPL = logical(electrode_table.stimulation_lateral);
            if strcmp(stimulation_group,'IPL')
                same_group = is_IPL;
            else
                same_group = ~is_IPL;
            end
            electrode_table = electrode_table(same_group,:);
            n_subjects = length(unique(electrode_table.subject_ID));
            n_electrodes = length(unique(electrode_table.electrode_ID));
            [statistics_table,tstat_box] = get_anatomical_variation_analyses(plot_file,electrode_table);
            save(mat_file,'statistics_table','n_subjects','n_electrodes');
%             do_brain_plots(resources_directory,plot_file,electrode_table,electrode_list);
%             do_box_plots(plot_file,tstat_box,statistics_table);

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

function electrode_table = check_previous_exclusions(analysis_directory,electrode_table)

exclusion_directory = fullfile(analysis_directory,'exclusion_lists');

%%% Load subject, session, electrode lists and exclusion lists
load(fullfile(exclusion_directory,'excluded_electrodes.mat'),'excluded_electrodes')

bad_electrodes = excluded_electrodes.electrode_ID;
exclude = ismember(electrode_table.electrode_ID,bad_electrodes);
electrode_table(exclude,:) = [];

end

function [statistics_table,tstat_box] = get_anatomical_variation_analyses(plot_file,electrode_table)
n_rows = height(electrode_table);

electrode_table.subject_ID = categorical(electrode_table.subject_ID);
electrode_table.electrode_ID = categorical(electrode_table.electrode_ID);
electrode_table.stimulation_left = logical(electrode_table.stimulation_left);

electrode_table.stimulation_power = mean(double(vertcat(electrode_table.slow_theta_stimulation_power{:}))/1000,2);
electrode_table.control_power = mean(double(vertcat(electrode_table.slow_theta_control_power{:}))/1000,2);
electrode_table.t_statistic = mean(double(vertcat(electrode_table.slow_theta_tstats{:}))/1000,2);

stimulation_effect_table = repmat(electrode_table,2,1);
stimulation_effect_table.condition = [true(n_rows,1);false(n_rows,1)];
stimulation_effect_table.power = [electrode_table.stimulation_power;electrode_table.control_power];

is_anterior = logical(electrode_table.anterior);
is_posterior = ~is_anterior;
is_left = logical(electrode_table.left);
is_right = ~is_left;
n_box_rows = max([sum(is_anterior),sum(is_posterior),sum(is_left),sum(is_right)]);
tstat_box = NaN(n_box_rows,4);

region = {'anterior';'posterior';'right';'left';'anterior-posterior';'left-right'};
n_regions = length(region);
t_statistic_lme = NaN(n_regions,1);
p_value_lme = NaN(n_regions,1);
n_subjects = NaN(n_regions,1);
n_electrodes = NaN(n_regions,1);
degrees_of_freedom = NaN(n_regions,1);
F_stat = NaN(n_regions,1);
p_F_stat = NaN(n_regions,1);
DF_numerator = NaN(n_regions,1);
DF_denominator = NaN(n_regions,1);
f_positive = NaN(n_regions,1);
p_binomial = NaN(n_regions,1);

for idx = 1:n_regions
    this_region = region{idx};
    
    switch this_region
        case 'anterior'
            this_electrode_table = electrode_table(is_anterior,:);
            this_effect_table = stimulation_effect_table(repmat(is_anterior,2,1),:);
        case 'posterior'
            this_electrode_table = electrode_table(is_posterior,:);
            this_effect_table = stimulation_effect_table(repmat(is_posterior,2,1),:);
        case 'left'
            this_electrode_table = electrode_table(is_left,:);
            this_effect_table = stimulation_effect_table(repmat(is_left,2,1),:);
        case 'right'
            this_electrode_table = electrode_table(is_right,:);
            this_effect_table = stimulation_effect_table(repmat(is_right,2,1),:);
        case {'anterior-posterior','left-right'}
            this_electrode_table = electrode_table;
            this_effect_table = stimulation_effect_table;
    end
    
    random_effect = '1 + condition';
    switch this_region
        case {'anterior','posterior','left','right'}
            formula = sprintf('power ~ condition + (%s|subject_ID) + (%s|subject_ID:electrode_ID) + (%s|subject_ID:stimulation_left) + (%s|stimulation_lateral:stimulation_left)',...
                random_effect,random_effect,random_effect,random_effect);
            statistics_index = 2;
        case {'anterior-posterior'}
            formula = sprintf('power ~ condition*anterior + (%s|subject_ID) + (%s|subject_ID:electrode_ID) + (%s|subject_ID:stimulation_left) + (%s|stimulation_lateral:stimulation_left)',...
                random_effect,random_effect,random_effect,random_effect);
            statistics_index = 4;
        case {'left-right'}
            formula = sprintf('power ~ condition*left + (%s|subject_ID) + (%s|subject_ID:electrode_ID) + (%s|subject_ID:stimulation_left) + (%s|stimulation_lateral:stimulation_left)',...
                random_effect,random_effect,random_effect,random_effect);
            statistics_index = 4;
    end
    
    n_subjects(idx) = length(unique(this_electrode_table.subject_ID));
    n_electrodes(idx) = length(unique(this_electrode_table.electrode_ID));
    n_rows = height(this_electrode_table);
    
    model = fitlme(this_effect_table,formula);
    [betas,beta_names,stats] = fixedEffects(model,'DFMethod','satterthwaite','alpha',0.05);
    H = ones(1,height(beta_names));
    H(1) = 0;
    C = H*betas;
    
    t_statistic_lme(idx) = stats.tStat(statistics_index);
    p_value_lme(idx) = stats.pValue(statistics_index);
    degrees_of_freedom(idx) = stats.DF(statistics_index);
    [p_F_stat(idx),F_stat(idx),DF_numerator(idx),DF_denominator(idx)] = coefTest(model,H,C,'DFMethod','satterthwaite');
    
    these_t_statistics = this_electrode_table.t_statistic;
    if idx < 5
        tstat_box(1:n_rows,idx) = these_t_statistics;
    end
    
%     [f_positive(idx),p_binomial(idx)] = do_binomial_and_piechart(plot_file,these_t_statistics,this_region);
 
end

statistics_table = table(region,t_statistic_lme,p_value_lme,degrees_of_freedom,n_subjects,n_electrodes,f_positive,p_binomial);

end

function [f_positive,p_binomial] = do_binomial_and_piechart(plot_file,t_statistics,region)

switch region
    case {'anterior','posterior','left','right'}
        figure_width = 20;
        figure_height = 20;
    case {'anterior-posterior','left-right'}
        figure_width = 36;
        figure_height = 36;
end

piechart_position = [0 0 figure_width figure_height];

probability_if_random = 0.5;

this_plot_file = [plot_file,sprintf('_%s_piechart',region)];
n_electrodes = length(t_statistics);
n_positive = sum(t_statistics > 0);
f_positive = n_positive/n_electrodes;
p_binomial = myBinomTest(n_positive,n_electrodes,probability_if_random,'two');

pie_input = [n_positive;n_electrodes-n_positive];
explode = [1;0];
labels = {'',''};

figure('Units','pixels','Visible','off','Position',piechart_position)
axes_handle = axes;
set(axes_handle,'Position',[0 0 1 1])
pie(pie_input,explode,labels);
pause(2)
print(this_plot_file,'-dsvg')
close all

end

function do_brain_plots(resources_directory,plot_file,electrode_table,electrode_list)
this_plot_file = [plot_file '_brain_plot'];
figure_width = 200;
figure_height = 200;
plot_position = [0 0 figure_width figure_height];
marker_size = 40;

mapx = makecolormap_EF('uniform1');

load(fullfile(resources_directory,'EC.mat'),'EC');
load(fullfile(resources_directory,'SURF.mat'),'surf');

list_indices = arrayfun(@(x) find(electrode_list.electrode_ID == x,1,'first'),electrode_table.electrode_ID);
electrode_table.coordinates = electrode_list.coordinates(list_indices);
no_coordinates = cellfun(@isempty,electrode_table.coordinates);
electrode_table(no_coordinates,:) = [];
n_electrodes = height(electrode_table);

t_statistics = mean(double(vertcat(electrode_table.slow_theta_tstats{:}))/1000,2);
[~,sorting_indices] = sortrows(abs(t_statistics),'ascend');
t_statistics = t_statistics(sorting_indices);
electrode_table = electrode_table(sorting_indices,:);

color_indices = round((t_statistics*1000)+3000);
under1 = color_indices < 1;
over6000 = color_indices > 6000;
color_indices(under1) = ones(sum(under1),1);
color_indices(over6000) = ones(sum(over6000),1)*6000;
colors = cell2mat(arrayfun(@(x) mapx(x,:),color_indices,'UniformOutput',false));

coordinates = vertcat(electrode_table.coordinates{:});
x = coordinates(:,1);
y = coordinates(:,2);
z = linspace(-100,-101,n_electrodes)';

figure('Units','pixels','Visible','off','Position',plot_position)

%%% 3D plot of a sample brain
hold on
axes_brain = trisurf(surf.tri,surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),'EdgeColor','none');
whitebg(gcf,EC.bak.color);
set(gcf,'Color',EC.bak.color,'InvertHardcopy','off');
eval(['material ',EC.glb.material,';']);
eval(['shading ',EC.glb.shading,';']);
set(axes_brain,'FaceColor',[0.95, 0.95, 0.95]);
set(axes_brain,'FaceAlpha',0.15);
daspect([1 1 1])
hold off

%%% 3D scatter plot of electrodes with colors corresponding to t-statistics
hold on
scatter3(x,y,z,marker_size,colors,'filled','Marker','o')
rotate3d on
axis tight 
axis vis3d off
view(180,-90)
lighting phong
camlight RIGHT
hold off

print(this_plot_file,'-dsvg')
close all
end

function do_box_plots(plot_file,tstat_box,statistics_table)
%%% Plot parameters
this_plot_file = [plot_file '_box_plot'];
figure_width = 164;
figure_height = 129;
box_position = [0 0 figure_width figure_height];
n_boxes = size(tstat_box,2);
x_limits = [0.5 n_boxes+0.5];
x_ticks = 1:n_boxes;
y_limits = [-3.5 3.5];
y_ticks = -3:1:3;
significance_y = 2.5;
point_colors = [[247,89,10];...
    [247,197,10];...
    [13,166,255];...
    [122,13,255]]/255;
box_colors = repmat([0 0 0],n_boxes,1);
marker_size = 7;

p_values = statistics_table.p_value_lme;

fig = figure('Units','pixels','Visible','off','Position',box_position);
axes('Parent',fig,'Units','pixels','Position',box_position);
hold on

for idx = 1:n_boxes+2
    
    switch num2str(idx)
        case {'1','2','3','4'}
            not_nan = ~isnan(tstat_box(:,idx));
            x_axis_scatter = normrnd(idx,0.05,sum(not_nan),1);
            t_statistics = tstat_box(not_nan,idx);
            scatter(x_axis_scatter,t_statistics,marker_size,point_colors(idx,:),'filled','Marker','o');
            significance_x = idx;
        case '5'
            significance_x = 1.5;
        case '6'
            significance_x = 3.5;
    end
    
    p_value = p_values(idx);
    
    if p_value < 0.01
        plot(significance_x-0.1,significance_y,'*k')
        plot(significance_x+0.1,significance_y,'*k')
    elseif p_value < 0.05
        plot(significance_x,significance_y,'*k')
    elseif p_value < 0.1
        plot(significance_x,significance_y,'+k')
    end
    
end

box_handle = boxplot(tstat_box,x_ticks,'Colors',box_colors,'Symbol','','DataLim',[-100 100]);
set(box_handle,{'linew'},{1.5})
plot(x_limits,[0 0],'-k')

xlim(x_limits); xticks(x_ticks); xticklabels([]);
ylim(y_limits); yticks(y_ticks); yticklabels([]);

hold off

print(this_plot_file,'-dsvg')
close all
end
