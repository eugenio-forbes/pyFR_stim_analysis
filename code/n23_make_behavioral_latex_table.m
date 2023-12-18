function n23_make_behavioral_latex_table(varargin)
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
statistics_directory = fullfile(analysis_directory,'statistics');
latex_directory = fullfile(analysis_directory,'latex');
if ~isfolder(latex_directory)
    mkdir(latex_directory);
end

%%% Load table with behavioral within group stimulation effects and between
%%% group interactions
load(fullfile(statistics_directory,'behavioral_stimulation_effects.mat'),'behavioral_stimulation_effects');
load(fullfile(statistics_directory,'behavioral_interactions.mat'),'behavioral_interactions');

%%% Define test variables, periods, groups, hemispheres to be printed
stimulation_hemispheres = {'either','left','right'};
stimulation_groups = {'IPL','retrosplenial'};
stimulation_periods = {'encoding','retrieval'};
behavioral_metrics = {'n_recalls','primacy','total_intrusions','SCF','TCF'};

%%% Make combinations of above in such a way that rows can be printed in
%%% neat hierarchical fashion
n_stimulation_hemispheres = numel(stimulation_hemispheres);
n_stimulation_groups = numel(stimulation_groups);
n_stimulation_periods = numel(stimulation_periods);
n_behavioral_metrics = numel(behavioral_metrics);

[Ax,Bx,Cx,Dx] = ndgrid(1:1:n_stimulation_hemispheres,1:n_stimulation_groups,1:n_stimulation_periods,1:n_behavioral_metrics);
n_combos = length(Ax(:));
behavioral_metric = behavioral_metrics(Dx(:));
stimulation_period = stimulation_periods(Cx(:));
stimulation_group = stimulation_groups(Bx(:));
stimulation_hemisphere = stimulation_hemispheres(Ax(:));

%%% Determine hierarchy branchpoint indices
metric_indices = 1:n_stimulation_hemispheres*n_stimulation_groups*n_stimulation_periods:n_combos;
period_indices = 1:n_stimulation_hemispheres*n_stimulation_groups:n_combos;
group_indices = 1:n_stimulation_hemispheres:n_combos;

%%% Open text file to print latex format table. First few lines to create
%%% table and headers
file_name = fullfile(latex_directory,'stimulation_effects_latex_table.txt');
file_id = fopen(file_name,'w');
fprintf(file_id,'\\begin{table}[h!]\n');
fprintf(file_id,'\\begin{tabular}{|c|c|c|c|c|c|c|}\n');
fprintf(file_id,'\\hline\n');
fprintf(file_id,'Metric & Stimulation Period & Region & Hemisphere & t-stat & p-value & N subjects \\\\ \\hline\n'); % Table headers

%%% Loop through combinations of results
for idx = 1:n_combos
    %%% Get combo info
    this_stimulation_hemisphere = stimulation_hemisphere{idx};
    this_stimulation_group = stimulation_group{idx};
    this_stimulation_period = stimulation_period{idx};
    this_behavioral_metric = behavioral_metric{idx};
    
    %%% Find index of this combo in table
    same_hemisphere = strcmp(behavioral_stimulation_effects.stimulation_hemisphere,this_stimulation_hemisphere);
    same_group = strcmp(behavioral_stimulation_effects.stimulation_group,this_stimulation_group);
    same_period = strcmp(behavioral_stimulation_effects.stimulation_period,this_stimulation_period);
    same_metric = strcmp(behavioral_stimulation_effects.behavioral_metric,this_behavioral_metric);
    same_row = same_hemisphere & same_group & same_period & same_metric;
    
    %%% Get t-stat, p-value, n-subjects
    tstat = behavioral_stimulation_effects.tstat(same_row);
    pvalue = behavioral_stimulation_effects.pvalue(same_row);
    n_subjects = behavioral_stimulation_effects.n_subjects(same_row);
    
    %%% Format combo parameteres based on whether idx is hierarchy
    %%% branchpoint
    if ismember(idx,metric_indices)
        switch this_behavioral_metric
            case 'n_recalls'
                metric = '\# recalls';
            case 'total_intrusions'
                metric = 'intrusions';
            otherwise
                metric = this_behavioral_metric;
        end
    else
        metric = '';
    end
    if ismember(idx,period_indices)
        period = this_stimulation_period;
    else
        period = '';
    end
    if ismember(idx,group_indices)
        if strcmp(this_stimulation_group,'IPL')
            group = 'AG';
        else
            group = 'PCC';
        end
    else
        group = '';
    end
    hemisphere = this_stimulation_hemisphere;
    
    %%% Format table row values in a line
    fprintf(file_id,'%s & %s & %s & %s & %.2f & %.2f & %d \\\\ \\hline\n',...
        metric,period,group,hemisphere,tstat,pvalue,n_subjects);
end
%%% Closing lines
fprintf(file_id, '\\end{tabular}\n');
fprintf(file_id, '\\end{table}\n');
fclose(file_id);


%%% Loop through combinations of periods and metrics to print interactions
%%% result table in latex format

[Ax,Bx] = ndgrid(1:n_stimulation_periods,1:n_behavioral_metrics);
n_combos = length(Ax(:));
behavioral_metric = behavioral_metrics(Bx(:));
stimulation_period = stimulation_periods(Ax(:));

%%% Get hierarchy branchpoints to print rows in hierarchical fashion
metric_indices = 1:n_stimulation_periods:n_combos;

%%% Opening lines  to create table and headers
file_name = fullfile(latex_directory,'interactions_latex_table.txt');
file_id = fopen(file_name,'w');
fprintf(file_id,'\\begin{table}[h!]\n');
fprintf(file_id,'\\begin{tabular}{|c|c|c|c|c|c|}\n');
fprintf(file_id,'\\hline\n');
fprintf(file_id,'Metric & Stimulation Period & t-stat & p-value & N AG & N PCC\\\\ \\hline\n'); % Table headers

%%% Loop through combinations
for idx = 1:n_combos
    %%% Get combo information
    this_stimulation_period = stimulation_period{idx};
    this_behavioral_metric = behavioral_metric{idx};
    
    %%% Get index with results for this combo
    same_period = strcmp(behavioral_interactions.stimulation_period,this_stimulation_period);
    same_metric = strcmp(behavioral_interactions.behavioral_metric,this_behavioral_metric);
    same_row = same_period & same_metric;
    
    %%% Get t-stat, p-value, N group 1, N group 2
    tstat = behavioral_interactions.tstat(same_row);
    pvalue = behavioral_interactions.pvalue(same_row);
    n_IPL = behavioral_interactions.n_subjects_IPL(same_row);
    n_RS = behavioral_interactions.n_subjects_RS(same_row);
    
    %%% Format based on row index and metric
    if ismember(idx,metric_indices)
        switch this_behavioral_metric
            case 'n_recalls'
                metric = '\# recalls';
            case 'total_intrusions'
                metric = 'intrusions';
            otherwise
                metric = this_behavioral_metric;
        end
    else
        metric = '';
    end
    period = this_stimulation_period;
    
    %%% Format row values in a line
    fprintf(file_id,'%s & %s & %.2f & %.2f & %d & %d \\\\ \\hline\n',...
        metric,period,tstat,pvalue,n_IPL,n_RS);
end
%%% Closing lines
fprintf(file_id, '\\end{tabular}\n');
fprintf(file_id, '\\end{table}\n');
fclose(file_id);

%%% Extra table comparing basline performance between groups with the same
%%% number of combinations as loop above.
%%% Opening line to create table and header.
file_name = fullfile(latex_directory,'baseline_latex_table.txt');
file_id = fopen(file_name,'w');
fprintf(file_id,'\\begin{table}[h!]\n');
fprintf(file_id,'\\begin{tabular}{|c|c|c|c|c|c|}\n');
fprintf(file_id,'\\hline\n');
fprintf(file_id,'Metric & Stimulation Period & t-stat & p-value & N AG & N PCC\\\\ \\hline\n'); % Table headers

%%% Loop through combinations
for idx = 1:n_combos
    %%% Get combo information
    this_stimulation_period = stimulation_period{idx};
    this_behavioral_metric = behavioral_metric{idx};
    
    %%% Get index with results for this combo
    same_period = strcmp(behavioral_interactions.stimulation_period,this_stimulation_period);
    same_metric = strcmp(behavioral_interactions.behavioral_metric,this_behavioral_metric);
    same_row = same_period & same_metric;
    
    %%% Get t-stat, p-value, N group 1, and N group 2
    tstat = behavioral_interactions.tstat_control(same_row);
    pvalue = behavioral_interactions.pvalue_control(same_row);
    n_IPL = behavioral_interactions.n_subjects_IPL(same_row);
    n_RS = behavioral_interactions.n_subjects_RS(same_row);
    
    %%% Format based on metric and index
    if ismember(idx,metric_indices)
        switch this_behavioral_metric
            case 'n_recalls'
                metric = '\# recalls';
            case 'total_intrusions'
                metric = 'intrusions';
            otherwise
                metric = this_behavioral_metric;
        end
    else
        metric = '';
    end
    period = this_stimulation_period;
    
    %%% Format row values in a line
    fprintf(file_id,'%s & %s & %.2f & %.2f & %d & %d \\\\ \\hline\n',...
        metric,period,tstat,pvalue,n_IPL,n_RS);
end
%%% Closing lines
fprintf(file_id, '\\end{tabular}\n');
fprintf(file_id, '\\end{table}\n');
fclose(file_id);
end