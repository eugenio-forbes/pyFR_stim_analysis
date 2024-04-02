function n10_behavioral_analysis(varargin)
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
resources_directory = fullfile(root_directory,'resources');
data_directory = fullfile(root_directory,'data');
statistics_directory = fullfile(root_directory,'statistics');
plots_directory = fullfile(root_directory,'plots/behavioral');
if ~isfolder(statistics_directory)
    mkdir(statistics_directory);
end
if ~isfolder(plots_directory)
    mkdir(plots_directory);
end

%%% Load subject list
load(fullfile(list_directory,'subject_list.mat'),'subject_list');

%%% Load session list
load(fullfile(list_directory,'session_list.mat'),'session_list');

%%% Analysis is limited to sessions that included a control condition (no
%%% stimulation). Exclude sessions without control
has_no_stimulation = session_list.has_no_stimulation;
session_list = session_list(has_no_stimulation,:);
subject_list = subject_list(ismember(subject_list.subject_ID,session_list.subject_ID),:);

%%% Load latent semantic analysis matrix
load(fullfile(resources_directory,'LSAmat.mat'),'LSAmat');
LSAmat = LSAmat * 1;

%%% Create table to loop through all combinations of conditions and
%%% session_IDs
conditions = [1;3;5];
session_IDs = session_list.session_ID;
[Ax,Bx] = ndgrid(1:numel(session_IDs),1:numel(conditions));
condition = conditions(Bx(:));
session_ID = session_IDs(Ax(:));
n_combos = length(Ax(:));

%%% Initialize cell arrays to hold table rows for each combination
behavioral_array = cell(n_combos,1);
n_recalls_array = cell(n_combos,1);
primacy_array = cell(n_combos,1);

%%% Initialize parpool
pool_object = gcp('nocreate');
if isempty(pool_object)
    parpool(parpool_n);
end

%%% Loop through sessions to gather behavioral data from events
for idx = 1:n_combos
    %%% Get combination session_ID and condition. Get row corresponding to
    %%% session
    this_session_ID = session_ID(idx);
    this_condition = condition(idx);
    this_session = session_list(session_list.session_ID == this_session_ID,:);
    
    %%% Check whether the session had corresponding condition before
    %%% continuing to extract behavioral metrics
    has_condition = false;
    if this_condition == 1 && this_session.has_encoding_stimulation
        has_condition = true;
    elseif this_condition == 3 && this_session.has_retrieval_stimulation
        has_condition = true;
    elseif this_condition == 5
        has_condition = true;
    end
    
    if has_condition
        [behavioral_array{idx},n_recalls_array{idx},primacy_array{idx}] = ...
            get_behavioral_metrics(data_directory,this_session,this_condition,LSAmat);
    end 
end
clear session_list combo LSAmat

%%% Concatenate individual session tables to be used in LME analysis
session_behavioral_table = vertcat(behavioral_array{:}); clear behavioral_array
n_recalls_table = vertcat(n_recalls_array{:}); clear n_recalls_array
primacy_table = vertcat(primacy_array{:}); clear primacy_array

%%% Save tables for later access
save(fullfile(statistics_directory,'session_behavioral_table.mat'),'session_behavioral_table');
save(fullfile(statistics_directory,'n_recalls_table.mat'),'n_recalls_table');
save(fullfile(statistics_directory,'primacy_table.mat'),'primacy_table');

%%% Get behavioral data for each subject averaged across sessions (weighted
%%% by number of lists in each session).
subject_behavioral_table = average_sessions(session_behavioral_table,subject_list,conditions);
save(fullfile(statistics_directory,'subject_behavioral_table.mat'),'subject_behavioral_table');

%%% Get table averaging encoding and retrieval stimulation data (weighted
%%% by number of lists in each condition)
any_stimulation_table = average_stimulation_conditions(subject_behavioral_table,subject_list);
save(fullfile(statistics_directory,'any_stimulation_table.mat'),'any_stimulation_table');

%%% Declare stimulation periods, regions, hemispheres; and behavioral
%%% metrics for LME analyses
stimulation_periods = {'encoding','retrieval','either'};
hemispheres = {'either','left','right'};
stimulation_groups = {'IPL','retrosplenial'};
variables = {'total_intrusions','TCF','SCF','overall_percent_recall','total_PLIs','total_XLIs',...
    'n_recalls','primacy'};

%%% Create columns for all combinations of parameters and preallocate
%%% cell arrays to hold statistics of stimulation effects on behavior
[Ax,Bx,Cx,Dx] = ndgrid(1:numel(stimulation_periods),1:numel(hemispheres),1:numel(stimulation_groups),1:numel(variables));
n_combos = length(Ax(:));
behavioral_metric = variables(Dx(:))';
stimulation_group = stimulation_groups(Cx(:))';
stimulation_hemisphere = hemispheres(Bx(:))';
stimulation_period = stimulation_periods(Ax(:))';
n_subjects = NaN(n_combos,1);
tstat = NaN(n_combos,1);
pvalue = NaN(n_combos,1);
df = NaN(n_combos,1);
F_stat = NaN(n_combos,1);
F_pval = NaN(n_combos,1);

%%% Loop through all combinations
for idx = 1:n_combos
    this_behavioral_metric = behavioral_metric{idx};
    %%% Based on behavioral metric select appropriate table, formula and
    %%% index to retrieve statistics from LME models.
    switch this_behavioral_metric
        case {'overall_percent_recall','total_PLIs','total_XLIs','total_intrusions','TCF','SCF'}
            analysis_table = session_behavioral_table;
            random_effect = '1 + condition';
            formula = sprintf('%s ~ condition + (%s|subject_ID) + (%s|subject_ID:session_ID)',...
                this_behavioral_metric,random_effect,random_effect);
            statistics_index = 2;
        case {'n_recalls'}
            analysis_table = n_recalls_table;
            random_effect = '1 + condition';
            formula = sprintf('%s ~ condition + (%s|subject_ID) + (%s|subject_ID:session_ID)',...
                this_behavioral_metric,random_effect,random_effect);
            statistics_index = 2;
        case {'primacy'}
            analysis_table = primacy_table;
            random_effect = '1 + condition + primacy + condition*primacy';          
            formula = sprintf('f_recall ~ condition*primacy + (%s|subject_ID) + (%s|subject_ID:session_ID)',...
                random_effect,random_effect);
            statistics_index = 4;
    end
    
    %%% Select columns of analysis table that have the stimulation region
    %%% and hemisphere of interest.
    this_stimulation_group = stimulation_group{idx};
    same_region = strcmp(analysis_table.stimulation_group,this_stimulation_group);
    analysis_table = analysis_table(same_region,:);
    
    %%% If analysis pools together both stimulation hemispheres for a given
    %%% region, include random effect for stimulation hemisphere
    this_stimulation_hemisphere = stimulation_hemisphere{idx};
    switch this_stimulation_hemisphere
        case {'left','right'}
            same_hemisphere = strcmp(analysis_table.stimulation_hemisphere,this_stimulation_hemisphere);
            analysis_table = analysis_table(same_hemisphere,:);
            final_formula = formula;
        case {'either'}
            extra_random_effect = sprintf(' + (%s|subject_ID:stimulation_hemisphere)',random_effect);
            final_formula = [formula extra_random_effect];
    end
    
    %%% If analysis pools together stimulation conditions use all data and
    %%% convert condition so that conditions for stimulation (1,3) are true
    %%% and false for control (5). If the analysis focuses on one type of
    %%% stimulation only include subjects with specific condition and
    %%% control and exclude other stimulation condition data.
    this_stimulation_period = stimulation_period{idx};
    switch this_stimulation_period
        case {'either'}
            new_condition = ismember(analysis_table.condition,[1,3]);
            analysis_table.condition = new_condition;
        case {'encoding'}
            excluded_condition = analysis_table.condition == 3;
            analysis_table(excluded_condition,:) = [];
            has_condition = logical(analysis_table.has_encoding_stimulation);
            analysis_table = analysis_table(has_condition,:);
            new_condition = analysis_table.condition == 1;
            analysis_table.condition = new_condition;
        case {'retrieval'}
            excluded_condition = analysis_table.condition == 1;
            analysis_table(excluded_condition,:) = [];
            has_condition = logical(analysis_table.has_retrieval_stimulation);
            analysis_table = analysis_table(has_condition,:);
            new_condition = analysis_table.condition == 3;
            analysis_table.condition = new_condition;
    end
    
    %%% Convert random effect variables to categorical type. Condition has
    %%% already been converted to logical.
    n_subjects(idx) = length(unique(analysis_table.subject_ID));
    analysis_table.subject_ID = categorical(analysis_table.subject_ID);
    analysis_table.session_ID = categorical(analysis_table.session_ID);
    analysis_table.stimulation_hemisphere = categorical(analysis_table.stimulation_hemisphere);
    
    %%% Perform linear mixed effects analysis with fully processed table
    %%% and formula. Gather tstatistic, p-value, and degrees of freedom
    %%% from effect.
    if ~isempty(analysis_table)
        lme_model = fitlme(analysis_table,final_formula);
        [betas,beta_names,stats] = fixedEffects(lme_model,'DFMethod','satterthwaite','alpha',0.05);
        H = ones(1,height(beta_names));
        H(1) = 0;
        C = H*betas;
        tstat(idx) = stats.tStat(statistics_index);
        pvalue(idx) = stats.pValue(statistics_index);
        df(idx) = stats.DF(statistics_index);
        [F_pval(idx),F_stat(idx),~,~] = coefTest(lme_model,H,C,'DFMethod','satterthwaite');
    end
end

%%% Make table for stimulation effects with preallocated columns and save
behavioral_stimulation_effects = table(behavioral_metric,stimulation_group,...
    stimulation_hemisphere,stimulation_period,n_subjects,tstat,pvalue,df,F_stat,F_pval);
save(fullfile(statistics_directory,'behavioral_stimulation_effects.mat'),'behavioral_stimulation_effects');
clear behavioral metric stimulation_region stimulation_hemisphere 
clear stimulation_period n_subjects tstat pvalue df

%%% Create columns for combinations of stimulation period and behavioral
%%% metrics to get interaction of stimulation group (pooling stimulation
%%% hemispheres in all cases)
[Ax,Bx] = ndgrid(1:numel(stimulation_periods),1:numel(variables));
n_combos = length(Ax(:));
behavioral_metric = variables(Bx(:))';
stimulation_period = stimulation_periods(Ax(:))';
n_subjects_IPL = NaN(n_combos,1);
n_subjects_RS = NaN(n_combos,1);
tstat = NaN(n_combos,1);
pvalue = NaN(n_combos,1);
df = NaN(n_combos,1);
tstat_control = NaN(n_combos,1);
pvalue_control = NaN(n_combos,1);
df_control = NaN(n_combos,1);
F_stat_control = NaN(n_combos,1);
F_pval_control = NaN(n_combos,1);
control_index = 2;

%%% Loop through all combinations
parfor idx = 1:n_combos
    this_behavioral_metric = behavioral_metric{idx};
    %%% Based on behavioral metric select appropriate table, formula and
    %%% index to retrieve statistics from LME models.
    switch this_behavioral_metric
        case {'overall_percent_recall','total_PLIs','total_XLIs','total_intrusions','TCF','SCF'}
            interaction_analysis_table = session_behavioral_table;
            interaction_random_effect = '1 + condition';
            interaction_formula = sprintf('%s ~ condition*stimulation_group + (%s|subject_ID) + (%s|subject_ID:session_ID)',...
                this_behavioral_metric,interaction_random_effect,interaction_random_effect);
            interaction_statistics_index = 4;
        case {'n_recalls'}
            interaction_analysis_table = n_recalls_table;
            interaction_random_effect = '1 + condition';
            interaction_formula = sprintf('%s ~ condition*stimulation_group + (%s|subject_ID) + (%s|subject_ID:session_ID)',...
                this_behavioral_metric,interaction_random_effect,interaction_random_effect);
            interaction_statistics_index = 4;
        case {'primacy'}
            interaction_analysis_table = primacy_table;
            interaction_random_effect = '1 + condition + primacy + condition*primacy';          
            interaction_formula = sprintf('f_recall ~ condition*primacy*stimulation_group + (%s|subject_ID) + (%s|subject_ID:session_ID)',...
                interaction_random_effect,interaction_random_effect);
            interaction_statistics_index = 8;
    end
    
    %%% If analysis pools together stimulation conditions use all data and
    %%% convert condition so that conditions for stimulation (1,3) are true
    %%% and false for control (5). If the analysis focuses on one type of
    %%% stimulation only include subjects with specific condition and
    %%% control and exclude other stimulation condition data.
    this_stimulation_period = stimulation_period{idx};
    switch this_stimulation_period
        case {'either'}
            new_condition = ismember(interaction_analysis_table.condition,[1,3]);
            interaction_analysis_table.condition = new_condition;
        case {'encoding'}
            excluded_condition = interaction_analysis_table.condition == 3;
            interaction_analysis_table(excluded_condition,:) = [];
            has_condition = logical(interaction_analysis_table.has_encoding_stimulation);
            interaction_analysis_table = interaction_analysis_table(has_condition,:);
            new_condition = interaction_analysis_table.condition == 1;
            interaction_analysis_table.condition = new_condition;
        case {'retrieval'}
            excluded_condition = interaction_analysis_table.condition == 1;
            interaction_analysis_table(excluded_condition,:) = [];
            has_condition = logical(interaction_analysis_table.has_retrieval_stimulation);
            interaction_analysis_table = interaction_analysis_table(has_condition,:);
            new_condition = interaction_analysis_table.condition == 3;
            interaction_analysis_table.condition = new_condition;
    end
    
    %%% Convert stimulation_group to logical where IPL is 1 and
    %%% retrosplenial is 0 so that positive statistic indicates that the
    %%% effect was more positive for IPL
    is_IPL = strcmp(interaction_analysis_table.stimulation_group,'IPL');
    interaction_analysis_table.stimulation_group = is_IPL;
    
    %%% Add random effect for stimulation hemisphere
    extra_random_effect = sprintf(' + (%s|subject_ID:stimulation_hemisphere) + (%s|stimulation_group:stimulation_hemisphere)',interaction_random_effect,interaction_random_effect);
    final_formula = [interaction_formula extra_random_effect];
    
    n_subjects_IPL(idx) = length(unique(interaction_analysis_table.subject_ID(is_IPL)));
    n_subjects_RS(idx) = length(unique(interaction_analysis_table.subject_ID(~is_IPL)));
    
    %%% Convert random effect variables to categorical type. Condition and
    %%% stimulation_group have already been converted to logical.
    interaction_analysis_table.subject_ID = categorical(interaction_analysis_table.subject_ID);
    interaction_analysis_table.session_ID = categorical(interaction_analysis_table.session_ID);
    interaction_analysis_table.stimulation_hemisphere = categorical(interaction_analysis_table.stimulation_hemisphere);
    
    %%% Perform linear mixed effects analysis with fully processed table
    %%% and formula. Gather tstatistic, p-value, and degrees of freedom
    %%% from effect.
    if ~isempty(interaction_analysis_table)
        lme_model = fitlme(interaction_analysis_table,final_formula);
        [betas,beta_names,stats] = fixedEffects(lme_model,'DFMethod','satterthwaite','alpha',0.05);
        H = ones(1,height(beta_names));
        H(1) = 0;
        C = H*betas;
        tstat(idx) = stats.tStat(interaction_statistics_index);
        pvalue(idx) = stats.pValue(interaction_statistics_index);
        df(idx) = stats.DF(interaction_statistics_index);
        tstat_control(idx) = stats.tStat(control_index);
        pvalue_control(idx) = stats.pValue(control_index);
        df_control(idx) = stats.DF(control_index);
        [F_pval_control(idx),F_stat_control(idx),~,~] = coefTest(lme_model,H,C,'DFMethod','satterthwaite');
    end
end

%%% Make table for stimulation intereactions with preallocated columns and save
behavioral_interactions = table(behavioral_metric,stimulation_period,...
    n_subjects_IPL,n_subjects_RS,tstat,pvalue,df,tstat_control,pvalue_control,df_control,F_stat_control,F_pval_control);
save(fullfile(statistics_directory,'behavioral_interactions.mat'),'behavioral_interactions');
clear behavioral_metric stimulation_period n_subjects_IPL n_subjects_RS tstat pvalue df
clear n_recalls_table primacy_table

%%% Declare behavioral metrics plotted with boxplots
box_plot_metrics = {'overall_percent_recall','total_PLIs','total_XLIs',...
    'total_intrusions','TCF','SCF','tstat_recalls','zvalue_recalls'};
plotting_options = {'merged','separate'};
[Ax,Bx] = ndgrid(1:numel(box_plot_metrics),1:numel(plotting_options));
n_combos = length(Ax(:));
plotting_option = plotting_options(Bx(:));
box_plot_metric = box_plot_metrics(Ax(:));
parfor idx = 1:n_combos
    this_plotting_option = plotting_option{idx};
    this_box_plot_metric = box_plot_metric{idx};
    switch this_plotting_option
        case 'separate'
            behavioral_boxplot(plots_directory,this_box_plot_metric,this_plotting_option,subject_behavioral_table,behavioral_stimulation_effects,behavioral_interactions)
        case 'merged'
            behavioral_boxplot(plots_directory,this_box_plot_metric,this_plotting_option,any_stimulation_table,behavioral_stimulation_effects,behavioral_interactions)
    end
end

%%% Plot serial position and lag-CRP curves
[Ax,Bx] = ndgrid(1:numel(stimulation_groups),1:numel(plotting_options));
n_combos = length(Ax(:));
plotting_option = plotting_options(Bx(:));
stimulation_group = stimulation_groups(Ax(:));
parfor idx = 1:n_combos
    this_plotting_option = plotting_option{idx};
    this_stimulation_group = stimulation_group{idx};
    switch this_plotting_option
        case 'separate'
            serial_position_curves(plots_directory,this_stimulation_group,this_plotting_option,subject_behavioral_table,behavioral_stimulation_effects)
            lag_crp_curves(plots_directory,this_stimulation_group,this_plotting_option,subject_behavioral_table)
        case 'merged'
            serial_position_curves(plots_directory,this_stimulation_group,this_plotting_option,any_stimulation_table,behavioral_stimulation_effects)
            lag_crp_curves(plots_directory,this_stimulation_group,this_plotting_option,any_stimulation_table)
    end
end


end

function [behavioral_table,n_recalls_table,primacy_table] = ...
        get_behavioral_metrics(data_directory,session_info,condition,LSAmat)
%%% Define serial positions
n_list_items = 10;
serial_positions = 1:n_list_items;

%%% Get session info to load events file
subject = session_info.subject{:};
task = session_info.task{:};
session = session_info.session{:};

events_file = fullfile(data_directory,subject,task,session,'events.mat');
load(events_file,'events');

events = struct2table(events);

%%% Get encoding and retrieval events of the corresponding condition.
%%% Get list numbers of corresponding lists and number of lists.
this_condition = events.stimcode == condition;    
condition_events = events(this_condition,:);
encoding_type = strcmp(condition_events.type,'WORD');
retrieval_type = strcmp(condition_events.type,'REC_WORD') & ~strcmp(condition_events.item,'<>');
encoding_events = condition_events(encoding_type,:);
retrieval_events = condition_events(retrieval_type,:);
condition_lists = unique(encoding_events.list);
n_lists = length(condition_lists);

%%% Extract behavioral metrics from encoding events:

%%% -Overall percent recall aggregating lists of the condition
overall_percent_recall = sum(encoding_events.recalled)/height(encoding_events);
%%% -Number of recalls made on each list of the condition
n_recalls_by_list = arrayfun(@(x) sum(encoding_events.recalled(encoding_events.list == x)),condition_lists);
%%% -Percent of recalls made for each encoding serial position across lists
f_recall_by_serialpos = arrayfun(@(x) sum(encoding_events.recalled(encoding_events.serialpos == x))/sum(encoding_events.serialpos == x),serial_positions);
%%% -Total prior-list intrusions across lists of the condition
total_PLIs = sum(retrieval_events.intrusion > 0);
%%% -Total extra-list intrusions across lists of the condition
total_XLIs = sum(retrieval_events.intrusion == -1);
%%% -Total intrusions of any kind across lists of the condition
total_intrusions = total_PLIs + total_XLIs;

%%% If the condition is an experimental condition get t-test and
%%% Mann-Whitney U Test of the number of recalls made on each list.
if condition < 5
    %%% Get number of recalls from control encoding events
    control_condition = events.stimcode == 5;
    control_events = events(control_condition,:);
    encoding_type = strcmp(control_events.type,'WORD');
    control_events = control_events(encoding_type,:);
    control_lists = unique(control_events.list);
    control_n_recalls_by_list = arrayfun(@(x) sum(control_events.recalled(control_events.list==x)),control_lists);
    
    %%% Two sample t-test with unknown unequal variance
    [~,p_ttest_recalls,~,statistics] = ttest2(n_recalls_by_list,control_n_recalls_by_list,'Vartype','unequal');
    tstat_recalls = statistics.tstat;
    
    %%% Mann-Whitney U test
    [p_ranksum_recalls,~,statistics] = ranksum(n_recalls_by_list,control_n_recalls_by_list);
    try
        zval_recalls = statistics.zval;
    catch
        zval_recalls = NaN;
    end
else
    %%% Otherwise will fill row values with NaN
    p_ttest_recalls = NaN;
    tstat_recalls = NaN;
    p_ranksum_recalls = NaN;
    zval_recalls = NaN;    
end

%%% Preallocate matrices for the performance of temporal/semantic
%%% clustering and lag-conditional response probabilities measurements:

%%% -Matrix indicating which encoded items (column) for every list (row)
%%% were recalled or not (ones/zeros);
recalls_matrix = zeros(n_lists,n_list_items);
%%% -Matrix listing serial positions of recalled encoded items in the order
%%% in which they were recalled, the rest of positions filled with zeros.
recalls_serial_positions = zeros(n_lists,n_list_items);
%%% -Matrix listing item numbers of recalled encoded items in the order in
%%% which they were recalled, the rest of the positions filled with zeros.
recalls_item_numbers = zeros(n_lists,n_list_items);
%%% -Matrix listing item numbers of all encoded items in the order in which
%%% they were presented
encoding_item_numbers = zeros(n_lists,n_list_items);

%%% Loop through all lists to fill in the matrices
for idx = 1:n_lists
    list_number = condition_lists(idx);
    
    is_this_list = encoding_events.list == list_number;
    this_list_encoding = encoding_events(is_this_list,:);
    is_this_list = retrieval_events.list == list_number;
    not_intrusion = retrieval_events.intrusion == 0;
    good_number = retrieval_events.itemno > 0;
    this_list_retrieval = retrieval_events(is_this_list&not_intrusion&good_number,:);
    
    recalled = this_list_encoding.recalled;
    
    recalls_matrix(idx,:) = recalled;
    encoding_items = this_list_encoding.itemno;
    change = ~ismember(encoding_items,1:308);
    changed_items = encoding_items(change);
    changed_to = ceil(rand(sum(change),1)*307)+1;
    encoding_items(change) = changed_to;
    encoding_item_numbers(idx,:) = encoding_items;
    this_list_encoding.itemno = encoding_items;
    
    retrieval_items = this_list_retrieval.itemno;
    if any(ismember(retrieval_items,changed_items))
        for jdx = 1:length(retrieval_items)
            if ismember(retrieval_items(jdx),changed_items)
                retrieval_items(jdx) = changed_to(ismember(changed_items,retrieval_items(jdx)));
            end
        end
    end
    this_list_retrieval.itemno = retrieval_items;
    
    recalled_items = this_list_encoding(logical(recalled),:);
    n_recalled_items = height(this_list_retrieval);
    indices = arrayfun(@(x) find(recalled_items.itemno==x),this_list_retrieval.itemno);
    
    recalls_serial_positions(idx,1:n_recalled_items) = recalled_items.serialpos(indices);
    recalls_items = recalled_items.itemno(indices);
    recalls_item_numbers(idx,1:n_recalled_items) = recalls_items;
    
end


%%% Use above matrices to get temporal/semantic clustering factors and
%%% lag-CRP using Penn behavioral toolbox functions. Functions
%%% automatically mask any cells remaining in zeros.
subject_indices = ones(n_lists,1); %%% Input meant for processing multiple subject simultaneously. Ones since only doing one at a time.
temporal_clustering_factor = temp_fact(recalls_serial_positions,subject_indices,n_list_items);
semantic_clustering_factor = dist_fact(recalls_item_numbers,encoding_item_numbers,subject_indices,LSAmat);
lag_CRP_vector = lag_crp(recalls_serial_positions,subject_indices,n_list_items);

%%% Make behavioral table row
behavioral_table = session_info;
behavioral_table.condition = condition;
behavioral_table.n_lists = n_lists;
behavioral_table.overall_percent_recall = overall_percent_recall;
behavioral_table.tstat_recalls = tstat_recalls;
behavioral_table.p_ttest_recalls = p_ttest_recalls;
behavioral_table.zvalue_recalls = zval_recalls;
behavioral_table.p_ranksum_recalls = p_ranksum_recalls;
behavioral_table.f_recall_by_serialpos = {f_recall_by_serialpos};
behavioral_table.total_PLIs = total_PLIs;
behavioral_table.total_XLIs = total_XLIs;
behavioral_table.total_intrusions = total_intrusions;
behavioral_table.TCF = temporal_clustering_factor;
behavioral_table.SCF = semantic_clustering_factor;
behavioral_table.lag_CRP = {lag_CRP_vector(6:14)}; %%%Originally probabilities from -9 to 9 transitions, keeping only -4 to 4.

%%% Make n_recalls_by_list table
n_recalls_table = repmat(session_info,n_lists,1);
n_recalls_table.condition = repmat(condition,n_lists,1);
n_recalls_table.n_recalls = n_recalls_by_list;

%%% Make primacy table
primacy_table = repmat(session_info,n_list_items,1);
primacy_table.condition = repmat(condition,n_list_items,1);
primacy_table.primacy = [true;false(n_list_items-1,1)]; %%% Labelling only first serial position as primacy
primacy_table.f_recall = f_recall_by_serialpos';
end

function subject_behavioral_table = average_sessions(session_behavioral_table,subject_list,conditions)
%%% Variables to be averaged for subject behavioral table
variable_names = {'n_lists','overall_percent_recall','tstat_recalls','zvalue_recalls',...
    'f_recall_by_serialpos','total_PLIs','total_XLIs','total_intrusions',...
    'TCF','SCF','lag_CRP'};

subject_IDs = subject_list.subject_ID;
[Ax,Bx] = ndgrid(1:numel(subject_IDs),1:numel(conditions));
combo = table;
combo.subject_ID = subject_IDs(Ax(:));
combo.condition = conditions(Bx(:));
n_combos = height(combo);

subject_array = cell(n_combos,1);

for idx = 1:n_combos
    subject_ID = combo.subject_ID(idx);
    condition = combo.condition(idx);
    subject_info = subject_list(subject_list.subject_ID == subject_ID,:);
    same_subject = session_behavioral_table.subject_ID == subject_ID;
    same_condition = session_behavioral_table.condition == condition;
    if sum(same_subject & same_condition) > 0
        subject_data = session_behavioral_table(same_subject & same_condition,:);
        n_sessions = height(subject_data);
        
        temp_table = subject_info;
        temp_table.condition = condition;
        
        if n_sessions == 1
            subject_table = [temp_table,subject_data(:,variable_names)];
        else
            subject_table = temp_table;
            n_lists = subject_data.n_lists;
            total_lists = sum(n_lists);
            for jdx = 1:length(variable_names)
                variable_name = variable_names{jdx};
                switch variable_name
                    case 'n_lists'
                        subject_table.n_lists = sum(subject_data.n_lists);
                    case {'f_recall_by_serialpos','lag_CRP'}
                        data = vertcat(subject_data.(variable_name){:});
                        subject_table.(variable_name) = {sum(data.*n_lists,1)/total_lists};
                    otherwise
                        data = subject_data.(variable_name);
                        subject_table.(variable_name) = sum(data.*n_lists)/total_lists;
                end         
            end
        end
        subject_array{idx} = subject_table; clear subject_table;
    end
end

subject_behavioral_table = vertcat(subject_array{:});

end

function any_stimulation_table = average_stimulation_conditions(subject_behavioral_table,subject_list)
%%% Variables to be averaged for stimulation average behavioral table
variable_names = {'n_lists','overall_percent_recall','tstat_recalls','zvalue_recalls',...
    'f_recall_by_serialpos','total_PLIs','total_XLIs','total_intrusions',...
    'TCF','SCF','lag_CRP'};

conditions = [0;1]; %%%0 for no stimulation, 1 for any stimulation

subject_IDs = subject_list.subject_ID;
[Ax,Bx] = ndgrid(1:numel(subject_IDs),1:numel(conditions));
combo = table;
combo.subject_ID = subject_IDs(Ax(:));
combo.condition = conditions(Bx(:));
n_combos = height(combo);

subject_array = cell(n_combos,1);

for idx = 1:n_combos
    subject_ID = combo.subject_ID(idx);
    condition = combo.condition(idx);
    if condition == 1
        corresponding_conditions = [1,3];
    else
        corresponding_conditions = 5;
    end
    subject_info = subject_list(subject_list.subject_ID == subject_ID,:);
    same_subject = subject_behavioral_table.subject_ID == subject_ID;
    same_condition = ismember(subject_behavioral_table.condition,corresponding_conditions);
    
    if sum(same_subject & same_condition) > 0
        subject_data = subject_behavioral_table(same_subject & same_condition,:);
        n_rows = height(subject_data);
        
        temp_table = subject_info;
        temp_table.condition = condition;
        
        if n_rows == 1
            subject_table = [temp_table,subject_data(:,variable_names)];
        else
            subject_table = temp_table;
            n_lists = subject_data.n_lists;
            total_lists = sum(n_lists);
            for jdx = 1:length(variable_names)
                variable_name = variable_names{jdx};
                switch variable_name
                    case 'n_lists'
                        subject_table.n_lists = sum(subject_data.n_lists);
                    case {'f_recall_by_serialpos','lag_CRP'}
                        data = vertcat(subject_data.(variable_name){:});
                        subject_table.(variable_name) = {sum(data.*n_lists,1)/total_lists};
                    otherwise
                        data = subject_data.(variable_name);
                        subject_table.(variable_name) = sum(data.*n_lists)/total_lists;
                end         
            end
        end
        subject_array{idx} = subject_table; clear subject_table;
    end
end

any_stimulation_table = vertcat(subject_array{:});

end

function behavioral_boxplot(plots_directory,behavioral_metric,plot_type,behavioral_table,effects,interactions)
%%% Plot paramateres
figure_width = 170;
figure_height = 180;
plot_file = fullfile(plots_directory,[behavioral_metric '_' plot_type]);
behavioral_observations = behavioral_table.(behavioral_metric);
n_rows = length(behavioral_observations);
conditions = behavioral_table.condition;
subject_IDs = behavioral_table.subject_ID;
stimulation_groups = behavioral_table.stimulation_group;
is_IPL = strcmp(stimulation_groups,'IPL');
colors = [[255 7 202];[252 127 30];[59 124 255];[167 255 5]]/255;

if ismember(behavioral_metric,{'tstat_recalls','zvalue_recalls'})
    relevant_metric = 'n_recalls';
else
    relevant_metric = behavioral_metric;
end

relevant_effects = strcmp(effects.behavioral_metric,relevant_metric) & strcmp(effects.stimulation_hemisphere,'either');
relevant_interactions = strcmp(interactions.behavioral_metric,relevant_metric);
effects(~relevant_effects,:) = [];
interactions(~relevant_interactions,:) = [];

switch plot_type
    case 'separate'
        encoding_stimulation = conditions == 1;
        retrieval_stimulation = conditions == 3;
        no_stimulation = conditions == 5;
        
        encoding_match = arrayfun(@(x) find(subject_IDs == x & no_stimulation),subject_IDs(encoding_stimulation));
        retrieval_match = arrayfun(@(x) find(subject_IDs == x & no_stimulation),subject_IDs(retrieval_stimulation));
        
        delta_behavior = NaN(n_rows,1);
        if ~ismember(behavioral_metric,{'tstat_recalls','zvalue_recalls'})
            delta_behavior(encoding_stimulation) = behavioral_observations(encoding_stimulation) - behavioral_observations(encoding_match);
            delta_behavior(retrieval_stimulation) = behavioral_observations(retrieval_stimulation) - behavioral_observations(retrieval_match);
        else
            delta_behavior(encoding_stimulation) = behavioral_observations(encoding_stimulation);
            delta_behavior(retrieval_stimulation) = behavioral_observations(retrieval_stimulation);
        end
        
        encoding_IPL = encoding_stimulation & is_IPL;
        retrieval_IPL = retrieval_stimulation & is_IPL;
        encoding_RS = encoding_stimulation & ~is_IPL;
        retrieval_RS = retrieval_stimulation &  ~is_IPL;
        n_largest = max([sum(encoding_IPL),sum(retrieval_IPL),sum(encoding_RS),sum(retrieval_RS)]);
        
        x_limits = [0.4,5.6];
        x_ticks = [1,2,4,5];
        
        boxplot_box = NaN(n_largest,5);
        boxplot_box(1:sum(encoding_IPL),1) = delta_behavior(encoding_IPL);
        boxplot_box(1:sum(retrieval_IPL),2) = delta_behavior(retrieval_IPL);
        boxplot_box(1:sum(encoding_RS),4) = delta_behavior(encoding_RS);
        boxplot_box(1:sum(retrieval_RS),5) = delta_behavior(retrieval_RS);
        boxplot_colors = zeros(5,3);
        boxplot_positions = 1:5;
        
        scatter_x_positions = NaN(n_rows,1);
        scatter_colors = NaN(n_rows,3);
        scatter_x_positions(encoding_IPL) = normrnd(1,0.03,sum(encoding_IPL),1);
        scatter_colors(encoding_IPL,:) = repmat(colors(1,:),sum(encoding_IPL),1);
        scatter_x_positions(retrieval_IPL) = normrnd(2,0.03,sum(retrieval_IPL),1);
        scatter_colors(retrieval_IPL,:) = repmat(colors(2,:),sum(retrieval_IPL),1);
        scatter_x_positions(encoding_RS) = normrnd(4,0.03,sum(encoding_RS),1);
        scatter_colors(encoding_RS,:) = repmat(colors(3,:),sum(encoding_RS),1);
        scatter_x_positions(retrieval_RS) = normrnd(5,0.03,sum(retrieval_RS),1);
        scatter_colors(retrieval_RS,:) = repmat(colors(4,:),sum(retrieval_RS),1);
        
        p_values_effect = NaN(4,1);
        p_values_effect(1) = effects.pvalue(strcmp(effects.stimulation_group,'IPL')&strcmp(effects.stimulation_period,'encoding'));
        p_values_effect(2) = effects.pvalue(strcmp(effects.stimulation_group,'IPL')&strcmp(effects.stimulation_period,'retrieval'));
        p_values_effect(3) = effects.pvalue(strcmp(effects.stimulation_group,'retrosplenial')&strcmp(effects.stimulation_period,'encoding'));
        p_values_effect(4) = effects.pvalue(strcmp(effects.stimulation_group,'retrosplenial')&strcmp(effects.stimulation_period,'retrieval'));
        p_values_interaction = NaN(2,1);
        p_values_interaction(1) = interactions.pvalue(strcmp(interactions.stimulation_period,'encoding'));
        p_values_interaction(2) = interactions.pvalue(strcmp(interactions.stimulation_period,'retrieval'));
        effect_x_positions = boxplot_positions;
        interaction_x_positions = [2.5 3.5];
        
    case 'merged'
        stimulation = logical(conditions);
        stimulation_match = arrayfun(@(x) find(subject_IDs == x & ~stimulation),subject_IDs(stimulation));
        
        if ~ismember(behavioral_metric,{'tstat_recalls','zvalue_recalls'})
            delta_behavior = NaN(n_rows,1);
            delta_behavior(stimulation) = behavioral_observations(stimulation) - behavioral_observations(stimulation_match);
        else
            delta_behavior = NaN(n_rows,1);
            delta_behavior(stimulation) = behavioral_observations(stimulation);
        end
        
        IPL = stimulation & is_IPL;
        retrosplenial = stimulation &  ~is_IPL;
        
        x_limits = [0.4,2.6];
        x_ticks = [1,2];        
        
        n_largest = max([sum(IPL),sum(retrosplenial)]);
        boxplot_box = NaN(n_largest,2);
        boxplot_box(1:sum(IPL),1) = delta_behavior(IPL);
        boxplot_box(1:sum(retrosplenial),2) = delta_behavior(retrosplenial);
        boxplot_colors = zeros(2,3);
        boxplot_positions = 1:2;

        scatter_x_positions = NaN(n_rows,1);
        scatter_colors = NaN(n_rows,3);
        scatter_x_positions(IPL) = normrnd(1,0.03,sum(IPL),1);
        scatter_colors(IPL,:) = repmat(colors(1,:),sum(IPL),1);
        scatter_x_positions(retrosplenial) = normrnd(2,0.03,sum(retrosplenial),1);
        scatter_colors(retrosplenial,:) = repmat(colors(3,:),sum(retrosplenial),1);
        
        p_values_effect = NaN(2,1);
        p_values_effect(1) = effects.pvalue(strcmp(effects.stimulation_group,'IPL')&strcmp(effects.stimulation_period,'either'));
        p_values_effect(2) = effects.pvalue(strcmp(effects.stimulation_group,'retrosplenial')&strcmp(effects.stimulation_period,'either'));
        p_values_interaction = interactions.pvalue(strcmp(interactions.stimulation_period,'either'));
        effect_x_positions = boxplot_positions;
        interaction_x_positions = 1.5;
end
bad_rows = isnan(delta_behavior);
delta_behavior(bad_rows) = [];
scatter_colors(bad_rows,:) = [];
scatter_x_positions(bad_rows,:) = [];
max_absolute = max(abs(delta_behavior));
y_limits = [-1.05,1.05]*max_absolute;
y_ticks = (-1:0.25:1)*max_absolute;
y_tick_position = x_limits(1)+(diff(x_limits)*0.05);
effect_y_position = y_limits(2)-(diff(y_limits)*0.05);
interaction_y_position = y_limits(1)+(diff(y_limits)*0.05);

fig = figure('Units','pixels','Visible','off','Position',[0 0 figure_width figure_height]);
axes('Parent',fig,'Units','pixels','Position',[0 0 figure_width figure_height]);
hold on
scatter(scatter_x_positions,delta_behavior,[],scatter_colors,'filled');
boxplot_handle = boxplot(boxplot_box,boxplot_positions,'Colors',boxplot_colors,'Symbol','','DataLim',[-100 100]);
plot(x_limits,[0 0],'-k')
set(boxplot_handle,{'linew'},{1.5});
xlim(x_limits);xticks(x_ticks);xticklabels([]);
ylim(y_limits);yticks(y_ticks);yticklabels([]);
text(y_tick_position,max_absolute,sprintf('%.2f',max_absolute));
text(y_tick_position,0,'0');
text(y_tick_position,-1*max_absolute,sprintf('%.2f',-1*max_absolute));
for idx = 1:length(p_values_effect)
    if p_values_effect(idx)<0.01
        plot([-0.1,0.1]+effect_x_positions(idx),repmat(effect_y_position,1,2),'*k');
    elseif p_values_effect(idx)<0.05
        plot(effect_x_positions(idx),effect_y_position,'*k');
    elseif p_values_effect(idx)<0.1
        plot(effect_x_positions(idx),effect_y_position,'+k');
    end
end
for idx = 1:length(p_values_interaction)
    if p_values_interaction(idx)<0.01
        plot([-0.1,0.1]+interaction_x_positions(idx),repmat(interaction_y_position,1,2),'*k');
    elseif p_values_interaction(idx)<0.05
        plot(interaction_x_positions(idx),interaction_y_position,'*k');
    elseif p_values_interaction(idx)<0.1
        plot(interaction_x_positions(idx),interaction_y_position,'+k');
    end
end
hold off

print(plot_file,'-dsvg')
        
end

function serial_position_curves(plots_directory,stimulation_group,plot_type,behavioral_table,effects)
%%% Plot paramateres
figure_width = 176;
figure_height = 161;
plot_file = fullfile(plots_directory,sprintf('serialpos_%s_%s',stimulation_group,plot_type));
stimulation_groups = behavioral_table.stimulation_group;
same_group = strcmp(stimulation_groups,stimulation_group);
behavioral_table(~same_group,:) = [];
behavioral_observations = vertcat(behavioral_table.f_recall_by_serialpos{:});
conditions = behavioral_table.condition;
n_words = size(behavioral_observations,2);
colors = [[255 7 202];[252 127 30];[59 124 255];[167 255 5]]/255;
x_positions = 1:n_words;
x_limits = [0 n_words] + 0.5;
x_ticks = x_positions;
y_limits = [0 0.7];
y_ticks = y_limits(1):0.1:y_limits(2);

relevant_effects = strcmp(effects.behavioral_metric,'primacy') & strcmp(effects.stimulation_group,stimulation_group) & strcmp(effects.stimulation_hemisphere,'either');
effects(~relevant_effects,:) = [];

switch plot_type
    case 'separate'
        encoding_stimulation = conditions == 1;
        retrieval_stimulation = conditions == 3;
        no_stimulation = conditions == 5;
        n_curves = 3;
        curve_values = cell(n_curves,1);
        curve_values{1} = behavioral_observations(no_stimulation,:);
        curve_values{2} = behavioral_observations(encoding_stimulation,:);
        curve_values{3} = behavioral_observations(retrieval_stimulation,:);
        if strcmp(stimulation_group,'IPL')
            curve_colors = [[0.35 0.35 0.35];colors(1,:);colors(2,:)];
        else
            curve_colors = [[0.35 0.35 0.35];colors(3,:);colors(4,:)];
        end
        p_values = NaN(2,1);
        p_values(1) = effects.pvalue(strcmp(effects.stimulation_period,'encoding'));
        p_values(2) = effects.pvalue(strcmp(effects.stimulation_period,'retrieval'));
        effect_x_positions = [2,3];
    case 'merged'
        stimulation = logical(behavioral_table.condition);
        n_curves = 2;
        curve_values = cell(n_curves,1);
        curve_values{1} = behavioral_observations(~stimulation,:);
        curve_values{2} = behavioral_observations(stimulation,:);
        if strcmp(stimulation_group,'IPL')
            curve_colors = [colors(2,:);colors(1,:)];
        else
            curve_colors = [colors(4,:);colors(3,:)];
        end
        p_values = effects.pvalue(strcmp(effects.stimulation_period,'retrieval'));
        effect_x_positions = 2;
end

effect_y_position = y_limits(2)-(diff(y_limits)*0.05);

fig = figure('Units','pixels','Visible','off','Position',[0 0 figure_width figure_height]);
axes('Parent',fig,'Units','pixels','Position',[0 0 figure_width figure_height]);
hold on
for idx = 1:n_curves
    stdshade(curve_values{idx},0.4,curve_colors(idx,:));
end
xlim(x_limits);xticks(x_ticks);xticklabels([]);
ylim(y_limits);yticks(y_ticks);yticklabels([]);
for idx = 1:length(p_values)
    if p_values(idx)<0.01
        plot([-0.1,0.1]+effect_x_positions(idx),repmat(effect_y_position,1,2),'*k');
    elseif p_values(idx)<0.05
        plot(effect_x_positions(idx),effect_y_position,'*k');
    elseif p_values(idx)<0.1
        plot(effect_x_positions(idx),effect_y_position,'+k');
    end
end
hold off

print(plot_file,'-dsvg')
end

function lag_crp_curves(plots_directory,stimulation_group,plot_type,behavioral_table)
%%% Plot paramateres
figure_width = 176;
figure_height = 161;
plot_file = fullfile(plots_directory,sprintf('lag-crp_%s_%s',stimulation_group,plot_type));
stimulation_groups = behavioral_table.stimulation_group;
same_group = strcmp(stimulation_groups,stimulation_group);
behavioral_table(~same_group,:) = [];
behavioral_observations = vertcat(behavioral_table.lag_CRP{:});
conditions = behavioral_table.condition;
n_transitions = (size(behavioral_observations,2)-1)/2;
portion1 = 1:n_transitions;
portion2 = n_transitions+2:(n_transitions*2)+1;
colors = [[255 7 202];[252 127 30];[59 124 255];[167 255 5]]/255;
x_limits = [0 portion2(end)] + 0.5;
x_ticks = [portion1,portion2];
y_limits = [0 0.7];
y_ticks = y_limits(1):0.1:y_limits(2);

switch plot_type
    case 'separate'
        encoding_stimulation = conditions == 1;
        retrieval_stimulation = conditions == 3;
        no_stimulation = conditions == 5;
        n_curves = 3;
        curve_values = cell(n_curves,1);
        curve_values{1} = behavioral_observations(no_stimulation,:);
        curve_values{2} = behavioral_observations(encoding_stimulation,:);
        curve_values{3} = behavioral_observations(retrieval_stimulation,:);
        if strcmp(stimulation_group,'IPL')
            curve_colors = [[0.35 0.35 0.35];colors(1,:);colors(2,:)];
        else
            curve_colors = [[0.35 0.35 0.35];colors(3,:);colors(4,:)];
        end     
    case 'merged'
        stimulation = logical(behavioral_table.condition);
        n_curves = 2;
        curve_values = cell(n_curves,1);
        curve_values{1} = behavioral_observations(~stimulation,:);
        curve_values{2} = behavioral_observations(stimulation,:);
        if strcmp(stimulation_group,'IPL')
            curve_colors = [colors(2,:);colors(1,:)];
        else
            curve_colors = [colors(4,:);colors(3,:)];
        end     
end

fig = figure('Units','pixels','Visible','off','Position',[0 0 figure_width figure_height]);
axes('Parent',fig,'Units','pixels','Position',[0 0 figure_width figure_height]);
hold on
for idx = 1:n_curves
    stdshade(curve_values{idx}(:,portion1),0.4,curve_colors(idx,:),portion1);
    stdshade(curve_values{idx}(:,portion2),0.4,curve_colors(idx,:),portion2);
end
xlim(x_limits);xticks(x_ticks);xticklabels([]);
ylim(y_limits);yticks(y_ticks);yticklabels([]);
hold off

print(plot_file,'-dsvg')
end
