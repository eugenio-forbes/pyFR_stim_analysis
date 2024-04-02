function n08_group_matching_tests(varargin)
if isempty(varargin)
    %%% Directory information
    root_directory = '/directory/to/pyFR_stim_analysis';
else
    root_directory = varargin{1};
end

%%% List directories
list_directory = fullfile(root_directory,'lists');
statistics_directory = fullfile(root_directory,'statistics/group_matching');
exclusions_directory = fullfile(root_directory,'exclusion_lists');
plots_directory = fullfile(root_directory,'plots');

%%% Create directories if needed
if ~isfolder(statistics_directory)
    mkdir(statistics_directory);
end
if ~isfolder(plots_directory)
    mkdir(plots_directory);
end

%%% Load subject list
load(fullfile(list_directory,'subject_list.mat'),'subject_list');

%%% Calculate missing columns from available columns
subject_list.Epilepsy_Years = subject_list.Age - subject_list.Onset;
subject_list.p_ictal = subject_list.n_ictal./subject_list.n_electrodes;
subject_list.p_interictal = subject_list.n_interictal./subject_list.n_electrodes;

%%% Define variable to be tested for group matching
binary_variables = {'Female','White','Black','Hispanic','Right_Handed',...
    'Left_Dominant','Heterotopia','Infarcts','Infectious','Febrile','Dysplasia',...
    'Encephalomalacia','Cavernoma','Cystic','AVM','Atrophy','Edema','Gliosis',...               
    'Autoimmune','MTS','Developmental','TBI','PSHx','Reimplant','Resection','VNS'};
numerical_variables = {'Age','Onset','Epilepsy_Years','Education','n_ictal','n_interictal','n_electrodes','p_ictal','p_interictal','stimulation_site_ictal','stimulation_site_interictal'};   
neuropsychometrics = {'FSIQ','VCI','PRI','VLT_Short','VLT_Long','VLT_Total','Depression'};

%%% To convert neuropsychometrics to numeric rank for Mann-Whitney U Test
neuropsychometrics_score_ranges = {'Impaired','Borderline','Low_Avg','Average','High_Avg','Superior','Very_Sup'};
depression_ranges = {'Minimal','Mild','Moderate','Severe'};

%%% Get Fisher's Exact Test for binarized categorical variables to test for
%%% matching of subject groups (retrosplenial and IPL stimulation groups)                                     
%%% Initialize table columns to hold statistics

n_variables = length(binary_variables);
test_variable = cell(n_variables,1);
group1 = repelem({'retrosplenial'},n_variables,1);
group2 = repelem({'IPL'},n_variables,1);
n_group1 = NaN(n_variables,1);
n_group2 = NaN(n_variables,1);
f_group1 = NaN(n_variables,1);
f_group2 = NaN(n_variables,1);
odds_ratio = NaN(n_variables,1);
confidence_interval = cell(n_variables,1);
p_value = NaN(n_variables,1);

for idx = 1:n_variables
    
    this_variable = binary_variables{idx};
    test_variable{idx} = this_variable;
    test_exclusions = subject_list.(this_variable) < 0;
    temp_subject_list = subject_list(~test_exclusions,:);
    retrosplenial_subjects = strcmp(temp_subject_list.stimulation_group,'retrosplenial');
    IPL_subjects = strcmp(temp_subject_list.stimulation_group,'IPL');
    retrosplenial_values = temp_subject_list.(this_variable)(retrosplenial_subjects);
    IPL_values = temp_subject_list.(this_variable)(IPL_subjects);
    retrosplenial_row = [sum(retrosplenial_values),sum(~retrosplenial_values)];
    IPL_row = [sum(IPL_values),sum(~IPL_values)];
    contingency_table = [retrosplenial_row;IPL_row];
    
    [~,p_val,stats] = fishertest(contingency_table);
    
    p_value(idx) = round(p_val*10000)/10000;
    odds_ratio(idx) = round(stats.OddsRatio*100)/100;
    confidence_interval{idx} = round(stats.ConfidenceInterval*100)/100;
    n_group1(idx) = sum(retrosplenial_row);
    n_group2(idx) = sum(IPL_row);
    f_group1(idx) = round((retrosplenial_row(1)/sum(retrosplenial_row))*1000)/10;
    f_group2(idx) = round((IPL_row(1)/sum(IPL_row))*1000)/10;
        
    clear temp_subject_list retrosplenial_values IPL_values retrosplenial_row IPL_row contingency_table p_val stats this_variable
end
binary_statistics = table(test_variable,group1,group2,n_group1,n_group2,f_group1,f_group2,odds_ratio,confidence_interval,p_value);
clear test_variable group1 group2 n_group1 n_group2 f_group1 f_group2 odds_ratio confidence_interval p_value

%%% Get Mann-Whitney U-test for numberical variables given that the
%%% distributions are not normal or close to normal. Initialize table
%%% columns to hold statistics.

n_variables = length(numerical_variables);
test_variable = cell(n_variables,1);
group1 = repelem({'retrosplenial'},n_variables,1);
group2 = repelem({'IPL'},n_variables,1);
n_group1 = NaN(n_variables,1);
n_group2 = NaN(n_variables,1);
m_group1 = NaN(n_variables,1);
m_group2 = NaN(n_variables,1);
sd_group1 = NaN(n_variables,1);
sd_group2 = NaN(n_variables,1);
z_value= NaN(n_variables,1);
p_value = NaN(n_variables,1);

for idx = 1:n_variables
    this_variable = numerical_variables{idx};
    test_variable{idx} = this_variable;
    
    retrosplenial_subjects = strcmp(subject_list.stimulation_group,'retrosplenial');
    IPL_subjects = strcmp(subject_list.stimulation_group,'IPL');
    
    retrosplenial_values = subject_list.(this_variable)(retrosplenial_subjects);
    IPL_values = subject_list.(this_variable)(IPL_subjects);
    
    [p_val,~,stats] = ranksum(retrosplenial_values,IPL_values);
    
    p_value(idx) = round(p_val*10000)/10000;
    z_value(idx) = round(stats.zval*100)/100;  
    n_group1(idx) = length(retrosplenial_values);
    n_group2(idx) = length(IPL_values);
    m_group1(idx) = round(mean(retrosplenial_values)*100)/100;
    m_group2(idx) = round(mean(IPL_values)*100)/100;
    sd_group1(idx) = round(std(retrosplenial_values)*100)/100;
    sd_group2(idx) = round(std(IPL_values)*100)/100;
    
    clear retrosplenial_values IPL_values p_val stats this_variable
end

numerical_statistics = table(test_variable,group1,group2,n_group1,n_group2,m_group1,m_group2,sd_group1,sd_group2,z_value,p_value);
clear test_variable group1 group2 n_group1 n_group2 m_group1 m_group2 sd_group1 sd_group2 z_value p_value

%%% Get Mann-Whitney U-test for neuropsychometrics given that the
%%% variables are ordinal categorical. Initialize table columns to hold 
%%% statistics.
n_variables = length(neuropsychometrics);
test_variable = cell(n_variables,1);
group1 = repelem({'retrosplenial'},n_variables,1);
group2 = repelem({'IPL'},n_variables,1);
n_group1 = NaN(n_variables,1);
n_group2 = NaN(n_variables,1);
m_group1 = NaN(n_variables,1);
m_group2 = NaN(n_variables,1);
z_value= NaN(n_variables,1);
p_value = NaN(n_variables,1);

for idx = 1:n_variables
    this_variable = neuropsychometrics{idx};
    test_variable{idx} = this_variable;
    
    temp_subject_list = subject_list(:,{'stimulation_group',this_variable});
    temp_subject_list = temp_subject_list(~strcmp(temp_subject_list.(this_variable),'N'),:);
    
    if strcmp(this_variable,'Depression')
        temp_subject_list.values = cellfun(@(x) find(strcmp(depression_ranges,x)),temp_subject_list.(this_variable));
    else
        temp_subject_list.values = cellfun(@(x) find(strcmp(neuropsychometrics_score_ranges,x)),temp_subject_list.(this_variable));
    end
    
    retrosplenial_subjects = strcmp(temp_subject_list.stimulation_group,'retrosplenial');
    IPL_subjects = strcmp(temp_subject_list.stimulation_group,'IPL');
    
    retrosplenial_values = temp_subject_list.values(retrosplenial_subjects);
    IPL_values = temp_subject_list.values(IPL_subjects);
    
    [p_val,~,stats] = ranksum(retrosplenial_values,IPL_values);
    
    p_value(idx) = round(p_val*10000)/10000;
    z_value(idx) = round(stats.zval*100)/100;  
    n_group1(idx) = length(retrosplenial_values);
    n_group2(idx) = length(IPL_values);
    m_group1(idx) = round(mean(retrosplenial_values)*100)/100;
    m_group2(idx) = round(mean(IPL_values)*100)/100;
    
    clear temp_subject_list retrosplenial_values IPL_values p_val stats this_variable
end

neuropsychometrics_statistics = table(test_variable,group1,group2,n_group1,n_group2,m_group1,m_group2,z_value,p_value);
clear test_variable group1 group2 n_group1 n_group2 m_group1 m_group2 z_value p_value

%%% Save statistics
save(fullfile(statistics_directory,'binary_statistics.mat'),'binary_statistics');
save(fullfile(statistics_directory,'numerical_statistics.mat'),'numerical_statistics');
save(fullfile(statistics_directory,'neuropsychometrics_statistics.mat'),'neuropsychometrics_statistics');
end    
