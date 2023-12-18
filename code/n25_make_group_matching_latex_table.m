function n25_make_group_matching_latex_table(varargin)
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
statistics_directory = fullfile(analysis_directory,'statistics/group_matching');
latex_directory = fullfile(analysis_directory,'latex');
if ~isfolder(latex_directory)
    mkdir(latex_directory);
end

%%% Load files with statistical results of group matching tests
load(fullfile(statistics_directory,'binary_statistics.mat'),'binary_statistics');
load(fullfile(statistics_directory,'numerical_statistics.mat'),'numerical_statistics');
load(fullfile(statistics_directory,'neuropsychometrics_statistics.mat'),'neuropsychometrics_statistics');

%%% Open file to print table with results comparing demographic proportions
%%% between groups. First few lines will contain latex commands to start a
%%% table and table headers
file_name = fullfile(latex_directory,'binary_statistics_latex_table.txt');
file_id = fopen(file_name,'w');
fprintf(file_id,'\\begin{table}[h!]\n');
fprintf(file_id,'\\begin{tabular}{|c|c|c|c|c|}\n');
fprintf(file_id,'\\hline\n');
fprintf(file_id,'Variable & PCC & AG & Odds Ratio & p-value \\\\ \\hline\n'); % Table headers

%%% Remove rows with bad statistical tests (variables with insufficient data)
bad_rows = cellfun(@(x) any(isinf(x)),binary_statistics.confidence_interval);
binary_statistics(bad_rows,:) = [];
n_rows = height(binary_statistics);

%%% Loop through table to rows to print results in predetermined format
for idx = 1:n_rows
    %%% Name of tested variable
    variable = binary_statistics.test_variable{idx};
    %%% Number of subjects in group 1
    retrosplenial = [sprintf('%.1f',binary_statistics.f_group1(idx)),'\%, N=',num2str(binary_statistics.n_group1(idx))];
    %%% Number of subjecs in group 2
    IPL = [sprintf('%.1f',binary_statistics.f_group2(idx)),'\%, N=',num2str(binary_statistics.n_group2(idx))];
    %%% Odds ratio of Fisher's Exact test
    odds_ratio = binary_statistics.odds_ratio(idx);
    %%% P value of Fisher's Exact test
    pvalue = binary_statistics.p_value(idx);
    
    %%% Print line with these values divided by columns
    fprintf(file_id,'%s & %s & %s & %.2f & %.2f \\\\ \\hline\n',...
        variable,retrosplenial,IPL,odds_ratio,pvalue);
end
%%% Last few lines contain commands to close latex table
fprintf(file_id, '\\end{tabular}\n');
fprintf(file_id, '\\end{table}\n');
fclose(file_id);

%%% Follow these same steps for between group comparisons of numerical variables
%%% Opening lines
file_name = fullfile(latex_directory,'numerical_statistics_latex_table.txt');
file_id = fopen(file_name,'w');
fprintf(file_id,'\\begin{table}[h!]\n');
fprintf(file_id,'\\begin{tabular}{|c|c|c|c|c|}\n');
fprintf(file_id,'\\hline\n');
fprintf(file_id,'Variable & PCC & AG & z-value & p-value \\\\ \\hline\n'); % Table headers

%%% Row exclusions
bad_rows = isnan(numerical_statistics.p_value);
numerical_statistics(bad_rows,:) = [];
n_rows = height(numerical_statistics);

%%% Row printing
for idx = 1:n_rows
    %%% Name of variable tested
    variable = numerical_statistics.test_variable{idx};
    %%% Number of subject in group 1
    retrosplenial = [sprintf('%.1f(%.1f)',numerical_statistics.m_group1(idx),numerical_statistics.sd_group1(idx)),', N=',num2str(numerical_statistics.n_group1(idx))];
    %%% Number of subjects in group 2
    IPL = [sprintf('%.1f(%.1f)',numerical_statistics.m_group2(idx),numerical_statistics.sd_group2(idx)),', N=',num2str(numerical_statistics.n_group2(idx))];
    %%% Z-value of Mann Whitney U test
    zvalue = numerical_statistics.z_value(idx);
    %%% P-value of Mann Whitney U test
    pvalue = numerical_statistics.p_value(idx);
    
    %%% Print values divided by column
    fprintf(file_id,'%s & %s & %s & %.2f & %.2f \\\\ \\hline\n',...
        variable,retrosplenial,IPL,zvalue,pvalue);
end

%%% Closing lines
fprintf(file_id, '\\end{tabular}\n');
fprintf(file_id, '\\end{table}\n');
fclose(file_id);

%%% Follow these same steps for between group comparisons of neuropsychometrics
%%% Opening lines
file_name = fullfile(latex_directory,'neuropsychometrics_latex_table.txt');
file_id = fopen(file_name,'w');

fprintf(file_id,'\\begin{table}[h!]\n');
fprintf(file_id,'\\begin{tabular}{|c|c|c|c|c|}\n');
fprintf(file_id,'\\hline\n');
fprintf(file_id,'Variable & PCC & AG & z-value & p-value \\\\ \\hline\n'); % Table headers

%%% Exclusions of rows
bad_rows = isnan(neuropsychometrics_statistics.p_value);
neuropsychometrics_statistics(bad_rows,:) = [];
n_rows = height(neuropsychometrics_statistics);

%%% Row printing
for idx = 1:n_rows
    %%% Name of variable tested
    variable = neuropsychometrics_statistics.test_variable{idx};
    %%% Number of subjects in group 1
    retrosplenial = [sprintf('%.1f',neuropsychometrics_statistics.m_group1(idx)),', N=',num2str(neuropsychometrics_statistics.n_group1(idx))];
    %%% Number of subjects in group 2
    IPL = [sprintf('%.1f',neuropsychometrics_statistics.m_group2(idx)),', N=',num2str(neuropsychometrics_statistics.n_group2(idx))];
    %%% Z-value of Mann Whitney U test
    zvalue = neuropsychometrics_statistics.z_value(idx);
    %%% P-value of Mann Whintey U test
    pvalue = neuropsychometrics_statistics.p_value(idx);
    
    %%% Print values divided by columns
    fprintf(file_id,'%s & %s & %s & %.2f & %.2f \\\\ \\hline\n',...
        variable,retrosplenial,IPL,zvalue,pvalue);
end

%%% Closing lines
fprintf(file_id, '\\end{tabular}\n');
fprintf(file_id, '\\end{table}\n');
fclose(file_id);

end