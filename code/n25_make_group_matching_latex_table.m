function n25_make_group_matching_latex_table(varargin)
if isempty(varargin)
    %%% Directory information
    root_directory = '/directory/to/pyFR_stim_analysis';
else
    root_directory = varargin{1};
end

%%% List directories
statistics_directory = fullfile(root_directory,'statistics/group_matching');
latex_directory = fullfile(root_directory,'latex');
if ~isfolder(latex_directory)
    mkdir(latex_directory);
end

load(fullfile(statistics_directory,'binary_statistics.mat'),'binary_statistics');
load(fullfile(statistics_directory,'numerical_statistics.mat'),'numerical_statistics');
load(fullfile(statistics_directory,'neuropsychometrics_statistics.mat'),'neuropsychometrics_statistics');

file_name = fullfile(latex_directory,'binary_statistics_latex_table.txt');
file_id = fopen(file_name,'w');
fprintf(file_id,'\\begin{table}[h!]\n');
fprintf(file_id,'\\begin{tabular}{|c|c|c|c|c|}\n');
fprintf(file_id,'\\hline\n');
fprintf(file_id,'Variable & PCC & AG & Odds Ratio & p-value \\\\ \\hline\n'); % Table headers

bad_rows = cellfun(@(x) any(isinf(x)),binary_statistics.confidence_interval);
binary_statistics(bad_rows,:) = [];
n_rows = height(binary_statistics);

for idx = 1:n_rows
    variable = binary_statistics.test_variable{idx};
    retrosplenial = [sprintf('%.1f',binary_statistics.f_group1(idx)),'\%, N=',num2str(binary_statistics.n_group1(idx))];
    IPL = [sprintf('%.1f',binary_statistics.f_group2(idx)),'\%, N=',num2str(binary_statistics.n_group2(idx))];
    odds_ratio = binary_statistics.odds_ratio(idx);
    pvalue = binary_statistics.p_value(idx);
    fprintf(file_id,'%s & %s & %s & %.2f & %.2f \\\\ \\hline\n',...
        variable,retrosplenial,IPL,odds_ratio,pvalue);
end
fprintf(file_id, '\\end{tabular}\n');
fprintf(file_id, '\\end{table}\n');
fclose(file_id);

file_name = fullfile(latex_directory,'numerical_statistics_latex_table.txt');
file_id = fopen(file_name,'w');
fprintf(file_id,'\\begin{table}[h!]\n');
fprintf(file_id,'\\begin{tabular}{|c|c|c|c|c|}\n');
fprintf(file_id,'\\hline\n');
fprintf(file_id,'Variable & PCC & AG & z-value & p-value \\\\ \\hline\n'); % Table headers

bad_rows = isnan(numerical_statistics.p_value);
numerical_statistics(bad_rows,:) = [];
n_rows = height(numerical_statistics);

for idx = 1:n_rows
    variable = numerical_statistics.test_variable{idx};
    retrosplenial = [sprintf('%.1f(%.1f)',numerical_statistics.m_group1(idx),numerical_statistics.sd_group1(idx)),', N=',num2str(numerical_statistics.n_group1(idx))];
    IPL = [sprintf('%.1f(%.1f)',numerical_statistics.m_group2(idx),numerical_statistics.sd_group2(idx)),', N=',num2str(numerical_statistics.n_group2(idx))];
    zvalue = numerical_statistics.z_value(idx);
    pvalue = numerical_statistics.p_value(idx);
    fprintf(file_id,'%s & %s & %s & %.2f & %.2f \\\\ \\hline\n',...
        variable,retrosplenial,IPL,zvalue,pvalue);
end

fprintf(file_id, '\\end{tabular}\n');
fprintf(file_id, '\\end{table}\n');
fclose(file_id);

file_name = fullfile(latex_directory,'neuropsychometrics_latex_table.txt');
file_id = fopen(file_name,'w');

fprintf(file_id,'\\begin{table}[h!]\n');
fprintf(file_id,'\\begin{tabular}{|c|c|c|c|c|}\n');
fprintf(file_id,'\\hline\n');
fprintf(file_id,'Variable & PCC & AG & z-value & p-value \\\\ \\hline\n'); % Table headers

bad_rows = isnan(neuropsychometrics_statistics.p_value);
neuropsychometrics_statistics(bad_rows,:) = [];
n_rows = height(neuropsychometrics_statistics);

for idx = 1:n_rows
    variable = neuropsychometrics_statistics.test_variable{idx};
    retrosplenial = [sprintf('%.1f',neuropsychometrics_statistics.m_group1(idx)),', N=',num2str(neuropsychometrics_statistics.n_group1(idx))];
    IPL = [sprintf('%.1f',neuropsychometrics_statistics.m_group2(idx)),', N=',num2str(neuropsychometrics_statistics.n_group2(idx))];
    zvalue = neuropsychometrics_statistics.z_value(idx);
    pvalue = neuropsychometrics_statistics.p_value(idx);
    fprintf(file_id,'%s & %s & %s & %.2f & %.2f \\\\ \\hline\n',...
        variable,retrosplenial,IPL,zvalue,pvalue);
end

fprintf(file_id, '\\end{tabular}\n');
fprintf(file_id, '\\end{table}\n');
fclose(file_id);

end
