function n24_make_demographics_latex_table(varargin)
if isempty(varargin)
    %%% Directory information
    root_directory = '/directory/to/pyFR_stim_analysis';
else
    root_directory = varargin{1};
end

%%% List directories
list_directory = fullfile(root_directory,'lists');
latex_directory = fullfile(root_directory,'latex');
if ~isfolder(latex_directory)
    mkdir(latex_directory);
end

load(fullfile(list_directory,'subject_list.mat'),'subject_list');
subject_list = sortrows(subject_list,{'stimulation_group','stimulation_hemisphere','Age','Right_Handed'});
n_subjects = height(subject_list);

file_name = fullfile(latex_directory,'demographics_latex_table.txt');
file_id = fopen(file_name,'w');
fprintf(file_id,'\\begin{table}[h!]\n');
fprintf(file_id,'\\begin{tabular}{|c|c|c|c|c|c|c|c|c|}\n');
fprintf(file_id,'\\hline\n');
fprintf(file_id,'Group & Hemisphere & Subject & Age & Years of Epilepsy & Race & Gender & Handedness & N Electrodes \\\\ \\hline\n'); % Table headers

group_indices = cellfun(@(x) find(strcmp(subject_list.stimulation_group,x),1,'first'),unique(subject_list.stimulation_group));
hemisphere_change = find(diff(strcmp(subject_list.stimulation_hemisphere,'right'))~=0)+1;
hemisphere_indices = [1;hemisphere_change];

for idx = 1:n_subjects
    
    if ismember(idx,group_indices)
        group = subject_list.stimulation_group{idx};
    else
        group = '';
    end
    
    if ismember(idx,hemisphere_indices)
        hemisphere = subject_list.stimulation_hemisphere{idx};
    else
        hemisphere = '';
    end
    
    subject = sprintf('%02d',idx);
    
    age = subject_list.Age(idx);
    onset = subject_list.Onset(idx);
    years_of_epilepsy = age-onset;
    
    is_female = subject_list.Female(idx);
    if is_female
        gender = 'F';
    else
        gender = 'M';
    end
    
    is_white = subject_list.White(idx);
    is_black = subject_list.Black(idx);
    is_hispanic = subject_list.Hispanic(idx);
    if is_white
        race = 'W';
    elseif is_black
        race = 'B';
    elseif is_hispanic
        race = 'H';
    end
    
    is_right_handed = subject_list.Right_Handed(idx);
    if is_right_handed
        handedness = 'R';
    else
        handedness = 'L';
    end
    
    n_electrodes = subject_list.n_electrodes(idx);
    
    fprintf(file_id,'%s & %s & %s & %.2f & %.2f & %s & %s & %s & %d \\\\ \\hline\n',...
        group,hemisphere,subject,age,years_of_epilepsy,race,gender,handedness,n_electrodes);
end
fprintf(file_id, '\\end{tabular}\n');
fprintf(file_id, '\\end{table}\n');
fclose(file_id);
end
