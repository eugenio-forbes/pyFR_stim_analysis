function n04_group_subjects(varargin)
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
list_directory = fullfile(analysis_directory,'lists');

%%% Load subject list file
load(fullfile(list_directory,'subject_list.mat'),'subject_list');

%%% Load session list
load(fullfile(list_directory,'session_list.mat'),'session_list');

%%% Grouping subjects based on whether the electrode is located in the
%%% retrosplenial region or whether it is in the inferior parietal lobule
n_sessions = height(session_list);
stimulation_hemisphere = cell(n_sessions,1);
stimulation_group = cell(n_sessions,1);

%%% Loop through sessions
for idx = 1:n_sessions
    stimulation_electrode_label = session_list.stimulation_electrode_label{idx};
    if strcmp(stimulation_electrode_label(1),'R')
        stimulation_hemisphere{idx} = 'right';
    else
        stimulation_hemisphere{idx} = 'left';
    end
    
    stimulation_neurologist_label = session_list(idx,:).stimulation_neurologist_label{:};
    if contains(stimulation_neurologist_label,{'post','cing'},'IgnoreCase',true)
        stimulation_group{idx} = 'retrosplenial';
    else
        stimulation_group{idx} = 'IPL';
    end
end
session_list.stimulation_group = stimulation_group;
session_list.stimulation_hemisphere = stimulation_hemisphere;

%%% Transfer information to subject list
[~,subject_indices,~] = unique(session_list.subject);
subject_list.stimulation_group = stimulation_group(subject_indices);
subject_list.stimulation_hemisphere = stimulation_hemisphere(subject_indices);
subject_list.stimulation_electrode_label = session_list.stimulation_electrode_label(subject_indices);
subject_list.stimulation_electrode_channel = session_list.stimulation_electrode_channel(subject_indices);
subject_list.stimulation_neurologist_label = session_list.stimulation_neurologist_label(subject_indices);
subject_list.stimulation_automatic_label = session_list.stimulation_automatic_label(subject_indices);
subject_list.stimulation_references = session_list.stimulation_references(subject_indices);
subject_list.stimulation_reference_type = session_list.stimulation_reference_type(subject_indices);
subject_list.stimulation_coordinates = session_list.stimulation_coordinates(subject_indices);

%%% Save lists
save(fullfile(list_directory,'subject_list.mat'),'subject_list');
save(fullfile(list_directory,'session_list.mat'),'session_list');
end