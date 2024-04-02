function n06_add_subject_data(varargin)
if isempty(varargin)
    %%% Directory information
    root_directory = '/directory/to/pyFR_stim_analysis';
else
    root_directory = varargin{1};
end

%%% List directories
list_directory = fullfile(root_directory,'lists');
resource_directory = fullfile(root_directory,'resources');

%%% Specify import options for excel sheet with subject data
excel_file_options = spreadsheetImportOptions('NumVariables', 42);
excel_file_options.Sheet = 'pyFR_stim_subject_data_binary';
excel_file_options.DataRange = 'A2:AP81';
excel_file_options.VariableNames = ...
    {'Subject', 'Age', 'Onset', 'Female', ...                             %%% UT and Penn subject ID's, age at SEEG implant, age of seizure onset, whether female (logical)
    'White','Black', 'Hispanic', 'Right_Handed', 'Left_Dominant',...               %%% logicals for race groups, handedness, and hemisphere with language dominance
    'Heterotopia','Infarcts', 'Infectious', 'Febrile', 'Dysplasia',...             %%% logicals for etiologies
    'Encephalomalacia','Cavernoma', 'Cystic', 'AVM', 'Atrophy',...                 %%% logicals for etiologies
    'Edema', 'Gliosis', 'Autoimmune','MTS', 'Developmental',...                    %%% logicals for etiologies
    'TBI', 'PSHx', 'Reimplant', 'Resection', 'VNS',...                             %%% logicals for incidence of TBI and PSHx
    'MRI_abnormalities', 'MRI_lateralization', 'MRI_distribution', 'Education',... %%% specific MRI abnormalities, their lateralization, and specific lobe distribution, years of education
    'FSIQ', 'VCI', 'PRI', 'VLT_Short', 'VLT_Long', 'VLT_Total',...                 %%% range descriptors for full-scale IQ, verbal comprehension index, perceptual reasoning index, verbal learning test (short/long delay and total)
    'Depression','Ictal', 'Interictal'};                                           %%% range descriptor for symptoms of depression, channel labels of ictal onset, channel labels of interictal spikes (only selective parietal and mesial temporal included)
excel_file_options.VariableTypes = ...
    {'char', 'double', 'double', 'logical',...
    'logical', 'logical', 'logical', 'logical', 'double',...                       %%% Saving 'Left_Dominant' as double to exclude the few cases with no language dominance testing done (fMRI,WADA)
    'logical', 'logical', 'logical', 'logical', 'logical',...
    'logical', 'logical', 'logical', 'logical', 'logical',...
    'logical', 'logical', 'logical', 'logical', 'logical',...
    'logical', 'logical', 'logical', 'logical', 'logical',...
    'char', 'char', 'char', 'double',...
    'char', 'char', 'char', 'char', 'char', 'char',...
    'char', 'char', 'char'};
excel_file_options = setvaropts(excel_file_options, 'Subject', 'WhitespaceRule', 'preserve');
excel_file_options = setvaropts(excel_file_options, {'Subject', 'MRI_abnormalities', 'MRI_lateralization', 'MRI_distribution', 'FSIQ', 'VCI', 'PRI', 'VLT_Short', 'VLT_Long', 'VLT_Total', 'Depression', 'Ictal', 'Interictal'}, 'EmptyFieldRule', 'auto');

%%% Import subject data from excel file
subject_data_file = fullfile(resource_directory,'pyFR_stim_subject_data_binary.xlsx');
subject_data = readtable(subject_data_file, excel_file_options, 'UseExcel', false);

%%% Load subject list
load(fullfile(list_directory,'subject_list.mat'),'subject_list');

%%% Load session list
load(fullfile(list_directory,'session_list.mat'),'session_list');

%%% Load electrode list
load(fullfile(list_directory,'electrode_list.mat'),'electrode_list');

%%% Get matching indices to filter and resort table, then merge with
%%% subject list.
subject_indices = cellfun(@(x) find(strcmp(subject_data.Subject,x(1:5))),subject_list.subject);
subject_data = subject_data(subject_indices,:);
subject_data.Subject = [];
subject_list = [subject_list,subject_data];

%%% Classify electrodes as ictal onset zone or interictal spiking zone
%%% based on subject ictal and interictal labels
n_electrodes = height(electrode_list);
ictal = false(n_electrodes,1);
interictal = false(n_electrodes,1);
stimulation_site_ictal = false(n_electrodes,1);
stimulation_site_interictal = false(n_electrodes,1);
for idx = 1:n_electrodes
    
    this_subject = electrode_list.subject{idx};
    this_label = electrode_list.label{idx};
    this_stimulation_electrode_label = electrode_list.stimulation_electrode_label(idx);
    
    subject_list_idx = strcmp(subject_list.subject,this_subject);
    
    subject_ictal_labels = subject_list.Ictal{subject_list_idx};
    subject_interictal_labels = subject_list.Interictal{subject_list_idx};
    
    if ~strcmp(subject_ictal_labels,'N')
        subject_ictal_labels = strsplit(subject_ictal_labels, ',');
        if ismember(this_label,subject_ictal_labels)
            ictal(idx) = true;
        end
        if ismember(this_stimulation_electrode_label,subject_ictal_labels)
            stimulation_site_ictal(idx) = true;
        end
    end
    
    if ~strcmp(subject_interictal_labels,'N')
        subject_interictal_labels = strsplit(subject_interictal_labels,',');
        if ismember(this_label,subject_interictal_labels)
            interictal(idx) = true;
        end
        if ismember(this_stimulation_electrode_label,subject_interictal_labels)
            stimulation_site_interictal(idx) = true;
        end
    end
end
electrode_list.ictal = ictal;
electrode_list.interictal = interictal;
electrode_list.stimulation_site_ictal = stimulation_site_ictal;
electrode_list.stimulation_site_interictal = stimulation_site_interictal;

%%% Save subject list
save(fullfile(list_directory,'subject_list.mat'),'subject_list');

%%% Save session list
save(fullfile(list_directory,'session_list.mat'),'session_list');

%%% Save electrode list
save(fullfile(list_directory,'electrode_list.mat'),'electrode_list');
end
