function n07_make_exclusions(varargin)
if isempty(varargin)
    %%% Directory information
    root_directory = '/directory/to/pyFR_stim_analysis';
else
    root_directory = varargin{1};
end

%%% List directories
list_directory = fullfile(root_directory,'lists');
exclusion_directory = fullfile(root_directory,'exclusion_lists');

%%% Load subject list and excluded subjects so far
load(fullfile(list_directory,'subject_list.mat'),'subject_list');
load(fullfile(exclusion_directory,'excluded_subjects.mat'),'excluded_subjects');

%%% Load session list and excluded sessions so far
load(fullfile(list_directory,'session_list.mat'),'session_list');
load(fullfile(exclusion_directory,'excluded_sessions.mat'),'excluded_sessions');

%%% Load electrode list and excluded electrodes so far
load(fullfile(list_directory,'electrode_list.mat'),'electrode_list');
load(fullfile(exclusion_directory,'excluded_electrodes.mat'),'excluded_electrodes');

%%% Perform exclusions of subjects, sessions, electrodes, based on whether:

%%% Whether subject has previously had neursurgical procedure with resection
subject_exclusion_criterion = subject_list.Resection; 

new_excluded_subjects = subject_list(subject_exclusion_criterion,{'subject','subject_ID'});
reason_for_exclusion = repelem({'Previous resection'},height(new_excluded_subjects),1);
new_excluded_subjects.reason_for_exclusion = reason_for_exclusion;

new_excluded_sessions = session_list(ismember(session_list.subject_ID,new_excluded_subjects.subject_ID),{'subject','task','session','subject_ID','session_ID'});
reason_for_exclusion = repelem({'Previous resection'},height(new_excluded_sessions),1);
new_excluded_sessions.reason_for_exclusion = reason_for_exclusion;

excluded_electrodes = electrode_list(ismember(electrode_list.subject_ID,new_excluded_subjects.subject_ID),{'subject','task','session','channel_number','subject_ID','session_ID','electrode_ID'});
reason_for_exclusion = repelem({'Previous resection'},height(excluded_electrodes),1);

excluded_electrodes.reason_for_exclusion = reason_for_exclusion;
excluded_subjects = [excluded_subjects;new_excluded_subjects]; clear new_excluded_subjects
excluded_sessions = [excluded_sessions;new_excluded_sessions]; clear new_excluded_sessions

%%% Whether subject had previously undergone depth/strip/grid electrode implantation prior to current sEEG monitoring.
subject_exclusion_criterion = subject_list.Reimplant; 

new_excluded_subjects = subject_list(subject_exclusion_criterion,{'subject','subject_ID'});
reason_for_exclusion = repelem({'Reimplantation'},height(new_excluded_subjects),1);
new_excluded_subjects.reason_for_exclusion = reason_for_exclusion;

new_excluded_sessions = session_list(ismember(session_list.subject_ID,new_excluded_subjects.subject_ID),{'subject','task','session','subject_ID','session_ID'});
reason_for_exclusion = repelem({'Reimplantation'},height(new_excluded_sessions),1);
new_excluded_sessions.reason_for_exclusion = reason_for_exclusion;

new_excluded_electrodes = electrode_list(ismember(electrode_list.subject_ID,new_excluded_subjects.subject_ID),{'subject','task','session','channel_number','subject_ID','session_ID','electrode_ID'});
reason_for_exclusion = repelem({'Reimplantation'},height(new_excluded_electrodes),1);
new_excluded_electrodes.reason_for_exclusion = reason_for_exclusion;

excluded_subjects = [excluded_subjects;new_excluded_subjects]; clear new_excluded_subjects
excluded_sessions = [excluded_sessions;new_excluded_sessions]; clear new_excluded_sessions
excluded_electrodes = [excluded_electrodes;new_excluded_electrodes]; clear new_excluded_electrodes;

%%% Whether stimulation site was partially contained within white matter,
%%% or partially extruded into sulcus, lateral ventricle, fissure, or resection.
subject_exclusion_criterion = contains(subject_list.stimulation_neurologist_label,{'WM','PCS','STS','LV','sulcus','fissure','resection','encephalomalacia'},'IgnoreCase',true);

new_excluded_subjects = subject_list(subject_exclusion_criterion,{'subject','subject_ID'});
reason_for_exclusion = repelem({'Outlying stimulation site'},height(new_excluded_subjects),1);
new_excluded_subjects.reason_for_exclusion = reason_for_exclusion;

new_excluded_sessions = session_list(ismember(session_list.subject_ID,new_excluded_subjects.subject_ID),{'subject','task','session','subject_ID','session_ID'});
reason_for_exclusion = repelem({'Outlying stimulation site'},height(new_excluded_sessions),1);
new_excluded_sessions.reason_for_exclusion = reason_for_exclusion;

new_excluded_electrodes = electrode_list(ismember(electrode_list.subject_ID,new_excluded_subjects.subject_ID),{'subject','task','session','channel_number','subject_ID','session_ID','electrode_ID'});
reason_for_exclusion = repelem({'Outlying stimulation site'},height(new_excluded_electrodes),1);
new_excluded_electrodes.reason_for_exclusion = reason_for_exclusion;

excluded_subjects = [excluded_subjects;new_excluded_subjects]; clear new_excluded_subjects
excluded_sessions = [excluded_sessions;new_excluded_sessions]; clear new_excluded_sessions
excluded_electrodes = [excluded_electrodes;new_excluded_electrodes]; clear new_excluded_electrodes;

%%% Whether stimulation site was ictal or interictal
subject_exclusion_criterion = subject_list.stimulation_site_ictal | subject_list.stimulation_site_interictal;

new_excluded_subjects = subject_list(subject_exclusion_criterion,{'subject','subject_ID'});
reason_for_exclusion = repelem({'Stimulation site ictal or interictal'},height(new_excluded_subjects),1);
new_excluded_subjects.reason_for_exclusion = reason_for_exclusion;

new_excluded_sessions = session_list(ismember(session_list.subject_ID,new_excluded_subjects.subject_ID),{'subject','task','session','subject_ID','session_ID'});
reason_for_exclusion = repelem({'Stimulation site ictal or interictal'},height(new_excluded_sessions),1);
new_excluded_sessions.reason_for_exclusion = reason_for_exclusion;

new_excluded_electrodes = electrode_list(ismember(electrode_list.subject_ID,new_excluded_subjects.subject_ID),{'subject','task','session','channel_number','subject_ID','session_ID','electrode_ID'});
reason_for_exclusion = repelem({'Stimulation site ictal or interictal'},height(new_excluded_electrodes),1);
new_excluded_electrodes.reason_for_exclusion = reason_for_exclusion;

excluded_subjects = [excluded_subjects;new_excluded_subjects]; clear new_excluded_subjects
excluded_sessions = [excluded_sessions;new_excluded_sessions]; clear new_excluded_sessions
excluded_electrodes = [excluded_electrodes;new_excluded_electrodes]; clear new_excluded_electrodes;

%%% Whether session was labeled with error
session_exclusion_criterion = contains(session_list.session,{'error'},'IgnoreCase',true);

new_excluded_sessions = session_list(session_exclusion_criterion,{'subject','task','session','subject_ID','session_ID'});
reason_for_exclusion = repelem({'session error'},height(new_excluded_sessions),1);
new_excluded_sessions.reason_for_exclusion = reason_for_exclusion;

new_excluded_subjects = subject_list(ismember(subject_list.subject_ID,new_excluded_sessions.subject_ID),{'subject','subject_ID'});
reason_for_exclusion = repelem({'session error'},height(new_excluded_subjects),1);
new_excluded_subjects.reason_for_exclusion = reason_for_exclusion;

new_excluded_electrodes = electrode_list(ismember(electrode_list.session_ID,new_excluded_sessions.session_ID),{'subject','task','session','channel_number','subject_ID','session_ID','electrode_ID'});
reason_for_exclusion = repelem({'session error'},height(new_excluded_electrodes),1);
new_excluded_electrodes.reason_for_exclusion = reason_for_exclusion;

excluded_subjects = [excluded_subjects;new_excluded_subjects]; clear new_excluded_subjects
excluded_sessions = [excluded_sessions;new_excluded_sessions]; clear new_excluded_sessions
excluded_electrodes = [excluded_electrodes;new_excluded_electrodes]; clear new_excluded_electrodes;

%%% Whether electrode was partially within amygdala, entorhinal cortex,
%%% parahippocampal gyrus, white matter, lateral ventricle according to
%%% neurologist label
electrode_exclusion_criterion = contains(electrode_list.label,{'WM','LV','ent','PHG','amyg'},'IgnoreCase',true);
new_excluded_electrodes = electrode_list(electrode_exclusion_criterion,{'subject','task','session','channel_number','subject_ID','session_ID','electrode_ID'});
reason_for_exclusion = repelem({'Outlying hippocampal contacts'},height(new_excluded_electrodes),1);
new_excluded_electrodes.reason_for_exclusion = reason_for_exclusion;
excluded_electrodes = [excluded_electrodes;new_excluded_electrodes];

%%% Remove rows of lists based on exclusion IDs
subject_list = subject_list(~ismember(subject_list.subject_ID,excluded_subjects.subject_ID),:);
session_list = session_list(~ismember(session_list.session_ID,excluded_sessions.session_ID),:);
electrode_list = electrode_list(~ismember(electrode_list.subject_ID,excluded_subjects.subject_ID),:);
electrode_list = electrode_list(~ismember(electrode_list.session_ID,excluded_sessions.session_ID),:);
electrode_list = electrode_list(~ismember(electrode_list.electrode_ID,new_excluded_electrodes.electrode_ID),:); clear new_excluded_electrodes

%%% Loop through subjects to see whether there where any remaining
%%% hippocampal contacts after exclusions
n_subjects = height(subject_list);
no_hippocampal_contacts = false(n_subjects,1);
for idx = 1:n_subjects
    this_subject = subject_list.subject{idx};
    these_electrodes = electrode_list(strcmp(electrode_list.subject,this_subject),:);
    if isempty(these_electrodes)
        no_hippocampal_contacts(idx) = true;
    end
end
if sum(no_hippocampal_contacts) > 0
    new_excluded_subjects = subject_list(no_hippocampal_contacts,{'subject','subject_ID'});
    reason_for_exclusion = repelem({'Excluded all contacts'},height(new_excluded_subjects),1);
    new_excluded_subjects.reason_for_exclusion = reason_for_exclusion;
    
    new_excluded_sessions = session_list(ismember(session_list.subject_ID,new_excluded_subjects.subject_ID),{'subject','task','session','subject_ID','session_ID'});
    reason_for_exclusion = repelem({'Excluded all contacts'},height(new_excluded_sessions),1);
    new_excluded_sessions.reason_for_exclusion = reason_for_exclusion;
    
    excluded_subjects = [excluded_subjects;new_excluded_subjects]; clear new_excluded_subjects
    excluded_sessions = [excluded_sessions;new_excluded_sessions]; clear new_excluded_sessions
    subject_list = subject_list(~ismember(subject_list.subject_ID,excluded_subjects.subject_ID),:);
    session_list = session_list(~ismember(session_list.session_ID,excluded_sessions.session_ID),:);
end

%%% Loop through subjects and count how many hippocampal electrodes were ictal or
%%% interictal and whether stimulation site was ictal or interictal
n_subjects = height(subject_list);
n_ictal = zeros(n_subjects,1);
n_interictal = zeros(n_subjects,1);
n_electrodes = zeros(n_subjects,1);
stimulation_site_ictal = false(n_subjects,1);
stimulation_site_interictal = false(n_subjects,1);

for idx = 1:n_subjects
    this_subject = subject_list.subject{idx};
    
    subject_electrode_indices = strcmp(electrode_list.subject,this_subject);
    these_electrodes = electrode_list(subject_electrode_indices,:);
    
    if ~isempty(these_electrodes)
        [~,unique_idx,~] = unique(these_electrodes.electrode_ID);
        these_electrodes = these_electrodes(unique_idx,:);
        n_ictal(idx) = sum(these_electrodes.ictal);
        n_interictal(idx) = sum(these_electrodes.interictal);
        n_electrodes(idx) = height(these_electrodes);
        stimulation_site_ictal(idx) = these_electrodes.stimulation_site_ictal(1);
        stimulation_site_interictal(idx) = these_electrodes.stimulation_site_interictal(1);
    else
        subject_ictal_labels = subject_list.Ictal{idx};
        subject_interictal_labels = subject_list.Interictal{idx};
        this_stimulation_electrode_label = subject_list.stimulation_electrode_label(idx);
        if ~strcmp(subject_ictal_labels,'N')
            subject_ictal_labels = strsplit(subject_ictal_labels, ',');
            if ismember(this_stimulation_electrode_label,subject_ictal_labels)
                stimulation_site_ictal(idx) = true;
            end
        end
        if ~strcmp(subject_interictal_labels,'N')
            subject_interictal_labels = strsplit(subject_interictal_labels,',');
            if ismember(this_stimulation_electrode_label,subject_interictal_labels)
                stimulation_site_interictal(idx) = true;
            end
        end
    end
end
subject_list.n_ictal = n_ictal;
subject_list.n_interictal = n_interictal;
subject_list.n_electrodes = n_electrodes;
subject_list.stimulation_site_ictal = stimulation_site_ictal;
subject_list.stimulation_site_interictal = stimulation_site_interictal;

%%% Loop through sessions to transfer information
n_sessions = height(session_list);
n_ictal = zeros(n_sessions,1);
n_interictal = zeros(n_sessions,1);
n_electrodes = zeros(n_sessions,1);
stimulation_site_ictal = false(n_sessions,1);
stimulation_site_interictal = false(n_sessions,1);

for idx = 1:n_sessions
    this_subject = session_list.subject{idx};
    
    subject_list_idx = strcmp(subject_list.subject,this_subject);
    
    n_ictal(idx) = subject_list.n_ictal(subject_list_idx);
    n_interictal(idx) = subject_list.n_interictal(subject_list_idx);
    n_electrodes(idx) = subject_list.n_electrodes(subject_list_idx);
    stimulation_site_ictal(idx) = subject_list.stimulation_site_ictal(subject_list_idx);
    stimulation_site_interictal(idx) = subject_list.stimulation_site_interictal(subject_list_idx);  
end

session_list.n_ictal = n_ictal;
session_list.n_interictal = n_interictal;
session_list.n_electrodes = n_electrodes;
session_list.stimulation_site_ictal = stimulation_site_ictal;
session_list.stimulation_site_interictal = stimulation_site_interictal;

%%%Save all lists
save(fullfile(list_directory,'subject_list.mat'),'subject_list');
save(fullfile(list_directory,'session_list.mat'),'session_list');
save(fullfile(list_directory,'electrode_list.mat'),'electrode_list');
save(fullfile(exclusion_directory,'excluded_subjects.mat'),'excluded_subjects');
save(fullfile(exclusion_directory,'excluded_sessions.mat'),'excluded_sessions');
save(fullfile(exclusion_directory,'excluded_electrodes.mat'),'excluded_electrodes');
end
