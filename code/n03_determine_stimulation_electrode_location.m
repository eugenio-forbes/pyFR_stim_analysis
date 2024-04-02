function n03_determine_stimulation_electrode_location(varargin)
if isempty(varargin)
    %%% Directory information
    root_directory = '/directory/to/pyFR_stim_analysis';
else
    root_directory = varargin{1};
end

%%% List directories
localization_directory = fullfile(root_directory,'iEEGxfMRI/Pipeline/5crossref');
coordinates_directory = fullfile(root_directory,'iEEGxfMRI/Pipeline/6finalize');
list_directory = fullfile(root_directory,'lists');

%%% Option to read MNI coordinates .csv file
coordinate_opts = delimitedTextImportOptions("NumVariables", 5);
coordinate_opts.DataLines = [2, Inf];
coordinate_opts.Delimiter = ",";
coordinate_opts.VariableNames = ["ElecNumber", "X", "Y", "Z", "Regions"];
coordinate_opts.VariableTypes = ["double", "double", "double", "double", "categorical"];
coordinate_opts.ExtraColumnsRule = "ignore";
coordinate_opts.EmptyLineRule = "read";
coordinate_opts = setvaropts(coordinate_opts, "Regions", "EmptyFieldRule", "auto");

%%% Load session list
load(fullfile(list_directory,'session_list.mat'),'session_list');

%%% Loop through session list and based on stimulation electrode label
%%% get channel number, neurologist localization, automatic localization,
%%% and Glasser label. Also get MNI RAS x,y,z coordinates of electrode.

n_sessions = height(session_list);

stimulation_electrode_channel = NaN(n_sessions,1);
stimulation_neurologist_label = cell(n_sessions,1);
stimulation_automatic_label = cell(n_sessions,1);
stimulation_glasser_label = cell(n_sessions,1);
stimulation_coordinates = cell(n_sessions,1);
stimulation_references = cell(n_sessions,1);
stimulation_reference_type = cell(n_sessions,1);

for idx = 1:n_sessions
    
    subject = session_list.subject{idx};
    subject = subject(1:5); %some subjects have an extra letter at the end
    stimulation_electrode_label = session_list.stimulation_electrode_label{idx};
    
    %%% Neurologist localization labels found in depth electrode information file
    depth_electrode_information_file = fullfile(localization_directory,sprintf('sub-%s/%s_depth_el_info.mat',subject,subject));
    %%% Automatic localization/Glasser labels found in local atlas information
    local_atlas_information_file = fullfile(localization_directory,sprintf('sub-%s/%s_local_atlas_info.mat',subject,subject));
    %%% MNI RAS coordinates file
    MNI_coordinates_file = fullfile(coordinates_directory,sprintf('sub-%s/sub-%s_MNIRAS.csv',subject,subject));
    
    if isfile(depth_electrode_information_file)
        
        load(depth_electrode_information_file,'depth_el_info')
        
        if isstruct(depth_el_info) %%% Some save as cell matrix, and some saved as struct
            depth_el_info = struct2table(depth_el_info);
            channel_numbers = depth_el_info.elec;
            labels = depth_el_info.contact;
            locations = depth_el_info.recon_label;
        elseif iscell(depth_el_info)
            channel_numbers = [depth_el_info{:,1}];
            labels = depth_el_info(:,2);
            locations = depth_el_info(:,3);
        end
        
        %%% Get row index based on labels and get channel number (channel
        %%% numbers skip)
        stimulation_electrode_idx = find(strcmp(labels,stimulation_electrode_label));
        
        if ~isempty(stimulation_electrode_idx)
            stimulation_electrode_channel(idx) = channel_numbers(stimulation_electrode_idx);
            stimulation_neurologist_label{idx} = locations{stimulation_electrode_idx};
        
            %%% Get indices of depth electrode info file of channels that
            %%% have the save label stem as stimulation site
            same_depth_electrode_indices = find(contains(labels,stimulation_electrode_label(1:2)));
            
            if ~isempty(same_depth_electrode_indices)
                same_depth_channel_numbers = channel_numbers(same_depth_electrode_indices);
                same_depth_locations = locations(same_depth_electrode_indices);
                same_depth_locations = strrep(same_depth_locations,' ','');
                same_depth_labels = labels(same_depth_electrode_indices);
                
                n_channels_depth_electrode = length(same_depth_labels);
                within_depth_electrode_idx = find(strcmp(same_depth_labels,stimulation_electrode_label));
                
                if any(strcmp(same_depth_locations,'WM')) %%% WM contact of same depth electrode if available
                    stimulation_reference_type{idx} = 'white-matter';
                    white_matter_indices = strcmp(same_depth_locations,'WM');
                    stimulation_references{idx} = same_depth_channel_numbers(white_matter_indices);
                else
                    stimulation_reference_type{idx} = 'bipolar'; %%% Otherwise get channel number of next adjacent channel
                    if within_depth_electrode_idx < n_channels_depth_electrode %Meaning it is not the most outermost channel
                        next_adjacent_idx = within_depth_electrode_idx + 1;
                    else
                        next_adjacent_idx = within_depth_electrode_idx - 1;
                    end
                    stimulation_references{idx} = same_depth_channel_numbers(next_adjacent_idx);
                end
                clear same_depth_channel_numbers same_depth_locations same_depth_labels
            end
        end
    end
    
    if isfile(local_atlas_information_file) && ~isnan(stimulation_electrode_channel(idx))
        load(local_atlas_information_file,'local_atlas_info')
        local_atlas_info = struct2table(local_atlas_info);
        
        channel_numbers = local_atlas_info.elec;
        automatic_labels = local_atlas_info.AAL_label;
        glasser_labels = local_atlas_info.Glasser_label;
        
        stimulation_electrode_idx = find(channel_numbers == stimulation_electrode_channel(idx));
        if ~isempty(stimulation_electrode_idx)
            stimulation_automatic_label{idx} = automatic_labels{stimulation_electrode_idx};
            stimulation_glasser_label{idx} = glasser_labels{stimulation_electrode_idx};
        end
    end
    
    if isfile(MNI_coordinates_file) && ~isnan(stimulation_electrode_channel(idx))
        MNI_coordinates_table = readtable(MNI_coordinates_file,coordinate_opts);
        channel_numbers = MNI_coordinates_table.ElecNumber;
        coordinates = [MNI_coordinates_table.X,MNI_coordinates_table.Y,MNI_coordinates_table.Z];
        
        stimulation_electrode_idx = channel_numbers == stimulation_electrode_channel(idx);
        if ~isempty(stimulation_electrode_idx)
            stimulation_coordinates{idx} = coordinates(stimulation_electrode_idx,:);
        end
    end
end

session_list.stimulation_electrode_channel = stimulation_electrode_channel;
session_list.stimulation_neurologist_label = stimulation_neurologist_label;
session_list.stimulation_automatic_label = stimulation_automatic_label;
session_list.stimulation_glasser_label = stimulation_glasser_label;
session_list.stimulation_coordinates = stimulation_coordinates;
session_list.stimulation_references = stimulation_references;
session_list.stimulation_reference_type = stimulation_reference_type;

%%% Some empty fields need to be manually added by looking at subject
%%% files, because of differences in labeling in the past or errors in localization pipeline.

%%% Save updated session list
save(fullfile(list_directory,'session_list.mat'),'session_list');
end
