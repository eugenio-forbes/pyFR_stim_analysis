function n05_make_electrode_list(varargin)
if isempty(varargin)
    %%% Directory information
    root_directory = '/directory/to/pyFR_stim_analysis';
else
    root_directory = varargin{1};
end

%%% List directories
list_directory = fullfile(root_directory,'lists');
localization_directory = fullfile(root_directory,'iEEGxfMRI/Pipeline/5crossref');
coordinates_directory = fullfile(root_directory,'iEEGxfMRI/Pipeline/6finalize');
subject_directory = fullfile(root_directory,'subject_files');

%%% Option to read MNI coordinates .csv file
coordinate_opts = delimitedTextImportOptions("NumVariables", 5);
coordinate_opts.DataLines = [2, Inf];
coordinate_opts.Delimiter = ",";
coordinate_opts.VariableNames = ["ElecNumber", "X", "Y", "Z", "Regions"];
coordinate_opts.VariableTypes = ["double", "double", "double", "double", "categorical"];
coordinate_opts.ExtraColumnsRule = "ignore";
coordinate_opts.EmptyLineRule = "read";
coordinate_opts = setvaropts(coordinate_opts, "Regions", "EmptyFieldRule", "auto");

%%% Load subject list file
load(fullfile(list_directory,'subject_list.mat'),'subject_list');

%%% Load session list
load(fullfile(list_directory,'session_list.mat'),'session_list');

%%% Loop through sessions and gather hippocampal channels (only those where
%%% neurologist and automatic localization match).
session_electrodes = cell(height(session_list),1);

for idx = 1:height(session_list)
    
    this_session = session_list(idx,:);
    subject = this_session.subject{:};
    
    %%% Load this subjects depth electrode information (neurologist
    %%% localization) and local atlas information (automatic).
    depth_electrode_information_file = fullfile(localization_directory,sprintf('sub-%s/%s_depth_el_info.mat',subject,subject));
    local_atlas_information_file = fullfile(localization_directory,sprintf('sub-%s/%s_local_atlas_info.mat',subject,subject));
    
    if isfile(depth_electrode_information_file)
        
        load(depth_electrode_information_file,'depth_el_info')
        
        %%% Get channel numbers, labels and locations
        if isstruct(depth_el_info) %%% Some save as cell matrix, and some saved as struct    
            depth_el_info = struct2table(depth_el_info);
            depth_channel_numbers = depth_el_info.elec;
            depth_labels = depth_el_info.contact;
            depth_locations = depth_el_info.recon_label;
        elseif iscell(depth_el_info)
            depth_channel_numbers = [depth_el_info{:,1}]';
            depth_labels = depth_el_info(:,2);
            depth_locations = depth_el_info(:,3);
        end
 
    else %%% Only one case were neither file was present in localization dir, so only getting neurologist localization
        %%% Would find text file version in subject directory
        depth_electrode_information_file = dir(fullfile(subject_directory,subject,'docs/*el_inf*.txt'));
   
        if ~isempty(depth_electrode_information_file)
            %%% Search for file given variable naming
            depth_file_name = {depth_electrode_information_file.name};
            depth_electrode_information_file = depth_electrode_information_file(~contains(depth_file_name,'._')); %%% MacOS acces makes these annoying copy files that start with ._
            
            if ~isempty(depth_electrode_information_file)
                depth_electrode_information_file = fullfile({depth_electrode_information_file.folder},{depth_electrode_information_file.name});
                depth_electrode_information_file = depth_electrode_information_file{1};
                
                %%% Read text file as table
                opts = detectImportOptions(depth_electrode_information_file,'Delimiter',' ');
                depth_el_info = readtable(depth_electrode_information_file,opts);
                
                depth_channel_numbers = depth_el_info{:,1};
                depth_labels = depth_el_info{:,2};
                depth_locations = arrayfun(@(x) strjoin(depth_el_info{x,3:end}),(1:height(depth_el_info))','UniformOutput',false); %%% Labels in this case will be merged because file is file is space delimited
            end
        end
        
    end
    clear depth_electrode_information_file
    
    %%% Find location matches for hippocampus
    hippocampus_hits = contains(depth_locations,'hip','IgnoreCase',true) & ~(contains(depth_locations,'para','IgnoreCase',true) & ~contains(depth_locations,{'vs','/'},'IgnoreCase',true));
    
    %%% Use automatic location for initial matches
    if isfile(local_atlas_information_file) %Missing for a few subjects
        load(local_atlas_information_file,'local_atlas_info');
        local_atlas_info = struct2table(local_atlas_info);
        automatic_channel_numbers = local_atlas_info.elec;
        automatic_locations = local_atlas_info.AAL_label;
        
        automatic_hits = contains(automatic_locations,'hip','IgnoreCase',true) & ~contains(automatic_locations,'para','IgnoreCase',true);
        
        automatic_channel_numbers = automatic_channel_numbers(automatic_hits);
        
        hippocampus_hits = hippocampus_hits | ismember(depth_channel_numbers,automatic_channel_numbers);
        
        clear automatic_locations local_atlas_info
    end
    clear local_atlas_information_file
    
    %%% Filter channel numbers, labels, and locations
    channel_number = depth_channel_numbers(hippocampus_hits);
    label = depth_labels(hippocampus_hits);
    location = depth_locations(hippocampus_hits);
    
    %%% If there were matches, add information about longitudinal location,
    %%% hemisphere, reference channel number and type, and MNI coordinates
    
    if ~isempty(channel_number)
        
        longitudinal_location = cell(length(channel_number),1);
        hemisphere = cell(length(channel_number),1);
        reference_channels = cell(length(channel_number),1);
        reference_type = cell(length(channel_number),1);
        coordinates = cell(length(channel_number),1);
        
        %%% MNI RAS coordinates file
        MNI_coordinates_file = fullfile(coordinates_directory,sprintf('sub-%s/sub-%s_MNIRAS.csv',subject,subject));
        if isfile(MNI_coordinates_file)
            MNI_coordinates_table = readtable(MNI_coordinates_file,coordinate_opts);
            MNI_channel_numbers = MNI_coordinates_table.ElecNumber;
            MNI_coordinates = [MNI_coordinates_table.X,MNI_coordinates_table.Y,MNI_coordinates_table.Z];
        end
        
        %%% Loop through channels to get information specific for each
        for jdx = 1:length(channel_number)
            
            this_label = label{jdx};
            this_location = location{jdx};
            this_channel_number = channel_number(jdx);
            
            if contains(this_location,{'ant','head'},'IgnoreCase',true)
                longitudinal_location{jdx} = 'anterior';
            elseif contains(this_location,{'pos','tail'},'IgnoreCase',true)
                longitudinal_location{jdx} = 'posterior';
            else
                longitudinal_location{jdx} = 'check';
            end
            
            if strcmp(this_label(1),'L')
                hemisphere{jdx} = 'left';
            elseif strcmp(this_label(1),'R')
                hemisphere{jdx} = 'right';
            else
                hemisphere{jdx} = 'check';
            end
            
            %%% Get indices of depth electrode info file of channels that
            %%% have the save label stem as this channel.
            same_depth_electrode_indices = find(contains(depth_labels,this_label(1:2)));
            
            if ~isempty(same_depth_electrode_indices)
                same_depth_channel_numbers = depth_channel_numbers(same_depth_electrode_indices);
                same_depth_locations = depth_locations(same_depth_electrode_indices);
                same_depth_locations = strrep(same_depth_locations,' ','');
                same_depth_labels = depth_labels(same_depth_electrode_indices);
                
                n_channels_depth_electrode = length(same_depth_labels);
                within_depth_electrode_idx = find(strcmp(same_depth_labels,this_label));
                
                if any(strcmp(same_depth_locations,'WM')) %%% Get channel number of deepest WM contact of same depth electrode if available
                    reference_type{jdx} = 'white-matter';
                    white_matter_indices = strcmp(same_depth_locations,'WM');
                    reference_channels{jdx} = same_depth_channel_numbers(white_matter_indices);
                else
                    reference_type{jdx} = 'bipolar'; %%% Otherwise get channel number of next adjacent channel
                    if within_depth_electrode_idx < n_channels_depth_electrode %Meaning it is not the most outermost channel
                        next_adjacent_idx = within_depth_electrode_idx + 1;
                    else
                        next_adjacent_idx = within_depth_electrode_idx - 1;
                    end
                    reference_channels{jdx} = same_depth_channel_numbers(next_adjacent_idx);
                end
                clear same_depth_channel_numbers same_depth_locations same_depth_labels
            end
            
            %%% Get coordinates
            if isfile(MNI_coordinates_file)
                this_MNI_idx = MNI_channel_numbers == this_channel_number;
                coordinates{jdx} = MNI_coordinates(this_MNI_idx,:);
            end
        end
        clear MNI_coordinates_file depth_channel_numbers depth_locations depth_labels
        
        %%% Make table based on session information and add information
        these_electrodes = repmat(this_session,length(channel_number),1);
        these_electrodes.channel_number = channel_number;
        these_electrodes.label = label;
        these_electrodes.location = location;
        these_electrodes.longitudinal_location = longitudinal_location;
        these_electrodes.hemisphere = hemisphere;
        these_electrodes.reference_channels = reference_channels;
        these_electrodes.reference_type = reference_type;
        these_electrodes.coordinates = coordinates;
        
        session_electrodes{idx} = these_electrodes;
    end
    clear these_electrodes subject channel_number label location longitudinal_location hemisphere reference_channels reference_type coordinates
end

%%% Concatenate all session electrode tables
electrode_list = vertcat(session_electrodes{:});

%%% Give each electrode a unique low memory ID for making tables,
%%% for use as categorical variable in random effects of LME, an for easier
%%% identification and filtering of exclusions.
[~,~,electrode_ID] = unique(strcat(electrode_list.subject,electrode_list.label));
electrode_list.electrode_ID = int16(electrode_ID);

%%% Save lists
save(fullfile(list_directory,'subject_list.mat'),'subject_list');
save(fullfile(list_directory,'session_list.mat'),'session_list');
save(fullfile(list_directory,'electrode_list.mat'),'electrode_list');
end
