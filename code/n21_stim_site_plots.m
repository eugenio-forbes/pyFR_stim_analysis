function n21_stim_site_plots(varargin)
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
table_directory = fullfile(analysis_directory,'tables');
plots_directory = fullfile(analysis_directory,'plots');
lock_directory = fullfile(analysis_directory,'locks');
lme_directory = fullfile(analysis_directory,'lme_results');
resources_directory = fullfile(analysis_directory,'resources');

%%% Declare regions of interest and periods of interest
stimulation_groups = {'IPL','retrosplenial'};
periods = {'encoding';'retrieval'};
[Ax,Bx] = ndgrid(1:numel(periods),1:numel(stimulation_groups));
n_combos = length(Ax(:));
stimulation_groups = stimulation_groups(Bx(:));
periods = periods(Ax(:));

%%%Initialize parpool
pool_object = gcp('nocreate');
if isempty(pool_object)
    parpool(24)
end

parfor idx = 1:n_combos
    period = periods{idx};
    stimulation_group = stimulation_groups{idx};
    
    folder_name = 'stim_site_plots';
    this_lme_directory = fullfile(lme_directory,folder_name);
    this_lock_directory = fullfile(lock_directory,folder_name);
    this_plot_directory = fullfile(plots_directory,folder_name);
    if ~isfolder(this_lme_directory)
        mkdir(this_lme_directory);
    end
    if ~isfolder(this_lock_directory)
        mkdir(this_lock_directory);
    end
    if ~isfolder(this_plot_directory)
        mkdir(this_plot_directory);
    end
    
    file_name = sprintf('%s_%s',period,stimulation_group);
    lock_file = fullfile(this_lock_directory,[file_name '_lock.txt']);
    done_file = fullfile(this_lock_directory,[file_name '_done.txt']);
    error_file = fullfile(this_lock_directory,[file_name '_error.txt']);
    plot_file = fullfile(this_plot_directory,file_name);
    
    pause(rand*5)
    
    if isfile(lock_file)
        fid = fopen(lock_file,'w'); fclose(fid);
        try
            electrode_table_file = sprintf('%s_stim_site_table',period);
            electrode_table = load(fullfile(table_directory,[electrode_table_file '.mat']));
            electrode_table = electrode_table.(electrode_table_file);
            
%             electrode_table = check_bad_subjects(analysis_directory,electrode_table);
            is_IPL = logical(electrode_table.stimulation_lateral);
            if strcmp(stimulation_group,'IPL')
                same_group = is_IPL;
            else
                same_group = ~is_IPL;
            end
            electrode_table = electrode_table(same_group,:);
            
            empty_cells = cellfun(@isempty,electrode_table.stimulation_coordinates);
            electrode_table(empty_cells,:) = [];
            
            is_left = logical(electrode_table.stimulation_left);
            left_electrodes = electrode_table(is_left,:);
            right_electrodes = electrode_table(~is_left,:);
            
            if sum(is_left)>0
                plot_stim_sites(resources_directory,plot_file,left_electrodes,stimulation_group,'left');
            end
            if sum(~is_left)>0
                plot_stim_sites(resources_directory,plot_file,right_electrodes,stimulation_group,'right')
            end
            
            fid = fopen(done_file,'w'); fclose(fid);
        catch this_error
            error_message = getReport(this_error, 'extended', 'hyperlinks', 'off');
            fid = fopen(error_file,'w');
            fprintf(fid,error_message);
            fclose(fid);
        end
    end
end

end

function plot_stim_sites(resources_directory,plot_file,electrode_table,stimulation_group,hemisphere)
this_plot_file = [plot_file '_' hemisphere];
figure_width = 200;
figure_height = 200;
plot_position = [0 0 figure_width figure_height];
marker_size = 40;
n_electrodes = height(electrode_table);

mapx = makecolormap_EF('uniform1');

load(fullfile(resources_directory,'EC.mat'),'EC');
load(fullfile(resources_directory,'SURF.mat'),'surf');

z_scores = double(electrode_table.slow_theta_power)/1000;
[~,sorting_indices] = sortrows(abs(z_scores),'ascend');
z_scores = z_scores(sorting_indices);
electrode_table = electrode_table(sorting_indices,:);

color_indices = round((z_scores*2000)+3000);
under1 = color_indices < 1;
over6000 = color_indices > 6000;
color_indices(under1) = ones(sum(under1),1);
color_indices(over6000) = ones(sum(over6000),1)*6000;
colors = cell2mat(arrayfun(@(x) mapx(x,:),color_indices,'UniformOutput',false));

coordinates = vertcat(electrode_table.stimulation_coordinates{:});
left_IPL = strcmp(hemisphere,'left') && strcmp(stimulation_group,'IPL');
right_retrosplenial = strcmp(hemisphere,'right') && strcmp(stimulation_group,'retrosplenial');
if left_IPL || right_retrosplenial
    x = linspace(-60,-61,n_electrodes);
    azimuth = -90;
else
    x = linspace(60,61,n_electrodes);
    azimuth = 90;
end
y = coordinates(:,2);
z = coordinates(:,3);
elevation = 0;

figure('Units','pixels','Visible','off','Position',plot_position)
%%% 3D plot of a sample brain
hold on
axes_brain = trisurf(surf.tri,surf.coord(1,:),...
    surf.coord(2,:),surf.coord(3,:),'EdgeColor','none');
whitebg(gcf,EC.bak.color);
set(gcf,'Color',EC.bak.color,'InvertHardcopy','off');
eval(['material ',EC.glb.material,';']);
eval(['shading ',EC.glb.shading,';']);
set(axes_brain,'FaceColor',[0.95, 0.95, 0.95]);
set(axes_brain,'FaceAlpha',0.15);
daspect([1 1 1])
hold off

%%% 3D scatter plot of electrodes with colors corresponding to t-statistics
hold on
scatter3(x,y,z,marker_size,colors,'filled','Marker','o')
rotate3d on
axis tight 
axis vis3d off
view(azimuth,elevation)
lighting phong
camlight(hemisphere)
hold off

print(this_plot_file,'-dsvg')
close all
end