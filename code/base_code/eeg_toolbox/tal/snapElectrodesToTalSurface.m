function bipolarSubjElecs = snapElectrodesToTalSurface(subject, doplot)

if ~exist('doplot','var')
    doplot = false;
end

% Load in the tal brain
picpath = fileparts(which('tal3d'));
picfile = fullfile(picpath,'mni_cortical_surface.mat');
load(picfile);  

% separate into left and right vertices
left_vertices = v(v(:,1)<0,:);
right_vertices = v(v(:,1)>0,:);

% This is for plotting purposes
indices = (1:size(v));
indices(v(:,1)<0)=[];
f_mask = any(ismember(f, indices),2);
f_left = f(~f_mask,:);
f_right = f(f_mask,:);
  
% Get the electrodes for the subject
bipolarSubjElecs = getBipolarSubjElecs(subject,true);
subj_xyz = [[bipolarSubjElecs.x]',[bipolarSubjElecs.y]',[bipolarSubjElecs.z]']; 

% Split the electrodes into left and right halves
left_locs = cellfun(@(x)~isempty(strfind(x,'Left')), {bipolarSubjElecs.Loc1});
right_locs = cellfun(@(x)~isempty(strfind(x,'Right')), {bipolarSubjElecs.Loc1});

% for plotting
left_xyz = subj_xyz(left_locs,:);
right_xyz = subj_xyz(right_locs,:);
    
% Snap electrodes on left half to left surface
% and right half to right surface
snapped_elecLocs_left = dsearchn(double(left_vertices), double(subj_xyz));
snapped_elecLocs_left = left_vertices(snapped_elecLocs_left,:);
% This field contains all electrodes, but only the left ones are snapped
% correctly. 
left_snapped_all_elecs = snapped_elecLocs_left;

snapped_elecLocs_right = dsearchn(double(right_vertices), double(subj_xyz));
snapped_elecLocs_right = right_vertices(snapped_elecLocs_right,:);
% This contains all electrodes, but only the right ones are snapped
% correctly.
right_snapped_all_elecs = snapped_elecLocs_right;


all_elecs_snapped = left_snapped_all_elecs;
all_elecs_snapped(right_locs) = right_snapped_all_elecs(right_locs);
x_snapped = all_elecs_snapped(:,1);
y_snapped = all_elecs_snapped(:,2);
z_snapped = all_elecs_snapped(:,3);
[bipolarSubjElecs.x_snapped] = deal(nan);
[bipolarSubjElecs.y_snapped] = deal(nan);
[bipolarSubjElecs.z_snapped] = deal(nan);
for i=1:length(x_snapped)
    bipolarSubjElecs(i).x_snapped = x_snapped(i);
    bipolarSubjElecs(i).y_snapped = y_snapped(i);
    bipolarSubjElecs(i).z_snapped = z_snapped(i);
end

if doplot
    plot3_wrapper(all_elecs_snapped,30,'g');
    hold all
    
    plot3_wrapper(left_xyz, 15, 'b');
    plot3_wrapper(right_xyz, 15, 'b');
    v_colors = repmat([.5,.5,.5],size(v,1),1);
    patch('faces',f_left,'vertices',v,'edgecolor','none','FaceColor',...
        'interp','FaceVertexCData',v_colors);
    patch('faces',f_right,'vertices',v,'edgecolor','none','FaceColor',...
        'interp','FaceVertexCData',v_colors);
end
