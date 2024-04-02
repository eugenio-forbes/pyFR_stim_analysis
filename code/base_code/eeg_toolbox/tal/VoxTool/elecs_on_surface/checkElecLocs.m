function checkElecLocs(subjList, startAt)
%CHECKELECLOCS(subjList, [startAt])
% Checks the electrode locations for a subject, showing diagnostic plots
%   subjList ----- the subject(s) to check
%   startAt  ----- The subject in subjList where it should start

close all

%which subjects
if ~exist('subjList','var') || isempty(subjList)
    get_subs_all;
elseif ~isa(subjList,'cell')
    subjList = {subjList};
end


% start at
if ~exist('startAt','var')
    currVal = 1;
else
    currVal = find(strcmp(startAt, subjList));
end
startVal = currVal;

% flags
%processing
%convMotherFlag = 0;
%snapFlag = 0;

%plot
plotFlag =1;
plotCompFlag = 0;
plotRawFlag = 1;
plotIndivFlag = 1;
plotAvgFlag = 1;
plotCompTalOld = 1;
plotCompTalAvg = 1;

if plotRawFlag
    rawPlot = figure('position',[992        1045         762         457]);
end
if plotIndivFlag
    indivPlot = figure('position', [2070        1027         770         461]);
end
if plotAvgFlag
    avgPlot = figure('position',[2851        1026         775         462]);
end
if plotCompFlag
    compPlot = figure;
end
if plotCompTalOld
    talOldPlot = figure('position', [823   132   634   450]);
    old_events = load('/data/eeg/tal/allTalLocs_GM.mat');
    old_events = old_events.events;

end
if plotCompTalAvg
    talAvgPlot = figure('position',[69   231   637   469]);
end

plotBrainFlag = 1;
checkDistsFlag = 1;

if checkDistsFlag
    distPlot = figure('position',[1814         252         728         549]);
end

% surf paths
surfDir = '/data/eeg/freesurfer/subjects/average/surf';
surfL = fullfile(surfDir,'lh.pial');[vL_avg,fL_avg] = read_surf_wrapper(surfL);
surfR = fullfile(surfDir,'rh.pial');[vR_avg,fR_avg] = read_surf_wrapper(surfR);
for s = startVal: length(subjList)
    fprintf(['\n' subjList{s}])
    try
%        if convMotherFlag
%            convertMotherToVoxCoords(subjList{s},[],0,1);
%        end
%         if snapFlag
%             updateTalStructWithSurfLocs(subjList{s},1);    
%         end

        % Get the bipolar tal structure
        talStruct = getBipolarSubjElecs(subjList{s},0);
        raw_xyz = cat(2,[talStruct.x]',[talStruct.y]',[talStruct.z]');
       
        % Get the indivSurf substructure
        tal_indiv = [talStruct.indivSurf];
        xyz_indiv = cat(2,[tal_indiv.x]',[tal_indiv.y]',[tal_indiv.z]');
        xyz_indiv_snap = cat(2,[tal_indiv.x_snap]',[tal_indiv.y_snap]',[tal_indiv.z_snap]');

        % get the avgSurf substructure
        tal_avg = [talStruct.avgSurf];
        xyz_avg = cat(2,[tal_avg.x]',[tal_avg.y]',[tal_avg.z]');
        xyz_avg_snap = cat(2,[tal_avg.x_eSnap]',[tal_avg.y_eSnap]',[tal_avg.z_eSnap]');

        if plotFlag 
            
            % plot compare flag
            if plotCompFlag
               figure(compPlot);clf; hold all
               plot3_wrapper(raw_xyz,10,'r');
               plot3_wrapper(xyz_avg,10,'b');
               plot3_wrapper(xyz_indiv,10,'g');
                
            end
            
            % plot raw
            if plotRawFlag && ~isnan(raw_xyz(1,1))
                figure(rawPlot);clf(rawPlot);hold off
                if plotBrainFlag
                    subplot(1,2,1)
                   tal3d([[1:size(raw_xyz,1)]' raw_xyz],[-90 0],[1 0 0],'num',[6], 0,0,1);
                   view([80 20]);
                   subplot(1,2,2);
                   tal3d([[1:size(raw_xyz,1)]' raw_xyz],[-90 0],[1 0 0],'num',[6], 0,0,1);
                   view([-80 -20]);
                else
                   plot3_wrapper(raw_xyz,20,'r');
                end
                title(subjList{s});
            elseif plotRawFlag && isnan(raw_xyz(1,1)) 
                figure(rawPlot);clf;hold off;
                title('NO RAW ELECS')
            end
            
            % plot indiv
            if plotIndivFlag && ~isnan(xyz_indiv(1,1))
               figure(indivPlot);hold off; clf;
               subplot(1,2,1)
               if plotBrainFlag 
                  [vL,fL] = read_surf_wrapper(tal_indiv(1).path2surfL);plotsurf_wrapper(vL,fL,[],false);
                  hold all;
                  [vR,fR] = read_surf_wrapper(tal_indiv(1).path2surfR);plotsurf_wrapper(vR,fR,[],false);
                  axis('off');
               end
               plot3_wrapper(xyz_indiv,20,'r');
               hold all;
               plot3_wrapper(xyz_indiv_snap,20,'b');
               title([subjList{s} ' INDIV']);
               view([80 20]);
               camlight
               
               subplot(1,2,2)
               if plotBrainFlag 
                  [vL,fL] = read_surf_wrapper(tal_indiv(1).path2surfL);plotsurf_wrapper(vL,fL,[],false);
                  hold all;
                  [vR,fR] = read_surf_wrapper(tal_indiv(1).path2surfR);plotsurf_wrapper(vR,fR,[],false);
                  axis('off');
               end
               plot3_wrapper(xyz_indiv,20,'r');
               hold all;
               plot3_wrapper(xyz_indiv_snap,20,'b');
               title([subjList{s} ' INDIV']);
               view([-80 -20]);
               camlight;
            elseif plotIndivFlag && isnan(xyz_indiv(1,1))
                figure(indivPlot); hold off;
                clf(indivPlot);
                title('NO INDIV ELECS');
            end
            
            % plot avg
            if plotAvgFlag && ~isnan(xyz_avg(1,1))
               figure(avgPlot);clf;hold off
               subplot(1,2,1)

               if plotBrainFlag
                  plotsurf_wrapper(vL_avg,fL_avg,[],false);
                  hold all;
                  plotsurf_wrapper(vR_avg,fR_avg,[],false);
                  axis('off');
               end
               plot3_wrapper(xyz_avg,20,'r');
               hold all
               plot3_wrapper(xyz_avg_snap,20,'b');
               title([subjList{s} ' AVG']);
               view([80 20]);
               camlight;
               
               subplot(1,2,2)
               if plotBrainFlag
                  plotsurf_wrapper(vL_avg,fL_avg,[],false);
                  hold all;
                  plotsurf_wrapper(vR_avg,fR_avg,[],false);
                  axis('off');
               end
               plot3_wrapper(xyz_avg,20,'r');
               hold all;
               plot3_wrapper(xyz_avg_snap,20,'b');
               title([subjList{s} ' AVG']);
               view([-80 -20]);
               camlight;
            elseif plotAvgFlag && isnan(xyz_avg(1,1))
                figure(avgPlot); hold off;
                clf(avgPlot);
                title('NO AVG ELECS');
            end
            
            if plotCompTalOld
                figure(talOldPlot); clf; hold off;
                compareTalToOld_local(subjList{s}, old_events, talStruct, false);
            end
            
            if plotCompTalAvg
                figure(talAvgPlot); clf; hold off;
                compareTalToAvg_local(talStruct );
            end
            
            % check distances for bipolar pairs
            if checkDistsFlag
                figure(distPlot); hold off; clf(distPlot)
                bp_distances = check_bp_distances(subjList{s});
                if ~any(abs(bp_distances-10)>10)

                    title(['NO BAD ELECS: ' subjList{s}]);
                else
                    title(subjList{s})
                end
            end
            pause;
        end

    catch e
       disp('**********ERROR*********')
       keyboard;
    end
    
    currVal = s; 
    disp(subjList{s})
end


function compareTalToAvg_local(talStruct)
avgSurf = [talStruct.avgSurf];
avg_xyz = [[avgSurf.x]; [avgSurf.y]; [avgSurf.z]]';
tal_xyz = [[talStruct.x]; [talStruct.y]; [talStruct.z]]';

  hold off;
 plot3_wrapper(avg_xyz,20,'blue','.');
 hold all
 plot3_wrapper(tal_xyz,5,'red','o');

 legend('Avg','Tal');
 
 distances = sqrt(sum((avg_xyz-tal_xyz).^2,2));
 stdDists = std(distances);
meanDists = mean(distances);

for i=1:length(distances)
    if distances(i)>meanDists+stdDists.*2
        
        plot3([avg_xyz(i,1);tal_xyz(i,1)], ...
            [avg_xyz(i,2);tal_xyz(i,2)], ...
            [avg_xyz(i,3);tal_xyz(i,3)], 'linestyle','-','color',[1 0 0]);
        text(avg_xyz(i,1),avg_xyz(i,2),avg_xyz(i,3), talStruct(i).tagName)
    else
        plot3([avg_xyz(i,1);tal_xyz(i,1)], ...
            [avg_xyz(i,2);tal_xyz(i,2)], ...
            [avg_xyz(i,3);tal_xyz(i,3)], 'linestyle','-','color',[.8 .8 .8]);
    end
end 

function compareTalToOld_local(subj, events, talStruct, bipol_flag)
 


 if bipol_flag
     try
         subj_events = load(['/data/eeg/tal/bipol/' subj '_BP_talairach_info.mat']);
         subj_events = subj_events.events;
     catch e
         title(sprintf('%s has no old coords',subj))
         return
     end
 else
     if strcmp(subj,'UP003')
         this_subj = 'UP003a';
     else
         this_subj = subj;
     end
     subj_events = events(strcmp({events.subject},this_subj));
     if isempty(subj_events)
         title(sprintf('%s has no old coords',subj));
         return
     end
 end
 
 if bipol_flag
     talEnames =  {talStruct.eNames};
 else
     talEnames = {talStruct.eName};
 end
 
 
 avgSurf = [talStruct.avgSurf];
 old_xyz = [[subj_events.x]; [subj_events.y]; [subj_events.z]]';
 new_xyz = [[talStruct.x]; [talStruct.y]; [talStruct.z]]';
% new_xyz = [[avgSurf.x]; [avgSurf.y]; [avgSurf.z]]';
 
   hold off;
 plot3_wrapper(old_xyz,20,'blue','.');
 hold all
 plot3_wrapper(new_xyz,8,'red','o');
 distances = nan(length(subj_events),1);
for i=1:length(subj_events)
    if bipol_flag
        elec_num = sprintf('%s-%s',subj_events(i).elec1, subj_events(i).elec2);
    else
        elec_num = num2str(subj_events(i).channel);
    end
    
    tal_i = (strcmp(talEnames, elec_num));
    if ~any(tal_i)
        fprintf('NONMATCHING ELEC: %s\n',elec_num)
        plot3_wrapper(old_xyz(i,:),20,'green');
        keyboard
        continue
    end
    distances(i) = sqrt(sum((old_xyz(i,:)-new_xyz(tal_i,:)).^2));
    

    plot3([old_xyz(i,1);new_xyz(tal_i,1)], ...
          [old_xyz(i,2);new_xyz(tal_i,2)], ...
          [old_xyz(i,3);new_xyz(tal_i,3)], 'linestyle','-','color',[.8 .8 .8]); 
end 

stdDists = std(distances);
meanDists = mean(distances);
bigDists = find((distances-meanDists)>stdDists*2);
%bigDists = find(distances>10);

 
 legend('Old','New');
 
for i=bigDists'
    if bipol_flag
        elec_num = sprintf('%s-%s',subj_events(i).elec1, subj_events(i).elec2);
    else
        elec_num = num2str(subj_events(i).channel);
    end
    
    tal_i = (strcmp(talEnames, elec_num));
    distances(i) = sqrt(sum((old_xyz(i,:)-new_xyz(tal_i,:)).^2));
    
    plot3([old_xyz(i,1);new_xyz(tal_i,1)], ...
          [old_xyz(i,2);new_xyz(tal_i,2)], ...
          [old_xyz(i,3);new_xyz(tal_i,3)], 'r-'); 
      
    text(old_xyz(i,1),old_xyz(i,2),old_xyz(i,3),elec_num)
end 
 


function bp_distances = check_bp_distances(subj)
% CHECK_BP_DISTANCES
% creates a plot of electrodes, with electrodes further than 16mm, or
% closer than 4mm, connected in red

talStruct = getBipolarSubjElecs(subj, false);
bpTalStruct = getBipolarSubjElecs(subj, true);

% bp_manual_file = fullfile('/data/eeg',subj,'docs/bp_manual.txt');
% fid = fopen(bp_manual_file);
% if fid==-1
%     fprintf('\nERROR: %s: no bp_manual file\n',subj)
%     bp_distances=[];
%     return
% end
% bp_manual = fscanf(fid,'%d %d %*s',[2,inf])';
% fclose(fid);

% try
%     bpPairs = load(fullfile('/data/eeg',subj,'tal/bpPairs.mat'));
%     bp_manual = bpPairs.allpairs;
% catch e
%     fprintf('\nERROR: %s: no bp_manual file\n',subj)
%     bp_distances=[];
%     return
% end

pairs = {bpTalStruct.channel};
bp_manual = cell2mat(pairs');
avgSurf = [talStruct.avgSurf];
xyz_all = [[avgSurf.x]', [avgSurf.y]', [avgSurf.z]'];
plot3_wrapper(xyz_all,10,'blue');
bp_avgSurf = [bpTalStruct.avgSurf];
xyz_bp = [[bp_avgSurf.x]', [bp_avgSurf.y]', [bp_avgSurf.z]'];
hold all;
plot3_wrapper(xyz_bp,2,'black','+');
[bp_distances, bad_tags] = getDistances(bp_manual, talStruct, xyz_all);


function [distances, bad_tags] = getDistances(bp_manual, talStruct, xyz)
distances = [];
bad_tags = {};
bp_manual = bp_manual';
for i=1:length(talStruct)
    [I,~] = find(bp_manual==talStruct(i).channel);
    index = find(bp_manual==talStruct(i).channel);
    index(I==1) = index(I==1)+1;
    index(I==2) = index(I==2)-1;
    relevant_xyz = xyz(ismember([talStruct.channel], bp_manual(index)),:);
    rep_xyz = repmat(xyz(i,:),size(relevant_xyz,1),1);
    xyz_diff = relevant_xyz-rep_xyz;
    this_distances = sqrt(sum(xyz_diff.^2,2));
    distances(end+1:end+size(xyz_diff,1),:) = this_distances;
    for j=find(abs(this_distances-10)<=6)'
        hold all
        plot3([relevant_xyz(j,1), repmat(xyz(i,1),length(j),1)],...
            [relevant_xyz(j,2), repmat(xyz(i,2),length(j),1)],...
            [relevant_xyz(j,3), repmat(xyz(i,3),length(j),1)], 'o-','color',[.8 .8 .8])
    end
    if any(abs(this_distances-10)>6)
        bad_tags{end+1} = talStruct(i).tagName;
         for j=find(abs(this_distances-10)>6)'
             hold all
             plot3([relevant_xyz(j,1), repmat(xyz(i,1),length(j),1)],...
                 [relevant_xyz(j,2), repmat(xyz(i,2),length(j),1)],...
                 [relevant_xyz(j,3), repmat(xyz(i,3),length(j),1)], 'ro-')
             text(xyz(i,1),xyz(i,2),xyz(i,3),talStruct(i).tagName);
         end
         drawnow;
    end
end 