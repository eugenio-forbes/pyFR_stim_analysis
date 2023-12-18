function newTimePlots_stim(type,band,stimReg,hipROI)
rootDir = '/project/TIBIR/Lega_lab/s427026/AGPCC';
norm = 'nostim';
var = 'normPow';
lat = 'ipscon';
times = [0,2];  

VarNames1 = {'sub','t','p'};
VarNames2 = {'sub','pow'};

switch band
    case 'slow_theta'
        freq_ind = 1;
        scale = 1;
    case 'fast_theta'
        freq_ind = 2;
        scale = 1;
end

if contains(type,'encoding')
    conditions = {'ES','NS'};
else
    conditions = {'RS','NS'};
end
load('hippoStruct.mat')
sesscodes = [hippoStruct(:).sesscode];
sesscodes = eval(sprintf('unique(sesscodes([hippoStruct(:).has%s]&[hippoStruct(:).has%s]))',conditions{1},conditions{2}));
getmean = ~contains(type,'median');

switch stimReg
    case 'AG'
        ROIs = {'AGL','AGR'};
    case 'PCC'
        ROIs = {'PCL','PCR'};
end

saveDir = fullfile(rootDir,'New_Time_Plots_trials2',type,lat,hipROI);
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end
plotName1 = fullfile(saveDir,sprintf('seppow_%s_%s.svg',stimReg,band));
plotName2 = fullfile(saveDir,sprintf('difpow_%s_%s.svg',stimReg,band));
plotName3 = fullfile(saveDir,sprintf('ttime_%s_%s.svg',stimReg,band));
plotName4 = fullfile(saveDir,sprintf('tcon_%s_%s.svg',stimReg,band));

dataDir = fullfile(rootDir,sprintf('Tables_WM_%s/%s/%s/%d',type,norm,var,freq_ind));
dataFile = fullfile(dataDir,'t.mat');
load(dataFile)
t = t(ismember(t.stimReg,ROIs),:);
t = t(ismember(t.sess,sesscodes),:);
t = t(ismember(t.time,times),:);
t = t(ismember(t.condition,conditions),:);

if strcmp(lat,'ipscon')
    Rt = t(contains(t.stimReg,'R'),:);
    Lt = t(contains(t.stimReg,'L'),:);
    clear t
    L_AR = Lt(ismember(Lt.hipReg,'AR'),:);
    L_AL = Lt(ismember(Lt.hipReg,'AL'),:);
    L_PR = Lt(ismember(Lt.hipReg,'PR'),:);
    L_PL = Lt(ismember(Lt.hipReg,'PL'),:);
    clear Lt
    L_AR.hipReg = strrep(L_AR.hipReg,'R','L');
    L_AL.hipReg = strrep(L_AL.hipReg,'L','R');
    L_PR.hipReg = strrep(L_PR.hipReg,'R','L');
    L_PL.hipReg = strrep(L_PL.hipReg,'L','R');
    Lt = [L_AR;L_AL;L_PR;L_PL];
    clear L_AR L_AL L_PR L_PL
    t = [Rt;Lt];
    clear Rt Lt
end
t = t(strcmp(t.hipReg,hipROI),:);
subcodes = unique(t.sub);
con1_time0 = [];
con1_time2 = [];
con2_time0 = [];
con2_time2 = [];
t_con1 = [];
t_con2 = [];
t_time0 = [];
t_time2 = [];
for idx = 1:length(subcodes)
    thiscode = subcodes(idx);
    sub = t(t.sub == thiscode,:);
    con1_time0_temp = sub([strcmp(sub.condition,conditions{1})&sub.time == 0],:).pow;
    con1_time2_temp = sub([strcmp(sub.condition,conditions{1})&sub.time == 2],:).pow;
    con2_time0_temp = sub([strcmp(sub.condition,conditions{2})&sub.time == 0],:).pow;
    con2_time2_temp = sub([strcmp(sub.condition,conditions{2})&sub.time == 2],:).pow;
    if getmean
        [~,p1,~,t1] = ttest([con1_time2_temp],[con1_time0_temp]);
        t1 = t1.tstat;
        [~,p2,~,t2] = ttest([con2_time2_temp],[con2_time0_temp]);
        t2 = t2.tstat;
        [~,p3,~,t3] = ttest2([con1_time0_temp],[con2_time0_temp]);
        t3 = t3.tstat;
        [~,p4,~,t4] = ttest2([con1_time2_temp],[con2_time2_temp]);
        t4 = t4.tstat;
        t_con1_temp = table(thiscode,t1,p1,'VariableNames',VarNames1);
        t_con1 = [t_con1;t_con1_temp]; clear t_con1_temp;
        t_con2_temp = table(thiscode,t2,p2,'VariableNames',VarNames1);
        t_con2 = [t_con2;t_con2_temp]; clear t_con2_temp;
        t_time0_temp = table(thiscode,t3,p3,'VariableNames',VarNames1);
        t_time0 = [t_time0;t_time0_temp]; clear t_time0_temp;
        t_time2_temp = table(thiscode,t4,p4,'VariableNames',VarNames1);
        t_time2 = [t_time2;t_time2_temp]; clear t_time2_temp;
        clear t1 t2 t3 t4 p1 p2 p3 p4
    end
    con1_time0_temp = table(thiscode,median(con1_time0_temp),'VariableNames',VarNames2);
    con1_time0 = [con1_time0;con1_time0_temp]; clear con1_time0_temp;
    con1_time2_temp = table(thiscode,median(con1_time2_temp),'VariableNames',VarNames2);
    con1_time2 = [con1_time2;con1_time2_temp]; clear con1_time2_temp;
    con2_time0_temp = table(thiscode,median(con2_time0_temp),'VariableNames',VarNames2);
    con2_time0 = [con2_time0;con2_time0_temp]; clear con2_time0_temp;
    con2_time2_temp = table(thiscode,median(con2_time2_temp),'VariableNames',VarNames2);
    con2_time2 = [con2_time2;con2_time2_temp]; clear con2_time2_temp;
    clear sub thiscode
end
clear t

sigDir = fullfile(rootDir,'LMEM_time_new',type,norm,var);
load(fullfile(sigDir,sprintf('Time_%s_0-2_%sB_slope/data.mat',conditions{1},stimReg)));
p_con1 = eval(sprintf('%s_c',hipROI)); clear AL_c AR_c PL_c PR_c
p_con1 = p_con1(freq_ind,2);
load(fullfile(sigDir,sprintf('Time_%s_0-2_%sB_slope/data.mat',conditions{2},stimReg)));
p_con2 = eval(sprintf('%s_c',hipROI)); clear AL_c AR_c PL_c PR_c
p_con2 = p_con2(freq_ind,2);
load(fullfile(sigDir,sprintf('Con_0_%sNS_%sB_slope/data.mat',conditions{1},stimReg)));
p_time0 = eval(sprintf('%s_c',hipROI)); clear AL_c AR_c PL_c PR_c
p_time0 = p_time0(freq_ind,2);
load(fullfile(sigDir,sprintf('Con_2_%sNS_%sB_slope/data.mat',conditions{1},stimReg)));
p_time2 = eval(sprintf('%s_c',hipROI)); clear AL_c AR_c PL_c PR_c
p_time2 = p_time2(freq_ind,2);
load(fullfile(sigDir,sprintf('TimexCon_%sNS_0-2_%sB_slope/data.mat',conditions{1},stimReg)));
p_interaction = eval(sprintf('%s_c',hipROI)); clear AL_c AR_c PL_c PR_c
p_interaction = p_interaction(freq_ind,6);

%%% separate pows
figure;
hold on
set(gcf,'Position',[0 0 80 100])
delta_box = [con1_time0.pow,con1_time2.pow,NaN(height(con1_time0),1),con2_time0.pow,con2_time2.pow];
for idx = 1:height(con1_time0)
    if getmean
        if t_con1(idx,:).p < 0.05
            plot([1 2],[delta_box(idx,1) delta_box(idx,2)],'-k','LineWidth',0.75);
        else
            plot([1 2],[delta_box(idx,1) delta_box(idx,2)],'-','Color',[0.5 0.5 0.5],'LineWidth',0.5);
        end
        if t_con2(idx,:).p < 0.05
            plot([4 5],[delta_box(idx,4) delta_box(idx,5)],'-k','LineWidth',0.75);
        else
            plot([4 5],[delta_box(idx,4) delta_box(idx,5)],'-','Color',[0.5 0.5 0.5],'LineWidth',0.5);
        end
    else
         plot([1 2],[delta_box(idx,1) delta_box(idx,2)],'-','Color',[0.5 0.5 0.5],'LineWidth',0.5);
         plot([4 5],[delta_box(idx,4) delta_box(idx,5)],'-','Color',[0.5 0.5 0.5],'LineWidth',0.5);
    end
end

ppos = 2.7 * scale;
if p_con1<0.01
    plot(1.4,ppos,'*k','MarkerSize',1.5)
    plot(1.6,ppos,'*k','MarkerSize',1.5)
elseif p_con1<0.05
    plot(1.5,ppos,'*k','MarkerSize',1.5)
elseif p_con1<0.1
    plot(1.5,ppos,'+k','MarkerSize',1.5)
end
if p_con2<0.01
    plot(4.4,ppos,'*k','MarkerSize',1.5)
    plot(4.6,ppos,'*k','MarkerSize',1.5)
elseif p_con2<0.05
    plot(4.5,ppos,'*k','MarkerSize',1.5)
elseif p_con2<0.1
    plot(4.5,ppos,'+k','MarkerSize',1.5)
end
if p_interaction<0.01
    plot(2.9,ppos,'*k','MarkerSize',1.5)
    plot(3.1,ppos,'*k','MarkerSize',1.5)
elseif p_interaction<0.05
    plot(3,ppos,'*k','MarkerSize',1.5)
elseif p_interaction<0.1
    plot(3,ppos,'+k','MarkerSize',1.5)
end
for plot_idx = [1,2,4,5]
    if plot_idx <3
        thiscolor = [0 0.6 1];
        leg_idx = 1;
    else
        thiscolor = [0.5 0.5 0.5];
        leg_idx = 2;
    end
    deltas = delta_box(:,plot_idx);
    for idx = 1:length(deltas)
        if idx == 1
            a(leg_idx) = plot(plot_idx,deltas(idx),'.','Color',thiscolor,'MarkerSize',3);
        else
            plot(plot_idx,deltas(idx),'.','Color',thiscolor,'MarkerSize',3);
        end
    end
end
colors = [[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0]];
plot([-100:100],zeros(201),'-k')
h = boxplot(delta_box,[1:5],'Colors',colors,'Symbol','','Widths',0.75);
set(h,{'linew'},{0.5})
ylim([-1 3]*scale); yticks([-1,0,3]*scale);yticklabels([])  %yticklabels({'-100','0','300'});
xlim([0 6]); xticks([1,2,4,5]);xticklabels([])  %xticklabels({'E','L','E','L'});
set(gcf,'Position',[0 0 80 100])
% legend(a,'Stim','No Stim','Orientation','vertical','Location','best')
hold off
saveas(gcf,plotName1)
close all

%%% pow differences
figure;
hold on
set(gcf,'Position',[0 0 80 100])
delta_box = [[con1_time0.pow-con1_time2.pow],[con2_time0.pow-con2_time2.pow]];
for idx = 1:height(con1_time0)
         plot([1 2],[delta_box(idx,1) delta_box(idx,2)],'-','Color',[0.5 0.5 0.5],'LineWidth',0.5);
end
for plot_idx = [1,2]
    if plot_idx < 1.5
        thiscolor = [0 0.6 1];
        leg_idx = 1;
    else
        thiscolor = [1 0 0.6];
        leg_idx = 2;
    end
    deltas = delta_box(:,plot_idx)
    
    for idx = 1:length(deltas)
        if idx == 1
            a(leg_idx) = plot(plot_idx,deltas(idx),'.','Color',thiscolor,'MarkerSize',3);
        else
            plot(plot_idx,deltas(idx),'.','Color',thiscolor,'MarkerSize',3);
        end
    end
end
ppos = 1.8*scale
if p_con1<0.01
    plot(0.9,ppos,'*k','MarkerSize',1.5)
    plot(1.1,ppos,'*k','MarkerSize',1.5)
elseif p_con1<0.05
    plot(1,ppos,'*k','MarkerSize',1.5)
elseif p_con1<0.1
    plot(1,ppos,'+k','MarkerSize',1.5)
end
if p_con2<0.01
    plot(1.9,ppos,'*k','MarkerSize',1.5)
    plot(2.1,ppos,'*k','MarkerSize',1.5)
elseif p_con2<0.05
    plot(2,ppos,'*k','MarkerSize',1.5)
elseif p_con2<0.1
    plot(2,ppos,'+k','MarkerSize',1.5)
end
if p_interaction<0.01
    plot(1.4,ppos,'*k','MarkerSize',1.5)
    plot(1.6,ppos,'*k','MarkerSize',1.5)
elseif p_interaction<0.05
    plot(1.5,ppos,'*k','MarkerSize',1.5)
elseif p_interaction<0.1
    plot(1.5,ppos,'+k','MarkerSize',1.5)
end
colors = [[0 0 0];[0 0 0]];
plot([-100:100],zeros(201),'-k')
h = boxplot(delta_box,[1:2],'Colors',colors,'Symbol','','Widths',0.75);
set(h,{'linew'},{0.5})
set(gcf,'Position',[0 0 80 100])

ylim([-2 2]); yticks([-2:1:2]);yticklabels([])  %yticklabels({'-200','','0','','200'});
xlim([0 3]); xticks([1,2]);xticklabels([])  %xticklabels({'Stim','No Stim'});
hold off
saveas(gcf,plotName2)
close all

%%% tstats per time
if getmean
    figure;
    hold on
    set(gcf,'Position',[0 0 80 100])
    delta_box = [t_con1.t,t_con2.t];
    for idx = 1:height(con1_time0)
        plot([1 2],[delta_box(idx,1) delta_box(idx,2)],'-','Color',[0.5 0.5 0.5],'LineWidth',0.5);
    end
    for plot_idx = [1,2]
        if plot_idx < 1.5
            thiscolor = [0 0.6 1];
            leg_idx = 1;
        else
            thiscolor = [0.5 0.5 0.5];
            leg_idx = 2;
        end
        deltas = delta_box(:,plot_idx);
        for idx = 1:length(deltas)
            if idx == 1
                a(leg_idx) = plot(plot_idx,deltas(idx),'.','Color',thiscolor,'MarkerSize',3);
            else
                plot(plot_idx,deltas(idx),'.','Color',thiscolor,'MarkerSize',3);
            end
        end
    end
    ppos = 4.5 * scale;
    if p_con1<0.01
        plot(0.9,ppos,'*k','MarkerSize',1.5)
        plot(1.1,ppos,'*k','MarkerSize',1.5)
    elseif p_con1<0.05
        plot(1,ppos,'*k','MarkerSize',1.5)
    elseif p_con1<0.1
        plot(1,ppos,'+k','MarkerSize',1.5)
    end
    if p_con2<0.01
        plot(1.9,ppos,'*k','MarkerSize',1.5)
        plot(2.1,ppos,'*k','MarkerSize',1.5)
    elseif p_con2<0.05
        plot(2,ppos,'*k','MarkerSize',1.5)
    elseif p_con2<0.1
        plot(2,ppos,'+k','MarkerSize',1.5)
    end
    if p_interaction<0.01
        plot(1.4,ppos,'*k','MarkerSize',1.5)
        plot(1.6,ppos,'*k','MarkerSize',1.5)
    elseif p_interaction<0.05
        plot(1.5,ppos,'*k','MarkerSize',1.5)
    elseif p_interaction<0.1
        plot(1.5,ppos,'+k','MarkerSize',1.5)
    end
    colors = [[0 0 0];[0 0 0]];
    plot([-100:100],zeros(201),'-k')
    h = boxplot(delta_box,[1:2],'Colors',colors,'Symbol','','Widths',0.75);
    set(h,{'linew'},{0.5})
    ylim([-5 5]); yticks([-5:5:5]);yticklabels([])  %yticklabels({'-5','0','5'});
    xlim([0 3]); xticks([1,2]);xticklabels([])  %xticklabels({'Stim','No Stim'});
%     legend(a,'Stim','No Stim','Orientation','vertical','Location','best')
    set(gcf,'Position',[0 0 80 100])
    hold off
    saveas(gcf,plotName3)
    close all
end

%%% tstats per con
if getmean
    figure;
    hold on
    set(gcf,'Position',[0 0 80 100])
    delta_box = [t_time0.t,t_time2.t];
    for idx = 1:height(con1_time0)
        plot([1 2],[delta_box(idx,1) delta_box(idx,2)],'-','Color',[0.5 0.5 0.5],'LineWidth',0.5);
    end
    for plot_idx = [1,2]
        if plot_idx < 1.5
            thiscolor = [0 0.6 1];
            leg_idx = 1;
        else
            thiscolor = [1 0 0.6];
            leg_idx = 2;
        end
        deltas = delta_box(:,plot_idx);
        for idx = 1:length(deltas)
            if idx == 1
                a(leg_idx) = plot(plot_idx,deltas(idx),'.','Color',thiscolor,'MarkerSize',3);
            else
                plot(plot_idx,deltas(idx),'.','Color',thiscolor,'MarkerSize',3);
            end
        end
    end
    ppos = 4.5*scale;
    if p_time0<0.01
        plot(0.9,ppos,'*k','MarkerSize',1.5)
        plot(1.1,ppos,'*k','MarkerSize',1.5)
    elseif p_time0<0.05
        plot(1,ppos,'*k','MarkerSize',1.5)
    elseif p_time0<0.1
        plot(1,ppos,'+k','MarkerSize',1.5)
    end
    if p_time2<0.01
        plot(1.9,ppos,'*k','MarkerSize',1.5)
        plot(2.1,ppos,'*k','MarkerSize',1.5)
    elseif p_time2<0.05
        plot(2,ppos,'*k','MarkerSize',1.5)
    elseif p_time2<0.1
        plot(2,ppos,'+k','MarkerSize',1.5)
    end
    if p_interaction<0.01
        plot(1.4,ppos,'*k','MarkerSize',1.5)
        plot(1.6,ppos,'*k','MarkerSize',1.5)
    elseif p_interaction<0.05
        plot(1.5,ppos,'*k','MarkerSize',1.5)
    elseif p_interaction<0.1
        plot(1.5,ppos,'+k','MarkerSize',1.5)
    end
    colors = [[0 0 0];[0 0 0]];
    plot([-100:100],zeros(201),'-k')
    h = boxplot(delta_box,[1:2],'Colors',colors,'Symbol','','Widths',0.75);
    set(h,{'linew'},{0.5})
    set(gcf,'Position',[0 0 80 100])
    ylim([-5 5]); yticks([-5:5:5]);yticklabels([]) %yticklabels({'-5','0','5'});
    xlim([0 3]); xticks([1,2]);xticklabels([]) %xticklabels({'Early','Late'});
    hold off
    saveas(gcf,plotName4)
    close all
end

end


