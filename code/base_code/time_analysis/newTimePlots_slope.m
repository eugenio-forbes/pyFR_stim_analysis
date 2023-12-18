function newTimePlots_slope(type,band,stimReg,hipROI)
rootDir = '/project/TIBIR/Lega_lab/s427026/AGPCC';
norm = 'nostim';
var = 'normPow';
lat = 'ipscon';

VarNames1 = {'sub','t','p'};
VarNames2 = {'sub','t'};

switch band
    case 'slow_theta'
        freq_ind = 1;
    case 'fast_theta'
        freq_ind = 2;
end

if contains(type,'encoding')
    conditions = {'ES','NS'};
else
    conditions = {'RS','NS'};
end

load('hippoStruct.mat')
sesscodes = [hippoStruct(:).sesscode];
sesscodes = eval(sprintf('unique(sesscodes([hippoStruct(:).has%s]&[hippoStruct(:).has%s]))',conditions{1},conditions{2}));

switch stimReg
    case 'AG'
        ROIs = {'AGL','AGR'};
    case 'PCC'
        ROIs = {'PCL','PCR'};
end

saveDir = fullfile(rootDir,'New_Time_Plots_slope',type,lat,hipROI);
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end
plotName1 = fullfile(saveDir,sprintf('conts_%s_%s.svg',stimReg,band));
plotName2 = fullfile(saveDir,sprintf('tslopes_%s_%s.svg',stimReg,band));

dataDir = fullfile(rootDir,sprintf('Tables_WM_%s_slope/%s/%s/%d',type,norm,var,freq_ind));
dataFile = fullfile(dataDir,'t.mat');
load(dataFile)
t = t(ismember(t.stimReg,ROIs),:);
t = t(ismember(t.sess,sesscodes),:);
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
con1_slopes = [];
con2_slopes = [];
t_slopes = [];
for idx = 1:length(subcodes)
    thiscode = subcodes(idx);
    sub = t(t.sub == thiscode,:);
    con1_slopes_temp = sub([strcmp(sub.condition,conditions{1})],:).slope;
    con2_slopes_temp = sub([strcmp(sub.condition,conditions{2})],:).slope;
    [~,p1,~,t1] = ttest2([con1_slopes_temp],[con2_slopes_temp]);
    t1 = t1.tstat;
    t_slopes_temp = table(thiscode,t1,p1,'VariableNames',VarNames1);
    t_slopes = [t_slopes;t_slopes_temp]; clear t_slopes_temp;
    clear t1 p1
    [~,~,~,m1] = ttest(con1_slopes_temp);
    m1 = m1.tstat;
    [~,~,~,m2] = ttest(con2_slopes_temp);
    m2 = m2.tstat;
    con1_slopes_temp = table(thiscode,m1,'VariableNames',VarNames2);
    con1_slopes = [con1_slopes;con1_slopes_temp]; clear con1_slopes_temp;
    con2_slopes_temp = table(thiscode,m2,'VariableNames',VarNames2);
    con2_slopes = [con2_slopes;con2_slopes_temp]; clear con2_slopes_temp;
    clear sub thiscode
end
clear t

sigDir = fullfile(rootDir,'LMEM_time_new_slope',type,norm,var);
load(fullfile(sigDir,sprintf('Con_%sNS_%sB_slope/data.mat',conditions{1},stimReg)));
p_con = eval(sprintf('%s_c',hipROI)); clear AL_c AR_c PL_c PR_c
p_con = p_con(freq_ind,2);

figure;
hold on
set(gcf,'Position',[0 0 80 100])
delta_box = [con1_slopes.t,con2_slopes.t];
for idx = 1:height(t_slopes)
    if t_slopes(idx,:).p < 0.05
        plot([1 2],[delta_box(idx,1) delta_box(idx,2)],'-k','LineWidth',0.5);
    else
        plot([1 2],[delta_box(idx,1) delta_box(idx,2)],'-','Color',[0.5 0.5 0.5],'LineWidth',0.5);
    end
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
ppos = 4.5;
if p_con<0.01
    plot(1.4,ppos,'*k','MarkerSize',3)
    plot(1.6,ppos,'*k','MarkerSize',3)
elseif p_con<0.05
    plot(1.5,ppos,'*k','MarkerSize',3)
elseif p_con<0.1
    plot(1.5,ppos,'+k','MarkerSize',3)
end
colors = [[0 0 0];[0 0 0]];
plot([-100:100],zeros(201),'-k')
h = boxplot(delta_box,[1:2],'Colors',colors,'Symbol','','Widths',0.75);
set(h,{'linew'},{0.5})
ylim([-5 5]); yticks([-5:5:5]); yticklabels([]) %yticklabels({'-5','0','5'});
xlim([0 3]); xticks([1,2]); xticklabels([]) %xticklabels({'Stim','No Stim'});
% legend(a,'Stim','No Stim','Orientation','vertical','Location','best')
hold off
saveas(gcf,plotName1)
close all

figure;
hold on
set(gcf,'Position',[0 0 80 100])
delta_box = [t_slopes.t];
thiscolor = [0 0.6 1];
leg_idx = 1;
deltas = delta_box(:,1);
for idx = 1:length(deltas)
    if idx == 1
        a(leg_idx) = plot(1,deltas(idx),'.','Color',thiscolor,'MarkerSize',5);
    else
        plot(1,deltas(idx),'.','Color',thiscolor,'MarkerSize',5);
    end
end
ppos = 4.5;
if p_con<0.01
    plot(0.9,ppos,'*k','MarkerSize',3)
    plot(1.1,ppos,'*k','MarkerSize',3)
elseif p_con<0.05
    plot(1,ppos,'*k','MarkerSize',3)
elseif p_con<0.1
    plot(1,ppos,'+k','MarkerSize',3)
end
colors = [0 0 0];
plot([-100:100],zeros(201),'-k')
h = boxplot(delta_box,'Symbol','','Colors',colors,'Widths',0.5);
set(h,{'linew'},{0.5})
ylim([-5 5]); yticks([-5:5:5]);yticklabels([]) %yticklabels({'-5','0','5'});
xlim([0.5 1.5]); xticks([1]); xticklabels([]) %xticklabels({'Stim - No Stim'});
hold off
saveas(gcf,plotName2)
close all

end
