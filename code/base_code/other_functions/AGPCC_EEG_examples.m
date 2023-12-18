clearvars; clc; warning off;
load('taskStruct.mat')
load('hippoStruct.mat')

saveDir = '/Users/forbes/Desktop/NewPlots';
subject = 'UT011';
type = 'WM_words';
chan = 31;

nostim_Color = [0 0 0];
base_Color = [173, 216, 230]/256;
no_Color = [0.9, 0.9, 0.9];
xlab = 'time (ms)';
ylab = 'Voltage (\muV)';
yy = [-150 150];


rootDir = '/Volumes/project/TIBIR/Lega_lab/s427026/AGPCC';
buffer = 500; 
lineNoise = [58 62; 98 102; 118 122; 178 182; 198 202; 238 242; 298 302; 358 362; 398 402];
kurtSkew = 1; 
kThresh = 4;

samplingFreq = 1000; nyquistFreq1 = floor(samplingFreq/2);
frequencies = (2.^((8:56)/8));
vFreqInd = find(frequencies<500);
vFreqs = frequencies(vFreqInd);

thisSub = taskStruct(strcmp({taskStruct.subject},subject));
task = thisSub(1).task;
sess = thisSub(1).session;
eventPath = sprintf(fullfile(rootDir,'/Events/%s/%s/%s/events.mat'),subject,task,sess);
load(eventPath);
eegfile = {events.eegfile};
eegfile = strcat('/Volumes',eegfile);
for idx = 1:length(events)
    events(idx).eegfile = eegfile{idx};
end
    

thisElec = hippoStruct(strcmp({hippoStruct.sub},subject));
thisElec = thisElec(ismember([thisElec.hipNum],chan));
thisElec = thisElec(1);
ref = [thisElec.ref];
chan = [chan,ref];

switch type
    case 'WM_words'
        period = 'Encoding';
        events = events(strcmp({events.type},'WORD'));
        posi = [0 0 660 282];
        xx = [0 2100];
        linePlot = true;
        linePos = 500;
        xx2 = [0:250:2100];
        xx3 = {'-500','-250','0','250','500','750','1000','1250','1500'};
        conditions = {'ES','NS'};
        offset = -500;
        stimCode = 1;
        dur = 2100;
        stim_Color = [214,28,78]/256;
    case 'WM_rec'
        period = 'Retrieval';
        events = events(strcmp({events.type},'REC_WORD'));
        events = events(~strcmp({events.item},'<>'));
        posi = [0 0 450 282];
        xx = [0 1500];
        linePlot = false;
        linePos = [];
        xx2 = [0:250:1500];
        xx3 = {'-1500','-1250','-1000','-750','-500','-250','0'};
        conditions = {'RS','NS'};
        offset = -1500;
        stimCode = 3;
        dur = 1500;
        stim_Color = [39,123,192]/256;
end

plot1 = sprintf('rawEEG_stim_%s_%s_%d.svg',period,subject,chan(1)); 
plot2 = sprintf('rawEEG_nostim_%s_%s_%d.svg',period,subject,chan(1));
plot3 = sprintf('lowEEG_stim_%s_%s_%d.svg',period,subject,chan(1));
plot4 = sprintf('lowEEG_nostim_%s.svg',period,subject,chan(1));
plot5 = sprintf('rawEEG_both_%s_%s_%d.svg',period,subject,chan(1));
plot6 = sprintf('lowEEG_both_%s_%s_%d.svg',period,subject,chan(1));

[eeg,kInd] = gete_MP_EF(chan,events,dur,offset,buffer,'freqs',vFreqs,'filtfreq',lineNoise,'filtorder',1,'filttype','stop',...
                'kthresh',kThresh,'keepk');
stimEEG = eeg([events.stimcode] == stimCode & ~ismember([1:length(events)],kInd),:);
nostimEEG = eeg([events.stimcode] == 5 & ~ismember([1:length(events)],kInd),:);

lowpass_EEG = gete_MP_EF(chan,events,dur,offset,buffer,'freqs',vFreqs,'filtfreq',12,'filtorder',1,'filttype','low');
stim_lowpass = lowpass_EEG([events.stimcode] == stimCode & ~ismember([1:length(events)],kInd),:);
nostim_lowpass = lowpass_EEG([events.stimcode] == 5 & ~ismember([1:length(events)],kInd),:);

stim_i = 1;
a = 0;
while stim_i < (size(stimEEG,1) + 1) & a == 0
    f1 = figure(1);
    fill([0,300,300,0],[-10000, -10000, 10000, 10000],base_Color)
    hold on
    fill([300,500,500,300],[-10000, -10000, 10000, 10000],no_Color)
    thisEEG = stimEEG(stim_i,:);
    plot(thisEEG,'-','Color',stim_Color,'lineWidth',2)
    plot([-5000 5000],[0 0],'-k')
    set(gcf,'Position',posi)
    xlabel(xlab);xlim(xx); xticks(xx2); xticklabels(xx3);
    ylabel(ylab);ylim(yy);
    if linePlot
        plot([linePos,linePos],[-10000,10000],'-k','lineWidth',2);
    end
    hold off
    
    f2 = figure(2);
    fill([0,300,300,0],[-10000, -10000, 10000, 10000],base_Color)
    hold on
    fill([300,500,500,300],[-10000, -10000, 10000, 10000],no_Color)
    thisEEG = stim_lowpass(stim_i,:);
    plot(thisEEG,'-','Color',stim_Color,'lineWidth',2)
    plot([-5000 5000],[0 0],'-k')
    set(gcf,'Position',posi)
    xlabel(xlab);xlim(xx); xticks(xx2); xticklabels(xx3);
    ylabel(ylab);ylim(yy);
    if linePlot
        plot([linePos,linePos],[-10000,10000],'-k','lineWidth',2);
    end
    hold off
    
    user_input = input('Did you like it? y/n:','s');
    if strcmp(user_input,'y')
        a = 1;
        saveas(figure(1),fullfile(saveDir,plot1));
        saveas(figure(2),fullfile(saveDir,plot3));
    else
        stim_i = stim_i + 1;
    end
    close all
end

nostim_i = 1;
a = 0;
while nostim_i < (size(nostimEEG,1) + 1) & a == 0
    f1 = figure(1);
    fill([0,300,300,0],[-10000, -10000, 10000, 10000],base_Color)
    hold on
    fill([300,500,500,300],[-10000, -10000, 10000, 10000],no_Color)
    thisEEG = nostimEEG(nostim_i,:);
    plot(thisEEG,'-','Color',nostim_Color,'lineWidth',2)
    plot([-5000 5000],[0 0],'-k')
    set(gcf,'Position',posi)
    xlabel(xlab);xlim(xx); xticks(xx2); xticklabels(xx3);
    ylabel(ylab);ylim(yy);
    if linePlot
        plot([linePos,linePos],[-10000,10000],'-k','lineWidth',2);
    end
    hold off
    
    f2 = figure(2);
    fill([0,300,300,0],[-10000, -10000, 10000, 10000],base_Color)
    hold on
    fill([300,500,500,300],[-10000, -10000, 10000, 10000],no_Color)
    thisEEG = nostim_lowpass(nostim_i,:);
    plot(thisEEG,'-','Color',nostim_Color,'lineWidth',2)
    plot([-5000 5000],[0 0],'-k')
    set(gcf,'Position',posi)
    xlabel(xlab);xlim(xx); xticks(xx2); xticklabels(xx3);
    ylabel(ylab);ylim(yy);
    if linePlot
        plot([linePos,linePos],[-10000,10000],'-k','lineWidth',2);
    end
    hold off
    
    user_input = input('Did you like it? y/n:','s');
    if strcmp(user_input,'y')
        a = 1;
        saveas(figure(1),fullfile(saveDir,plot2));
        saveas(figure(2),fullfile(saveDir,plot4));
    else
        nostim_i = nostim_i + 1;
    end
    close all
end

stim_low = stim_lowpass(stim_i,:);
stim_raw = stimEEG(stim_i,:);
nostim_low = nostim_lowpass(nostim_i,:);
nostim_raw = nostimEEG(nostim_i,:);

figure;
fill([0,300,300,0],[-10000, -10000, 10000, 10000],base_Color)
hold on
fill([300,500,500,300],[-10000, -10000, 10000, 10000],no_Color)
plot(stim_raw,'-','Color',stim_Color,'lineWidth',2)
plot(nostim_raw,'-','Color',nostim_Color,'lineWidth',2)
plot([-5000 5000],[0 0],'-k')
set(gcf,'Position',posi)
xlabel(xlab);xlim(xx); xticks(xx2); xticklabels(xx3);
ylabel(ylab);ylim(yy);
if linePlot
    plot([linePos,linePos],[-10000,10000],'-k','lineWidth',2);
end
hold off
saveas(gcf,fullfile(saveDir,plot5))
close all

figure;
fill([0,300,300,0],[-10000, -10000, 10000, 10000],base_Color)
hold on
fill([300,500,500,300],[-10000, -10000, 10000, 10000],no_Color)
plot(stim_low,'-','Color',stim_Color,'lineWidth',2)
plot(nostim_low,'-','Color',nostim_Color,'lineWidth',2)
plot([-5000 5000],[0 0],'-k')
set(gcf,'Position',posi)
xlabel(xlab);xlim(xx); xticks(xx2); xticklabels(xx3);
ylabel(ylab);ylim(yy);
if linePlot
    plot([linePos,linePos],[-10000,10000],'-k','lineWidth',2);
end
hold off
saveas(gcf,fullfile(saveDir,plot6))
close all

