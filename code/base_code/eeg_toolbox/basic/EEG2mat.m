
% filename ='/Users/JL/ExperimentalData/subjFiles/UT053/eeg.reref/UT053_pyFR_stim_encode_0_15Jun17_1424.006';
% filename = '/Users/JL/ExperimentalData/subjFiles/UT017/eeg.reref/UT017_29Oct15_1026.038';
filename = '/Volumes/project/TIBIR/Lega_lab/shared/lega_ansir/subjFiles/UT052/eeg.reref/UT052_FR1_0_09Jun17_2218.001';
[dat,samplerate] = loadData(filename);

fileLength = num2str((length(dat)/samplerate)/60);

figure(2)
plot(1:length(dat),dat)
set(gcf,'position',[200 600 2000 750]);

function [dat,samplerate] = loadData(filename)
%LOADDATA - Load EEG data from a file
%

% get the data format
[samplerate,nBytes,dataformat,gain] = GetRateAndFormat(filename);

% Open and load the file
fid = fopen(filename, 'r','l');
dat =  fread(fid, inf, dataformat);
fclose(fid);

% apply the gain
dat = dat.*gain;
end

% figure(1)
% set(gcf,'position',[200 600 2000 750]);
% plot(1:length(all),all)
% hold on
% plot(1:length(HL),HL,'r')