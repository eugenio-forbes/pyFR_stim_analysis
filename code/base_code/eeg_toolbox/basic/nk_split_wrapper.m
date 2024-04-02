% USER INPUTS
subj       = 'UT030';
rawDirName = '6_21_2016 CCEP';
badTags    = {'X' 'BAD' 'REF' 'ECG' 'DC'}; %all non-brain channels (to exclude from reference)

%---------------------------------------
%---------------------------------------
%---------------------------------------
dataDir  = '/Users/brad/ExperimentData/subjFiles';
subjDir  = fullfile(dataDir,subj);
outDir   = fullfile(subjDir,'eeg.noreref');
rerefDir = fullfile(subjDir,'eeg.reref');
rawDir   = fullfile(subjDir,'raw',rawDirName);
docsDir  = fullfile(subjDir,'docs');
talDir   = fullfile(subjDir,'tal');

tagNameFile = fullfile(docsDir,'tagnameorder_script.m');

run(tagNameFile)

nk_split(subj,rawDir,outDir,tagNameOrder);

keyboard

% reref it
[grids fileStems]=nk_process(subj,rawDir,tagNameOrder,badTags);
tagNameFile = fullfile(docsDir,'electrodes.m');
fprintf('\n\n')


fprintf('If sync channel is not EKG, you must now delete the sync channels from leads.txt and good_leads.txt.\n')

keyboard

reref(fileStems,grids,rerefDir,talDir);
