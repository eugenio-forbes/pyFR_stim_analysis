%Written by JFB
% Updates by AGR: 
% (1) automatically make good_leads
% (2) re-write jacksheet so it matches up with vox tool
% USER INPUTS
subj       = 'UT004';
rawDirName = 'EEG2100';
badTags    = {'EKG' 'DC' 'BAD' 'EDF' 'E'}; % exclude from re-ref
%exclude these from jacksheet and re-ref
badLabels = {'100','101','102','103','104','105','106','107',...
    '108','109','110','111','112','113'}; 
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
fprintf('\n')
eegFile  = getEDF_eeg_files(rawDir);
for k=1:length(eegFile)
  fprintf(' Splitting %d of %d files: \n\t%s',k,length(eegFile),eegFile{k})  
  edf_split(eegFile{k},subj,outDir,badLabels)
end
fprintf('\n\n')

% copy the jacksheet to the docs directory
jackMotherFile = fullfile(docsDir,'jacksheet.txt');
if ~exist(jackMotherFile,'file')
  fprintf('Making jacksheet in docs...')
  jackDaughter = fullfile(outDir,'jacksheet.txt');
  if ~exist(jackDaughter,'file');
    error(sprintf('\n\n\tNo jacksheet.txt in %s...not right\n\n',outDir))
  end
  system(['cp ' jackDaughter ' ' jackMotherFile]);
  pause(.1)
  fprintf('done\n')
else
  fprintf('Jacksheet exists in docs\n')
end

% make the electrode.m file
electrodesFile = fullfile(docsDir,'electrodes.m');
if ~exist(electrodesFile,'file')
  fprintf('Making electrodes.m in docs...')
  fid_jack = fopen(jackMotherFile,'r');
  X        = textscan(fid_jack,'%d%s');
  fclose(fid_jack);
  elecNum             = X{1};
  elecNam             = X{2};
  elecNam_stripped    = regexprep(elecNam,'\d','');
  un_elecNam_stripped = unique(elecNam_stripped);
  minMaxEls           = [];
  for k=1:length(un_elecNam_stripped)
    thisTag    = un_elecNam_stripped{k};
    indThisTag = strcmp(elecNam_stripped,thisTag);
    elsThisTag = elecNum(indThisTag);
    minMaxEls  = cat(1,minMaxEls,[min(elsThisTag) max(elsThisTag)]);
  end
  
  % now sort them
  [~,sortInd] = sort(minMaxEls(:,1));
  tagToWriteInElectrodesDotM = un_elecNam_stripped(sortInd);
  numToWriteInElectrodesDotM = minMaxEls(sortInd,:);
  
  % now write the electrodes.m file
  fid_elecs = fopen(electrodesFile,'w');
  %fprintf(fid_elecs,'%% electrodes.m file\n');
  %fprintf(fid_elecs,'%% Made within edf_split_wrapper.m for %s\n',subj);
  %fprintf(fid_elecs,'%%  ''''  using information from raw directory %s\n',rawDirName);
  %[~,user]=system('whoami');user=regexprep(user,'\n','');
  %fprintf(fid_elecs,'%%  ''''  by %s on %s\n',user,datestr(now));
  %fprintf(fid_elecs,'%%\n%%\n\n');
  fprintf(fid_elecs,'r=[\n');
  for k=1:size(numToWriteInElectrodesDotM,1)
    if all(cellfun(@isempty,regexp(tagToWriteInElectrodesDotM{k},badTags)))
    %if isempty(strfind(tagToWriteInElectrodesDotM{k},'EKG'));
      fprintf('  including %s in reference\n',tagToWriteInElectrodesDotM{k})
      fprintf(fid_elecs,'%3.0d , %3.0d;  %%%s  (OPTIONAL: fill in full name here)\n',...
	      numToWriteInElectrodesDotM(k,:),tagToWriteInElectrodesDotM{k});
    else
      fprintf('  **NOT** including %s in reference\n',tagToWriteInElectrodesDotM{k})
    end
  end
  fprintf(fid_elecs,'];\n\n');
  fclose(fid_elecs);  
  pause(.1)
  fprintf('done\n')

% no need to do this because edf_split has been updated - agr 6.29.2014
%   %Here, we will re-write the jacksheet from electrodes.m in the docs directory so that 
%   %it has electrode names that match up with vox tool
%   cd(docsDir)
%   !mv jacksheet.txt jacksheet_orig.txt
%   electrodes2jacksheet(subj)

 else
   fprintf('electrodes.m exists in docs\n')
 end


% make the leads.txt file
leadsFile = fullfile(talDir,'leads.txt');
if ~exist(leadsFile,'file')
  fprintf('Making leads.txt in docs...')
  fid_jack  = fopen(fullfile(docsDir,'jacksheet.txt'),'r');
  fid_leads = fopen(leadsFile,'w');
  X         = textscan(fid_jack,'%d%s');
  fclose(fid_jack);
  elecNum   = X{1};
  for e=1:length(elecNum)
    fprintf(fid_leads,'%d\n',elecNum(e));
  end  
  fclose(fid_leads);
  fprintf('done.\n')
  
  %write out a good_leads (can be edited later)
  cd(talDir)
  !cp leads.txt good_leads.txt

else
  fprintf('leads.txt exists in tal\n')
end

% reref it
fprintf('\nreref me: \n')
for k=1:length(eegFile)  
  thisFileExt = edf_split(eegFile{k},subj,outDir);
  keyboard
  fprintf('fix electrodes .m file \n')
  
  tagNameFile = fullfile(docsDir,'electrodes.m');
  run(tagNameFile)
  reref({thisFileExt},r,rerefDir,talDir);
end

% mk the bipolar thing
%mkGridStructForBipolarFromJacksheet(subj,badTags)
