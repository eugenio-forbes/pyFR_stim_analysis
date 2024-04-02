function prep_egi_data2(subject,sessdir,varargin)
%PREP_EGI_DATA2   Process one session of a pyEPL experiment with EGI recordings.
%   PREP_EGI_DATA2(SUBJECT,SESSDIR,VARARGIN) post-processes data in SESSDIR.
%   SUBJECT is the string identifier of the subject the data were recorded from.
%
%   VARARGIN is a set of parameter-value pairs that can be used to change
%   options from their default values. NOTE: all paths can be either absolute
%   or relative to SESSDIR.
%
%   Parameters:
%     'steps_to_run'  Specifies which parts of post-processing to run. 
%                     Cell array which may contain: 'split', 'reref',
%                     and/or 'align'. Default: {}
%     'agethresh'     If steps_to_run is empty, each step will be
%                     run if relevant files have been modified in the
%                     last [agethresh] days. Default: 1
%     'eventfiles'    Cell array of paths to .mat files containing
%                     events structures. Default: {'events.mat'}
%     'badchanfiles'  Cell array with paths to .txt files with one channel 
%                     number per row, indicating which channels were "bad" in
%                     this session. Default: {'bad_chan.txt'}
%     'norerefdir'    Specifies where split, but not rereferenced,
%                     EEG files should be saved.  Default: 'eeg/eeg.noreref'
%     'rerefdir'      Specifies where rereferenced EEG files should
%                     be saved. Default: 'eeg/eeg.reref'
%     'captype'       String indicating the type of cap used. Determines
%                     which channels recorded EOG. 'HCGSN' (default) or
%                     'GSN200'
%     'channels'      Vector of channel numbers indicating which channels
%                     should be included in rereferencing. Default:
%                     1:129
%     'ziprawfile'    Logical value; if true, the raw file is
%                     bzipped after splitting.  Default: 1.
%     'unixclean'     UNIX command to be run after rereferencing to
%                     delete or move unneeded files. Default: deletes
%                     all channel files in eeg/eeg.noreref

% 11/19/2008: RWL: added support for multiple bad channel files and shell-style
%             comments in bad channel files
% 12/9/2008:  RWL: added support for saving split files elsewhere

% set up directory structure
eegdir = fullfile(sessdir,'eeg');
cd(sessdir);

% process varargin
def.eventfiles = {'events.mat'};
def.badchanfiles = {'bad_chan.txt'};
def.norerefdir = fullfile(eegdir,'eeg.noreref');
def.rerefdir = fullfile(eegdir,'eeg.reref');
def.norerefdir = fullfile(eegdir,'eeg.noreref');
def.captype = 'HCGSN';
def.channels = 1:129;
def.steps_to_run = {};
def.agethresh = .8;
def.ziprawfile = 1;
def.unixclean = ['rm ' fullfile(def.norerefdir, '*.[0-9][0-9][0-9]')];

[eid,emsg,eventfiles,badchanfiles,norerefdir,rerefdir,captype,channels,steps_to_run,agethresh,ziprawfile,unixclean] = getargs(fieldnames(def),struct2cell(def),varargin{:});

% if manual directions specfied, read them in
steps = {'split','reref','align'};
if ~iscell(steps_to_run)
  steps_to_run = {steps_to_run};
end

for i=1:length(steps)
  if ismember(steps{i},steps_to_run)
    val = 1;
    else
    val = 0;
  end
  run.(steps{i}) = val;
end

% prepare input files
if ~iscell(badchanfiles)
  badchanfiles = {badchanfiles};
end
for i=1:length(badchanfiles)
  badchanfiles{i} = fullfile(eegdir,badchanfiles{i});
end;
raw = dir(fullfile(eegdir,'*.raw'));
raw = [raw; dir(fullfile(eegdir,'*.raw.bz2'))];

if ~(run.split | run.reref | run.align)
  % if none manually specified, go to automatic;
  % which parts get run depends on which files have been modified
  if splitcheck(raw,agethresh); % modified .raw file
    % run everything
    run.split = 1;
    run.reref = 1;
    run.align = 1;
  elseif filecheck(badchanfiles,agethresh); % modified badchan files
    if length(dir(fullfile(norerefdir,'*.001')))==0
      % need to split to get noreref files
      run.split = 1;
    end
    run.reref = 1;
  end
  
  if filecheck(eventfiles,agethresh); % modified events
    run.align = 1;
    if length(dir(fullfile(rerefdir,'*.001')))==0
      % need to split, reref
      run.split = 1;
      run.reref = 1;
    end      
  end
end

% split EEG data into channels
if run.split
  if length(raw)==0
    error('No .raw file found in %s.', eegdir)
  end
  
  for i=1:length(raw)
    % unzip if necessary
    rawfile = fullfile(eegdir,raw(i).name);
    fixed = strrep(rawfile, ' ', '\ ');
    if strcmp(rawfile(end-3:end),'.bz2')
      unix(['bunzip2 -v ' fixed]);
      rawfile = rawfile(1:end-4);
      fixed = fixed(1:end-4);
      ziprawfile = 1; % we want to re-zip after splitting             
    end
    
    % split this .raw file
    egi_split(rawfile,subject,norerefdir);
    
    if ziprawfile
      % zip
      unix(['bzip2 -v ' fixed]);
    end
  end
end

% rereferencing
if run.reref
  try
    % read the bad channel files if they exist; concatenate all bad
    % channels into one set
    badchan = [];
    for i=1:length(badchanfiles)
      fname = badchanfiles{i};
      fid = fopen(fname);
      % shell style comments are allowed in bad channel files
      tmp = textscan(fid,'%d', 'CommentStyle', '#');

      % sanity check
      if size(tmp{1}',1)~=1
	error('Bad file: %s. There should be only one channel per row.', ...
	      fname)
      end

      badchan = [badchan tmp{1}'];
    end
    
  catch
    % continue without excluding bad channels from rereferencing
    fprintf('Warning: unable to read all bad channels; excluding none.\n');
    badchan = [];
  end

  % save out all bad channels (possibly none) to new file in reref directory
  if ~exist(rerefdir,'dir')
    mkdir(rerefdir);
  end
  fprintf(fopen(fullfile(rerefdir, 'bad_chan.txt'), 'w'), '%d\n', badchan);

  % setup's done...time to do the actual re-referencing
  goodchan = setdiff(channels,badchan);

  % get the basename of each .raw file that has been split
  chan = dir(fullfile(norerefdir,'*.001'));
  if length(chan)==0
    error('No channel files found in %s.', norerefdir)
  end
  for i=1:length(chan)
    basename = chan(i).name(1:end-4);
    noreref_files{i} = fullfile(norerefdir,basename);
  end

  chan_endpts = [min(channels) max(channels)];
  % rereference using "good" channels to calculate average
  reref(noreref_files,chan_endpts,rerefdir,{channels,goodchan});
  
  % if specified, run clean script
  if ~isempty(unixclean)
    unix(unixclean);
  end
end

% alignment and artifact detection
if run.align
  % get the basename of each .raw file that has been rereferenced
  chan = dir(fullfile(rerefdir,'*.001'));
  if length(chan)==0
    error('No channel files found in %s.', rerefdir)
  end
  for i=1:length(chan)
    basename = chan(i).name(1:end-4);
    reref_files{i} = fullfile(rerefdir,[basename '.001']);
    
    % get the EEG sync pulse file
    sync = dir(fullfile(norerefdir,[basename '.DIN*']));
    if length(sync)==0
      error('Pulse file not found in %s.', norerefdir)
    end
    eegsyncfile{i} = fullfile(norerefdir, sync(1).name);
  end

  % there should be only one behavioral sync pulse file
  behsyncfile = {fullfile(sessdir,'eeg.eeglog.up')};
  if ~exist(behsyncfile{1},'file')
    % if we haven't already, extract the UP pulses
    fixEEGLog(fullfile(sessdir,'eeg.eeglog'), behsyncfile{1});
  end
  
  % get the samplerate
  samplerate = GetRateAndFormat(norerefdir);
  
  % run the alignment
  runAlign(samplerate,behsyncfile,eegsyncfile,reref_files,eventfiles,'mstime',0,1);
  % add artifact info
  [eog, perif] = getcapinfo(captype);
  for i=1:length(eventfiles)
    addArtifacts(eventfiles{i}, eog, 100, 0);
  end
end


function update = filecheck(files,agethresh)
  %FILECHECK   See if any files have been recently modified.
  %   UPDATE is false if none of the files in cell array
  %   FILES have been modified in the past AGETHRESH days.
  %
  
  update = 0;
  if length(files)==0
    error('files must be a cell array containing paths')
  end
  
  for i=1:length(files)
    d = dir(files{i});
    if length(d)==0
      return
    end
    age = now - datenum(d.date);
    if age<agethresh
      update = 1;
      return
    end
  end
%endfunction

function update = splitcheck(raw,agethresh)
  %EEGUPDATE   Check if EEG data need to be split.
  update = 0;

  % look for a .raw file in eegdir
  if length(raw)==0
    error('No .raw file found.')
  end

  % see if any .raw files are recently modified
  for i=1:length(raw)
    age = now - datenum(raw(i).date);
    if age<agethresh
      update = 1;
      return
    end
  end
%endfunction

function fixEEGLog(infile,outfile)
  %FIXEEGLOG - Fix pyepl EEG Logfile leaving only UP pulses.
  %
  %
  % FUNCTION:
  %   fixEEGLog(infile,outfile)
  %
  % INPUT ARGS:
  %   infile = 'eeg.eeglog';
  %   outfile = 'eeg.eeglog.up';
  %

  % read in the logfile
  [mstime, maxoffset, type] = textread(infile,'%s%n%s%*[^\n]','emptyvalue',NaN);

  % write out new file
  fid = fopen(outfile,'w');
  for i = 1:length(type)
    if strcmp(type{i},'UP')
      % save it to file
      fprintf(fid,'%s\t%d\t%s\n',mstime{i},maxoffset(i),type{i});
    end
  end
  fclose(fid);
%endfunction
