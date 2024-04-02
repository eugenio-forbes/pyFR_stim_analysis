function split_ncs(data_dir,out_dir)
%process_ncs - process ncs files and place in the specified directory.

%2-22-10 JRM,JJ     wrote it.
%2-22-10 JRM        added sync pulse file (sync.txt)
%2-23-10 JRM        change sync pulse file format to plain text
%2-23-10 JRM        timestamps are each associated with 512 samples

SAMPLES_PER_BLOCK = 512; %how many samples per timestamp?

if ~exist(out_dir,'dir')
    mkdir(out_dir)
end

%create channels files
csc_files = dir(fullfile(data_dir,'*.ncs'));
for fname = {csc_files.name}
    chan = regexp_grab_token(fname{1},'(\d+)');
    outfile = fullfile(out_dir,['chan.',sprintf('%03d',chan)]);
    
    if exist(outfile,'file')
        fprintf('Skipping pre-existing file: %s\n',outfile);
        continue;
    end
    
    try
        [data,info,timestamps] = load_ncs(fullfile(data_dir,fname{1}));
    catch
        continue;
    end
    if ~exist('gain','var')
        gain = info.bitVolts;
    elseif gain ~= info.bitVolts
        error('voltage range does not agree for all channels');
    end
       
    fd = fopen(outfile,'w','l');
    nsamples = fwrite(fd,data,'short');
    fclose(fd);
    assert(nsamples == length(data),'some sync pulses were not written to sync.txt!');
end

%create params file
fname = fullfile(out_dir,'params.txt');
fd = fopen(fname,'w');
fprintf(fd,'samplerate %0.50f\n',info.actualSampleRate);
fprintf(fd,'dataformat ''short''\n');
fprintf(fd,'gain %0.20f\n',info.bitVolts);
fclose(fd);

%create sync.txt
fname = fullfile(out_dir,'sync.txt');
fd = fopen(fname,'w');
[upstrokes,downstrokes] = load_nev(fullfile(data_dir,'Events.nev'));
if (upstrokes(2) - upstrokes(1)) > 1e5
    upstrokes = upstrokes(2:end);
end
if (downstrokes(2) - downstrokes(1)) > 1e5
    downstrokes = downstrokes(2:end);
end
if upstrokes(1) < downstrokes(1)
    pulse_times = upstrokes;
else
    pulse_times = downstrokes;
end
sync_samples = round(interp1(timestamps,1+(SAMPLES_PER_BLOCK*(0:(length(timestamps)-1))),pulse_times)); %JRM NOTE: THIS LINE NEEDS DEBUGGING!!
fprintf(fd,'%d\n',sync_samples);
fclose(fd);




function[n] =regexp_grab_token(string,pat)
%takes a regexp that captures one token, and returns it as a number (if it
%can). 
%find_num_in_str(string,pat)
[s,f,t]=regexp(string,pat);

if length(t)==0
  n=[];
  return;
end

token=string(t{1}(1):t{1}(2));

n=str2num(token);

if isempty(n) %if we couldn't return it as a number
  n=token; %then return a string. better than nothing, i guess.
end
