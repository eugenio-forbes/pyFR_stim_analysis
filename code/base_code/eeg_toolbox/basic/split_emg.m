function[] = split_emg(emg_path)
%
% 3-3-10   JRM, MN   Wrote it.

lfpdir = fullfile(emg_path,'microeeg');
mkdir(lfpdir);

emg_files = dir(fullfile(emg_path,'*.raw.emg'));

fprintf('splitting %d files',length(emg_files));
for i = 1:length(emg_files)   
    [~,raw_file] = fileparts(emg_files(i).name);
    [~,chan_filename] = fileparts(raw_file);        
    data = emgpool({fullfile(emg_path,emg_files(i).name)});    
    for chan = data.channels
        outfile = fullfile(lfpdir,sprintf('%s.%03d',chan_filename,chan));
        if ~exist(outfile,'file')
            chan_data = data.read(chan,data.samples_sum,1);
            fd = fopen(outfile,'w','l');
            nsamples = fwrite(fd,chan_data,'single');
            fclose(fd);
            assert(nsamples == length(chan_data),'some sync pulses were not written to sync.txt!');
        end
        fprintf('.');
    end
    fprintf('-');
end
fprintf('\n');

fname = fullfile(lfpdir,'params.txt');
fd = fopen(fname,'w');
fprintf(fd,'samplerate %0.5f\n',19996.2); %JRM: FIX THIS!!
fprintf(fd,'dataformat ''single''\n');
fprintf(fd,'gain %f\n',1.0);
fclose(fd);