function emg_to_ncs(emg_files, csc_dir)
% Converts EMG files recorded at Freiburg to Neuralynx Ncs files.
% 
% INPUT:
% emg_files: Cell array of EMG file names; the data will be concatenated
%            in the given order.
%
% csc_dir:   Directory to drop CSC* files into.

RECORD_SIZE = 512;

for fno = 1:length(emg_files)

  % Log last timestamp to fix time offset below.
  if fno > 1
    last_timestamp = timestamps(end);
  end

  % Read file.
  clear timestamps data_samples; % We can run out of memory
                                 % if these are loaded while reading
                                 % the next file.
  [timestamps,data_samples,sample_freq] = read_emg(emg_files{fno});
  
  % Fix time offset -- we're assuming this is a contiguous
  % recording, but each file starts at time 0.
  if fno > 1
    timestamps = last_timestamp + mean(diff(timestamps)) + timestamps;
  end

  % File info. 
  num_channels = size(data_samples,1);
  num_samples  = size(data_samples,2);
  num_records  = length(timestamps);

  % Scale to volts.
  data_samples = data_samples .* 1e6;

  % Write each channel to a separate file.
  for chno = 1:num_channels
    f = fopen(fullfile(csc_dir, ['CSC' num2str(chno) '.Ncs']), 'a', 'l');

    % Write header the first time around,
    % otherwise append data to end.
    if fno == 1
      % Wave_clus relies on a very particular (hard coded) format 
      % for this, which we replicate here although some of the 
      % fields carry no extra information.
      headerSize = 16384;
      fprintf(f,'######## Neuralynx Data File Header\n');
      fprintf(f,'## File Name: CSC%d.Ncs\n', chno);
      D = dir(emg_files{1});
      fprintf(f,'## Time Opened: (m/d/y): %s\tAt Time: %s\n',...
                datestr(D.datenum,'mm/dd/yyyy'),...
                datestr(D.datenum,'HH:MM:SS.FFF'));
      fprintf(f,'-CheetahRev NaN\n');
      fprintf(f,'-NLX_Base_Class_Name CSC%d\n', chno);
      fprintf(f,'-NLX_Base_Class_Type NaN\n');
      fprintf(f,'-RecordSize 1044\n');
      fprintf(f,'-ADChannel NaN\n');
      fprintf(f,'-ADGain NaN\n');
      fprintf(f,'-AmpGain NaN\n');
      fprintf(f,'-AmpLowCut NaN\n');
      fprintf(f,'-AmpHiCut NaN\n');
      fprintf(f,'-SubSamplingInterleave NaN\n');
      fprintf(f,'-SamplingFrequency %d\n', sample_freq);
      fprintf(f,'-ADBitVolts %f\n',  1/1e6);
      fprintf(f,'-ADMaxValue NaN\n');
      fwrite(f,zeros(headerSize-ftell(f),1));
    else
      fseek(f,0,'eof');
    end

    for rno = 1:num_records
      % Timestamp.
      fwrite(f,timestamps(rno),'uint64');

      % Channel.
      fwrite(f,chno,'uint32');

      % Sampling frequency.
      fwrite(f,sample_freq,'uint32');

      % Num. valid samples.
      if rno < num_records
        num_valid_samples = RECORD_SIZE;
      else
        num_valid_samples = rem(num_samples, RECORD_SIZE);
      end
      fwrite(f,num_valid_samples,'uint32');

      % Finally, the data, making sure to fill the record
      % to the end even if there's partial data.
      begin_sample = (rno-1)*RECORD_SIZE+1;
      end_sample   = min(begin_sample+RECORD_SIZE-1,num_samples);
      fwrite(f,data_samples(chno,begin_sample:end_sample),'int16');
      padding = RECORD_SIZE-rem(end_sample,RECORD_SIZE);
      if padding > 0 && padding < RECORD_SIZE
        fwrite(f,zeros(padding,1),'int16');
      end
    end

    fclose(f);
  end
end

function [timestamps,data_samples,sample_freq] = read_emg(emg_file)

HEADERSIZE = 56788;
CHANNEL_NO_OFS = 42;     % int16
SAMPLING_RATE_OFS = 76;  % float32

BYTES_PER_SAMPLE = 4;

SAMPLES_PER_BLOCK = 512;

fid = fopen(emg_file, 'r');

% read header information
fseek(fid, CHANNEL_NO_OFS, 'bof');
channelNo = fread(fid, 1, 'int16');
fseek(fid, SAMPLING_RATE_OFS, 'bof');
sample_freq = fread(fid, 1, 'float32');

% read data
fseek(fid, HEADERSIZE, 'bof');
data_samples = fread(fid, [channelNo, Inf], 'float32');
timestamps = ((0:SAMPLES_PER_BLOCK:size(data_samples,2)) / sample_freq) * 1e6;

fclose(fid);

