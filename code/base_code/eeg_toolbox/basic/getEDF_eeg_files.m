function EEG_file = getEDF_eeg_files(d_in)
%
%
%
%
%
%

d=dir([d_in '/*.edf']);  


if length(d)<1;
  error('Expected 1 .EEG file, found none')
end


if length(d)>1;
  fprintf('\n')
  fprintf('  WARNING: Found %d .EEG files: Splitting Both.\n',length(d))
  for k=1:length(d)
    fprintf('             - %s\n',d(k).name)
  end
  fprintf('\n')
end
fprintf('\n')

EEG_file={};
for k=1:length(d)
  EEG_file{k} = fullfile(d_in,d(k).name);
end