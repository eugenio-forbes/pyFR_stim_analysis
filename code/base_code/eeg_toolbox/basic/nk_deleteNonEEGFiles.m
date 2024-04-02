function[] = nk_deleteNonEEGFiles(nk_dir)
%This function deletes all non-eeg files within an nk_directory.
%Particularly files that dont have .21E and .EEG 

%Inputs
%nk_dir         the directory where the nk data lives

%written by ashwin g ramayya (12-16-2013)

cd(nk_dir);

d = dir;

for i = 1:length(d)
   if ~d(i).isdir
      c = textscan(d(i).name,'%s'); c = c{1}{1};
      if strcmp(c((end-4):end),'.21E') || strcmp(c((end-4):end),'.EEG')
          continue
      else
          status = system(['rm ' c],'-echo');
      end
   end
end