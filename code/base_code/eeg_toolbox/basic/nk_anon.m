function addresses = nk_anon(nk_dir);
% This function anonymizes a Nihon Kohden file.  It does three things:
% 1. removes the patient's name from the EEG file
% 2. it deletes the CMT file (i don't know what this file is for anyhow)
% 3. it deletes the BFT file (we don't konw what that does either...)
% updated 12-16-2013 to allow for de-id multiple eeg files within a directory - ashwin g ramayya

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% .EEG: print and write over name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=dir([nk_dir '/*.EEG']); 
fprintf('\nFound %d .EEG file(s)',length(d)
for i = 1:length(d)
    EEG_file=fullfile(nk_dir,d(i).name);
    %assert(length(d)==1,'Expected 1 .EEG file, found %d',length(d));
    fid=fopen(EEG_file,'r+');
    fseek(fid,79,'bof'); 
    name=fread(fid,32,'*char')';
    XX=repmat('X',1,32);
    fseek(fid,79,'bof');
    fwrite(fid,XX,'char'); %overwrite the name in the device block
    fclose(fid);
    fprintf('\n%s: found NAME=%s...REMOVING.\n',EEG_file,name);
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % .PNT: print and write over name
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d=dir([nk_dir '/*.pnt']);  PNT_file=fullfile(nk_dir,d.name);
% assert(length(d)==1,'Expected 1 .PNT file, found %d',length(d));
% fid=fopen(PNT_file,'r+');
% fseek(fid,79,'bof'); %overwrite the name in the device block
% name=fread(fid,32,'*char')';
% XX=repmat('X',1,32);
% fseek(fid,79,'bof');
% fwrite(fid,XX,'char');
% fclose(fid);
% fprintf('%s: found NAME=%s...REMOVING.\n',PNT_file,name);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % .LOG: print and write over name
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d=dir([nk_dir '/*.LOG']);  LOG_file=fullfile(nk_dir,d.name);
% assert(length(d)==1,'Expected 1 .LOG file, found %d',length(d));
% fid=fopen(LOG_file,'r+');
% fseek(fid,79,'bof'); 
% name=fread(fid,32,'*char')';
% XX=repmat('X',1,32);
% fseek(fid,79,'bof');
% fwrite(fid,XX,'char'); %overwrite the name in the device block
% fclose(fid);
% fprintf('%s: found NAME=%s...REMOVING.\n',LOG_file,name);

%%%%%%%%%%%%%%%%%%
%Completion status
%%%%%%%%%%%%%%%%%%
d=dir([nk_dir '/*.21E']);  REF_file=fullfile(nk_dir,d.name);
fprintf('\nDone anonymizing. Run again to make sure names are gone.\n')
fprintf('OK to copy %s and %s to RHINO.\n\n',EEG_file,REF_file);

% %remove the CMT and BFT files
% d=dir([nk_dir '/*.CMT']);  
% if length(d)>0
%   CMT_file=fullfile(nk_dir,d.name);
%   fprintf('deleting %s\n',CMT_file);
%   delete(CMT_file);
% end

% d=dir([nk_dir '/*.BFT']);  
% if length(d)>0
%   BFT_file=fullfile(nk_dir,d.name);
%   fprintf('deleting %s\n',BFT_file);
%   delete(BFT_file);
% end


