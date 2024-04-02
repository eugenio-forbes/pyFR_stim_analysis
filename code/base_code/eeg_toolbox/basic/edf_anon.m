function edf_anon(edf_dir);
%
% DESCRIPTION:
%  This function anonymizes .edf files.  
%
% FUNCTION:
%  function addresses = edf_anon(edf_dir);
%
% INPUT:
%  edf_dir: The directory where the edf file(s) is located. 
%
% OUTPUT:
%  The function replaces the name string with a seris of "XXXX" and
%  re-writes the file. 
%
% NOTES: 
%  (1) written by jfburke (4/14/11)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% .edf: the first 192*8 bits represent 192 'uchar' variables (8-bit
% unsigned character).  The name filed begins at the 9th 'uchar'
% (i.e. 72 bit) and is 80 uchar's long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d1=dir([edf_dir '/*.edf']);
d2=dir([edf_dir '/*.EDF']);
d = [d1 d2];

%assert(length(d)==1,'Expected 1 .edf file, found %d',length(d));
fprintf('\n\nFOUND %d .EDF FILES(S)\n',length(d))
for k=1:length(d)
  EEG_file=fullfile(edf_dir,d(k).name);
  anonThisFile_local(EEG_file)
end

function anonThisFile_local(EEG_file)
  fid=fopen(EEG_file,'r+');

  % seek to the 8th uchar.. the nest ting you read in will be the 9th
  % uchar
  fseek(fid,8,'bof'); 
  
  % read in the name
  PID_raw = fread(fid,80,'uchar');
  PID=deblank(char(PID_raw)');
  
  % read in other variables
  record_ID=char(fread(fid,80,'uchar'))';
  DATE=char(fread(fid,8,'uchar'))';
  TIME=char(fread(fid,5,'uchar'))';
  fseek(fid,3,'cof'); 
  HeadLen   = char(fread(fid,8,'uchar'))';
  reserved1 = char(fread(fid,44,'uchar'))';
  NRec      = str2double(char(fread(fid,8,'uchar'))');     % 8 Bytes  # of data records
  Dur       = str2double(char(fread(fid,8,'uchar'))');     % 8 Bytes  # duration of data record in sec
  NS        = str2double(char(fread(fid,4,'uchar'))');     % 4 Bytes  # of signals
  
  % replace patient name with this
  XX=repmat('X',1,80);
  
  % print a display
  fprintf('\n  Anonymizing %s....\n',EEG_file)
  fprintf('\n')
  fprintf('\tFound patient name: %s\n',PID)
  fprintf('\n')
  fprintf('\tSession details:\n',PID)
  fprintf('\t   date.......%s\n',DATE) 
			   %
  fprintf('\t   time.......%s\n',TIME) 
			   %
  fprintf('\t   length.....%2.2f min\n',Dur*NRec/60)
			   %
  fprintf('\t   num Elecs..%d\n',NS) 
			   %
  fprintf('\n')    
  %fprintf('\tReplacing patient name with: %s\n',XX)    
  fprintf('\tReplacing patient name...')
  fseek(fid,8,'bof');
  fwrite(fid,XX,'char'); %overwrite the name in the device block
  fclose(fid);
  fprintf('DONE.\n')
  fprintf('\tRe-run to verify correct anonymization.\n\n')
  

