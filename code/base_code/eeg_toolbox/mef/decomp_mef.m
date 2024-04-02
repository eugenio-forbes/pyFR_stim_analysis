function eeg_data = decomp_mef(file_name, start_index, stop_index, password)
%DECOMP_MEF Retrieves a section of EEG data from a compressed MEF file
%   eeg_data = decomp_mef(file_name, start_index, stop_index, password)
%
% INPUTS:
%   file_name: path (absolute or relative) of the target MEF file
%   start_index: sample number of the start of the section to extract
%   end_index: sample number of the end of the section to extract
%   password: password (subject or session) of the file to open
%
% OUTPUTS:
%   eeg_data: a Dx1 int32 vector of EEG data, where D is the number of
%       samples between start_index and end_index (inclusive)
