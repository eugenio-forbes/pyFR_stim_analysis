function [eeg_data,eeg_lengths] = decomp_mef_events(file_name, start_index_array, duration, password)
%DECOMP_MEF_EVENTS Retrieves multiple events from a compressed MEF file
%   [eeg_data,eeg_lengths] = decomp_mef_events(file_name, start_index_array, 
%       duration, password)
%
% INPUTS:
%   file_name: path (absolute or relative) of the target MEF file
%   start_index_array: vector of the start indices (in samples) of all 
%       events to be extracted
%   duration: scalar representing the length in samples of extracted events 
%       (i.e., all events must be the same length)
%   password: password (subject or session) of the file to open
%
% OUTPUTS:
%   eeg_data: a DxN int32 matrix of EEG data (D = duration, N = number of
%       events)
%   eeg_lengths: a Nx1 vector of the durations (in samples) of the
%       extracted events. Should always be equal to the input duration,
%       unless the start_index + duration indicates a point past the end of
%       the file