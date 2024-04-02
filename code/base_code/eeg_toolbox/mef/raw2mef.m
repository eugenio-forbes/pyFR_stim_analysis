function raw2mef(inData, samplingRate, outFileName, subjectPassword, sessionPassword, block_interval)
%RAW2MEF Converts raw EEG data to a compressed MEF file
%   raw2mef(inData, samplingRate, outFileName, subjectPassword, sessionPassword, block_interval)
%
% INPUTS:
%   inData: vector of input data (in int32 format)
%   samplingRate: sampling rate of the data
%   outFileName: path of the output file to save the compressed data to
%   subjectPassword: password for everything (can open everything,
%       including personal patient data)
%   sessionPassword: secondary password (can only open the EEG
%       data and non-confidential information about the session)
%   [block_interval = 1]: the length (in seconds) of the compression block
%       size. Longer blocks give higher compression ratios, but also take
%       longer to extract the data.