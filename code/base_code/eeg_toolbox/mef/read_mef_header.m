function header = read_mef_header(file_name, password)
%READ_MEF_HEADER Retrieves the header information from a MEF file
%	header = read_mef_header(file_name, password)
%
% INPUTS:
%   file_name: path (absolute or relative) of the target MEF file
%   password: session or subject password of the file
%
% OUTPUTS:
%   header: a structure containing all the header fields. If the session
%       password was used, confidential fields will appear as 'encrypted'.
