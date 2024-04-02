function setDir = getLastUsedDirectory( handles, numBack )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if ~exist('numBack','var')
    numBack = 1;
end
if isfield(handles, 'mriDir') || isfield(handles,'ctDir')
    if isfield(handles,'mriDir')
        setDir = handles.mriDir;
    else
        setDir = handles.ctDir;
    end
    for i=1:numBack
        directorySplit = regexp(setDir,'/','split');
        setDir = strrep(setDir,['/',directorySplit{length(directorySplit)}],'');
    end
end

