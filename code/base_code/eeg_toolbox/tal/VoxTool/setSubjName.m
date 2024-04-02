function handles = setSubjName(handles)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    handles.subjDir = getLastUsedDirectory(handles,2);
    handles.sourceDir = handles.subjDir;
    set(handles.sourceDirEdit,'string',handles.subjDir);
    set(handles.subjDirEdit,'string',handles.subjDir);
    subjDirSplit = regexp(handles.subjDir,'/','split');
    handles.subjName = subjDirSplit{length(subjDirSplit)};
    set(handles.subjNameEdit,'string',handles.subjName);
    handles.subjSourceName = handles.subjName;
    set(handles.subjSourceEdit,'string',handles.subjSourceName);
    if exist([handles.subjDir '/docs/jacksheet.txt'])
        handles.jacksheet = [handles.subjDir '/docs/jacksheet.txt'];
        set(handles.jacksheetEdit,'string',handles.jacksheet);
    end
end

