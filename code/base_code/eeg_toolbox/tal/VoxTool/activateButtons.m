function activateButtons( handles )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if isfield(handles,'mriDir') && isfield(handles,'ctDir') && ...
        isfield(handles,'subjName') && isfield(handles,'subjSourceName') &&...
        isfield(handles,'subjDir')
    set(handles.localizeButton,'enable','on')
    if isfield(handles,'voxMother') && isfield(handles,'username') &&...
            isfield(handles,'jacksheet')
        set(handles.registerButton,'enable','on')
        if handles.isRegistered
            set(handles.updateButton,'enable','on')
        else
            set(handles.updateButton,'enable','off')
        end
    else
        set(handles.registerButton,'enable','off')
        set(handles.updateButton,'enable','off')
    end
else
    set(handles.localizeButton,'enable','off');
    set(handles.registerButton,'enable','off');
    set(handles.updateButton,'enable','off');
end

