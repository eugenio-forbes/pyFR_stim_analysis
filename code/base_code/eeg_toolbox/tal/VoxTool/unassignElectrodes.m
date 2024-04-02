function handles = unassignElectrodes( handles, selections )
%UNASSIGNELECTRODES removes points from the assigned electrodes list, and
%deletes all points having to do with them being assigned from the CT

%If unassignedElectrodes already exists, find that field so that it can be
%added to
un= isfield(handles,'unassignedElectrodes');
if un
    unassigned = handles.unassignedElectrodes;
end

%Make copies of all of the fields in handles.setPoints so they can be
%removed from and then reassigned.
assigned_coords = handles.setPoints.coords;
assigned_names = handles.setPoints.names;
assigned_baseNames = handles.setPoints.baseNames;
assigned_groupNums = handles.setPoints.groupNum;
assigned_size = handles.setPoints.size;
assigned_type = handles.setPoints.type;

%The list of points to be deleted
to_delete = nan(length(selections),1);
for i=1:length(selections)
    s_i = selections(i);
    % If unassigned exists, add the newly removed points back to the
    % unassigned list
    if un
        unassigned(length(unassigned)+1,1) = assigned_names(s_i);
    end
    
    % Remove points tagged with this name from the graph
    tag = [assigned_names{selections(i)}];
    h = findobj(handles.mainAxes,'tag',tag);
    delete(h);
    
    to_delete(i) = s_i;
end

%remove all of the to-be-deleted points from each of the setPoints fields
assigned_coords(to_delete,:)=[];
assigned_names(to_delete)=[];
assigned_baseNames(to_delete) = [];
assigned_groupNums(to_delete) = [];
assigned_size(to_delete,:) = [];
assigned_type(to_delete) = [];

%Sort all points alphabetically
assigned_names_padded = cellfun(@(x)regexprep(x,'([^0-9])([0-9])$','$10$2'),assigned_names,'uniformoutput',false);
[~,ind]= sort(assigned_names_padded);
assigned_names = assigned_names(ind);
assigned_coords = assigned_coords(ind,:);
assigned_names = assigned_names(ind);
assigned_baseNames = assigned_baseNames(ind);
assigned_groupNums = assigned_groupNums(ind);
assigned_size= assigned_size(ind, :);
assigned_type = assigned_type(ind);

% If the unassigned list exists, reassign to it
% TODO: make sure this works. I think it doesn't.
if un
    [unassigned1,ind] = sort(unassigned(:,1));
    unassigned2 = unassigned(ind,2);
    unassigned = {unassigned1,unassigned2};
end

% reassign the new fields to handles and save to both figures
handles.setPoints.coords = assigned_coords;
handles.setPoints.names = assigned_names;
handles.setPoints.baseNames = assigned_baseNames;
handles.setPoints.groupNum = assigned_groupNums;
handles.setPoints.size = assigned_size;
handles.setPoints.type = assigned_type;
if un
    handles.unassignedElectrodes = unassigned;
end

handles.saved = false;

guidata(handles.controlFigure,handles);
guidata(handles.CT3D_fig, handles);

%make the unassigned and assigned lists.
if un
makeList(handles.unassignedList, unassigned{:,1})
end
makeList(handles.assignedList, assigned_names, assigned_coords, assigned_type, assigned_size);

end

