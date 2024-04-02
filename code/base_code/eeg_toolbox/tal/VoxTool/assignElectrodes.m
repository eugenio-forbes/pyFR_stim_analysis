function handles = assignElectrodes( handles, newElecNames )
%ASSIGNELECTRODES assigns an electrode or set of electrodes specified in
% newElecNames to the coordinates specified in 
% handles.selectedPoints.coords.
% Returns handles, a structure that includes the changes that are made to
% the assigned and unassigned electrode lists
    
%Gets currently assigned points
coords = handles.selectedPoints.coords;

%If unassignedElectrodes exists, attempt to pull from that list, otherwise
%ignore.
if isfield(handles,'unassignedElectrodes')
    oldUnass=handles.unassignedElectrodes;
else
    oldUnass={};
end
newUnass = oldUnass;
removed = [];

%Make copies of all of the handles.setPoints fields to be added to
newCoords = handles.setPoints.coords;
numExisting = length(handles.setPoints.names);
newNames = cell(numExisting+length(newElecNames), 1);    
newBaseNames = cell(numExisting+length(newElecNames),1);
newTypes = nan(numExisting+length(newElecNames),1);
newSizes = nan(numExisting+length(newElecNames),2);
newGroupNum = nan(numExisting+length(newElecNames),1);

%If assigned electrodes already exist, place those as the first elements in
%the list
if numExisting>0
    newNames(1:numExisting)=handles.setPoints.names;
    newBaseNames(1:numExisting)=handles.setPoints.baseNames;
    newTypes(1:numExisting) = handles.setPoints.type;
    disp(size(handles.setPoints.groupNum));
    newGroupNum(1:numExisting,:) = handles.setPoints.groupNum;
    newSizes(1:numExisting, :) = handles.setPoints.size;
end

%Get all of the unique basenames (basename=RSTB for electrode RSTB6, e.g.)
uniqueBaseNames = unique(handles.setPoints.baseNames);

%Loop through all of the assigned coordinates
for i=1:size(coords, 1)
    %Switch the coloring of the selected point to cyan
    locatePoint(handles, [coords(i,1), coords(i,2) , coords(i,3)], ...
        false, 'cyan', newElecNames{i});
    
    %Check if the new electrode name is in the unassigned electrode list.
    %If it is, add that index to a list of indices to be removed from the
    %unassigned list.
    index = strmatch(newElecNames{i}, oldUnass);
    removed=[removed, index];
    if ~isempty(index)
        newElecName = oldUnass(index,:);
    else
        newElecName = newElecNames{i};
    end

    %Get the basename and number
    splitName = regexp(newElecName,'\d+','split');
    baseName = splitName{1};
    number = regexp(newElecName,'\d+','match');
    if iscell(number)
        try
        number = str2double(number{1});
        catch
            number = 1;
        end
    end
    %If there is no number, assume it is the first
    if isempty(number)
        number = 1;
    end
    
    %Get the size and type of the electrode array
    sizes = cellstr(get(handles.elecSizePopup,'String'));
    this_size = sizes{get(handles.elecSizePopup,'value')};
    this_size = str2double(regexp(this_size,'x','split'));
    type = get(get(handles.electrodeTypePanel,'SelectedObject'),'String');
    type = type(1);
    
    %Checm to see if the basename of this electrode has been added before
    if any(ismember(handles.setPoints.baseNames,baseName))
        old_i = find(ismember(handles.setPoints.baseNames,baseName),1,'first');
        old_type = handles.setPoints.type(old_i);
        old_size = handles.setPoints.size(old_i,:);
        old_groupNum = handles.setPoints.groupNum(old_i);
        old_name = handles.setPoints.names{old_i};
        
        % Error checking. Too many elecs?
        if number>this_size(1)*this_size(2)
            button = questdlg(sprintf('Electrode %s has size %dx%d. Continue anyway?',...
                newElecName, this_size(1), this_size(2)),'Electrode size mismatch','Yes','No','No');
            if strcmp(button,'No')
                return
            end
        end
        %make sure the type of electrode did not change
        if type~=old_type 
            msgbox(sprintf('Cannot assign.\nElectrode %s has type %s, electrode %s has type %s',...
                old_name, old_type, newElecName, type))
            return 
        %make sure the size of the electrode array did not change
        elseif ~all(this_size==old_size)
            msgbox(sprintf('Cannot assign.\nElectrode %s has size %dx%d, electrode %s has size %dx%d',...
                old_name, old_size(1), old_size(2), newElecName, this_size(1), this_size(2)))
            return
        %If everything else is good, assigne the properties
        else
            newTypes(numExisting+i) = type;
            newBaseNames{numExisting+i} = baseName;
            newSizes(numExisting+i,:) = this_size;
            newGroupNum(numExisting+i) = old_groupNum;
            uniqueBaseNames{length(uniqueBaseNames)+1} = baseName;
        end
    else
        newTypes(numExisting+i) = type;
        newBaseNames{numExisting+i} = baseName;
        newSizes(numExisting+i,:) = this_size;
        if any(ismember(uniqueBaseNames,baseName))
            newGroupNum(numExisting+i) = find(ismember(uniqueBaseNames, baseName),1,'first');
        else
            newGroupNum(numExisting+i) = length(uniqueBaseNames)+1;
            uniqueBaseNames{length(uniqueBaseNames)+1} = baseName;
        end
    end
    newCoords = [newCoords;coords];
    newNames{length(handles.setPoints.names)+i} = newElecName;
end
handles.selectedPoints.coords = [];

%If the unassigned electrode field exists, remove the added electrodes from
%it
newUnass(removed,:)=[];
if isfield(handles,'unassignedElectrodes')
    newUnass = sort(newUnass);
    handles.unassignedElectrodes = newUnass;
end

%color the points that have been assigned based on their group number, such
%that all electrodes with the same basename will be colored the same.
for i=1:length(newElecNames)
    splitName = regexp(newElecNames{i},'\d+','split');
    baseName = splitName{1};
    if any(strcmp(handles.colormap_labels, baseName))
        colormap_i = find(strcmp(handles.colormap_labels, baseName));
    else
        handles.colormap_labels{...
            find(cellfun(@(x)isempty(x),handles.colormap_labels),1, 'first')} = ...
            baseName;
        colormap_i = find(strcmp(handles.colormap_labels, baseName));
    end
    %colormap_i = newGroupNum(numExisting+i);
    [~, cluster] = getPointCluster(handles, coords(i,:)');
    h = plot3(cluster(:,1), cluster(:,2), cluster(:,3), '.');
    set(h,'color',handles.colormap(colormap_i,:));
    set(h,'tag',newElecNames{i});
end

%sort all of the electrodes alphabetically
newNames_padded = cellfun(@(x)regexprep(x,'([^0-9])([0-9])$','$10$2'),newNames,'uniformoutput',false);
[~, ind] = sort(newNames_padded);
newNames = newNames(ind);
newCoords = newCoords(ind,:);
newGroupNum = newGroupNum(ind);
newBaseNames = newBaseNames(ind);
newTypes = newTypes(ind);
newSizes = newSizes(ind,:);





%assign back to handles structure
handles.setPoints.coords = newCoords;
handles.setPoints.names = newNames;
handles.setPoints.groupNum = newGroupNum;
handles.setPoints.baseNames = newBaseNames;
handles.setPoints.type = newTypes;
handles.setPoints.size = newSizes;
if isfield(handles,'unassignedElectrodes')
    makeList(handles.unassignedList, newUnass);
end

%make the list and save it
handles.assignedElectrodes = newNames;
makeList(handles.assignedList, newNames, newCoords, newTypes, newSizes);
handles.saved = false;
guidata(handles.assignButton, handles)


