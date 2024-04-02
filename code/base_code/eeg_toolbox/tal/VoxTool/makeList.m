function makeList( handleToList , names, coords, types, sizes)
%MAKELIST makes the cell string used to make the lists. Depending on the
%number of arguments passed in, it either makes a list with lal of them, or
%just with the names.
    newList = cell(length(names),1);
    for i=1:length(newList)

        if exist('coords','var')
            newList{i} = sprintf('%s %d %d %d %s %dx%d',...
                names{i}, coords(i,2), coords(i,1), coords(i,3), ...
                types(i), sizes(i,1), sizes(i,2));
        else
            newList{i} = [names{i}];
        end
        
    end

    set(handleToList,'string',newList);
    set(handleToList,'value',1);

end