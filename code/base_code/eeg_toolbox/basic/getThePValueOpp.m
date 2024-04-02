function [pos] = getThePValueOpp(tmpVect,tmpScal,nVallnDist)

tmpVect = sort(tmpVect);

tmp = find(tmpVect>tmpScal);

if~isempty(tmp)
    pos = ((tmp(1)-1) / nVallnDist);


else
    if tmpScal>=max(tmpVect)
        if tmpScal~=min(tmpVect)
        pos=1;
        else
        pos = NaN
        end
    else
    pos = NaN
    end
end  