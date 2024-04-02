function [] = writesurf_w_elecs(path2surf,elecLocs,eNames,radius,path2surf_w_elecs)
% [path2surf_w_elecs] = writesurf_w_elecs(path2surf,elecLocs,eNames,radius,path2surf_w_elecs);
%   This function writes labels the faces associated with each electrodes
%   and saves it as a .mat file in the subject's tal directory. It does
%   this by appending a fourth column to the faces matrix.
%
% INPUTS:
%   path2surf -         The path to the surface file for this subject
%                       (lh.pial)
%   electrodeLocs   -   %matrix of XYZ coords (nElecs x 3)
%                        -must be in surface RAS space (snapped or not)
%                        -must only include locations from a single hemisphere
%   eNames  -    unique identifier of the electrode (e.g., '101-102', 'eNames' in talStruct)
%   radius  -    mm to use when assigning faces to a particular electrode
%   path2surf_w_elecs - path where the new surface will be saved.
% OUTPUTS:
% created 11/13 by Ashwin Ramayya (ashwinramayya@gmail.com)
%% error check
if length(eNames) ~= size(elecLocs,1)
    error ('eNames must correspond to elecLocs')
end

%% load surf
[v,f] = read_surf(path2surf);

%% update faces so it is one indexed
f = f + 1;

%% update faces with a fourth column
% append a fourth column to faces, this will be used to color the surface
f = [f, zeros(size(f,1),1)];

% get idx assoc with vertices
v_idx = 1:length(v);

for i=1:size(elecLocs,1)
    point = elecLocs(i,:);
    pointRep = repmat(point, size(v,1),1);
    distances = sqrt(sum((pointRep-v).^2,2));
    relevant_idx = v_idx(distances<radius);

    % update 4th column of faces to reflect electrode group
    f(any(ismember(f, relevant_idx),2),4) = i;

end

%% save structure
save(path2surf_w_elecs,'v','f','eNames')
