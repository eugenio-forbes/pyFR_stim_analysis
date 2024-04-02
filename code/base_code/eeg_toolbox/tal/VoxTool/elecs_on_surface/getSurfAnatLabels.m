function [anatLabels] = getSurfAnatLabels(elecLocs,path2surf,path2annot, c)
%getSurfAnatLabels(electrodeLocs, path2Annot)
%   This function gets surface anatomical labels associated with a
%   particular electrode
%
% INPUTS:
%   electrodeLocs   -   %matrix of XYZ coords (nElecs x 3)
%                        -must be in surface RAS space (snapped or not)
%                        -must only include locations from a single hemisphere
%
%   path2surf -         The path to the surface file for this subject
%                       (lh.pial)
%
%   path2annot -  The path to the annotation file for the same
%                       hemisphere (lh.aparc.annot)
%
% OR:
%   electrodeLocs   -   %matrix of XYZ coords (nElecs x 3)
%                        -must be in surface RAS space (snapped or not)
%                        -must only include locations from a single hemisphere
%
%   V      -         The vertices read from the surface file for this subj
%                       (lh.pial)
%
%   labels -         The labels read from the annot file for this subj
%                       hemisphere (lh.aparc.annot)
%
%   colortable -     The colortable read from the annot for this subj
%
% OUTPUTS:
%   anatLabels - cell array of anatomical labels associated with the elctrodes
% created 11/13 by Ashwin Ramayya (ashwinramayya@gmail.com)
% Default uses aparc.annot
%http://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation
%     'unknown'
%     'bankssts'
%     'caudalanteriorcingulate'
%     'caudalmiddlefrontal'
%     'corpuscallosum'
%     'cuneus'
%     'entorhinal'
%     'fusiform'
%     'inferiorparietal'
%     'inferiortemporal'
%     'isthmuscingulate'
%     'lateraloccipital'
%     'lateralorbitofrontal'
%     'lingual'
%     'medialorbitofrontal'
%     'middletemporal'
%     'parahippocampal'
%     'paracentral'
%     'parsopercularis'
%     'parsorbitalis'
%     'parstriangularis'
%     'pericalcarine'
%     'postcentral'
%     'posteriorcingulate'
%     'precentral'
%     'precuneus'
%     'rostralanteriorcingulate'
%     'rostralmiddlefrontal'
%     'superiorfrontal'
%     'superiorparietal'
%     'superiortemporal'
%     'supramarginal'
%     'frontalpole'
%     'temporalpole'
%     'transversetemporal'
%     'insula'

if ischar(path2surf)
    %% load surf
    [v_surf] = read_surf(path2surf);
    %% load annot
    [~,l,c] = read_annotation(path2annot);
else
    v_surf = path2surf;
    l = path2annot;
end
%% id vertices assoc with elec locs
elecIdx = dsearchn(double(v_surf), double(elecLocs));

%% get anat labels assoc with each elec
anatLabels = {};
for i = 1:length(elecIdx)
    try
        anatLabels{i} = c.struct_names{l(elecIdx(i))==c.table(:,5)};
    catch
        anatlabels{i} = 'unknown';
    end
end