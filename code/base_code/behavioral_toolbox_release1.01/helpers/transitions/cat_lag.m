function cl = cat_lag(from_pt, to_pt, params)
% CAT_LAG
% function cl = cat_lag(from_pt, to_pt, params)
%
% INPUTS
% params: must have the following fields
%       pres_catlabels
%       cat_type

from_label = params.pres_catlabels(from_pt);
if params.cat_type == 1
  cat_sps = find(params.pres_catlabels==from_label);
else
  % from_pt is one of the list, so alter its label, anything
  % positive will keep it valid
  params.pres_catlabels(from_pt) = params.pres_catlabels(from_pt) + 1;
  cat_sps = find(params.pres_catlabels~=from_label);
end

cl = find(to_pt==cat_sps) - find(from_pt==cat_sps);