function [transits_array] = cat_transitions(serial_position, ...
				            prior_recalls, transition, ...
				            params)
% CAT_TRANSITIONS  Returns the possible lags from a given serial
%                  position, excluding lags to serial positions which
%                  have already been recalled, conditional on those
%                  lags being of the same category (or a different
%                  category) from the just-recalled item
%
%  [transits_array] = possible_transitions(serial_position, prior_recalls, ...
%                                          transition, params)
%
%  NOTE: this function is meant to be passed as a condition function to
%  conditional_transitions(); its arguments are dictated by the requirements
%  of that function.
%
%  INPUTS:
%  serial_position:  the serial position from which possible transitions should
%                    be calculated.  If this value is less than 1 (i.e., the
%                    'recall' was an intrusion or empty cell), an empty
%                    array is returned.
%
%    prior_recalls:  a row vector of serial positions which have already been
%                    recalled; transitions to these serial positions are
%                    excluded from the output
%
%       transition:  the current transition value (accepted here to meet the
%                    requirements of conditional_transitions(), but not used)
%
%           params:  a structure containing a field 'list_length' which
%                    specifies the length of the list, a field
%                    'cat_type', which specifies whether possible
%                    same-category, or possible other-category
%                    transitions are returned, and a field
%                    'pres_catlabels' which specifies the category
%                    labels of all of the list items.
%
%  OUTPUTS:
%   transits_array:  a row vector of possible transitions from the current
%                    serial position (excluding transitions to
%                    previously-recalled serial positions).
%                    transitions are returned as positions relative
%                    to the originating item, ignoring invalid
%                    intervening items.
%
%  EXAMPLES:
%  >> sp = 4;
%  >> prior_recalls = [3 2];
%  >> params.list_length = 6;
%  >> params.cat_type = 1;
%  >> params.pres_catlabels = [1 0 1 1 0 1];
%  
%  % with no prior recalls 
%  >> cat_transitions(sp, [], -1, params)
%  ans =
%     -2    -1    1
%     
%  % with two prior recalls, transitions to those positions are excluded:
%  >> cat_transitions(sp, prior_recalls, -1, params)
%  ans = 
%     -2     1
%  

% sanity checks
if ~exist('serial_position', 'var')
  error('You must pass a serial position')
elseif ~exist('prior_recalls', 'var')
  error('You must pass a prior_recalls vector')
elseif ~exist('params', 'var')
  error('You must pass a params struct')
elseif ~isfield(params, 'list_length')
  error('params must have a list_length field')
elseif ~isfield(params, 'pres_catlabels')
  error('params must have a pres_catlabels field')
elseif ~isfield(params, 'cat_type')
  error('params must have a cat_type field')
end

list_length = params.list_length;
pres_catlabels = params.pres_catlabels;
cat_type = params.cat_type;

if ~isfield(params, 'to_mask_pres')
  params.to_mask_pres = ones(1,list_length);
end
to_mask_pres = params.to_mask_pres;

% calculate the cat_position of this serial_position
from_label = pres_catlabels(serial_position);
if cat_type == 1
  cat_sps = find(pres_catlabels==from_label);
  cat_pos = find(serial_position==cat_sps);
  total_poss = sum(pres_catlabels==from_label);
else
  pres_catlabels(serial_position) = pres_catlabels(serial_position) + 1;
  cat_sps = find(pres_catlabels~=from_label);
  cat_pos = find(serial_position==cat_sps);
  total_poss = sum(pres_catlabels~=from_label);
end

% Transitions from an intrusion or empty cell are never allowed
if serial_position < 1
  transits_array = [];
  return
end

% Nor are transitions from a repeated word 
if any(ismember(prior_recalls, cat_pos))
  transits_array = [];  
  return
end

% Generally speaking, all transitions of magnitude greater than 0
% are possible, from (-serial_position + 1) to 
% (list_length - serial_position), inclusive

transits_array = [-cat_pos + 1 : -1, 1 : total_poss - cat_pos];

% Remove transitions to previously recalled serial positions
% translate into conditional category space
disallowed_positions = [find(to_mask_pres==0) prior_recalls];

disallowed_cat_pos = [];
for i = 1:length(disallowed_positions)
  if ~isempty(find(disallowed_positions(i)==cat_sps))
    disallowed_cat_pos(i) = find(disallowed_positions(i)==cat_sps);
  end
end

disallowed_transits = unique(disallowed_cat_pos - cat_pos);

transits_array = setdiff(transits_array, disallowed_transits);