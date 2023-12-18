function temp_facts = temp_fact(recalls_matrix, subjects, list_length,...
    pres_mask_from,pres_mask_to,rec_mask_from,rec_mask_to)
%  TEMP_FACT  Computes a lag-based temporal clustering factor.
%
%
%  temp_facts = temp_fact(recalls_matrix, subjects, list_length, ...
%                         mask, pres_mask)
%
%  INPUTS:
%  recalls_matrix:  a matrix whose elements are serial positions of recalled
%                   items.  The rows of this matrix should represent recalls
%                   made by a single subject on a single trial.
%
%        subjects:  a column vector which indexes the rows of recalls
%                   with a subject number (or other identifier).  That is, 
%                   the recall trials of subject S should be located in
%                   recalls_matrix(find(subjects==S), :)
%
%     list_length:  a scalar indicating the number of serial positions in the
%                   presented lists.  serial positions are assumed to run 
%                   from 1:list_length.
%
%  pres_mask_from:  if given, a logical matrix of the same shape as the
%                   item presentation matrix (i.e., num_trials x
%                   list_length).  This mask should be true at
%                   positions (t, sp) where an item in the condition
%                   of interest was presented at serial position sp on
%                   trial t; and false everywhere else.
%
%   pres_mask_to:   a logical matrix with the same structure as
%                   pres_mask_to. pres_mask_from excludes transitions that
%                   terminate on item (t,sp+1).  
%
%   rec_mask_from:  if given, a logical matrix of the same shape as
%                   recalls_matrix, which is false at positions (i, j)
%                   where the value at recalls_matrix(i, j) should be
%                   excluded from the calculation of the probability
%                   of recall.  If NOT given, a standard clean recalls
%                   mask is used, which excludes repeats, intrusions
%                   and empty cells.  from_mask excludes transitions
%                   that initiate with item j.
%
%     rec_mask_to:  a logical matrix with the same properties as
%                   rec_mask_from. to_mask excludes transitions that 
%                   terminate on item j+1.
% 
%
%  OUTPUTS:
%      temp_facts:  a vector of temporal clustering factors, one
%                   for each subject.
%
%  NOTE:
%      This version aggregates by trial, then takes the set of
%      trial means.  Prior versions aggregated all transitions per
%      subject and then took the mean.
%

% sanity checks
if ~exist('recalls_matrix', 'var')
  error('You must pass a recalls matrix.')
elseif ~exist('subjects', 'var')
  error('You must pass a subjects vector.')
elseif ~exist('list_length', 'var')
  error('You must pass a list length.') 
elseif size(recalls_matrix, 1) ~= length(subjects)
  error('recalls matrix must have the same number of rows as subjects.')
end

if ~exist('pres_mask_from', 'var')
  pres_mask_from = ones(length(subjects),list_length);
end
if ~exist('pres_mask_to', 'var')
  pres_mask_to = pres_mask_from;
end
if ~exist('rec_mask_from','var')
    rec_mask_from = make_clean_recalls_mask2d(recalls_matrix);
end
if ~exist('rec_mask_to','var')
    rec_mask_to = rec_mask_from;
end

% this information will get passed to general_temp_fact, which expects
% each presented item has an associated category. so, set the category
% labels as if all items were presented with the same category.
pres_catlabels = ones(length(subjects),list_length);
cat_type = 1;

% general_temp_fact will do the work:
temp_facts = general_temp_fact(recalls_matrix, pres_catlabels, ...
			       subjects, cat_type, ...
			       rec_mask_from, rec_mask_to, ...
			       pres_mask_from, pres_mask_to, ...
                               false);
