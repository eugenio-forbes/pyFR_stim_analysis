function rec_mask = update_recalls_mask(recalls, rec_mask, pres_mask)
% UPDATE_RECALLS_MASK update the recalls mask (rec_mask) with invalid
% positions as specified in the presentation mask (pres_mask).
%
% function rec_mask = update_recalls_mask(recalls, rec_mask, pres_mask)
%
% INPUTS:
% recalls:   a matrix whose elements are serial positions of
%            recalled items.  The rows of this matrix should
%            represent recalls made by a single subject on a
%            single trial.
%
% rec_mask:  a logical matrix of the same shape as
%            recalls, which is false at positions (i,j)
%            where the value at recalls(i,j) should be
%            excluded from the calculation of the probability
%            of recall.
%
% pres_mask: a logical matrix of the same shape as the
%            item presentation matrix (i.e., num_trials x
%            list_length).  This mask should be true at
%            positions (t, sp) where an item in the condition
%            of interest was presented at serial position sp on
%            trial t, and false everywhere else.



if ~exist('rec_mask','var')
  error('You must provide a recalls mask.');
elseif ~exist('pres_mask','var')
  error('You must provide a presentation mask.');
elseif ~exist('recalls','var')
  error('You must provide a recalls matrix.');
elseif size(rec_mask,1) ~= size(pres_mask,1)
  error('The two masks must have the same number of rows.')
elseif size(recalls,1) ~= size(pres_mask,1)
  error(['The recalls matrix and the mask must have the same number of rows.']);
end

for i=1:size(rec_mask,1)
  mask_these = find(pres_mask(i,:) == 0);
  for j=1:length(mask_these)
    rec_mask(i,recalls(i,:)==mask_these(j)) = 0;
  end
end