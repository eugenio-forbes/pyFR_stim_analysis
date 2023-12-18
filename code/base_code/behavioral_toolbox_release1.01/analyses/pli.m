function [plis] = pli(intrusions_matrix,subjects,rec_mask)
% PLI the number of Prior List Intrusions made by each subject
%
% [plis] = pli(recalls_itemnos_matrix,pres_itemnos_matrix, subjects,rec_mask)
%
%  INPUTS: 
%  intrusions_matrix: for each trial (row) and output position (column),
%                     indicates the type of intrusion that was made:
%
%                     0: no intrusion
%
%                     -1: extra-list intrusion (XLI)
%
%                     positive integer: prior-list intrusion (PLI), 
%                     indicating the number of lists back from which the 
%                     PLI was recalled
%
%        subjects:    a column vector which indexes the rows of
%                     intrusions_matrix with a subject number (or other 
%                     identifier).  That is, the recall trials of subject S 
%                     should be located in:
%                     intrrusions_matrix(find(subjects==S), :)
%
%
%        rec_mask:    if given (e.g. to exclude an intrusion repeated
%                     within a trial), a logical matrix of the same shape as 
%                     intrusions_matrix, which is false at positions (i, j)
%                     where the value at intrusions_matrix(i, j) should be
%                     excluded from the calculation of the probability of
%                     recall. If NOT given, a standard blank mask is used,
%                     which is true at all output positions.
%                     
%
% OUTPUTS:
%        plis:        a vector of total number of PLIs. Its rows are
%                     indexed by subject.


% sanity checks:
if ~exist('intrusions_matrix', 'var')
  error('You must pass an intrusions_matrix.')
elseif ~exist('subjects', 'var')
  error('You must pass a subjects vector.')
elseif size(intrusions_matrix, 1) ~= length(subjects)
  error('intrusions_matrix must have the same number of rows as subjects.')
end

if ~exist('rec_mask', 'var')
  % create blank mask if none was given
  rec_mask = make_blank_mask(intrusions_matrix);
elseif size(rec_mask) ~= size(intrusions_matrix)
  error('intrusions_matrix and rec_mask must have the same shape.')
end

% apply_by_index applies pli_for_subj on each subjects recall data 
% separately
% As with other functions, pli_for_subj does most of the work here
plis = apply_by_index(@pli_for_subj, subjects, 1, {intrusions_matrix,rec_mask}); 

function pli_count = pli_for_subj(intrusions_matrix,rec_mask)
  % add up the total number of positive integers in a subject's intrusions
  % matrix
  pli_count = sum(intrusions_matrix(rec_mask)>0);
  