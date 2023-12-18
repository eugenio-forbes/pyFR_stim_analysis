function [xlis] = xli(intrusions_matrix,subjects,rec_mask)
% XLI the number of Extra List Intrusions made by each subject
%
% [xlis] = xli(intrusions_matrix,subjects,rec_mask)
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
%        xlis:        a vector of total number of XLIs. Its rows are indexed
%                     by subject.

% sanity checks:
if ~exist('intrusions_matrix', 'var')
  error('You must pass an intrusions matrix.')
elseif ~exist('subjects', 'var')
  error('You must pass a subjects vector.')
elseif size(intrusions_matrix, 1) ~= length(subjects)
  error('intrusions matrix must have the same number of rows as subjects.')
end

if ~exist('rec_mask', 'var')
  % create blank mask if none was given
  rec_mask = make_blank_mask(intrusions_matrix);
elseif size(rec_mask) ~= size(intrusions_matrix)
  error('intrusions_matrix and rec_mask must have the same shape.')
end

% apply_by_index applies xli_for_subj on each subjects recall data separately
% As with other functions, xli_for_subj does most of the work here
xlis = apply_by_index(@xli_for_subj_avg, subjects, 1, {intrusions_matrix,rec_mask}); 

function xli_count = xli_for_subj_avg(intrusions_matrix,rec_mask)
    % add up the total number of -1's in a subject's intrusions matrix
    xli_count = sum(intrusions_matrix(rec_mask)==-1);