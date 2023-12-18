function lag_crps = lag_crp(recalls_matrix, subjects, list_length, ...
			from_mask, to_mask, pres_mask)
%  LAG_CRP  Computes (lag) conditional response probabilities from a matrix
%       of recalled serial positions.
%
%  lag_crps = lag_crp(recalls_matrix, subjects, list_length, ...
%                 from_mask, to_mask, pres_mask)
%
%  INPUTS:
%  recalls_matrix:  a matrix whose elements are serial positions of recalled
%                   items.  The rows of this matrix should represent recalls
%                   made by a single subject on a single trial.
%
%        subjects:  a column vector which indexes the rows of recalls_matrix
%                   with a subject number (or other identifier).  That is, 
%                   the recall trials of subject S should be located in
%                   recalls_matrix(find(subjects==S), :)
%
%     list_length:  a scalar indicating the number of serial positions in the
%                   presented lists.  serial positions are assumed to run 
%                   from 1:list_length.
%
%       from_mask:  if given, a logical matrix of the same shape as 
%                   recalls_matrix, which is false at positions (i, j) where
%                   the transition FROM recalls_matrix(i, j) to
%                   recalls_matrix(i, j+1) should be excluded from
%                   the calculation of the CRP.  
%
%         to_mask:  if given, a logical matrix of the same shape as
%                   recalls_matrix, which is false at positions (i,
%                   j) where the transition from
%                   recalls_matrix(i, j-1) TO recalls_matrix(i, j)
%                   should be excluded from the calculation of the CRP.
% 
%                   If neither from_mask nor to_mask is given, a
%                   standard clean recalls mask is used, which
%                   excludes repeats, intrusions and empty cells.  If
%                   from_mask is given but to_mask is not, from_mask
%                   will be used for both masks (i.e., transitions
%                   both to and from masked out elements will be
%                   excluded).
% 
%       pres_mask:  if given, a logical matrix of the same shape as
%                   the presented items matrix, i.e., 
%                   (num_trials x list_length).  This mask should
%                   be true at positions (t, sp) where an item in
%                   the condition of interest was presented at
%                   serial position sp on trial t; and false
%                   everywhere else.  If NOT given, a blank mask of
%                   this shape is used.
%
%  OUTPUTS:
%        lag_crps:  a matrix of lag-CRP values.  Each row contains the values
%                   for one subject.  It has as many columns as there are
%                   possible transitions (i.e., the length of
%                   (-list_length + 1) : (list_length - 1) ).
%                   The center column, corresponding to the "transition of
%                   length 0," is guaranteed to be filled with NaNs.
%
%                   For example, if list_length == 4, a row in lag_crps
%                   has 7 columns, corresponding to the transitions from
%                   -3 to +3:
%                   lag-CRPs:     [ 0.1  0.2  0.3  NaN  0.3  0.1  0.0 ]
%                   transitions:    -3   -2    -1   0    +1   +2   +3
%

% sanity checks:
if ~exist('recalls_matrix', 'var')
  error('You must pass a recalls matrix.')
elseif ~exist('subjects', 'var')
  error('You must pass a subjects vector.')
elseif ~exist('list_length', 'var')
  error('You must pass a list length.') 
elseif size(recalls_matrix, 1) ~= length(subjects)
  error('recalls matrix must have the same number of rows as subjects.')
elseif ~exist('from_mask', 'var')
  % create standard clean recalls mask if none was given
  from_mask = make_clean_recalls_mask2d(recalls_matrix);
end

if ~exist('to_mask', 'var')
  % assume to_mask should be the same as from_mask (i.e.,
  % transitions to and from the same points should be excluded)
  to_mask = from_mask;
end

if size(from_mask) ~= size(recalls_matrix) | ...
   size(to_mask) ~= size(recalls_matrix)
  error('recalls_matrix and from and to masks must have the same shape.')
end

pres_mask_size = [size(recalls_matrix, 1), list_length];
if ~exist('pres_mask', 'var')
  % assume pres_mask should allow all presented items
  pres_mask = true(pres_mask_size);
elseif size(pres_mask) ~= pres_mask_size
  error(['pres_mask must have the same number of rows as' ...
	 ' recalls_matrix and list_length columns'])
end

% crp_for_subj does all the real work
lag_crps = apply_by_index(@crp_for_subj, ...
			  subjects, ...
			  1, ...
			  {recalls_matrix, from_mask, to_mask, pres_mask}, ...
			  list_length);
%endfunction

function subj_crp = crp_for_subj(recalls, from_mask, to_mask, ...
				 pres_mask, list_length)
  % Helper for crp:
  % calculates the lag-CRP for each possible transition for one subject's
  % recall trials; returns a row vector of CRPs, indexed by
  % (-list_length + index)
  
  % the golden goose: these will store actual and possible transitions for
  % this subject's trials
  actual_transitions = [];
  poss_transitions = [];

  % arguments for conditional_transitions:
  step = 1;   % for now, step of 1 is hard-coded; this may change   
  params =  struct('list_length', list_length);

  % conditional_transitions returns the actual and possible transitions for
  % each trial; we simply concatenate them all together and use collect() to
  % gather them up, since we don't care about distinctions between individual
  % trials within a subject
  num_trials = size(recalls, 1);
  for i = 1:num_trials
    % possible_transitions will exclude serial positions masked out by
    % params.to_mask_pres from trial_possibles
    params.to_mask_pres = pres_mask(i, :);
    [trial_actuals, trial_possibles] = conditional_transitions(...
					recalls(i, :), ...
    	                                from_mask(i, :), to_mask(i, :), ...
					@lag, @possible_transitions, ...
			                step, params);
    actual_transitions = [actual_transitions, trial_actuals];
    poss_transitions = [poss_transitions, catcell(trial_possibles)];
  end

  % the range of all possible transitions depends only on the list length;
  % we want to gather all the subject's actual and possible transitions in
  % this range
  all_possible_transitions = [-list_length + 1 : list_length - 1];
  actuals_counts = collect(actual_transitions, all_possible_transitions);
  possibles_counts = collect(poss_transitions, all_possible_transitions);

  % and that's it!
  subj_crp = actuals_counts ./ possibles_counts;

%endfunction

function d = lag(sp1, sp2, params)
  d = sp2 - sp1;
%endfunction