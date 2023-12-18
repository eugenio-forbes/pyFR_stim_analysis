function [trans_row, cond_cell] = conditional_transitions(...
    data_row, from_mask, to_mask, transit_func, condition, step, params)
% [trans_row, cond_cell] = conditional_transitions(...
%     data_row, from_mask, to_mask, transit_func, condition, step, params)
%
%  Returns transitions in data_row and a cell array of conditions.
%
%  This function behaves similarly to transitions(), with two
%  important differences: first, that the inclusion in the output of a
%  particular transition between two elements is conditional upon that
%  transition being in a vector returned by a condition function; and
%  second, that two values (one of actual transitions, one of
%  conditions) are returned.  Effectively, this function allows you to
%  dynamically determine which transitions to count, as opposed to the
%  static approach encouraged by transitions() (i.e., using a mask to
%  exclude certain positions).
%
%  condition should be a function handle which returns an array
%  representing possible transitions.  It should accept the current
%  element of data_row as a first argument, an array of the previously
%  encountered elements of data_row as a second argument, the
%  transition (as calculated by transit_func) from the current element
%  as a third argument, and any additional arguments passed to
%  condtional_transitions via params.
%
%  cond_cell is a cell array of the condition arrays returned by
%  the condition function.  The transition in trans_row(i) had the
%  condition in cond_cell{i}.
% 
%  For example, suppose you wanted to only consider transitions at
%  small lags:
%
%  >> data_row = [22, 19, 23, 4, 5, 6]
%  >> function cond_row = exclude_big_lags(element, priors, transit, params)
%  ..    % always allows lags in the -3..3 range
%  ..    cond_row =  [-3:3]
%  .. end
%  >> transitions_row, conditions_cell = conditional_transitions(
%  ..    data_row, logical(data_row), logical(data_row), 
%  ..    @lag, @exclude_big_lags, 1, struct())
%  >> transitions_row
%  [-3,  NaN, NaN, 1,  1]
%  >> conditions_cell 
%  {
%    [-3, -2, -1, 0, 1, 2, 3],
%    ...
%  }
%
%  A particular transition between two items, as computed by
%  transit_func, will be included in the output array of transitions
%  IF AND ONLY IF its value is found in the array returned by the
%  condition function.  Otherwise, a NaN will be placed transitions
%  array.  In the example above, the transitions between 22 and 19, 4
%  and 5, and 5 and 6 are included in the output because they each
%  fell in the range -3..3.  The other lags were excluded because they
%  fell outside this range; they therefore are represented by NaNs.
%
%  Note that, if no transitions are possible from a given element of
%  data_row, condition should return an empty vector, so that the
%  expression:
%
%  >> any(ismember(condition_row,transition))
%
%  evaluates to false, instead of throwing an error.
%
%  Masks take precedence over condition calculations.  If either
%  from_mask or to_mask excludes an element in data_row, no condition
%  calculation is performed for the transition to or from that
%  element:
%
%  >> kill_from1 = [false, true, true, true, true, true]
%  >> kill_to5 = [true, true, true, false, true, true]
%  >> transitions_row, conditions_cell = conditional_transitions(
%  ..    data_row, kill_from1, kill_to5, @distance, 
%  ..    @exclude_big_lags, 1, struct())
%  >> transitions_row
%  [NaN, 4, -19, NaN, 1]
%
%  All other arguments to conditional_transitions are the same as
%  would be expected by the transitions() function; consult its
%  docstring for more information.

% sanity checks
if ~islogical(from_mask) || ~islogical(to_mask)
  error('Masks must be logical variables');
end

if step < 1
  error('Non-positive steps are not supported');
end

if ~isa(condition,'function_handle')
  error('condition must be a function handle');
end

% cond_helper does all the work
[trans_row, cond_cell] = iterative_cond_helper(data_row, ...
				    from_mask, to_mask, step, transit_func, ...
				    condition, params);
%endfunction

function [trans_row, cond_cell] = iterative_cond_helper(data_row, ...
				     from_mask, to_mask, step, transit_func, ...
				     cond_func, params)

  % initialization: trans_row should have as many elements as
  % data_row - step, as should cond_cell, since that is the maximum
  % number of transitions and conditions we can calculate
  row_length = size(data_row, 2);  
  trans_row = NaN(1, row_length - step);
  cond_cell = cell(1, row_length - step);
  priors = [];

  for i = 1:row_length
    if i + step > row_length
      % return when we fall off the end of the list, i.e., the element at
      % i + step does not exist. 
      return
    end

    from_pt = data_row(i);
    to_pt = data_row(i + step);   
    from_pt_mask = from_mask(i);
    to_pt_mask = to_mask(i + step);
    if ~from_pt_mask || ~to_pt_mask
      % if either the from point or the to point is masked out, continue
      % without updating trans_row or cond_cell
      priors = [priors, from_pt];

      % trans_row(i) will remain NaN and cond_cell{i} will remain empty
      continue
    end

    % otherwise, calculate the current transition and condition and append them
    % to trans_row and cond_cell 
    transition = transit_func(from_pt, to_pt, params);
    condition = cond_func(from_pt, priors, transition, params);

    if any(ismember(condition, transition))
      % the transition meets the condition; include it in trans_row
      trans_row(i) = transition;
    end
    % whether or not the transition meets the condition, we update
    % cond_cell and priors
    cond_cell{i} = condition;
    priors = [priors, from_pt];
  end
%endfunction