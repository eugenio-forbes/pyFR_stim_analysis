function arr = catcell(c)
% CATCELL(C) Concatenates the row vectors in a cell array into a single
%            vector
%
% Internal function used by conditional_transitions; the explanations here
% are relatively simplistic since this function will probably only get used
% internally by other functions rather than anything that needs to get
% called directly by the user.
%
% INPUTS:
%     c:  a cell array containing only row vectors or scalars
% 
% OUTPUTS:
%   arr:  a row vector containing the values found in c

if ~exist('c','var')
   error('Must pass in a cell array, c, to catcell.') 
end
if ~iscell(c)
    error('c must be a cell array.')
end

  arr = [c{:}];
%endfunction
