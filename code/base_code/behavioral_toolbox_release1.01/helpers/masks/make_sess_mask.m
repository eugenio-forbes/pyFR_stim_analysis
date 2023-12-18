function sess_mask = make_sess_mask(recalls, session_index, sessions)
% MAKE_SESS_MASK  Creates a mask which excludes trials by session
% (or other index)
%
% function sess_mask = make_sess_mask(recalls, session_index, sessions);
%
% INPUTS:
%   recalls:	
%     a recalls matrix, or any matrix in which rows represent events
%     by a single subject on a single trial.
%  session_index:
%     a single column of the same height as the recalls
%     matrix. Each element is the session label of the
%     corresponding row in the recalls matrix
%  sessions:
%     a vector of session labels that should be included
%
%
% OUTPUTS:
%  sess_mask:	
%     a mask of the same size as recalls which is true in
%     rows i where session_index(i) is a member of sessions, and
%     false everywhere else
%
%
% EXAMPLE:
% >> recalls = [9 8 4 2 0; ...
%               8 9 7 3 1; ...
%               8 7 4 2 3; ...
%               7 8 9 6 1];
% >> session_index = [1; 2; 3; 3];
% >> sessions = [1 2]; % want to include only sessions 1 and 2
% >> sess_mask = make_sess_mask(recalls, session_index, sessions)
% >> sess_mask =
%
%     1     1     1     1     1
%     1     1     1     1     1
%     0     0     0     0     0
%     0     0     0     0     0
%          

% sanity checks:
if ~exist('sessions', 'var')
  error('You must pass a vector of sessions to include')
elseif ~exist('session_index', 'var')
  error('You must pass a session index vector')
elseif ~exist('recalls', 'var')
  error('You must pass a recalls matrix or similar')
elseif size(session_index, 1) ~= size(recalls, 1)
  error('session_index must have the same number of rows as recalls')
end

sess_mask = false(size(recalls));
rows_to_unmask = ismember(session_index, sessions);
sess_mask(rows_to_unmask, :) = true;

%endfunction
