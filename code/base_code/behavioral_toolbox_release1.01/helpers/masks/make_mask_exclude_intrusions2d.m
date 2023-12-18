function mask = make_mask_exclude_intrusions2d(recalls_matrix)
  %  mask = make_mask_exclude_intrusions2d(recalls_matrix)
  %  Makes a mask of the same shape as recalls_matrix which is false at 
  %  positions (i,j) such that recalls_matrix(i, j) is an intrusion.
  
  % sanity:
  if ~exist('recalls_matrix', 'var')
    error('You must pass a recalls matrix.')
  end

  mask = make_blank_mask(recalls_matrix);

  mask(recalls_matrix < 1) = false;

%endfunction
