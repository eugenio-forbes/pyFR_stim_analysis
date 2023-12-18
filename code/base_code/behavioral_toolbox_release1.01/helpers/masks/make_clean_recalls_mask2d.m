function mask = make_clean_recalls_mask2d(recalls_matrix)
  %  mask = make_clean_recalls_mask2d(recalls_matrix)
  %  Makes a mask of the same shape as recalls_matrix which is false at 
  %  positions (i,j) if recalls_matrix(i,j) is an intrusion, 
  %  repeat, or empty cell.
  
  % sanity:
  if ~exist('recalls_matrix', 'var')
    error('You must pass a recalls matrix.')
  end
  
    mask = and(make_mask_exclude_repeats2d(recalls_matrix), ...
	     make_mask_exclude_intrusions2d(recalls_matrix));
  end

%endfunction
