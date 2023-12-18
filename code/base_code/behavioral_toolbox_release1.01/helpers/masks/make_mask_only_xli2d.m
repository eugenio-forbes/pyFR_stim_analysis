function mask = make_mask_only_xli2d(intrusion_matrix)
%  mask = make_mask_only_xli2d(intrusion_row)
%  Makes a mask of the same shape as intrusion_matrix which is true at 
%  positions that are extra list intrusions
%  

% sanity:
if ~exist('intrusion_matrix', 'var')
  error('You must pass an intrusions matrix.')
end

mask = ~make_blank_mask(intrusion_matrix);
mask(intrusion_matrix==-1) = true;

%endfunction