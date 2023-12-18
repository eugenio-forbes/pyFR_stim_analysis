function mask = make_mask_only_pli2d(intrusion_matrix)
%  mask = make_mask_only_xli2d(intrusion_row)
%  Makes a mask of the same shape as intrusion_matrix which is true at 
%  positions that are prior list intrusions
%  

% sanity:
if ~exist('intrusion_matrix', 'var')
  error('You must pass an intrusions matrix.')
elseif ndims(intrusion_matrix) ~= 2
  error('intrusion_matrix must be two-dimensional.')
end

mask = ~make_blank_mask(intrusion_matrix);
mask(intrusion_matrix>0) = true;

%endfunction