function blank_mask = make_blank_mask(data_matrix)
 % blank_mask = make_blank_mask(data_matrix)
 % Returns an all-true (blank) mask of the same shape as data_matrix

 % a blank mask is just one that is true everywhere
 blank_mask = true(size(data_matrix));

%endfunction
