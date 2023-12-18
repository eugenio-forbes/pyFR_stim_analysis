function mask = make_mask_exclude_repeats2d(recalls_itemnos)
%  mask = make_clean_recalls_mask2d(recalls_itemnos)
%  Makes a mask of the same shape as recalls_itemnos which is false at
%  positions (i,j) such that recalls_itemnos(i, j) is a repeat.

% sanity:
if ~exist('recalls_itemnos', 'var')
    error('You must pass a recalls matrix.')
end

rows = size(recalls_itemnos, 1);
cols = size(recalls_itemnos, 2);
mask = make_blank_mask(recalls_itemnos);

for i = 1:rows
    recalls_row = recalls_itemnos(i, :);
    prev_recalls = [];
    for j = 1:cols
        recall = recalls_row(j);
        if find(prev_recalls == recall)
            % mask out this position if previously recalled
            mask(i,j) = false;
        else
            % leave the mask blank (true) at this position, but update prev_recalls
            prev_recalls = [prev_recalls, recall];
        end
    end
end

%endfunction
