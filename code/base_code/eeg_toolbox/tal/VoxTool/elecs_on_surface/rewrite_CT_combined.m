function rewrite_CT_combined( CT_combined, CT_combined_rewrite, MRI_combined,...
    CT2MRI)
%REWRITE_CT_COMBINED(CT_combined, CT_combined_rewrite, MRI_combined, CT2MRImat)
% writes the CT_combined file from a .img to .nii.gz (or whatever is
% specified). CT_combined_rewrite should have a different basename than
% CT_combined or freesurfer may run into problems.
%
% Also "flirt"s the file to MRI_combined, if it is passed in
%   CT_combined - location of the original CT_combined file
%   CT_combined_rewrite - location of the new CT_combined. Must end with a
%   .nii.gz
%   MRI_combined (optional) - location of the MRI_combined file to be used
%                             in flirt
%   CT2MRI - location where the CT2MRI.mat and CT2MRI.nii.gz files
%               should be written (without the extension)

disp('reading...')
ct_struct = MRIread(CT_combined);

disp('writing...')
MRIwrite(ct_struct, CT_combined_rewrite);
disp('rewriting')
%system(sprintf('mri_convert %s %s',CT_combined, CT_combined_rewrite));

disp('flirting...')
system(sprintf('flirt -in %s -ref %s -out %s.nii.gz -omat %s.mat',CT_combined_rewrite, MRI_combined, CT2MRI, CT2MRI))
%system(sprintf('bash run_flirt.sh %s %s %s',CT_combined_rewrite, MRI_combined, CT2MRI))
beep
end

