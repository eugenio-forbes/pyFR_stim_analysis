function [subjnum_c,ROI1_c,ROI1_Bi_c,ROI2_c,ROI2_Bi_c] = common_subj( subjnum1,ROI1,ROI1_Bi,subjnum2,ROI2,ROI2_Bi)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

subjnum_c = intersect(subjnum1,subjnum2);

for n = 1:length(subjnum_c)
    subj_bi1 = strcmp(subjnum1,subjnum_c{n});
    ind1(n) = find(subj_bi1);

    subj_bi2 = strcmp(subjnum2,subjnum_c{n});
    ind2(n) = find(subj_bi2);
end

ROI1_c = ROI1(ind1);
ROI2_c = ROI2(ind2);

ROI1_Bi_c = ROI1_Bi(ind1); 
ROI2_Bi_c = ROI2_Bi(ind2);