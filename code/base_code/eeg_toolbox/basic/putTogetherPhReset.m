function putTogetherPhReset(dateX,band,test,startStr,endStr)


subj_array = {'087','088','092','101','103'};

%subj_array = {'087','092','101','103','096'};

%subj_array = {'096'};

sum_zcorr = [];sum_zcorr_foo = [];


for m = subj_array

    m = m{1};
    subj = sprintf('%s',m);
    side = 'L';
    
    allpFoo = []; psum = [];
    
          loadPath = fullfile('~/DBS_data/kareem_analysis/sigtests/RedoSN/phaseReset/',dateX,[side band test subj],[startStr '-' endStr 'phReset31.mat']);

        load(loadPath);
        
        allpFoo = norminv(allpFoo);
        psum = norminv(psum);
        
         if isempty(sum_zcorr)
               sum_zcorr = psum;
               sum_zcorr_foo = allpFoo;
               
         else
             sum_zcorr = sum_zcorr + psum;
             sum_zcorr_foo = sum_zcorr_foo + allpFoo;
             
         end

end
         
     right_subj_array = {'087','092','094','100','103'};

     %right_subj_array = {'087','092','094','100','103','096'};
     
     %right_subj_array = {'096'};
     
     for n = right_subj_array
    
    n = n{1};
    subj = sprintf('%s',n);
    side = 'R';
      
        rightLoadPath= fullfile('~/DBS_data/kareem_analysis/sigtests/RedoSN/phaseReset/',dateX,[side band test subj],[startStr '-' endStr 'phReset31.mat']);

        load(rightLoadPath);
        
        allpFoo = norminv(allpFoo);
        psum = norminv(psum);
        
         if isempty(sum_zcorr)
               sum_zcorr = psum;
               sum_zcorr_foo = allpFoo;
               
         else
             sum_zcorr = sum_zcorr + psum;
             sum_zcorr_foo = sum_zcorr_foo + allpFoo;
             
         end
        
         
     end
     
     for yy = 1:size(sum_zcorr,2)
            cor_boot = getThePValueOpp(sum_zcorr_foo(:,yy),sum_zcorr(yy),size(sum_zcorr_foo,1));
            
              all_cor(yy) = cor_boot;
            %all_pred(yy,jj) = pred_boot;
            
            %all_err(yy) = err_boot;
            
        end

        
            all_cor(all_cor < .000001) = 1/(size(sum_zcorr_foo,3)); all_cor(all_cor > .9999) = .9999;
%all_pred(all_pred < .000001) = 1/(size(sum_zcorr_foo,3)); all_pred(all_pred > .9999) = .9999;
            %all_err(all_err < .000001) = 1/(size(sum_zerr_foo,3)); all_err(all_err > .9999) = .9999;


load('~/valid_freqs.mat');

keyboard  
     
             
             