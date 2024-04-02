function [zVal pVal] = ranksum_shuff(recVect, nonVect, shuffles)



		z_true = [];
        ind1 = recVect; ind2 = nonVect;
		
		
        [p1 h1 stats1] = ranksum_ci(ind1,ind2,.05,1);
        
        pv = stats1.ranksum; 
		
        z_true = p1;
        
        ind1_size = size(ind1,1);
        ind2_size = size(ind2,1);
        
        rn = [];
        rn = vertcat(ind1,ind2);
        z_foo = [];
            for r = 1:shuffles;
                
                randind = randperm(length(rn));
                
                ind1_foo_data = rn(randind(1:ind1_size));
                
                ind2_foo_data = rn(randind(ind1_size+1:end));
                
                                
                [p_foo h_foo stats1_foo] = ranksum_ci(ind1_foo_data, ind2_foo_data, .05,1);
                
                z_foo(r) = p_foo;
                
            end
        
        pVal = getThePValueOpp(z_foo,z_true,r);
      


pVal(pVal<.00001) = .00001;
pVal(pVal>.99999) = .99999;

zVal = norminv(pVal);
       
       
       
