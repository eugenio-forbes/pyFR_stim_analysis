for f = 1:45
    for t= 1:900
        [pval,~,stat] = ranksum_ci(zMatrixRecAll(:,f,t),zMatrixNonAll(:,f,t),.05,-1);
        pmat(f,t) = pval;
        zmat(f,t) = stat.zval;
    end
end


