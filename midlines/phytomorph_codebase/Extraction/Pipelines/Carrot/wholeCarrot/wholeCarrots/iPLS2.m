function [] = iPLS2(X1,Y1,X2,Y2,L1,L2,WIDTH,WIDTH2)
    
    X1_u = bsxfun(@times,X1,WIDTH.^-1);
    for e = 1:100
       

        [Xloadings,Yloadings,Xscores,Yscores,betaCOUNT,pctVar,mse,stats,Weights] = plsregress(X1_u,Y1,L1);
        COUNT = [ones(size(X1_u,1),1) X1_u]*betaCOUNT;
        CV(e,1) = corr(COUNT,Y1);
        
    
        X2_u = bsxfun(@times,X2,WIDTH2.^-1);
        COUNT_2 = [ones(size(X2_u,1),1) X2_u]*betaCOUNT;
        X2_u = bsxfun(@times,X2,COUNT_2.^-1);
        [Xloadings,Yloadings,Xscores,Yscores,betaWIDTH,pctVar,mse,stats,Weights] = plsregress(X2_u,Y2,L2);
        WIDTH2 = [ones(size(X2_u,1),1) X2_u]*betaWIDTH;
        
        
        
        CV(e,2) = corr(WIDTH2,Y2);
        
        
        X1_u = bsxfun(@times,X1,COUNT.^-1);
        WIDTH = [ones(size(X1_u,1),1) X1_u]*betaWIDTH;
        X1_u = bsxfun(@times,X1,WIDTH.^-1);
        plot(CV)
        drawnow
    end
    

end