function [wFunc] = getPLSDR_func(X,Y,d)
    [Xloadings,Yloadings,Xscores,Yscores, beta,pctVar,mse,stats,Weights] = plsregress(X',Y',d);
    U = mean(X,2);
    wFunc = @(X,e0,e1)(((Weights')*bsxfun(@minus,X,U)));
end