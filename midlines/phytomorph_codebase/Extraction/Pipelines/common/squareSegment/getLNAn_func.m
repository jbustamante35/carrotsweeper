function [wFunc] = getLNAn_func(X,Y,d)
    Weights = mynLDA(X',Y',1,d);
    U = mean(X,2);
    wFunc = @(X,e0,e1)(((Weights')*bsxfun(@minus,X,U)));
end