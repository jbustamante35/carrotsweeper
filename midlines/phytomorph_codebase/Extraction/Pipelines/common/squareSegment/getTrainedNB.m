function [func] = getTrainedNB(X,Y)
    nb = fitcnb(X',Y);
    func = @(X,e0,e1)predict(nb,X');
end