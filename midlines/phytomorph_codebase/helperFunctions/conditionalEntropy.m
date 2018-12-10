function [E] = conditionalEntropy(X,Y)
    %uqX = unique(X);
    uqX = [0 1];
    for u = 1:numel(uqX)
        tmpBV = X==uqX(u);
        N(u) = sum(tmpBV);
        E(u) = shannonEntopy(Y(tmpBV));
    end
    P = N * size(X,1)^-1;
    E = P*E';
end