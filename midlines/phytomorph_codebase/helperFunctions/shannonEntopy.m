function [E] = shannonEntopy(X)
    %uqX = unique(X);
    uqX = [0 1];
    for u = 1:numel(uqX)
        N(u) = sum(X==uqX(u));
    end
    P = N * size(X,1)^-1;
    F = log(P);
    F(isinf(F)) = 0;
    E = -sum(P.*F);
end