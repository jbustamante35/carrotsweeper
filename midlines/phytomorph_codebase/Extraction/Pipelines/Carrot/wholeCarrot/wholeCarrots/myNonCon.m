function [C,Ceq] = myNonCon(X,N)
    Ceq = zeros(size(X,1),1);
    C = -(abs(X)-N);
end