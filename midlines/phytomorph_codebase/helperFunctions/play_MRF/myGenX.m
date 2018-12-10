function [S] = myGenX(nB,e,nb,n)
    N = n^2;
    %nB = floor(nb/64);
    %X = [];
    X = ones(1,nB*64);
    %ex = rem(e,64);
    ex = e;
    nbex = min(64,N-nB*64);
    Y = zeros(1,N - nB*64-nbex);
    X = [X bitget(uint64(ex),1:nbex) Y];
    S = double(reshape(X,[n n]));
end