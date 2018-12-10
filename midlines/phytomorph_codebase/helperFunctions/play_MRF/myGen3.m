function [S] = myGen3(e)
    n = 3;
    N = n^2;
    
    nB = floor(e/64);
    
    X = [];
    
    for e = 0:(nB - 1)
        X = [X ones(1,64)];
    end
    
    ex = rem(e,64);
    
    X = [X bitget(uint64(ex-1),1:N)];
    
    S = double(reshape(,[n n]));
end