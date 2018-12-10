function [M] = mfftm(L,n)
    k = 0:(L-1);
    W = exp(-1i*2*pi/L);
    for e = 1:n
        M(e,:) = W.^((k)*(e-1));
    end
end
