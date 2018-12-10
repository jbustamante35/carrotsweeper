function [V,D,alpha] = ker_eig(D,N,P)
    sz = size(D,2);
    D = {D};    
    alpha = size(D{1},1);
    K = ker1(D,P,alpha);
    [V,D] = eigs(K,N);
end
