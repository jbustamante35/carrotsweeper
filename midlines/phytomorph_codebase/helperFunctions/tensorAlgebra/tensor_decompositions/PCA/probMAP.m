function [C BV L LAM] = probMAP(M,COM)
[BV L] = eigs(M,COM);
LAM = L;
L = cumsum(diag(L))*sum(diag(L))^-1;
C = M*BV;


