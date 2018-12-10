function [S C U BV L sL ERR] = PCA_FIT_FULL_wW(M,COM)
U = mean(M,1);
for i = 1:size(M,1)
    M(i,:) = M(i,:) - U;
end
Std = std(M,1,1);
for i = 1:size(M,1)
    M(i,:) = M(i,:).*Std.^-1;
end
COV = cov(M);
[BV L] = eigs(COV,COM);
sL = L;
L = cumsum(diag(L))*sum(diag(L))^-1;
C = M*BV;
% create the simulated signal
S = (BV*C')';
for i = 1:size(S,1)
    S(i,:) = S(i,:) + U;
end
ERR = sum((S - M).^2,2).^.5;

