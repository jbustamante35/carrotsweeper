function [ct] = nCm(n,bt)
    % obtain location values    
    tlv = clv(bt);        
    % compute digits
    ct = cd(n,tlv);
end
%%%%%%%%%%%%%%%%%%
% compute location values
function [lv] = clv(b)
    lv = cumprod([1 b(1:end-1)]);
    lv = fliplr(lv);
end
%%%%%%%%%%%%%%%%%%
% compute digits
function [c] = cd(n,lv)
    for i = 1:numel(lv)
        c(i) = floor(n/lv(i));
        n = n - lv(i)*c(i);
    end    
end
%{
    ns = [2 7];
    bs = [10 10];
    bt = [3 4 6 8];
    nt = nC(ns,bs,bt);
%}