function [v] = myMax(T,para)
    [JUNK v] = max(T.T(T.fidx,para{1}),[],1);
    fidx = find(T.fidx);
    v = fidx(v);
end