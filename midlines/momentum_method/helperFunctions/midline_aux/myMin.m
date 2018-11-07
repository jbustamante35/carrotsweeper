function [v] = myMin(T,para)
    [JUNK v] = min(T.T(T.fidx,para{1}),[],1);
    fidx = find(T.fidx);
    v = fidx(v);
end