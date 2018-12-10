function [v] = myMean(T,para)
    v = mean(T.T(T.fidx,:),para{1});
    del = T.T - repmat(v,[size(T.T,1) 1]);
    [JUNK v] = min(sum(del.*del,2));
end