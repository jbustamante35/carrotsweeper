function [d] = normalizeData(d)
    init = mean(d(1:3,:),1);
    d = bsxfun(@minus,d,init);
end