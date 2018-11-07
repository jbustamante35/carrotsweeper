function [idx] = snapTo(X,P)
    idx = [];
    delta = bsxfun(@minus,X,P);
    delta = sum(delta.*delta,2);
    [~,idx] = min(delta);
end