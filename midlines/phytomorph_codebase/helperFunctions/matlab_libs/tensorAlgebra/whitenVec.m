function [d,u,s] = whitenVec(d)
    u = mean(d,1);
    s = std(d,1,1);
    s = s ^ -1;
    d = bsxfun(@minus,d,u);
    d = bsxfun(@times,d,s);
end