function [d] = hilbertNormalize(d)
    r = sum(d.*d,2).^-.5;
    d = bsxfun(@times,d,r);
end