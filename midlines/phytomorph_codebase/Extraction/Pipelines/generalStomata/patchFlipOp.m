function [d,do] = patchFlipOp(d)
    do = cat(3,d,flipdim(d,1),flipdim(d,2),flipdim(flipdim(d,1),2));
    d = reshape(d,[size(d,1)*size(d,2) size(d,3)]);
end