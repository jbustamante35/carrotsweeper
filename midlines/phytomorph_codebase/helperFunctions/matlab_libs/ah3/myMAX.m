function [bv] = myMAX(v,n)
    [j,sidx] = sort(v,'descend');
    bv = zeros(size(v));
    bv(sidx(1:n)) = 1;
end