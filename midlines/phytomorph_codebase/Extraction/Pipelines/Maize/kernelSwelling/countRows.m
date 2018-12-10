function [rowCount] = countRows(SZ,C,structSZ)
    M = zeros(SZ);
    C = round(C);
    for e = 1:size(C,1)
        M(C(e,1),C(e,2)) = 1;
    end
    M = imdilate(M,strel('disk',structSZ));
    M = any(M,2);
    R = regionprops(M);
    rowCount = numel(R);
end