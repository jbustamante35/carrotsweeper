function [M] = fillSmallHoles(M,T)
    iM = ~M;
    fM = bwareaopen(iM,T);
    fM = fM ~= iM;
    M = fM + M;
end