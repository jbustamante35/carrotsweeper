function [DEPTH SZ] = setinfo(SET,sep)
    fidx = strfind(SET{1},sep);
    DEPTH = size(fidx,2);
    SZ = numel(SET);
end
