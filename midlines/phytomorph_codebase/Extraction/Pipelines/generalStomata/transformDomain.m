function [D] = transformDomain(D,T)
    if size(D,2) ~= size(T,2)
        D = [D ones(size(D,1),1)];
    end
    D = mtimesx(T,D,'T')';
end