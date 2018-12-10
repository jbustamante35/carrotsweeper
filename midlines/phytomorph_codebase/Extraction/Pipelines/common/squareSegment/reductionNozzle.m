function [C] = reductionNozzle(D,U,E)
    D = bsxfun(@minus,D,U);
    C = E'*D;
end