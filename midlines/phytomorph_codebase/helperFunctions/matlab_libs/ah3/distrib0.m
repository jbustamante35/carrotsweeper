function [Y] = distrib0(D,P)
    Y = poisspdf(D,P(1));
    for e = 1:(numel(P)-1)/2
        Y = normpdf(P((e-1)*2+1:(e-1)*2+2));
    end
end