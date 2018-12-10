function [l] = length(S)
    dS = diff(S,1,1);
    l = sum(dS.*dS,2).^.5;
    l = cumsum([0;l],1);
end