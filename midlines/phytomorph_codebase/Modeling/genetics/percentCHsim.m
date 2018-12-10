function [p sites] = percentCHsim(C1,C2)
    nonEq = C1.ch == C2.ch;
    sites = find(nonEq);
    p = sum(nonEq);
    p = p * numel(C1.ch).^-1;
end