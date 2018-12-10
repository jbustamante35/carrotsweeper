function [Cn] = crossInv(C1,C2,CHp,xM1,xM2)
    if nargin == 3
        % generate xOver
        [xM1] = generateXover(CHp.xSites,CHp.cLv,CHp.exeProb);
        [xM2] = generateXover(CHp.xSites,CHp.cLv,CHp.exeProb);
    end


    l1 = {};
    % generate two games
    [G1] = miosis(xM1,C1);
    [l1{1},l1{2}] = selectGam(G1);
    % select gam from two
    Gs1 = l1{round(rand(1))+1};

    l2 = {};
    % generate two games
    [G2] = miosis(xM2,C2);
    [l2{1},l2{2}] = selectGam(G2);
    % select gam from two
    Gs2 = l2{round(rand(1))+1};


    % fuse games
    Cn = fuseGam(Gs1,Gs2);
end