function [xM] = generateXover(xSites,cLv,exeProb)
    TxM = {};
    for e = 1:numel(cLv)
        [TxM{e}] = generateXoverMax(xSites{e},cLv(e),exeProb{e},'sparse');
    end
    xM = blkdiag(TxM{:});
end