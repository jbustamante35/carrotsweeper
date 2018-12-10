function [CHp] = generateCHpara(chNRange,lengthRange)
    if nargin < 1
        chNRange = {1 5};
    end
    if nargin < 2
        lengthRange = {100 200};
    end
    debug = 0;
    %%%%%%%%%%%%%%%%%%%%%%
    % generate the number of CH
    cN = generateCHn('Uniform',chNRange,{1});
    if debug;cN = 2;end
    % generate the length of each chr
    [cLv] = generateCHlength('Uniform',lengthRange,cN);
    if debug;cLv = [100 100];end
    %%%%%%%%%%%%%%%%%%%%%%
    % init vars
    xSites = {};
    xVec = {};
    xP = {};
    exeProb = {};
    %%%%%%%%%%%%%%%%%%%%%%
    % generate  crossover matrix
    for e = 1:numel(cLv)
        %%%%%%%%%%%%%%%%%%%%%%
        % generate number of xOver sites
        xP{e} = generateNumberOfXoverSites('Uniform',{0 10});
        if debug;xP{e} = 2;end
        %%%%%%%%%%%%%%%%%%%%%%
        % generate the xover sites
        [xSites{e} xVec{e}] = generateXoverLocs('Uniform',{1 xP{e}},{1 cLv(e)});
        %%%%%%%%%%%%%%%%%%%%%%
        % generate the threshold for xovr to happen
        [exeProb{e}] = generateExecuteProb(numel(xSites{e}));
        if debug;exeProb{e} = zeros(size(exeProb{e}));end
    end
    CHp.cLv = cLv;
    CHp.cN = cN;
    CHp.xSites = xSites;
    CHp.xVec = xVec;
    CHp.exeProb = exeProb;
end