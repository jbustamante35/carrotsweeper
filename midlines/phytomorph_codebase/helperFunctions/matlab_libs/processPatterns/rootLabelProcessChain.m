function [H] = rootLabelProcessChain()
    %%%%%%%%%
    % root tip process chain
    p{1}.c(1).f = @cwtK;
    p{1}.c(1).p{1} = [17:23];
    p{1}.c(2).f = @watershedMax;
    p{1}.c(2).p{1} = 30;
    p{1}.i(1).f = @selectPattern;
    p{1}.i(1).p = [];
    p{1}.i(2).f = @endReturnPattern;
    p{1}.i(2).p = [];
    %%%%%%%%%
    % base process chain
    p{2}.c(1).f = @tensorEdge;
    p{2}.c(1).p{1} = 2;
    p{2}.c(2).f = @myMean;
    p{2}.c(2).p{1} = 1;
    p{2}.i(1).f = @selectPattern;
    p{2}.i(1).p = [];
    p{2}.i(2).f = @endReturnPattern;
    p{2}.i(2).p = [];
    %%%%%%%%%
    % lowerLeftcorner process chain
    p{3}.c(1).f = @tensorEdge;
    p{3}.c(1).p{1} = 2;
    p{3}.c(2).f = @myMax;
    p{3}.c(2).p{1} = 1;
    p{3}.i(1).f = @selectPattern;
    p{3}.i(1).p = [];
    p{3}.i(2).f = @endReturnPattern;
    p{4}.i(2).p = [];
    %%%%%%%%%
    % upperRightcorner process chain
    p{4}.c(1).f = @tensorEdge;
    p{4}.c(1).p{1} = 2;
    p{4}.c(2).f = @myMin;
    p{4}.c(2).p{1} = 1;
    p{4}.i(1).f = @selectPattern;
    p{4}.i(1).p = [];
    p{4}.i(2).f = @endReturnPattern;
    p{4}.i(2).p = [];
    %%%%%%%%%
    % parallel process chain
    F{1} = @(x)processChainPattern(p{1},x);
    F{2} = @(x)processChainPattern(p{2},x);
    F{3} = @(x)processChainPattern(p{3},x);
    F{4} = @(x)processChainPattern(p{4},x);
    
    H = @(x)ParallelProcessInteratorPattern(F,x);
    
end
