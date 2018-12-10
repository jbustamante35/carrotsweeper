function [Y] = ParallelProcessInteratorPattern(F,X)
    for f = 1:numel(F)
        Y{f} = F{f}(X);
    end
end