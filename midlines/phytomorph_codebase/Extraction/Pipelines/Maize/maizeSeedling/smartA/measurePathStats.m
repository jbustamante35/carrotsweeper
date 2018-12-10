function [measureBlocks,P] = measurePathStats(path,measureBlocks,toM)
    try
        %prob = NaN*ones(1,numel(measureBlocks));
        
        for m = 1:numel(measureBlocks)
            prob{m}(1) = 0;
            ptr = numel(path);
            basis = 8.^([1:ndims(measureBlocks(m).statBlock)]-1);
            N = 1;
            while (ptr - measureBlocks(m).history(end)) >= 1
                for e = 1:numel(measureBlocks(m).history)
                    ss(e) = path(ptr - measureBlocks(m).history(e));
                end
                idx = basis*(ss)'+1;
                measureBlocks(m).statBlock(idx) = measureBlocks(m).statBlock(idx) + 1;
                if ~isempty(toM)
                    if N == 1
                        prob{m}(N) = [];
                    end
                    prob{m} = [prob{m} toM(m).statBlock(idx)];
                    N = N + 1;
                end
                ptr = ptr - 1; 
            end
            if ~isempty(toM)
                tmp = sort(prob{m},'descend');
                P(m) = mean(tmp(1:10));
            end
        end
        
        %prob = sum(prob);
        P = sum(P);
    catch ME
        ME;
    end
end