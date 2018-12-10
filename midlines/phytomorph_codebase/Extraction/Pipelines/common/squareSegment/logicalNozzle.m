function [res] = logicalNozzle(X,func,dim,toT)
    % this nozzle will make the inputs logical
    if toT
        X = permute(X,[2 1 3]);
    end
    
    res = func(X,dim);
end