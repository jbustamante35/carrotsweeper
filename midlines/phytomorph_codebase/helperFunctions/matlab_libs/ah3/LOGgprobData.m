function [p] = LOGgprobData(data,model)
    
    %%%%%%%%%%%%%%%%%%%
    % look up prob of data
    p = LOGprobData(data,model);
    
    %%%%%%%%%%%%%%%%%%%
    % soft threshold
    p = normpdf(p,0,1);
    
    %%%%%%%%%%%%%%%%%%%
    % normalize
    MAX = normpdf(0,0,1);
    p = p / MAX;
end