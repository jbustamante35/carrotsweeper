function [p] = gprobData(data,model)
    %%%%%%%%%%%%%%%%%%%
    % alpha for cut off
    alpha = [10 0 100];
    %%%%%%%%%%%%%%%%%%%
    % look up prob of data
    p = nprobData(data,model);
    %%%%%%%%%%%%%%%%%%%
    % soft threshold
    p = glf(p,alpha);
end