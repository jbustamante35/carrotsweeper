function [vec] = normalizeVec(vec)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % euclidean normalization of vector
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vec = vec * norm(vec)^-1;
end
