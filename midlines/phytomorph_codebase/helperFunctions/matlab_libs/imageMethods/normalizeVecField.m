function [vec] = normalizeVecField(vec)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % euclidean normalization of vector
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vec = vec.*repmat(sum(vec.*vec,3).^-.5,[1 1 size(vec,3)]);
end