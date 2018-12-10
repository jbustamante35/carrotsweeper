function [newV] = countParticleBlobs(spray,v)
    if isempty(v)
        v = -inf;
    end
    R = regionprops(logical(spray));
    newV = max(v,numel(R));
end