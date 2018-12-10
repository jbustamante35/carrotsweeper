function [o] = buildAffine(T,P)
    o = [reshape(T,[numel(T)^.5 numel(T)^.5]) P'];
end