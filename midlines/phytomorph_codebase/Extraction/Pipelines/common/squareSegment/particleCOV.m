function [C] = particleCOV(particleSpray,M,C)
    if isempty(C)
        C = zeros(size(particleSpray,1));
    end
    tmp = bsxfun(@minus,particleSpray,M{2}*M{1}^-1);
    C = C + (tmp*tmp')*M{1}^-1;
end