function [cnt] = particleCounter(particleSpray,cnt)
    n = ndims(particleSpray);
    if isempty(cnt)
        cnt = 0;
    end
    cnt = cnt + size(particleSpray,n);
end