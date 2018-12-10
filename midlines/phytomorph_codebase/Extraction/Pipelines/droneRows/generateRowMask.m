function [mask] = generateRowMask(subI,f,phase)
    mask = zeros(size(subI));
    for r = 1:size(subI,1)
        mask(r,:) = cos(2*pi*(1:size(subI,2))*f + phase(r));
    end
    
end