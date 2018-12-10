function [F] = extractStrain(velP,smoothP)
    vPs = imfilter(velP,fspecial('average',smoothP),'replicate');
    sr = gradient(vPs')';
    F = cat(3,vPs,sr);
end