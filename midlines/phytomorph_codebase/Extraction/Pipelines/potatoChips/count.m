function [fidx sel uF] = count(feature)
    sr = [];
    for k = 1:numel(feature)
        RAT = round(feature.*feature(k).^-1);
        sr(k) = sum(RAT==1);    
    end

    MO = mode(sr);
    fidx = (sr == MO);
    uF = mean(feature(fidx));
    sel = round(feature.*uF^-1);    
    fidx = (sel==1);
end