function [fidx] = count2(feature,mag)
    sr = [];
    for k = 1:numel(feature)
        RAT = round(mag*feature.*feature(k).^-1);
        sr(k) = sum(RAT==1);    
    end
    MO = mode(sr);
    fidx = (sr == MO);
    uF = mean(feature(fidx));
    sel = round(mag*feature.*uF^-1);    
    %fidx = (sel==1);
    fidx = sel;
end