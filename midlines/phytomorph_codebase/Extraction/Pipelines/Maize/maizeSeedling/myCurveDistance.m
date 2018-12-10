function [d] = myCurveDistance(x,E,U,rawData)
    M = PCA_BKPROJ(x,E,U);
    idx = find(~isnan(rawData));
    d = nanmean((rawData - M).^2);
    %d = mean(log(normpdf(rawData(idx) - M(idx),0,1)));
    %s = nanstd(rawData-M);
    
    %plot(M);hold on;plot(rawData,'r')
    %drawnow
end