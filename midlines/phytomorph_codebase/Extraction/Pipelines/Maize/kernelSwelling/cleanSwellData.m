function [d jd] = cleanSwellData(d,t)
    dd = diff(d,1,1);
    rmidx = find(any(abs(dd) > t,1) | all(d==0,1) | any(isnan(d),1) | any(isinf(d),1) | sum(d<0,1) > 10);
    jd = d(:,rmidx);
    d(:,rmidx) = [];
end