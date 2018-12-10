function [sig var] = subBL(sig,fFrames)
    p = polyfit((1:fFrames)',sig(1:fFrames),1);
    [brob,stats] = robustfit((1:fFrames)',sig(1:fFrames));
    var = stats.robust_s;
    var = std(sig(1:fFrames));
    baseLine = polyval(fliplr(brob'),1:numel(sig));
    sig = sig - baseLine';
end