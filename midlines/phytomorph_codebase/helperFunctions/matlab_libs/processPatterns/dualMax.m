function [fidx] = dualMax(sig,para)
    %%% find point on the upper scale.
    %%% isolate region around it [2*PAD]
    %%% find max within pad
    sig = abs(sig);
    [JUNK fidx] = max(sig(2,:));
    sig(:,1:(fidx-para{1})) = 0;
    sig(:,(fidx+para{1}):end) = 0;
    [JUNK fidx] = max(sig(1,:));
end