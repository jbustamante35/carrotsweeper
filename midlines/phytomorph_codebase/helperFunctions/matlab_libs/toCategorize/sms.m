function [fidx] = dualMax(sig,PAD)
    %%% find point on the upper scale.
    %%% isolate region around it [2PAD]
    %%% find max within pad
    [JUNK fidx] = max(sig(1,:));
    sig(2,1:(fidx-PAD)) = 0;
    sig(2,(fidx+PAD):end) = 0;
    [JUNK fidx] = max(sig(2,:));
end