function pad = extendDimension(msk, val, THRESH)
%% extendDimension: extend dimension of mask by defined number of pixels
% The threshold parameter THERSH defines the minimum length of either dimension.
% If a dimension is lower than the threshold value, it is padded to match the
% threshold value with the value defined in the val parameter. This function
% automatically detects which dimension should be padded.
% THRESH = 300; % Original 200
chk    = size(msk) - THRESH;
idx    = chk < 0;

if idx
    pad = msk;
else
    pad = padarray(msk, abs(chk(idx)), val);
end

end