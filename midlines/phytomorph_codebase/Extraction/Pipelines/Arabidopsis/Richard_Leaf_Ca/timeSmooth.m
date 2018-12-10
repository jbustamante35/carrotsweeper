function [sig] = timeSmooth(sig,N)
    % reshape and smooth each pixel value over time
    sz = size(sig);
    sig = reshape(sig,[sz(1)*sz(2) sz(3)]);
    sig = imfilter(sig,fspecial('average',[1 N]),'replicate');
    sig = reshape(sig,sz);
end