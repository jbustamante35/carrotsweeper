function [ smoothed ] = smoothTassel( tassel, sigma, kernelDim, level)
%SMOOTHTASSEL Run a gaussian kernel over a binary tassel image and rethreshold
%the smoothed image.
%   tassel: binary image of tassel
%   sigma: positive value for standard deviation of gaussian kernel
%   kernelDim: width to extend the kernel to either side of 0.  Note this
%       is approx half the width of the final kernel, which is 2*kernelDim
%       +  1
%   level: used to determine threshold for re-binarizing the smoothed image

% Define kernel
K=inline('exp(-(x.^2+y.^2)/2/sig^2)');
[dx,dy]=meshgrid([-kernelDim:kernelDim]);
kern = K(sigma, dx, dy)/sum(sum(K(sigma, dx, dy)));

smoothed = conv2(double(tassel), kern, 'same');

level = graythresh(smoothed);
%smoothed = smoothed > level;
smoothed = smoothed > 0.5*level;

L = bwlabel(smoothed);
R = regionprops(logical(smoothed));
    
keep = find([R.Area] == max([R.Area]));
smoothed(L ~= keep) = 0;


end

