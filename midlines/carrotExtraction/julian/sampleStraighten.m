function [vec, vecM] = sampleStraighten(mline, msk, img)
%% sampleStraighten: midline-based straightener
% This function straightens an object in an image by extending from the midline
% coordinates and extracting the pixel values each coordinate corresponds to.
%
% Usage:
%   [vec, vecM] = sampleStraighten(midline, carrotMask, carrotImage)
%
% Input:
%   mline: midline coordinates
%   msk: binary mask image
%   img: grayscale or rgb image corresponding to binary mask image
%
% Output:
%   vec: pixel values from raw image
%   vecM: raw coordinates used for vec 
% 

%%
[dS, dG] = extendCarrotMidline(mline, [0 0], msk);
dsz      = size(dG);
imgSize  = size(img, 3);
vec      = [];

%% 
for k = 1 : imgSize
    vec(:, k) = ba_interp2(double(img(:,:,k)) / 255, dS(:,2), dS(:,1));
end

%% 
vecM = ba_interp2(double(msk) / 255, dS(:,2), dS(:,1));
vec  = reshape(vec,  [dsz(1) dsz(2) 3]);
vecM = reshape(vecM, [dsz(1) dsz(2)]);

end