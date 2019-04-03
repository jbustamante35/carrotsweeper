function [vecS, vecM] = sampleStraighten(mline, msk, img)
%% sampleStraighten: midline-based straightener
% This function straightens an object in an image by extending from the midline
% coordinates and extracting the pixel values each coordinate corresponds to.
%
% Usage:
%   [vecS, vecM] = sampleStraighten(mline, msk, img)
%
% Input:
%   mline: midline coordinates
%   msk: binary mask image
%   img: grayscale or rgb image corresponding to binary mask image
%
% Output:
%   vecS: pixel values from mask image
%   vecM: pixel values from grayscale or original image
%

%%
try
    [dS, dG] = extendCarrotMidline(mline, [0 0], msk);
    dSize    = size(dG);
    chnls    = size(img, 3);
    vecS     = zeros(length(dS), chnls);
    
    %%
    for chnl = 1 : chnls
        vecS(:,chnl) = ba_interp2(double(img(:,:,chnl)) / 255, dS(:,2), dS(:,1));
    end
    
    vecM = ba_interp2(double(msk) / 255, dS(:,2), dS(:,1));
    
    %%
    vecS = reshape(vecS, dSize(1:2));
    vecM = reshape(vecM, dSize(1:2));
    
catch e
    fprintf(2, 'Error straightening mask\n%s\n', e.getReport);
    vecS = [];
    vecM = [];
    
end

end