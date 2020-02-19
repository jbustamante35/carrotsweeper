function [smsk, sdata] = getStraightenedMask(crds, img, BNZ, SCL)
%% getStraightenedMask: straighten image from midline
% This is a modified version of the sampleStraighten function that straightens
% an object in an image by extending normal vectors along each coordinate of a
% midline.
%
% Usage:
%   [smsk, sdata] = getStraightenedMask(crds, img, BNZ, SCL)
%
% Input:
%   crds: x-/y-coordinates of curve to map
%   img: image to interpolate pixels from coordinates
%   BNZ: boolean to binarize the final mask [for bw objects]
%   SCL: scaler to extend normal to desired distance [in pixels]
%
% Output:
%   smsk: straightened image
%   sdata: extra data for visualization or debugging
%

%% Create envelope structure
% Set unit length vector to place outer boundary
if nargin < 3
    % Default binarization on and envelope size to half the width of the image
    BNZ = 1;
    SCL = ceil(size(img,1) / 2);
end

tng  = gradient(crds')';
d2e  = sum((tng .* tng), 2).^(-0.5);
ulng = bsxfun(@times, tng, d2e) * SCL;

% Compute distances from midline points to edge of envelope
bndsOut = [-getDim(ulng, 2) , getDim(ulng, 1)] + crds;
bndsInn = [getDim(ulng, 2) , -getDim(ulng, 1)] + crds;

%% Map curves to image
[envO, datO] = map2img(double(img), crds, bndsOut, SCL, BNZ);
[envI, datI] = map2img(double(img), crds, bndsInn, SCL, BNZ);

if BNZ
    %% For CarrotSweeper straightener
    smsk = handleFLIP([flipud(envO) ; envI], 3);
    
    % Extract largest object from binarized image
    prp                            = regionprops(smsk, 'Area', 'PixelIdxList');
    [~ , maxIdx]                   = max(cell2mat(arrayfun(@(x) x.Area, ...
        prp, 'UniformOutput', 0)));
    smsk                           = zeros(size(smsk));
    smsk(prp(maxIdx).PixelIdxList) = 1;
    
else
    %% For HypoQuantyl S-Patches
    smsk = [fliplr(envI') , envO'];
    
end

% Extra data for visualization or debugging
sdata = struct('OuterData', datI, 'InnerData', datO);

end

function [env, edata] = map2img(img, crds, ebnds, dscl, bnz)
%% map2img: interpolate pixel intensities from curve coordinates
% Create the envelope structure
[eCrds, eGrid] = generateFullEnvelope(crds, ebnds, dscl, 'cs');

% Map curves to image
sz     = [size(eGrid,1), length(crds)];
mapimg = ba_interp2(double(img), eCrds(:,1), eCrds(:,2)); % Don't stretch vector
% mapimg = interp2(double(img), eCrds(:,1), eCrds(:,2));

% Binarize for CarrotSweeper
if bnz
    env = logical(reshape(mapimg, sz));
    %     env = imbinarize(env, 'adaptive', 'Sensitivity', 1);
    %     mapimg = imbinarize(mapimg, 'adaptive', 'Sensitivity', 1);
else
    env = reshape(mapimg, sz);
    
    % Replace out-of-frame pixels with median background intensity
    msk             = img > graythresh(img / 255) * 255;
    bk              = median(img(msk(:)));
    env(isnan(env)) = bk;
    
end

% Extra data for visualization or debugging
edata = struct('eCrds', eCrds, 'GridSize', size(eGrid));

end



