function smsk = getStraightenedMask(mline, msk)
%% getStraightenedMask: straighten image from midline
% This is a modified version of the sampleStraighten function that straightens
% an object in an image by extending normal vectors along each coordinate of a
% midline. 
%
% 
%
% Usage:
%   smsk = getStraightenedMask(mline, msk)
%
% Input:
%   mline: midline coordinates
%   msk: binary mask image
%
% Output:
%   smsk: straightened mask
%

%% Create envelope structure
% Set unit length vector to place outer boundary
dscl = ceil(size(msk,1) / 2);
tng  = gradient(mline')';
d2e  = sum((tng .* tng), 2).^(-0.5);
ulng = bsxfun(@times, tng, d2e) * dscl;

% Compute distances from midline points to edge of envelope
eO = [-getDim(ulng, 2) , getDim(ulng, 1)] + mline;
eI = [getDim(ulng, 2) , -getDim(ulng, 1)] + mline;

% Create
[eOut, sOut] = generateFullEnvelope(mline, eO, dscl, 'cs');
[eInn, sInn] = generateFullEnvelope(mline, eI, dscl, 'cs');

%% Map curves to image
sz   = [size(sInn,1), length(mline)];

outI  = imbinarize(ba_interp2(double(msk), eInn(:,1), eInn(:,2)), ...
    'adaptive', 'Sensitivity', 1);
fullI = reshape(outI, sz);

outO  = imbinarize(ba_interp2(double(msk), eOut(:,1), eOut(:,2)), ...
    'adaptive', 'Sensitivity', 1);
fullO = reshape(outO, sz);

fullE = handleFLIP([flipud(fullI) ; fullO],3);

% Extract largest object and resize to specified dimensions
prp                           = regionprops(fullE, 'Area', 'PixelIdxList');
[~ , maxIdx]                  = max(cell2mat(arrayfun(@(x) x.Area, ...
    prp, 'UniformOutput', 0)));
smsk                           = zeros(size(fullE));
smsk(prp(maxIdx).PixelIdxList) = 1;

end



