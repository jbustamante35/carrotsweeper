function [pmsk, crv, mline, smsk, tcrd, dsts] = runStraighteningPipeline(img)
%% runStraighteningPipeline: core pipeline for straightening algorithm
% This is the core function that runs the straightening pipeline. It takes a
% binarized image, processes it to face left-to-right, then:
%   1) Extracts the contour from the binarized mask
%   2) Identifies the tip - the highest curvature coordinate
%   3) Generates the midline using Nathan's walking algorithm from the tip
%   4) Straightens the mask by sampling the image from the interpolated midline
%
% Usage:
%   [pmsk, crv, mline, smsk, tcrd, dsts] = runStraighteningPipeline(img, vis)
%
% Input:
%   img: binarized mask of the image
%   vis: boolean to show results
%
% Output:
%   pmsk: processed left-right-facing mask
%   crv: x-/y-coordinates of the mask's contour
%   mline: x-/y-coordinates of the mask's midline
%   smsk: straightened mask
%   tcrd: coordinate of carrot tip
%   dsts: distance transform of midline
%

%% Some constants to consider playing around with
THRESH = 300; % Minimum length to pad one or both dimensions of image
FACE   = 2;   % Direction to point straightened images (original 3)

%% Prepare mask for extraction functions
msk = extendDimension(img, 0, THRESH);
msk = double(imcomplement(msk));

%% Run processed mask through extraction functions
[pmsk, crv, mline, tcrd, dsts] = getContourAndMidline(msk, 0);
smsk                           = sampleStraighten(mline, flip(pmsk, 3), pmsk);

% Post-process data [Flip right-to-left and binarize]
% pmsk  = handleFLIP(pmsk, FACE);
% crv   = fliplr(crv);
% mline = fliplr(mline);
smsk  = handleFLIP(imbinarize(smsk), FACE);

end