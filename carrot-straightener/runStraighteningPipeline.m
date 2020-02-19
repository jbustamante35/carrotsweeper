function [pmsk, crv, mline, smsk, tcrd, dsts, nrms] = runStraighteningPipeline(img)
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
FACE   = 3;   % Direction to point straightened images (original 3)

%% Prepare mask for extraction functions
msk = extendDimension(img, 0, THRESH);
msk = double(imcomplement(msk));

%% Run processed mask through extraction functions
tic;
fprintf('\nGetting contour and midline...');
[pmsk, crv, mline, tcrd] = getContourAndMidline(msk, 0);
fprintf('Done...[%.02f sec]\n', toc);

try
    tic;
    fprintf('Running straightening...');
    % smsk = sampleStraighten(mline, pmsk);
    [smsk, nrms]  = getStraightenedMask(mline, pmsk);
    fprintf('Done...[%.02f sec] \n', toc);
    
    
    % Post-process data [Flip right-to-left and binarize]
    smsk = handleFLIP(smsk, FACE);
    dsts = sum(smsk);
catch 
    fprintf(2, '\nError running algorithm\n');
    smsk = [];
    nrms = [];
    dsts = [];
end

end

