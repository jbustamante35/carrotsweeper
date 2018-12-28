function gmm = constructGMM(fin, numComponents, sampleSize, maxIter, vis)
%% constructGMM: 
% 
%
% Usage:
%   gmm = constructGMM(fin, numComponents, sampleSize, maxIter, vis)
%
% Input:
%   fin: cell array of file paths to images
%   numComponents: number of clusters to extract for GMM
%   sampleSize: number of pixels to sample per image (unclear?)
%   maxIter: number of iterations to fit model
%   vis: boolean to visualize clusters overlaid on sample image
%
% Output:
%   gmm: gaussian mixture model
%

%% Prep image data for constructing GMM 
% Generate n x 3 matrix of RGB values and creates a random permutation of all
% images for unbiased model 
TS = {};
parfor e = 1:numel(fin)
    img   = imread(fin{e});
    img   = imresize(img, 0.1);
    img   = permute(img,[3 1 2]);
    sz    = size(img);
    img   = reshape(img, [sz(1) prod(sz(2:3))])';
    img   = img(randperm(size(img, 1)), :);
    TS{e} = img(1:sampleSize, :);
end

%% Set-up model's components
MS  = zeros(sampleSize * numel(fin), 3);
str = 1;
for e = 1:numel(TS)
    stp            = str + sampleSize - 1;
    MS(str:stp, :) = TS{e};
    str            = stp + 1;
end

%% Construct GMM
opts = statset('Display', 'iter', 'MaxIter', maxIter);
gmm  = fitgmdist(MS, numComponents, 'Options', opts);

%% Apply GMM and view overlay of clusters
if vis
    % Overlay clusters onto first image
    imgCluster = imread(fin{1});
    imgOrg     = imgCluster;
    orgSz      = size(imgCluster);
    imgCluster = permute(imgCluster, [3 1 2]);
    imgSz      = size(imgCluster);
    imgCluster = reshape(imgCluster, [imgSz(1) prod(imgSz(2:3))])';
    imgCluster = double(imgCluster);
    kidx       = gmm.cluster(imgCluster);
    kidx       = reshape(kidx, orgSz(1:2));
    lRGB       = label2rgb(kidx);
    imshow(lRGB,[]);

    % look at overlay of clusters
    for e = 1 : numComponents
        out = flattenMaskOverlay(imgOrg, kidx == e);
        imshow(out, []);
        ttl = sprintf('Cluster %d', e);
        title(ttl);
        drawnow;
        waitforbuttonpress;
    end
end
end
