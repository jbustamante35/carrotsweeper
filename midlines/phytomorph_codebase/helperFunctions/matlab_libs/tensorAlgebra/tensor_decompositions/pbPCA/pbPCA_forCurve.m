function [frameSet] = pbPCA_forCurve(point,imageStack,Domain,extraDim)
    if nargin == 3;extraDim = [];end
    tm = point(3);
    % get the image name
    imgF = imageStack{tm};
    % construct parameters for readin
    affine = eye(3);
    affine(1:2,3) = point(1:2)';
    % get the image name
    imgF = imageStack{tm};
    frameSet = pbPCA(imgF,affine,Domain.d,@lt);
    affineD = frameSet(1:end-1,end);
    affineD = point';
    % append dims
    if ~isempty(extraDim);frameSet = insertExtraDim(frameSet(1:end-1,1:end-1),extraDim);end;
    % replace affine
    frameSet = insertExtraDim(frameSet,1);
    frameSet(1:numel(affineD),end) = affineD;
    % return
    frameSet = phytoAaffine(frameSet);
end

function [frameSet] = insertExtraDim(frameSet,extraDim)
    Z1 = zeros(size(frameSet,1),size(extraDim,2));
    Z2 = zeros(size(extraDim,1),size(frameSet,2));
    frameSet = [frameSet Z1];
    extraDim = [Z2 extraDim];
    frameSet = [frameSet;extraDim];
end