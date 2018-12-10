function [basisStack] = transForm_extractToFeatureTensor(basisStack)
     amp = thawTensor(basisStack,1);
     ang = thawTensor(basisStack,2);
     cang = cumsum(cos(ang),1);
     sang = cumsum(sin(ang),1);
     basisStack = cat(4,amp,cang,sang);
     basisStack = permute(basisStack,[1 2 4 3]);
end