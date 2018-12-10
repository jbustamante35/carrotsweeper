function [FC,toFeatureFORCNN] = prepareForNetworks(T,FE,FU)
    opData1 = thawTensor(T,4);
    opData2 = thawTensor(T,5);
    toFeature = cat(4,opData1,cos(opData2),sin(opData2));
  
    toFeature = permute(toFeature,[2 1 4 3]);
    toFeatureFORCNN = toFeature;
    toFeature(:,:,2,:) = cumsum(toFeature(:,:,2,:),2);
    toFeature(:,:,3,:) = cumsum(toFeature(:,:,3,:),2);
    %toFeatureFORCNN = toFeature;
    
    toFeature = permute(toFeature,[1 2 4 3]);
    FC = [];
    for d = 1:size(toFeature,4)
        tmp = toFeature(:,:,:,d);
        sz = size(tmp);
        tmp = reshape(tmp,[prod(sz(1:2)) sz(3)]);
        tmpFC = PCA_REPROJ(tmp',FE{d},FU{d});
        FC = cat(2,FC,tmpFC);
    end
end