function layers = convBlock(filterSize,numFilters,numConvLayers)
    layers = [
        convolution2dLayer(filterSize,numFilters,'Padding',round((filterSize-1)/2))
        reluLayer];
    layers = repmat(layers,numConvLayers,1);
end