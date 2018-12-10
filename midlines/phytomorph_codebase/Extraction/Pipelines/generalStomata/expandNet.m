function [newLayers] = expandNet(net)
    newLayers = net.Layers(1);
    for e = 2:numel(net.Layers)
        if isa(net.Layers(e),'nnet.cnn.layer.Convolution2DLayer')
            filterSize = net.Layers(e).FilterSize+2;
            numLayers = net.Layers(e).NumFilters;
            layer = convolution2dLayer(filterSize,numLayers,'Padding',round((7-1)/2));
            layer.Weights = padarray(net.Layers(e).Weights,[1 1],0,'both');
            newLayers = [newLayers 
                        layer];
            
        else
            newLayers = [newLayers 
                        net.Layers(e)];
        end
    end


  
end