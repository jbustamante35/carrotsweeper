function [] = myNetworkBuilder(imageSize,initFilterSize,proposedDepth,numberOfFiltersPerLayer,strideSize,strideNumber)
    % number of filters per layer can increase,stable,decrease.
    % current setting is a consstant number of filters per layer
    numberOfFiltersPerLayer = repmat(numberOfFiltersPerLayer,proposedDepth,1);
    % constant number of blocks per layer
    numberBlocks = repmat(numberOfFiltersPerLayer,proposedDepth,1);
    % init with image input layers
    layers = imageInputLayer(imageSize);

     [layers] = buildLayers(initFilterSize,numberOfFiltersPerLayer,strideSize,strideNumber,numberBlocks,totalLayers)
end



function [value] = parameterRateFunction(para,depth)
    initValue = para(1);
    rate = para(2);
    value = round(initValue - rate*depth);
end


function [layers] = constructLayers(imageSize,layerInitSize,strideSize,proposedDepth)
    %[proposedDepth] = calcNumberOfLayers(layerInitSize,strideSize,proposedDepth);
    %layers = [];
    layers = [imageInputLayer(imageSize)];
end

function layers = buildBlockLayer(filterSize,numFilters,numBlocks)
    layers = [];
    for block = 1:numBlocks
        convolutionName = [baseName '_convolution_' num2str(block)];
        reluName = [baseName '_relu_' num2str(block)];
        layers = [layers 
                    convolution2dLayer(filterSize,numFilters,'Padding',(filterSize-1)/2,'Name',convolutionName)
                    reluLayer('Name',reluName)];
    end
end

function [depth] = calcNumberOfLayers(para,proposedDepth)
    initValue = para(1);
    rateValue = para(2);
    acutalDepth = 1;
    for e = 1:proposedDepth
        if floor(layerInitSize/strideSize) > 1
            acutalDepth = acutalDepth + 1;
        end
    end
end

function [layers] = buildLayers(para)
    layers = [];
    
    

    for d = 1:proposedDepth
        % variable filterSize,numberOfFilters,numberOfBlocks,strideSize,strideNumber with layer depth 
        filterSize = parameterRateFunction(para.filterSizeParameter,d);
        numberOfFilters = parameterRateFunction(para.filterNumberParameter,d);
        numberOfBlocks = parameterRateFunction(para.filterBlockNumber,d);
        strideSize = parameterRateFunction(para.strideSize,d);
        strideNumber = parameterRateFunction(para.strideNumber,d);
        % construct block layer
        toAdd = buildBlockLayer(filterSize,numberOfFilters,numberOfBlocks);
        % cap off block layer
        toAdd = [toAdd maxPooling2dLayer(strideSize,'Stride',strideNumber)];
        % stack block layer 
        layers = [layers toAdd];
    end

end



%{
imageSize = [61 61 1];
initFilterSize = 7;
proposedDepth = 2;
numberOfFiltersPerLayer = 4;
strideSize = 2;
strideNumber = 2;
myNetworkBuilder(imageSize,initFilterSize,proposedDepth,numberOfFiltersPerLayer,strideSize,strideNumber)
%}