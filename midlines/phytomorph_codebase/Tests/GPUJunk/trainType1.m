function [net] = trainType1(X,Y,layers,maxE)
    %imageAugmenter = imageDataAugmenter('RandRotation',[-90 90],'RandXScale',[.6 1.2],'RandYScale',[.6 1.2]);
    %imageSize = [21 21 1];
    %datasource = augmentedImageSource(imageSize,X,categorical(Y),'DataAugmentation',imageAugmenter);
    options = trainingOptions('sgdm','ExecutionEnvironment','multi-gpu','MaxEpochs',maxE);
    net = trainNetwork(X,Y,layers,options);
end