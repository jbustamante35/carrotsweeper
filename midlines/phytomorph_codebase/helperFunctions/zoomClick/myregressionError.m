function [err,net] = myregressionError(nTRAIN,trainX,trainY,testX,testY,exeEnvironment,OPTI)



    xSize = size(trainX);
    xSize(end) = [];
    ySize = size(trainY,2);
    layers = [imageInputLayer(xSize,'Normalization','None');
              convolution2dLayer([OPTI.nHood OPTI.nHood],OPTI.nKernels);
              reluLayer();
              maxPooling2dLayer(2,'Stride',2);
              fullyConnectedLayer(ySize);
              regressionLayer();];
    options = trainingOptions('sgdm',...
        'InitialLearnRate',OPTI.initLearnRate,...
        'MaxEpochs',nTRAIN,...
        'L2Regularization',OPTI.L2Regularization,...
        'Momentum',OPTI.Momentum,...
        'ExecutionEnvironment',exeEnvironment);
    
    
    
    net = trainNetwork(trainX,trainY,layers,options);
    
    preY = net.predict(testX);
    delta = preY - testY;
    
    
    err = sum(sum(delta.*delta,2).^.5,1);
end