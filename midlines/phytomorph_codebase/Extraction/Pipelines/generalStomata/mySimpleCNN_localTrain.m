function [trainedNetwork] = mySimpleCNN_localTrain(xData,yData,para,maxE)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make options for training
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    options = trainingOptions('sgdm',...
        'LearnRateSchedule',char(para.LearnRateSchedule),...
        'InitialLearnRate',para.InitialLearnRate,...
        'LearnRateDropFactor',para.LearnRateDropFactor,...
        'LearnRateDropPeriod',para.LearnRateDropPeriod,...
        'MaxEpochs',maxE,...
        'L2Regularization',para.L2Regularization,...
        'Momentum',para.Momentum,...
        'Plots','training-progress',...
        'Verbose',true,...
        'ExecutionEnvironment','parallel');
    
    middleLayers = createMiddle(para.layerSize,para.layerNumber,para.MiddleNumber);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make layers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    layers = [ ...
        imageInputLayer([size(xData,1) size(xData,2) 1])
        middleLayers
        fullyConnectedLayer(size(yData,2))
        regressionLayer];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % train network
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    trainedNetwork = trainNetwork(xData,yData,layers,options);

end


function [layers] = createMiddle(layerSize,layerNumber,numberBlocks)
    layers = [...
        convolution2dLayer([layerSize layerSize],layerNumber,'Padding','same')
        reluLayer];
    layers = repmat(layers,numberBlocks,1);
end