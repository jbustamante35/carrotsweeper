function [BayesObject,trainedNetwork] = mySimpleCnn(para,xTrain,yTrain,xTest,yTest,maxTIME,exeEnvironment,maxE)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % optimization of parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['***************************************\n']);
    fprintf(['start optimization of parameters\n']);
    fprintf(['***************************************\n']);
    % create objective function
    func = @(P)mySimpleObjectiveFunction(P,xTrain,yTrain,xTest,yTest,exeEnvironment,maxE);
    % optimize the objective function
    BayesObject = bayesopt(func,para,'Verbose',1,'UseParallel',false,'MaxTime',maxTIME);
    fprintf(['***************************************\n']);
    fprintf(['end optimization of parameters\n']);
    fprintf(['***************************************\n']);



    fprintf(['***************************************\n']);
    fprintf(['start training network\n']);
    fprintf(['***************************************\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % train network
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    para = BayesObject.XAtMinEstimatedObjective;
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
        'Plots','none',...
        'Verbose',true,...
        'ExecutionEnvironment','parallel');

    middleLayers = createMiddle(para.layerSize,para.layerNumber,para.MiddleNumber);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make layers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    layers = [ ...
        imageInputLayer([size(xTrain,1) size(xTrain,2) 1])
        middleLayers
        fullyConnectedLayer(size(yTrain,2))
        regressionLayer];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % train network
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    trainedNetwork = trainNetwork(cat(4,xTrain,xTest),cat(1,yTrain,yTest),layers,options);
    fprintf(['***************************************\n']);
    fprintf(['end training network\n']);
    fprintf(['***************************************\n']);

    
end


function [d] = mySimpleObjectiveFunction(para,xTrain,yTrain,xTest,yTest,exeEnvironment,maxE)
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
        'Plots','none',...
        'Verbose',true,...
        'ExecutionEnvironment',exeEnvironment);


    middleLayers = createMiddle(para.layerSize,para.layerNumber,para.MiddleNumber);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make layers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    layers = [ ...
        imageInputLayer([size(xTrain,1) size(xTrain,2) 1])
        middleLayers
        fullyConnectedLayer(size(yTrain,2))
        regressionLayer];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % train network
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    network = trainNetwork(xTrain,yTrain,layers,options);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % slow options
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    options = trainingOptions('sgdm',...
        'LearnRateSchedule','none',...
        'InitialLearnRate',para.InitialLearnRate/10,...
        'L2Regularization',para.L2Regularization,...
        'Momentum',para.Momentum,...
        'Plots','none',...
        'ExecutionEnvironment',exeEnvironment);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % slow train network
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    network = trainNetwork(xTrain,yTrain,network.Layers,options);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % predict and measure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yPredict = network.predict(xTest);
    d = norm(yPredict - yTest);
end

function [layers] = createMiddle(layerSize,layerNumber,numberBlocks)
    layers = [...
        convolution2dLayer([layerSize layerSize],layerNumber,'Padding','same')
        reluLayer];
    layers = repmat(layers,numberBlocks,1);
end