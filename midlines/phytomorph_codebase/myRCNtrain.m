function [detector] = myRCNtrain(trainTable,trainTableData,oPath)
    mkdir(oPath);
    for e = 1:size(trainTableData,4)
        trainTable.imageFileName{e} = [oPath num2str(e) '.tif'];
        imwrite(trainTableData(:,:,:,e),trainTable.imageFileName{e});
    end
    %% make RCN for meta box
    % Create image input layer.
    inputLayer = imageInputLayer([64 64 3]);
    % Define the convolutional layer parameters.
    filterSize = [6 6];
    numFilters = 20;
    % Create the middle layers.
    middleLayers = [
        convolution2dLayer(filterSize, numFilters, 'Padding', 1)
        reluLayer()
        convolution2dLayer(filterSize, numFilters, 'Padding', 1)
        reluLayer()
        maxPooling2dLayer(3, 'Stride',2)
        ];

    finalLayers = [

        % Add a fully connected layer with 64 output neurons. The output size
        % of this layer will be an array with a length of 64.
        fullyConnectedLayer(64)

        % Add a ReLU non-linearity.
        reluLayer()

        % Add the last fully connected layer. At this point, the network must
        % produce outputs that can be used to measure whether the input image
        % belongs to one of the object classes or background. This measurement
        % is made using the subsequent loss layers.
        fullyConnectedLayer(width(trainTable))

        % Add the softmax loss layer and classification layer.
        softmaxLayer()
        classificationLayer()
    ];
    layers = [
        inputLayer
        middleLayers
        finalLayers
        ];

    % Options for step 1.
    optionsStage1 = trainingOptions('sgdm', ...
        'MaxEpochs', 10, ...
        'InitialLearnRate', 1e-5, ...
        'CheckpointPath', tempdir,...
        'ExecutionEnvironment','auto','Plots','training-progress');

    % Options for step 2.
    optionsStage2 = trainingOptions('sgdm', ...
        'MaxEpochs', 10, ...
        'InitialLearnRate', 1e-5, ...
        'CheckpointPath', tempdir,...
        'ExecutionEnvironment','auto','Plots','training-progress');

    % Options for step 3.
    optionsStage3 = trainingOptions('sgdm', ...
        'MaxEpochs', 10, ...
        'InitialLearnRate', 1e-6, ...
        'CheckpointPath', tempdir,...
        'ExecutionEnvironment','auto','Plots','training-progress');

    % Options for step 4.
    optionsStage4 = trainingOptions('sgdm', ...
        'MaxEpochs', 10, ...
        'InitialLearnRate', 1e-6, ...
        'CheckpointPath', tempdir,...
        'ExecutionEnvironment','auto','Plots','training-progress');

    options = [
        optionsStage1
        optionsStage2
        optionsStage3
        optionsStage4
        ];
    %%
    detector = trainFasterRCNNObjectDetector(trainTable, layers, options, ...
            'NegativeOverlapRange', [0 0.3], ...
            'PositiveOverlapRange', [0.6 1], ...
            'BoxPyramidScale', 1.2);
        
end
%{
%%
        tableFile = '/mnt/tetra/nate/maizeSeedlingData/tableFile.mat';
        load(tableFile);

%}