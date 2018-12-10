function [valError,cons,fileName] = makeObjFcn(optVars,trainImages,trainLabels,valImages,valLabels,toSLOW,DISP,MaxE,exeEnvironment,numClasses)
            

            %%%%%%%%%%%%%%%%%%%%
            %
            %%%%%%%%%%%%%%%%%%%%
            options = trainingOptions('sgdm',...
                'InitialLearnRate',optVars.InitialLearnRate,...
                'Momentum',optVars.Momentum,...
                'MaxEpochs',MaxE,...
                'L2Regularization',optVars.L2Regularization,...
                'Shuffle','every-epoch',...
                'Verbose',true,...
                'ExecutionEnvironment',exeEnvironment,...
                'Plots',DISP);

            
            
            SZ = size(trainImages);
            imageSize = SZ(1:3);
            
            
            
            
            initialNumFilters = round(32/4/sqrt(optVars.NetworkDepth));
            initialNumFilters = 4;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % layers
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % init layers vector to inpout image
            layers = imageInputLayer(imageSize);
            % build up other layers
            for e = 1:optVars.NetworkDepth
            layers = [layers 
                    convBlock(optVars.filterSize,initialNumFilters*e,1)
                    maxPooling2dLayer(2,'Stride',2)];
            end
            % remove last max pool
            layers(end) = [];
            % make fully connected layer and other final layers
            layers = [layers 
                fullyConnectedLayer(numClasses)
                softmaxLayer
                classificationLayer];
            
            
            imageAugmenter = imageDataAugmenter('RandRotation',[-90 90],'RandXScale',[.6 1.2],'RandYScale',[.6 1.2]);
            datasource = augmentedImageSource(imageSize,trainImages,trainLabels,...
                'DataAugmentation',imageAugmenter);
            
            
            trainedNet = trainNetwork(datasource,layers,options);
            

            



            % if slow then get steady state ~ for stuff
            if toSLOW
                options = trainingOptions('sgdm',...
                'InitialLearnRate',optVars.InitialLearnRate/10,...
                'Momentum',optVars.Momentum,...
                'MaxEpochs',MaxE,...
                'L2Regularization',optVars.L2Regularization,...
                'Shuffle','every-epoch',...
                'Verbose',true,...
                'ExecutionEnvironment',exeEnvironment);
                trainedNet = trainNetwork(datasource,trainedNet.Layers,options);
            end

            if ~isempty(valImages)
                predictedLabels = classify(trainedNet,valImages);
                valAccuracy = mean(predictedLabels == valLabels);
                valError = 1 - valAccuracy;
            else
                valError = 0;
                valAccuracy = 1;
            end

            fileName = num2str(valError,10) + ".mat";
            save(fileName,'trainedNet','valError','options')
            cons = [];
end

       
