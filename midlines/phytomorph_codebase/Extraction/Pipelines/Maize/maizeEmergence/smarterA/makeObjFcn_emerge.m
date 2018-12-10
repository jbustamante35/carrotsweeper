function [valError,cons,trainedNet] = makeObjFcn_emerge(...
                                    optVars,...
                                    trainImages,...
                                    trainLabels,...
                                    valImages,...
                                    valLabels,...
                                    toSLOW,...
                                    DISP,...
                                    MaxE,...
                                    exeEnvironment,...
                                    nDIMS,...
                                    TYPE)

        numClasses = numel(unique(trainLabels));
                                
        if nargin == 10
            TYPE = 'sequence';
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create options for major train
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        options = trainingOptions('sgdm',...
            'InitialLearnRate',optVars.InitialLearnRate,...
            'Momentum',optVars.Momentum,...
            'MaxEpochs',MaxE,...
            'L2Regularization',optVars.L2Regularization,...
            'Shuffle','every-epoch',...
            'Verbose',true,...
            'ExecutionEnvironment',exeEnvironment,...
            'Plots',DISP);




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make layers
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        layers = [ ...
            sequenceInputLayer(nDIMS)
            lstmLayer(optVars.lstmStates,'OutputMode',TYPE)
            fullyConnectedLayer(numClasses)
            softmaxLayer
            classificationLayer];



        trainedNet = trainNetwork(trainImages,trainLabels,layers,options);



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create options for slow train
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            trainedNet = trainNetwork(valImages,valLabels,trainedNet.Layers,options);
        end

        if ~isempty(valImages)
            predictedLabels = classify(trainedNet,valImages);
            for e = 1:numel(predictedLabels)
                
                if iscell(valLabels)
                    tmp1 = find(double(valLabels{e})==2);
                else
                    tmp1 = find(double(valLabels(e))==2);
                end
                
                
                if iscell(predictedLabels)
                    tmp2 = find(double(predictedLabels{e})==2);
                else
                    tmp2 = find(double(predictedLabels(e))==2);
                end
                
                
                if ~isempty(tmp1) && ~isempty(tmp2)
                    delta(e) = abs(tmp1(1)-tmp2(1));
                elseif ~isempty(tmp1)
                    delta(e) = tmp1(1);
                elseif ~isempty(tmp2)
                    delta(e) = tmp2(1);
                end
            end
            valError = mean(delta);
        else
            valError = 0;

        end
        cons = [];

        
        %{
        fileName = num2str(valError,10) + ".mat";
        save(fileName,'trainedNet','valError','options')
        cons = [];
        %}
end

       
