function [BO,net] = gpuTrainRegression(nTrain,trainX,trainY,testX,testY,exeEnvironment,para)

    [func] = generateZoomNetworkLayers(nTrain(1),trainX,trainY,testX,testY,exeEnvironment);

    BO = bayesopt(func,para,'Verbose',2,...
        'AcquisitionFunctionName','expected-improvement-plus');
    
    
    
    [err,net] = myregressionError(nTrain(2),trainX,trainY,testX,testY,exeEnvironment,BO.XAtMinEstimatedObjective);
end