function [func] = generateZoomNetworkLayers(nTRAIN,trainX,trainY,testX,testY,exeEnvironment)
    


    func = @(X)myregressionError(nTRAIN,trainX,trainY,testX,testY,exeEnvironment,X);

    
    
end