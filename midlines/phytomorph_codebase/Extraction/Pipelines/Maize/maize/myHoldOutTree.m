function [delta] = myHoldOutTree(factors,xTrain,yTrain,xTest,yTest)
    deltaTest = [];
    deltaTrain = [];
   
    for e = 1:size(yTrain,2)
        %R = classregtree(xTrain,yTrain(:,e));
        %R = RegressionTree.fit(xTrain,yTrain(:,e));
        ens = fitensemble(xTrain,yTrain(:,e),'Bag',factors,'Tree','type','regression');
        deltaTest = [deltaTest ens.predict(xTest)];
        deltaTrain = [deltaTrain ens.predict(xTrain)];
        %deltaTest = [deltaTest predict(R,xTest)];
        %deltaTrain = [deltaTrain predict(R,xTrain)];
        %deltaTest = [deltaTest eval(R,xTest)];
        %deltaTrain = [deltaTrain eval(R,xTrain)];
    end
    
    
    deltaTest = deltaTest - yTest;
    deltaTest = sum(deltaTest.*deltaTest,2).^.5;
    
    deltaTrain = deltaTrain - yTrain;
    deltaTrain = sum(deltaTrain.*deltaTrain,2).^.5;
    
    
    
    
    
    labels = [zeros(size(deltaTrain));ones(size(deltaTest))];
    delta = [deltaTrain;deltaTest];
    delta = [labels delta];
    
end