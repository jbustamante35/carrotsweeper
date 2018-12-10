function [delta] = myFUN(factors,xTrain,yTrain,xTest,yTest)

    

    %{
    %net = fitnet(factors);
    %net = cascadeforwardnet(factors);
    net = feedforwardnet([factors 3 factors]);
    %net.trainParam.showWindow = false;
    net = train(net,xTrain',yTrain');
    
    deltaTest = net(xTest')';
    deltaTest = deltaTest - yTest;
    deltaTest = sum(deltaTest.*deltaTest,2).^.5;
    
    
    deltaTrain = net(xTrain')';  
    deltaTrain = deltaTrain - yTrain;
    deltaTrain = sum(deltaTrain.*deltaTrain,2).^.5;
    
    
    labels = [zeros(size(deltaTrain));ones(size(deltaTest))];
    delta = [deltaTrain;deltaTest];
    delta = [labels delta];
    %}
    %{
    %[A,B,r,U,V,stats] = canoncorr(xTrain,yTrain);
    muX = mean(xTrain,1);
    muY = mean(yTrain,1);
    xTrain = bsxfun(@minus,xTrain,muX);
    yTrain = bsxfun(@minus,yTrain,muY);
    [A,B,U,V,D] = myCCA(xTrain,yTrain,factors);
    
    
    U = bsxfun(@minus,xTest,muX)*A;
    deltaTest = (pinv(B')*U')';
    deltaTest = deltaTest - yTest;
    deltaTest = sum(deltaTest.*deltaTest,2).^.5;
    
    
    U = xTrain*A;
    deltaTrain = (pinv(B')*U')';    
    deltaTrain = deltaTrain - yTrain;
    deltaTrain = sum(deltaTrain.*deltaTrain,2).^.5;
    
    
    labels = [zeros(size(deltaTrain));ones(size(deltaTest))];
    delta = [deltaTrain;deltaTest];
    delta = [labels delta];
    %}
    
    
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(xTrain,yTrain,factors);    
    deltaTest = [ones(size(xTest,1),1),xTest]*BETA;
    deltaTest = deltaTest - yTest;
    deltaTest = sum(deltaTest.*deltaTest,2).^.5;
    
    deltaTrain = [ones(size(xTrain,1),1),xTrain]*BETA;
    deltaTrain = deltaTrain - yTrain;
    deltaTrain = sum(deltaTrain.*deltaTrain,2).^.5;
    
    
    labels = [zeros(size(deltaTrain));ones(size(deltaTest))];
    delta = [deltaTrain;deltaTest];
    delta = [labels delta];
    
     
    
end