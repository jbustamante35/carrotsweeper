function [net] = myGPUtest(X,Y,layers,EE,maxE)
%{
    layers = [imageInputLayer([size(X,1) size(X,2) 1],'Normalization','None')
          convolution2dLayer([11,11],20)
          reluLayer
          maxPooling2dLayer([2 2],'Stride',2)
          convolution2dLayer([7,7],3)
          fullyConnectedLayer(2)
          softmaxLayer
          classificationLayer()];
    %}
    %options = trainingOptions('sgdm','ExecutionEnvironment','parallel','MaxEpochs',3,'CheckpointPath','/mnt/snapper/nate/CP/','Plots','training-progress','ValidationFrequency',200,'ValidationData',VALData);
    options = trainingOptions('sgdm','ExecutionEnvironment',EE,'MaxEpochs',maxE);
    net = trainNetwork(X,Y,layers,options);
end

%{

    net = myGPUtest(rand(50,50,1,100),[ones(1,50) zeros(1,50)]');
    func = cFlow('myGPUtest');


    func.setMCRversion('v930');
    func.setGPU(4);
    func(rand(50,50,1,100),[ones(1,50),zeros(1,50)]');




    auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
    auth = auth{1};
    func.submitDag(auth,50,50);

%}