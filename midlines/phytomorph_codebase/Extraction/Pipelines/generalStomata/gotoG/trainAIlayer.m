function [AI_layer] = trainAIlayer(labeledTrainingPackage,fraction)

    %trainX = permute(labeledTrainingPackage.C,[1 3 2]);
    trainX = labeledTrainingPackage.C;
    trainY = labeledTrainingPackage.M;
  
    
    trainY(1,:,:) = [];
    trainY(:,1,:) = [];
    szY = size(trainY);
    trainY = reshape(trainY,[prod(szY(1:3)) 1])';
    trainY_FULL = trainY;
    
    szX = size(trainX);
    trainX = reshape(trainX,[prod(szX(1)) prod(szX(2:3))]);
    
    fidx1 = find(trainY == 1);
    fidx0 = find(trainY == 0);
    % permute the 0 class
    fidx0 = fidx0(randperm(numel(fidx0)));
    numZeros = min(round(fraction(1)*numel(fidx1)),numel(fidx0));
    
    
    
    trainY = [trainY(fidx0(1:numZeros)) trainY(fidx1)];
    trainX = [trainX(:,fidx0(1:numZeros)) trainX(:,fidx1)]; 
    
    
    cv = cvpartition(trainY','KFold',2);
    train1 = find(cv.training(1));
    train2 = find(cv.training(2));
    [IDX, Z] = rankfeatures(trainX, trainY);
    
    
    
    fprintf(['Starting training of PLS.\n']);
    for l = 1:20
        [~,~,~,~, beta,~,~,~,~] = plsregress(trainX(:,train1)',trainY(train1)',l);
        yPre = [ones(numel(train2),1) trainX(:,train2)']*beta;
        CV(l) = corr(yPre,trainY(train2)');
        fprintf(['Done fitting model (' num2str(l) ').\n'])
    end
    fprintf(['End training of PLS.\n']);
    [~,~,~,~, beta,~,~,~,~] = plsregress(trainX(:,train1)',trainY(train1)',5);
    
    fprintf(['Starting training of NN.\n']);
    CV = [];
    for l = 1:20
        NET = patternnet([l]);
        [NET,tr] = train(NET,trainX,full(ind2vec((trainY+1))),'useParallel','yes');
        CV(l) = min(tr.vperf);
        fprintf(['Done fitting model (' num2str(l) ').\n'])
    end
    fprintf(['End training of NN.\n']);
    NET = patternnet([15]);
    [NET,tr] = train(NET,trainX,full(ind2vec((trainY+1))),'useParallel','yes');
    
    
    
    % tree
    fprintf(['Starting training of TREE.\n']);
    tree = fitctree(trainX(IDX(1:3),:)',trainY(:)','ScoreTransform','logit');
    fprintf(['End training of TREE.\n']);
    
    % glm
    fprintf(['Starting training of GLM.\n']);
    GLM = fitglm(trainX(IDX(1:6),:)',trainY','quadratic','Distribution','binomial','link','logit','BinomialSize',1);
    fprintf(['End training of GLM.\n']);
    
    % step
    fprintf(['Starting training of STEP.\n']);
    STEP = stepwiseglm(trainX(IDX(1:6),:)',trainY','quadratic','Distribution','binomial','link','logit');
    fprintf(['End training of STEP.\n']);
    
    % start lasso
    fprintf(['Starting training of STEP.\n']);
    [B,FitInfo] = lassoglm(trainX(:,train1)',trainY(train1)','binomial','link','probit','CV',10);
    indx = FitInfo.Index1SE;
    B0 = B(:,indx);
    cnst = FitInfo.Intercept(indx);
    B1 = [cnst;B0];
    
    
  
    
    % start ridge
    cv = cvpartition(trainY','KFold',2);
    train1 = find(cv.training(1));
    [B,FitInfo] = lassoglm(trainX(:,train1)',trainY(train1)','binomial','link','probit','CV',10,'Alpha',.001);
    indx = FitInfo.Index1SE;
    B0 = B(:,indx);
    cnst = FitInfo.Intercept(indx);
    B2 = [cnst;B0];
    
    %preds = glmval(B1,trainX(:,train2)','logit');
    %C = confusionmat(logical(trainY(train2))',preds>.5); % plot residual
    AI_layer.LASSO = B1;
    AI_layer.RIDGE = B2;
    AI_layer.PLS = beta;
    AI_layer.TREE = tree;
    AI_layer.GLM = GLM;
    AI_layer.STEP = STEP;
    AI_layer.IDX = IDX;
    AI_layer.FDA = fitcdiscr(trainX',trainY','DiscrimType','quadratic');
    AI_layer.NN = NET;
    
    
    
    
    szT = size(labeledTrainingPackage.T);
    trainT = reshape(labeledTrainingPackage.T,[szT(1:3) prod(szT(4:5))]);
    %
    layers = [imageInputLayer(szT(1:3));
          convolution2dLayer(5,13);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(2);
          softmaxLayer();
          classificationLayer()];
    options = trainingOptions('sgdm','Plots','training-progress','MaxEpochs',20,'InitialLearnRate',0.0001,'ExecutionEnvironment','parallel');
    if fraction(2) ~= 1
        trainT = trainT(:,:,:,1:fraction(2):end);
        trainY_FULL = trainY_FULL(1:fraction(2):end);
    end
    AI_layer.CNN = trainNetwork(trainT,categorical(trainY_FULL'),layers,options);
    
      
      
      
    AI_layer.NN_func = @(X)current_sorghum_network(X);
    genFunction(AI_layer.NN,'/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/generalStomata/gotoG/current_sorghum_network.m');

      
      
      
      %{
    HOPE2 = AI_layer.FDA.predict(labeledTrainingPackage.C(:,:,1)');
    
    HOPE = AI_layer.NN(labeledTrainingPackage.C(:,:,1));
     %}
end