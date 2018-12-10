function [para,labeledTrainingPackage] = optimizeAIoutput(AI_layer,labeledTrainingPackage,initPara)
    
    C = labeledTrainingPackage.C;
    T = labeledTrainingPackage.T;
    parfor e = 1:size(labeledTrainingPackage.M,3)
        tmpData = struct('C',[],'T',[]);
        fprintf(['Start:apply AI layer to training package.\n']);
        tmpData.C = C(:,:,e);
        tmpData.T = T(:,:,:,:,e);
        probMap(:,:,:,e) = applyAIlayer(tmpData,AI_layer,'','',[108 108]);
        fprintf(['End:apply AI layer to training package.\n']);
    end
    
    
    %load('./faster.mat','probMap');

    labeledTrainingPackage.probMap = probMap;
    
    paraLabels = [2*ones(1,13) 3*ones(1,13)];
    paraLabels = [2*ones(1,5) 3*ones(1,5)];
    
    init = [.1 1*3/9 1*3/9 1*3/9 200 800 .12 .03 100 400 .12 .03 5];
    delta = [.05 .2 .2 .2 50 50 .1 .1 20 20 .1 .1 1]*.05^-1;
    
    f = logical([1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 1]);    

    initDelta{1} = [init delta];

    initDelta{1} = initDelta{1}(f(14:end));
    if nargin == 2
        delta_init{1} = delta;
    else
        delta_init{1} = initPara{1}(1:13);
    end
    delta_init{1} = delta_init{1}(f(1:13));
    
    
    
    
    
    
    
    mm{3} = @(X,Y)matthews_correlation(X,Y);
    mm{2} = @(X,Y)myMetric(X,Y);
    mm{1} = @(X,Y)positiveLikehoodRatio(X,Y);
    para = {};
    
    per = .2;
    

    
    for o = 1:numel(delta_init)%1%numel(mm)
        
        %func = @(X)opti2(probMap,labeledTrainingPackage.M(1:108,1:108,:),X,labeledTrainingPackage.oI(1:108,1:108,:),initDelta{o},paraLabels);
        func = @(X)integrateAndmeasureProbMaps(probMap,labeledTrainingPackage.M(1:108,1:108,:),X,labeledTrainingPackage.oI(1:108,1:108,:),initDelta{o},paraLabels)
        ops = optimoptions('particleswarm','Display','iter');
        para{o} = particleswarm(func,5,delta_init{o}-per*delta(f(1:13)),delta_init{o}+per*delta(f(1:13)),ops);
        %para{o} = fminsearch(func,delta_init{o}(1:end),ops);
        ops = optimset('Display','iter');
        para{o} = fminsearch(func,delta_init{o}(1:end),ops);
        %[~,T] = opti2(probMap,labeledTrainingPackage.M(1:108,1:108,:),para{o},labeledTrainingPackage.oI(1:108,1:108,:));
        para{o} = [para{o} initDelta{o}];
        [~,PT{o},PT2{o}] = integrateAndmeasureProbMaps(probMap,labeledTrainingPackage.M(1:108,1:108,:),para{o}(1:5),labeledTrainingPackage.oI(1:108,1:108,:),initDelta{o},paraLabels,[]);
    end
    

   
end