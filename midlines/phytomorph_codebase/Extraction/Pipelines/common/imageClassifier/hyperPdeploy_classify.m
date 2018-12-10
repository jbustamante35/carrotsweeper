function [BayesObject] = hyperPdeploy_classify(optimVars,X,Y,vX,vY,exeEnvironment,maxE,toSlow,maxEval,maxTime,nDIMS,TYPE)
    if nargin == 11
        TYPE = 'sequence';
    end

    
    
    ObjFcn = @(hyperParameters,trainX,trainY,valX,valY)makeObjFcn_classify(hyperParameters,trainX,trainY,valX,valY,toSlow,'None',maxE,exeEnvironment{1},nDIMS,TYPE);
    
    
    


    BayesObject = bayesopt(@(OPS)ObjFcn(OPS,X,Y,vX,vY),...
        optimVars,...
        'MaxObj',maxEval,...
        'MaxTime',maxTime,...
        'IsObjectiveDeterministic',false,...
        'UseParallel',exeEnvironment{2},...
        'PlotFcn',{});
end