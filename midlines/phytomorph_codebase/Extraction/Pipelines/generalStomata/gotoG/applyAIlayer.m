function [probMap] = applyAIlayer(data,published_AI_layer,basisU,basisE,SZ)
    fprintf(['Starting application of AI layer.\n']);tic
    % load from file if needed
    if ischar(data)
        data = myRemoteLoader(data,'T','rT');
    end
    
    
    % if C does not exist then compute
    if ~isfield(data,'C')
        data = prepareData(data,basisU,basisE);
    end
    
    
    % NN
    fprintf(['Starting NN subLayer. \n']);
    NN_out = published_AI_layer.NN_func(data.C);
    NN_out = NN_out(2,:)';
    % FDA
    fprintf(['Starting LDA subLayer. \n']);
    [~,FDA_out] = published_AI_layer.FDA.predict(data.C');
    FDA_out = FDA_out(:,2);
    % CNN
    fprintf(['Starting CNN subLayer. \n']);
    CNN_out = published_AI_layer.CNN.predict(data.T);
    CNN_out = CNN_out(:,2);
    % LASSO
    fprintf(['Starting LASSO subLayer. \n']);
    LASSO_out = glmval(published_AI_layer.LASSO,data.C','logit');
    % RIDGE
    fprintf(['Starting LASSO subLayer. \n']);
    RIDGE_out = glmval(published_AI_layer.RIDGE,data.C','logit');
    % RIDGE
    fprintf(['Starting LASSO subLayer. \n']);
    PLS_out = [ones(size(data.C',1),1) data.C']*published_AI_layer.PLS;
    TREE_out = published_AI_layer.TREE.predict(data.C(published_AI_layer.IDX(1:3),:)');
    STEP_out = published_AI_layer.STEP.predict(data.C(published_AI_layer.IDX(1:6),:)');
    
    GLM_out = published_AI_layer.GLM.predict(data.C(published_AI_layer.IDX(1:6),:)');
    
    probMap = [NN_out FDA_out CNN_out LASSO_out RIDGE_out PLS_out TREE_out STEP_out GLM_out];
    fprintf(['Starting Stack of AI subLayer output(s). \n']);
    probMap = reshape(probMap,[SZ size(probMap,2)]);
    fprintf(['Ending application of AI layer.' num2str(toc) '\n']);
end