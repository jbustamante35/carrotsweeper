function [model] = buildModel(data,model)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT:
    %       data:       data.domain is the tensor for the domain
    %                   data.codomain is the tensor for the codomain
    %       model:      the map from data.domain -> data.codomain
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %       outLabel:   output labels
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    reportGen = 0;
    log = 1;
    %%%%%%%%%%%%%%%%%%%%
    % report
    %%%%%%%%%%%%%%%%%%%%
    fprintf(['building model@init.......\n']);
    
    %%%%%%%%%%%%%%%%%%%%
    % construct group structure -  this is a meta-method
    % each method will store results in column vector in gS cell array
    %%%%%%%%%%%%%%%%%%%%
    fprintf(['grouping@data.......\n']);
    para = [];
    gS = buildGroups(data,para);
    
    
    
    %%%%%%%%%%%%%%%%%%%%
    % build linear model subgroup of data.domain
    %%%%%%%%%%%%%%%%%%%%
    fprintf(['model@data.......\n']);
    para.nComp = 3;
    %data_model = build_pwlDataModel(data,gS,para);
    data_model = build_dataModel(data,gS,para);
    
    
    %%%%%%%%%%%%%%%%%%%%
    % project all of data.domain into model -> latent
    %%%%%%%%%%%%%%%%%%%%
    fprintf(['latent@data.......\n']);
    para = [];
    latent_data = eval_pwlDataModel(data,data_model,para);    
    
    
    %%%%%%%%%%%%%%%%%%%%
    % project a special section
    %%%%%%%%%%%%%%%%%%%%
    map = [];
    fprintf(['learn@domain->codomain.......\n']);
    tdata.domain = latent_data{1}{2};
    tdata.codomain = data.codomain;
    
    %%%%%%%%%%%%%%%%%%%%
    % learn on latent parameters
    %%%%%%%%%%%%%%%%%%%%
    learningMap = buildLearningFunction(map,tdata);
    
    
    
    
    
    model.data_model = data_model;
    model.learningMap = learningMap;

    %{
    if log
        model.log.expV.g0.u0.value = u0;
        model.log.expV.g0.s0.value = s0;
        
        model.log.expV.g1.u1.value = u1;
        model.log.expV.g1.s1.value = s1;
        
        model.log.hist.data.domain.value = X;
        model.log.hist.data.codomain.value = [pdf0;pdf1];
    end
    
    
    
    if reportGen
        % generate data along covec
        plot(X,pdf0,'k');
        hold on
        plot(X,pdf1,'r');
    end
    %}
end