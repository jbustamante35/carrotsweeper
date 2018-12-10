function [grade] = evalModel(data,model)
    %%%%%%%%%%%%%%%%%%%%
    % project all of data.domain into model -> latent
    %%%%%%%%%%%%%%%%%%%%
    para = [];
    fprintf(['latent@data.......\n']);
    latent_data = eval_pwlDataModel(data,model.data_model,para);
    
    %%%%%%%%%%%%%%%%%%%%
    % project a special section - only the group==1
    %%%%%%%%%%%%%%%%%%%%
    map = [];
    fprintf(['learn@domain->codomain.......\n']);
    tdata.domain = latent_data{1}{2};
    tdata.codomain = data.codomain;
    
    %%%%%%%%%%%%%%%%%%%%
    % grade
    %%%%%%%%%%%%%%%%%%%%
    grade = eval_learningModel(tdata.domain,model.learningMap);
    
end