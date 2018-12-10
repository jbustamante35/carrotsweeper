function [model grade] = initEvalModel(data,model)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT:
    %       data:       data.domain is the tensor for the domain
    %                   data.codomain is the tensor for the codomain
    %       model:      the map from data.domain -> data.codomain
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %       outLabel:   output labels
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create grade for output
    grade = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % build model phase
    % if domain the domain is not empty and the model is null and the
    % codomain is not empty then build model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % notes: added(13.01.25) - there is not yet a
    %        grading structure to determine if a model is "good" or stable
    %        also the notion of testing. the idea of testing seems 
    %        to be "outside" the buildModel, however, the buildModel
    %        could have internal splits on the data.domain
    %        such that there is a testing and validation.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~(isempty(data.domain)) 
        % if the model is null and there there is a domain and codomain then build model
        if strcmp(model.type,'null') && ~isempty(data.codomain)     
            model = buildModel(data,model);
        end
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % eval model phase
    % if there is a model the eval the model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % grade the data via the model
    if ~strcmp(model,'null')
        grade = evalModel(data,model);
        grade = grade(1).binary;
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HACK
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(grade)
        grade = rand(1,size(data.domain,1));
        grade = grade > .6;
    end
end