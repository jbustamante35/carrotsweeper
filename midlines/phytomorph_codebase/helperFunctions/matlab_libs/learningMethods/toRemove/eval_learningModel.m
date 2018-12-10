function [grade] = eval_learningModel(data,model)
    % NOTE :    QUESTION : why not use full pdf. why project to ray which
    % best solves.  answer might be about noise. and overlap.
    % loop over each model type and eval
    for meth = 1:numel(model)        
        switch model(meth).type
            case 'lda'
                grade(meth).prob = data*model(meth).para.lambdaSplit.value;
                if ~isempty(model(meth).para.threshold.value)
                    grade(meth).binary = grade(meth).prob > model(meth).para.threshold.value;
                end
            case 'pdf'
                grade(meth).prob = nprobData(data,model(meth).para.mv.value);                
                if ~isempty(model(meth).para.threshold.value)
                    grade(meth).binary = grade(meth).prob > model(meth).para.threshold.value;
                end
        end
    end
end