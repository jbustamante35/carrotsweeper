function [model] = buildLearningFunction(model,data)    
    %{
    %%%%%%%%%%%%%%%%%
    % generate LDA vector
    lambdaSplit = myLDA(data.domain,data.codomain);    
    fprintf(['building model@map.normalize.......\n']);
    % flip to point towards ones
    grade = data.domain*lambdaSplit;
    if (grade'*data.codomain < 1)
        lambdaSplit = -lambdaSplit;
    end
    
    % generate threshold
    % regrade
    grade = data.domain*lambdaSplit;
    idx0 = data.codomain==0;
    idx1 = data.codomain==1;
    u1 = mean(grade(idx1));
    s1 = std(grade(idx1));
    u0 = mean(grade(idx0));
    s0 = std(grade(idx0));
    X = linspace(u0,u1,100);
    pdf0 = normpdf(X,u0,s0);
    pdf1 = normpdf(X,u1,s1);
    [JUNK vidx] = min(abs(pdf0 - pdf1));
    % store LDA vector
    model(1).para.lambdaSplit.value = lambdaSplit;
    model(1).para.threshold.value = X(vidx);
    model(1).name = 'lda split';
    model(1).type = 'lda';
    %}
    
    
    %%%%%%%%%%%%%%%%%
    % generate vars for pdf model
    idx = find(data.codomain);
    u = mean(data.domain(idx,:),1);
    s = std(data.domain(idx,:),1,1);
    model(end+1).para.mv.value = [u s];
    model(end).name = 'maxi-LLhood';
    model(end).type = 'pdf';
    model(end).para.max.value = probData(u,model(end).para.mv.value);
    model(end).para.threshold.value = [];
    %%%%%%%%%%%%%%%%%
    % generate threshold
    tempGrade = eval_learningModel(data.domain,model(1));    
    grade = tempGrade.prob;
    idx0 = data.codomain==0;
    idx1 = data.codomain==1;
    u1 = mean(grade(idx1));
    s1 = std(grade(idx1));
    u0 = mean(grade(idx0));
    s0 = std(grade(idx0));
    X = linspace(u0,u1,100);
    pdf0 = normpdf(X,u0,s0);
    pdf1 = normpdf(X,u1,s1);
    [JUNK vidx] = min(abs(pdf0 - pdf1));
    %%%%%%%%%%%%%%%%%
    % store theshold
    model(1).para.threshold.value = X(vidx);
    
end