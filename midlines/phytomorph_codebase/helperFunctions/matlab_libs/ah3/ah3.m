function [dist svec] = ah3(model,data,target,groups)
    TYPE  = 'prob';

    % for each model
    for mod = 1:size(model,1)
        switch TYPE
            case 'prob'
                % calculate the likehood of of the data
                p(:,mod) = gprobData(data,model(mod,:));
            case 'log'
                % calculate the -log likehood of of the data
                p(:,mod) = LOGgprobData(data,model(mod,:));
        end
        % integrate the number of objects in each group
        if (nargout == 1)
            attempt(mod,:) = grpIntegration(groups,p(:,mod),target);
        else
            [attempt(mod,:) g]  = grpIntegration(groups,p(:,mod),target);
        end
        % calculate the error in counting
        delta(mod,:) = abs(target-attempt(mod,:));
        % take the norm as measure of error
        dist(mod) = norm(delta(mod,:));
    end
    % find the best model
    [guess midx] = min(dist);
    fprintf(['Best guess is @' num2str(guess) '\n'])
end