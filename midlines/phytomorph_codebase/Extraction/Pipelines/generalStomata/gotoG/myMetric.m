function [m] = myMetric(targets,outputs)
    C = confusionmat(logical(targets(:,1)),logical(outputs(:,1)),'order',logical([0 1]));
    precision = C(1,1)*sum(C(:,1)).^-1;
    tpr = C(1,1)*sum(C(1,:)).^-1;
    m = 2*((precision*tpr)*(precision + tpr)^-1);
end