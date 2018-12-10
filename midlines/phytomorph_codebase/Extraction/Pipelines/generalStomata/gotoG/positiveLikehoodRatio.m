function [m] = positiveLikehoodRatio(targets,outputs)
    [tpr,fpr] = roc(double(targets(:,1))',double(outputs)');
    %C = confusionmat(double(targets(:,1))',double(outputs)');
    %tpr = C(1,1)*sum(C(1,:)).^-1;
    %fpr = C(2,1)*sum(C(2,:)).^-1;
    m = tpr(2)*fpr(2)^-1;
end