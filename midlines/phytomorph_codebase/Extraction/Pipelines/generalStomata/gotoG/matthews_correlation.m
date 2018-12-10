function [m] = matthews_correlation(targets,outputs)
    C = confusionmat(logical(targets(:,1)),logical(outputs(:,1)),'order',logical([0 1]));
    TOP = prod(diag(C)) - prod(diag(imrotate(C,90)));
    BOTTOM = prod([sum(C,1) sum(C,2)']).^.5;
    m = TOP.*BOTTOM.^-1;
end