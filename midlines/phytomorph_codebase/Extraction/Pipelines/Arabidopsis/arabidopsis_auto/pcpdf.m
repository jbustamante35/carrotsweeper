function [p] = pcpdf(PC,P)
    PC = bsxfun(@minus,PC,P);
    p = mvnpdf(PC,[0 0 0],100*diag(ones(1,3)));
end