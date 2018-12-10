function [Yp] = myNonPre(rX,P,beta)
    %tic;
    sz = numel(P)^.5;
    P = reshape(P,[sz sz]);
    pX = [];
    for p = 1:size(P,1)
        pX = [pX prod(bsxfun(@power,rX,P(p,:)),2)];
    end
    % predict
    Yp = [ones(size(pX,1),1) pX]*beta;
end