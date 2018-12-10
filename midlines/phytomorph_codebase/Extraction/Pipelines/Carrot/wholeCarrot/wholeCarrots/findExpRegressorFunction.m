function [C,beta] = findExpRegressorFunction(P,rX,rY)

   
  
    %{
    idx = find(abs(P) < 1);
    P(idx) = 1;
%}
    %tic;
    sz = numel(P)^.5;
    P = reshape(P,[sz sz]);
    pX = [];
    for p = 1:size(P,1)
        pX = [pX prod(bsxfun(@power,rX,P(p,:)),2)];
    end
    
    
    [beta,C] = findRegressor((pX),rY,min(5,size(pX,2)));
    C = -(C);
    
    
end


