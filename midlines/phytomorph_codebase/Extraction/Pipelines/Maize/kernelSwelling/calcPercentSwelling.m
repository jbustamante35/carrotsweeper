function [S] = calcPercentSwelling(A,numInit)
    S = zeros(size(A));
    if nargin == 1
        numInit = 1;
    end
    if numel(A) >= numInit
        initA = mean(A(:,1:numInit),2);
        A = bsxfun(@minus,A,initA);
        S = bsxfun(@times,A,initA.^-1);    
    end
end