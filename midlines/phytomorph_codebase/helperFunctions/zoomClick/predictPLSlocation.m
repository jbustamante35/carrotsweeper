function [DELTA] = predictPLSlocation(X,beta,SLIDE)
    X1 = [ones(size(X,1),1) X]*beta(:,:,1);
    X2 = [ones(size(X,1),1) X]*beta(:,:,2);
    fidx1 = find(X1 > .5);
    fidx2 = find(X2 > .5);
    if isempty(fidx1)
        fidx1 = 1;
    end
    if isempty(fidx2)
        fidx2 = 1;
    end
    DELTA = [SLIDE(fidx1(end)),SLIDE(fidx2(end))];
end