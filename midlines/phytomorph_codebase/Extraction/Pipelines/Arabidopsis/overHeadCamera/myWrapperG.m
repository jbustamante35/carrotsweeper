function [Y] = myWrapperG(data,mu)
    nd = numel(mu)/2;
    for k = 1:nd
        [dis] = myWrapperDistance(data,mu(k),1);
        Y(:,k) = normpdf(dis,mu(k),mu(k+nd));
    end
    Y = -sum(log(max(Y,[],2)));
end