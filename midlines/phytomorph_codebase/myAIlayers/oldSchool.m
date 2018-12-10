function [] = oldSchool(X,V,P)

    R = zeros(size(tmp,1),size(V,2)^2);
    for e = 1:size(tmp,1)
        w = tmp(e,:)'*tmp(e,:);
        R(e,:) = w(:)';
    end
    r = mtimesx(P,'T',R);
end