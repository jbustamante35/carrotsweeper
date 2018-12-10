function [X,T] = position_flf(X,p,dT,frac)
    v = flf(X,p);
    cnt = 1;
    T(cnt) = 0;
    while v < frac*max(p(:,1))
        X(:,cnt+1) = X(:,cnt) + dT*v;
        T(cnt+1) = T(cnt)+dT;
        v = flf(X(:,cnt+1),p);
        cnt = cnt + 1;
    end
    

end