function [wx,wy,sx,sy,PRE] = myPLS_EZ(X,Y,k)
    Sxy = X'*Y;
    Syx = Y'*X;
    MT = [zeros(size(Sxy,1),size(Syx,2)) Sxy];
    MB = [Syx zeros(size(Syx,1),size(Sxy,2))];
    M = [MT;MB];
    [V,D] = eigs(M,2*k);
    [J sidx] = sort(diag(D),'descend');
    V = V(:,sidx);
    
    V = V(:,1:k);
    wx = V(1:size(X,2),:);
    wy = V(size(X,2)+1:end,:);
    %wx = V(1:size(X,2),:);
    %wy = V(size(X,2)+1:end,:);
    nx = sum(wx.*wx,1).^.5;
    ny = sum(wy.*wy,1).^.5;
    wx = bsxfun(@times,wx,nx.^-1);
    wy = bsxfun(@times,wy,ny.^-1);
    sx = X*wx;
    sy = Y*wy;
    
    %{
    sgn = sign(diag(D));
    for i = 1:size(V,2)
        V(:,i) = -sgn(i)*V(:,i);
    end
    %}
    
    for e = 1:size(sx,2)
        p = std(sy(:,e))*std(sx(:,e)).^-1;
        PRE(e) = p(1);
    end
    %{
    for e = 1:size(sx,2)
        p =  polyfit(sx(:,e),sy(:,e),1);
        PRE(e) = p(1);
    end
    %}
end