function [wx,wy,pX,pY,D] = myCCA(X,Y,k)

    Sxy = X'*Y;
    Syx = Y'*X;
    Sxx = X'*X;
    Syy = Y'*Y;
    LHS = [[zeros(size(Sxy,1),size(Syx,2)),Sxy];[Syx,zeros(size(Syx,1),size(Sxy,2))]];
    RHS = [[Sxx,zeros(size(Sxx,1),size(Syy,2))];[zeros(size(Syy,1),size(Sxx,2)),Syy]];
    [V,D] = eigs(LHS,RHS,2*k);
    
    
    [J sidx] = sort(diag(D),'descend');
    %{
    sgn = sign(diag(D));
    for i = 1:size(V,2)
        V(:,i) = -sgn(i)*V(:,i);
    end
    %}
    
    V = V(:,sidx);
    
    
    V = V(:,1:k);
    wx = V(1:size(X,2),:);
    wy = V(size(X,2)+1:end,:);
    %{
    for e = 1:size(wx,2)
        wx(:,e) = wx(:,e)/norm(wx(:,e));
        wy(:,e) = wy(:,e)/norm(wy(:,e));
    end
    %}
    
    pX = X*wx;
    wx = bsxfun(@times,wx,std(pX,1,1).^-1); % scale for unit norm    
    pX = X*wx;
    pY = Y*wy;
    wy = bsxfun(@times,wy,std(pY,0,1).^-1);
    pY = Y*wy;
end