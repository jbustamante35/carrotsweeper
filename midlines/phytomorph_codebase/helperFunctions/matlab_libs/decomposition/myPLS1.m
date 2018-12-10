function [wx,wy] = myPLS1(X,Y,k)
    Sxy = X'*Y;
    [wx,S,wy] = svd(Sxy,0);
    
    pX = X*wx;
    wx = bsxfun(@times,wx,sum(pX.*pX,1).^-.5); % scale for unit norm
    
    
    pY = Y*wy;
    wy = bsxfun(@times,wy,sum(pY.*pY,1).^-.5);
end